"""Batch docking on GPU with Uni-Dock.

The CPU path starts one Vina process per ligand and pays process startup every time. Uni-Dock
reads a *list* of ligands and docks the whole list in one call, which is the entire point of
this module: the GPU stays busy across the batch instead of being set up and torn down per
peptide. With several GPUs the batch is split into one shard per device and the shards run at
the same time.

This module is GPU-only on purpose. There is no CPU fallback anywhere in it: a missing GPU, a
missing `unidock` binary, or a Uni-Dock call that exits non-zero raises. Falling back quietly
would mean a run that was asked for the GPU path silently produces CPU numbers hours later,
and inventing a score for a ligand whose output could not be read would put a number in the
ranking that no docking ever produced. Both failures look like data, which is the one kind of
failure a screening run cannot absorb.

    from unidock_runner import dock_batch
    scores = dock_batch(receptor_qt, ligand_paths, "dock_out", center, size)
    print(f"{len(scores)} scored, {len(scores.missing)} without a usable score")
"""

import os
import shlex
import subprocess
import sys
from collections import Counter
from concurrent.futures import ThreadPoolExecutor

# gpu_backends is imported as a module rather than with `from gpu_backends import ...` so every
# call below is a fresh attribute lookup: that is the one place the GPU gate lives, and the tests
# replace those attributes instead of a copy bound at import time.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gpu_backends  # noqa: E402

# Uni-Dock names each output file after its input ligand and drops it in --dir. Reading the
# scores back means rebuilding those names, so the rule is written down once here.
OUTPUT_SUFFIX = "_out.pdbqt"

# Both Vina and Uni-Dock write the pose affinity on this remark line.
RESULT_MARKER = "REMARK VINA RESULT:"

# Uni-Dock's three presets, slowest last. Rejecting anything else early matters because a typo
# here would otherwise surface as a failed call after the receptor has already been prepared.
SEARCH_MODES = ("fast", "balance", "detail")


class BatchScores(dict):
    """{ligand stem: best affinity in kcal/mol}, plus the ligands that produced no usable score.

    A plain dict would hide the interesting half of a bad batch: the caller gets 48 scores for
    50 ligands and nothing says which two vanished or why. So the mapping carries:

        missing   (ligand stem, reason) pairs for every ligand with no usable score, i.e. no
                  output pdbqt in the output dir, or one with no parseable REMARK VINA RESULT
        commands  the Uni-Dock command lines that were run, one per shard, for the run log
        complete  True when nothing was dropped

    `dock_batch` also prints every dropped ligand, so a silent drop is impossible: the caller
    has to ignore both a warning on stdout and a non-empty `missing` to miss it.
    """

    def __init__(self, scores=None, missing=None, commands=None):
        super().__init__(scores or {})
        self.missing = list(missing or [])
        self.commands = [list(command) for command in (commands or [])]

    @property
    def complete(self):
        """True when every submitted ligand came back with a score."""
        return not self.missing


def ligand_stem(ligand_path):
    """Basename without extension. This is the key used in the returned scores."""
    return os.path.splitext(os.path.basename(ligand_path))[0]


def output_path_for(ligand_path, output_dir):
    """Where Uni-Dock will write this ligand's docked poses."""
    return os.path.join(output_dir, ligand_stem(ligand_path) + OUTPUT_SUFFIX)


def write_ligand_index(ligand_paths, index_path):
    """Write Uni-Dock's --ligand_index file, one ligand path per line. Returns index_path.

    An empty ligand list raises instead of writing an empty index: Uni-Dock accepts an empty
    index, exits 0 and docks nothing, which reads downstream like a batch where every ligand
    failed rather than like the caller mistake it is.

    Paths are written exactly as given. Uni-Dock resolves them against its own working
    directory, so the pipeline should pass absolute paths -- the docking stages chdir between
    folders, and a relative path in here is read from wherever Uni-Dock happens to start.
    """
    ligand_paths = list(ligand_paths)
    if not ligand_paths:
        raise ValueError(
            "write_ligand_index needs at least one ligand path, got an empty list. Refusing to "
            "write an empty ligand index: Uni-Dock would exit 0 having docked nothing.")
    parent = os.path.dirname(os.path.abspath(index_path))
    os.makedirs(parent, exist_ok=True)
    with open(index_path, "w") as handle:
        handle.write("\n".join(ligand_paths) + "\n")
    return index_path


def build_command(receptor_qt, ligand_index, output_dir, center, size, *, binary,
                  device_index=None, search_mode="balance"):
    """Build the Uni-Dock command line for one batch of ligands.

    `center` and `size` are the gridbox the pipeline already computed from the pocket atoms and
    are passed through unchanged. The GPU path has to search exactly the box the CPU path
    searched; a box that moves makes the new scores incomparable with every result already on
    disk, and the discrepancy would be read as a difference between the docking engines.

    `device_index`, when given, becomes --device_id, which is how one shard is pinned to one GPU.
    """
    if search_mode not in SEARCH_MODES:
        raise ValueError(
            f"search_mode must be one of {', '.join(SEARCH_MODES)}, got {search_mode!r}")
    command = [
        binary,
        "--receptor", receptor_qt,
        "--ligand_index", ligand_index,
        "--dir", output_dir,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--search_mode", search_mode,
    ]
    if device_index is not None:
        command += ["--device_id", str(device_index)]
    return command


def parse_score_from_pdbqt(text):
    """Affinity of the first pose in a Uni-Dock output pdbqt, or None when there is none.

    The poses are written best first, so the first `REMARK VINA RESULT:` line is the number the
    ranking wants. If that first line is present but does not parse, this returns None rather
    than reading on to the next pose: the next pose is a worse one, and reporting it as the best
    would quietly understate the ligand. None sends the ligand to `BatchScores.missing`, where
    it can be seen, which is the better outcome.
    """
    for line in (text or "").splitlines():
        stripped = line.strip()
        if not stripped.startswith(RESULT_MARKER):
            continue
        fields = stripped[len(RESULT_MARKER):].split()
        if not fields:
            return None
        try:
            return float(fields[0])
        except ValueError:
            return None
    return None


def shard(items, buckets):
    """Split `items` into `buckets` contiguous, as-even-as-possible lists.

    Contiguous rather than round-robin so a shard's ligand index reads in the same order as the
    input, which makes a failing shard easy to reproduce by hand. The result always has exactly
    `buckets` lists, so a shard lines up with the device at the same position; the trailing ones
    are empty only when there are fewer items than buckets.
    """
    if buckets < 1:
        raise ValueError(f"buckets must be at least 1, got {buckets}")
    items = list(items)
    base, remainder = divmod(len(items), buckets)
    shards = []
    start = 0
    for position in range(buckets):
        # The first `remainder` buckets take one extra item, which keeps the largest and
        # smallest shard within one item of each other.
        length = base + (1 if position < remainder else 0)
        shards.append(items[start:start + length])
        start += length
    return shards


def _run_shards(commands, runner=None):
    """Run one Uni-Dock call per shard, concurrently, and raise on the first failure.

    Threads, not processes: each call is a blocking subprocess pinned to its own --device_id, so
    the GIL is released while it runs and the GPUs really do work at the same time. Sharding
    across devices and then running the shards one after another would leave every GPU but one
    idle, which would give up the reason for sharding in the first place.
    """
    if runner is None:
        runner = subprocess.run

    def call(command):
        try:
            return runner(command, check=False, capture_output=True, text=True)
        except OSError as exc:
            raise RuntimeError(
                f"Uni-Dock could not be executed ({exc}). Command:\n  {shlex.join(command)}"
            ) from exc

    with ThreadPoolExecutor(max_workers=len(commands)) as pool:
        futures = [pool.submit(call, command) for command in commands]
    # Leaving the pool has joined every shard, so no Uni-Dock process is still writing into the
    # output dir by the time a failure is raised below and the outputs are read.
    for command, future in zip(commands, futures):
        result = future.result()
        if result.returncode != 0:
            # Uni-Dock reports some failures (bad gridbox, unreadable receptor) on stdout, so
            # fall back to it: an error message with nothing in it is not actionable.
            detail = (result.stderr or "").strip() or (result.stdout or "").strip()
            raise RuntimeError(
                f"Uni-Dock exited {result.returncode}. Command:\n  {shlex.join(command)}\n"
                f"stderr:\n{detail or '(no output captured)'}")


def _read_text(path):
    """File contents, or None when it is not there or cannot be read."""
    try:
        with open(path) as handle:
            return handle.read()
    except OSError:
        return None


def dock_batch(receptor_qt, ligand_paths, output_dir, center, size, *, env=None, runner=None,
               devices=None, search_mode="balance"):
    """Dock a whole batch of ligands on GPU and return their best affinities.

    Requires a GPU and the `unidock` binary through gpu_backends; GpuUnavailableError propagates
    untouched, because this path has no CPU fallback by design. The ligands are sharded across
    the usable devices, one ligand index and one Uni-Dock call per shard, and each shard's
    output pdbqt is then read back from `output_dir`.

    Returns a `BatchScores` -- a dict of {ligand stem: best affinity} that also carries
    `missing`, the (stem, reason) pairs for ligands with no usable score, and `commands`, the
    calls that were made. Ligands in `missing` are absent from the mapping and are also printed
    as a warning, so a dropped ligand cannot pass unnoticed.

    Raises ValueError for an empty ligand list or ligands that share a basename, and RuntimeError
    when a Uni-Dock call fails, carrying the command and what it printed.

    `runner` replaces subprocess.run and `devices` the detected device list; both exist so this
    can be tested without a GPU.
    """
    ligand_paths = list(ligand_paths)
    if not ligand_paths:
        raise ValueError("dock_batch needs at least one ligand, got an empty list.")

    # Uni-Dock names its outputs after the ligand basename, so two ligands with the same
    # basename overwrite each other in output_dir and one score would be reported for both.
    duplicates = sorted(
        stem for stem, count in Counter(ligand_stem(path) for path in ligand_paths).items()
        if count > 1)
    if duplicates:
        raise ValueError(
            "Ligands must have distinct basenames, because Uni-Dock names its output files "
            f"after them and they all land in {output_dir}. Repeated: {', '.join(duplicates)}")

    # The GPU gate. Every call here raises GpuUnavailableError on a machine that cannot run this,
    # and that is deliberately not caught: the run stops instead of drifting onto the CPU path.
    if devices is None:
        # visible_device_indices already applies the require_gpus check and then narrows it to
        # what CUDA_VISIBLE_DEVICES allows, so asking it twice would only re-run nvidia-smi.
        devices = gpu_backends.visible_device_indices(env=env, runner=runner)
    else:
        # A caller-supplied device list still has to clear the gate: GPU-only means the hardware
        # is verified even when the caller names the devices itself.
        gpu_backends.require_gpus(runner=runner)
    binary = gpu_backends.require_binary("unidock", env=env)
    devices = list(devices)
    if not devices:
        raise ValueError("devices is an empty list; pass None to detect the usable GPUs.")

    os.makedirs(output_dir, exist_ok=True)

    # Drop any output left by an earlier run of these same ligands. Uni-Dock overwrites what it
    # does produce, but a ligand that fails this time would otherwise be scored from the file the
    # previous run left behind: a real number, from a docking nobody asked about, indistinguishable
    # from a fresh one.
    for ligand_path in ligand_paths:
        stale = output_path_for(ligand_path, output_dir)
        if os.path.exists(stale):
            os.remove(stale)

    # One shard per device, but never more shards than ligands: an extra shard would be an empty
    # ligand index, i.e. a Uni-Dock call that docks nothing.
    shards = shard(ligand_paths, min(len(devices), len(ligand_paths)))

    commands = []
    for position, shard_ligands in enumerate(shards):
        index_path = os.path.join(output_dir, f"ligand_index_shard{position}.txt")
        write_ligand_index(shard_ligands, index_path)
        commands.append(build_command(
            receptor_qt, index_path, output_dir, center, size, binary=binary,
            device_index=devices[position], search_mode=search_mode))

    _run_shards(commands, runner=runner)

    scores = {}
    missing = []
    for ligand_path in ligand_paths:
        stem = ligand_stem(ligand_path)
        out_path = output_path_for(ligand_path, output_dir)
        text = _read_text(out_path)
        if text is None:
            missing.append((stem, f"no output pdbqt at {out_path}"))
            continue
        score = parse_score_from_pdbqt(text)
        if score is None:
            missing.append((stem, f"no parseable '{RESULT_MARKER}' in {out_path}"))
            continue
        scores[stem] = score

    result = BatchScores(scores, missing=missing, commands=commands)
    if missing:
        # Printed as well as returned: a drop that only lives on an attribute is a drop nobody
        # notices until the ranking comes out short, long after the GPU is free.
        print(f"WARNING: Uni-Dock produced no usable score for {len(missing)} of "
              f"{len(ligand_paths)} ligands:")
        for stem, why in missing:
            print(f"  {stem}: {why}")
    return result
