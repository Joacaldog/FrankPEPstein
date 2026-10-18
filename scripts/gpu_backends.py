"""GPU backends for the docking and rescoring stages, opt-in with a CPU fallback.

The pipeline's default behaviour does not change: every backend here is off unless an
environment variable asks for it, and every backend degrades to the CPU path it replaces
instead of failing the run.

Two environment variables drive it:

    FRANKPEP_DOCK_BACKEND     vina (default) | unidock
    FRANKPEP_RESCORE_BACKEND  vina (default) | gnina | consensus

`consensus` keeps the Vina score and adds GNINA's CNN score as an extra column, which is the
conservative way to try a new scoring function: nothing is discarded on its word alone.

Hardware note, because it decides what is worth enabling on which machine: Uni-Dock and GNINA
run on Pascal cards (GTX 1080 / 1080 Ti, compute capability 6.1), so they need no new hardware.
Co-folding models such as Boltz-2 need tensor cores and bfloat16, which means compute capability
8.0 or newer; `describe_environment()` reports both so a run can say which path it can take.
"""

import os
import shutil
import subprocess

DEFAULT_DOCK_BACKEND = "vina"
DEFAULT_RESCORE_BACKEND = "vina"

DOCK_BACKENDS = ("vina", "unidock")
RESCORE_BACKENDS = ("vina", "gnina", "consensus")

# Compute capability at or above which bfloat16 and tensor cores are available, i.e. the floor
# for the co-folding models (Boltz-2, Chai-1). Pascal sits below it; Ampere and newer above.
MODERN_COMPUTE_CAPABILITY = 8.0

_NVIDIA_SMI_QUERY = "index,name,compute_cap,memory.total"


class Gpu:
    """One GPU as nvidia-smi reports it."""

    def __init__(self, index, name, compute_cap, memory_mib):
        self.index = index
        self.name = name
        self.compute_cap = compute_cap
        self.memory_mib = memory_mib

    @property
    def supports_modern_models(self):
        """True when this card can run bfloat16 transformer models such as Boltz-2."""
        return self.compute_cap is not None and self.compute_cap >= MODERN_COMPUTE_CAPABILITY

    def __repr__(self):
        cap = "unknown" if self.compute_cap is None else f"{self.compute_cap:g}"
        return f"Gpu(index={self.index}, name={self.name!r}, compute_cap={cap}, memory_mib={self.memory_mib})"


def _run_nvidia_smi(runner=None):
    """Return nvidia-smi's CSV output, or None when the tool is absent or fails."""
    if runner is None:
        # Only the real path needs the lookup; an injected runner stands in for the tool.
        if shutil.which("nvidia-smi") is None:
            return None
        runner = subprocess.run
    try:
        result = runner(
            ["nvidia-smi", f"--query-gpu={_NVIDIA_SMI_QUERY}",
             "--format=csv,noheader,nounits"],
            check=False, capture_output=True, text=True)
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout


def parse_nvidia_smi(output):
    """Parse the CSV rows of `nvidia-smi --query-gpu` into Gpu objects.

    Kept separate from the call so it can be tested without a GPU present. Rows that do not
    parse are skipped rather than raising: a partially readable nvidia-smi is still useful.
    """
    gpus = []
    for line in (output or "").splitlines():
        fields = [field.strip() for field in line.split(",")]
        if len(fields) < 4:
            continue
        try:
            index = int(fields[0])
        except ValueError:
            continue
        try:
            compute_cap = float(fields[2])
        except ValueError:
            compute_cap = None
        try:
            memory_mib = int(float(fields[3]))
        except ValueError:
            memory_mib = None
        gpus.append(Gpu(index=index, name=fields[1], compute_cap=compute_cap,
                        memory_mib=memory_mib))
    return gpus


def detect_gpus(runner=None):
    """List the GPUs visible to this process. Empty list when there are none."""
    return parse_nvidia_smi(_run_nvidia_smi(runner=runner))


def _binary_for(backend, env=None):
    """Path to a backend's executable, honouring an explicit override.

    FRANKPEP_UNIDOCK_PATH / FRANKPEP_GNINA_PATH let a machine point at a binary that is not on
    PATH, which is the normal case for a shared lab machine.
    """
    env = os.environ if env is None else env
    override = env.get(f"FRANKPEP_{backend.upper()}_PATH")
    if override:
        return override if os.path.exists(override) else None
    return shutil.which(backend)


def resolve_dock_backend(env=None, runner=None):
    """Decide which docking backend to use.

    Returns (name, reason). `name` is always a backend that can actually run here, so callers
    never have to handle a missing binary; `reason` says why, for the run log.
    """
    env = os.environ if env is None else env
    requested = (env.get("FRANKPEP_DOCK_BACKEND") or DEFAULT_DOCK_BACKEND).strip().lower()

    if requested == DEFAULT_DOCK_BACKEND:
        return DEFAULT_DOCK_BACKEND, "default CPU Vina"
    if requested not in DOCK_BACKENDS:
        return (DEFAULT_DOCK_BACKEND,
                f"unknown backend {requested!r}, falling back to {DEFAULT_DOCK_BACKEND}")
    if not detect_gpus(runner=runner):
        return DEFAULT_DOCK_BACKEND, f"{requested} requested but no GPU is visible"
    if _binary_for(requested, env=env) is None:
        return DEFAULT_DOCK_BACKEND, f"{requested} requested but its binary was not found"
    return requested, f"{requested} on GPU"


def resolve_rescore_backend(env=None, runner=None):
    """Decide which rescoring backend to use. Same contract as resolve_dock_backend."""
    env = os.environ if env is None else env
    requested = (env.get("FRANKPEP_RESCORE_BACKEND") or DEFAULT_RESCORE_BACKEND).strip().lower()

    if requested == DEFAULT_RESCORE_BACKEND:
        return DEFAULT_RESCORE_BACKEND, "default Vina scoring"
    if requested not in RESCORE_BACKENDS:
        return (DEFAULT_RESCORE_BACKEND,
                f"unknown backend {requested!r}, falling back to {DEFAULT_RESCORE_BACKEND}")
    if not detect_gpus(runner=runner):
        return DEFAULT_RESCORE_BACKEND, f"{requested} requested but no GPU is visible"
    if _binary_for("gnina", env=env) is None:
        return DEFAULT_RESCORE_BACKEND, f"{requested} requested but the gnina binary was not found"
    return requested, f"{requested} using gnina on GPU"


def gnina_score_command(receptor_qt, ligand, env=None):
    """Command that rescores an existing pose with GNINA's CNN, without re-docking it.

    --score_only mirrors what the Vina rescore does: the pose is already optimised, what is
    being asked is how good it looks to a different scoring function.
    """
    binary = _binary_for("gnina", env=env) or "gnina"
    return [binary, "--score_only", "--receptor", receptor_qt, "--ligand", ligand]


def parse_gnina_score(stdout):
    """Pull the CNN affinity out of GNINA's output. None when it is not there.

    GNINA prints both `CNNscore` (a pose-quality probability) and `CNNaffinity` (a pKd-like
    number). The affinity is the one comparable in spirit to Vina's estimate, so that is what
    the pipeline records; higher is better, which is the opposite of Vina's convention.
    """
    for line in (stdout or "").splitlines():
        stripped = line.strip()
        if stripped.startswith("CNNaffinity"):
            parts = stripped.replace(":", " ").split()
            for part in reversed(parts):
                try:
                    return float(part)
                except ValueError:
                    continue
    return None


def unidock_dock_command(receptor_qt, ligand_index_file, output_dir, center, size, env=None):
    """Command that docks a whole batch of ligands in one Uni-Dock call.

    The batch is the point: Vina pays process startup per ligand, Uni-Dock reads a list of
    ligands and keeps the GPU busy across all of them. `center` and `size` are the gridbox the
    pipeline already computes, passed through unchanged so the search space does not move.
    """
    binary = _binary_for("unidock", env=env) or "unidock"
    return [
        binary,
        "--receptor", receptor_qt,
        "--ligand_index", ligand_index_file,
        "--dir", output_dir,
        "--center_x", str(center[0]),
        "--center_y", str(center[1]),
        "--center_z", str(center[2]),
        "--size_x", str(size[0]),
        "--size_y", str(size[1]),
        "--size_z", str(size[2]),
        "--search_mode", "balance",
    ]


def describe_environment(env=None, runner=None):
    """Lines describing what this machine can run. Printed at the start of a stage."""
    gpus = detect_gpus(runner=runner)
    lines = []
    if not gpus:
        lines.append("GPU: none visible; every stage runs on CPU.")
    else:
        for gpu in gpus:
            cap = "unknown" if gpu.compute_cap is None else f"{gpu.compute_cap:g}"
            memory = "unknown" if gpu.memory_mib is None else f"{gpu.memory_mib} MiB"
            suffix = ("tensor cores available, can run co-folding models"
                      if gpu.supports_modern_models
                      else "no tensor cores, docking and rescoring only")
            lines.append(f"GPU {gpu.index}: {gpu.name}, compute {cap}, {memory} ({suffix})")

    dock, dock_reason = resolve_dock_backend(env=env, runner=runner)
    rescore, rescore_reason = resolve_rescore_backend(env=env, runner=runner)
    lines.append(f"Docking backend: {dock} ({dock_reason})")
    lines.append(f"Rescoring backend: {rescore} ({rescore_reason})")
    return lines
