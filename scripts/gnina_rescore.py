"""Rescore existing poses with GNINA's CNN scoring function, on GPU only.

This is a second opinion on poses another stage already produced, not a second dock: every
call is `--score_only`, so the pose is left exactly where the refinement put it and the only
question asked is how it looks to a different scoring function.

GPU only, with no CPU fallback. A missing GPU, a missing `gnina` binary, or a failed GNINA call
stops the stage: the GpuUnavailableError that `gpu_backends` raises is left to propagate, and a
failed call raises RuntimeError carrying the command and GNINA's stderr. Nothing here invents a
score or quietly reports a worse one, because a rescoring step that degrades silently is
indistinguishable from one that worked.

The three numbers GNINA prints disagree about which direction is good, which is the whole
reason this module parses them into named fields instead of a bare float:

    CNNaffinity  pKd-like, HIGHER is better (the opposite of Vina's kcal/mol)
    CNNscore     pose-quality probability in [0, 1], not an affinity at all
    Affinity     GNINA's own Vina-style energy in kcal/mol, LOWER is better

Usage:

    import gnina_rescore

    scores = gnina_rescore.rescore_poses("receptor.pdbqt", poses)
    for pose, pose_scores in scores.items():
        print(pose, pose_scores.cnn_affinity)

With `on_error="skip"` a pose whose GNINA call failed is left out of the returned mapping and
recorded in its `.skipped` attribute, so a caller can still see what was dropped and why:

    scores = gnina_rescore.rescore_poses("receptor.pdbqt", poses, on_error="skip")
    for dropped in scores.skipped:
        print(dropped.describe())
"""

import math
import os
import re
import shlex
import subprocess
import sys

# Same sibling import the other stage scripts use, so this module works whatever the cwd is.
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import gpu_backends  # noqa: E402

ON_ERROR_MODES = ("raise", "skip")

# One pattern per field, anchored at the start of the (stripped) line. The anchor is what keeps
# `CNNaffinity:` out of the `Affinity:` field, and the negative lookahead keeps GNINA's
# `CNNaffinity_variance:` line from being mistaken for the affinity itself.
_SCORE_PATTERNS = (
    ("cnn_affinity", re.compile(r"^CNNaffinity(?!\w)\s*[:=]?\s*(\S+)", re.IGNORECASE)),
    ("cnn_score", re.compile(r"^CNNscore(?!\w)\s*[:=]?\s*(\S+)", re.IGNORECASE)),
    ("affinity", re.compile(r"^Affinity(?!\w)\s*[:=]?\s*(\S+)", re.IGNORECASE)),
)


class PoseScores:
    """What GNINA reported for one pose. A field GNINA did not report stays None.

    None means "not reported", never zero: a 0.0 CNNaffinity is a real, very weak prediction,
    so filling a missing value in with it would read as a confident answer.
    """

    def __init__(self, cnn_affinity=None, cnn_score=None, affinity=None):
        # pKd-like CNN affinity: HIGHER IS BETTER, the opposite of the Vina columns.
        self.cnn_affinity = cnn_affinity
        # CNN pose-quality probability in [0, 1]: higher is a more plausible pose. NOT an affinity.
        self.cnn_score = cnn_score
        # GNINA's own Vina-style binding energy in kcal/mol: LOWER IS BETTER.
        self.affinity = affinity

    def __eq__(self, other):
        if not isinstance(other, PoseScores):
            return NotImplemented
        return (self.cnn_affinity, self.cnn_score, self.affinity) == (
            other.cnn_affinity, other.cnn_score, other.affinity)

    def __repr__(self):
        return (f"PoseScores(cnn_affinity={self.cnn_affinity}, cnn_score={self.cnn_score}, "
                f"affinity={self.affinity})")


class SkippedPose:
    """A pose whose GNINA call failed and that `on_error="skip"` dropped.

    `returncode` is None when GNINA could not be launched at all; `stderr` is whatever it wrote
    before giving up, which is where GNINA says things like which CUDA device it could not open.
    """

    def __init__(self, pose_path, command, returncode, stderr):
        self.pose_path = pose_path
        self.command = list(command)
        self.returncode = returncode
        self.stderr = stderr

    def describe(self):
        """One actionable message: which pose failed, how, what GNINA said, and the command."""
        status = "could not be launched" if self.returncode is None \
            else f"exited with code {self.returncode}"
        stderr = (self.stderr or "").strip() or "<no stderr captured>"
        return (f"GNINA rescoring failed for {self.pose_path}: {status}.\n"
                f"command: {shlex.join(self.command)}\n"
                f"stderr: {stderr}")

    def __repr__(self):
        return (f"SkippedPose(pose_path={self.pose_path!r}, returncode={self.returncode!r})")


class RescoreResults(dict):
    """`{pose_path: PoseScores}`, plus the poses that were dropped.

    A plain dict for every caller that only wants the scores, with one extra attribute so a
    dropped pose cannot pass unnoticed: `.skipped` holds a SkippedPose per failed call under
    `on_error="skip"`, and is empty under `on_error="raise"`, where a failure stops the run
    instead. Note that dict.copy() returns a plain dict and loses the attribute.
    """

    def __init__(self, *args, **kwargs):
        super().__init__(*args, **kwargs)
        self.skipped = []


def _finite_float(token):
    """Parse one token into a float, or None when it is not a usable number.

    "nan" and "inf" parse as floats but are not scores, and NaN in particular compares False
    against everything, so it would sort to an arbitrary end of a ranking instead of being
    reported as the missing value it is.
    """
    try:
        value = float(token)
    except (TypeError, ValueError):
        return None
    return value if math.isfinite(value) else None


def build_command(receptor, ligand, *, binary, device_index=None):
    """The GNINA command that scores `ligand` against `receptor` without moving it.

    `--score_only` is the whole point: the pose comes from a stage that already optimised it,
    so re-docking it here would answer a different question and discard that work.
    """
    command = [binary, "--score_only", "--receptor", receptor, "--ligand", ligand]
    if device_index is not None:
        # Device 0 is a real device, so this tests for None rather than for truthiness.
        command += ["--device", str(device_index)]
    return command


def parse_scores(stdout):
    """Pull CNNaffinity, CNNscore and Affinity out of GNINA's stdout.

    Tolerant by design: GNINA surrounds its numbers with a banner, warnings and timing lines,
    and the fields moved between releases. A field that is absent, non-numeric or NaN stays
    None rather than being guessed at.

    The first value reported for each field wins. A multi-model ligand file makes GNINA print
    one block per pose, and this module scores one pose per call, so the first block is the one
    being asked about; later blocks are ignored rather than silently overwriting it.
    """
    values = {name: None for name, _ in _SCORE_PATTERNS}
    seen = set()
    for raw_line in (stdout or "").splitlines():
        line = raw_line.strip()
        for name, pattern in _SCORE_PATTERNS:
            if name in seen:
                continue
            match = pattern.match(line)
            if match is None:
                continue
            # Marked as seen even when the value is junk, so a later pose's block cannot be
            # read as this pose's missing number.
            seen.add(name)
            values[name] = _finite_float(match.group(1))
            break
        if len(seen) == len(_SCORE_PATTERNS):
            break
    return PoseScores(**values)


def rescore_poses(receptor, pose_paths, *, env=None, runner=None, devices=None,
                  on_error="raise"):
    """Rescore every pose in `pose_paths` against `receptor` with GNINA, one call per pose.

    Returns a RescoreResults mapping each pose path to its PoseScores. GNINA is run once per
    pose, round-robin over the usable device indices.

    GPU only: the GpuUnavailableError raised by gpu_backends when there is no usable card or no
    `gnina` binary is deliberately not caught, so the stage stops instead of degrading.

    `env` is the mapping gpu_backends reads the FRANKPEP_GNINA_PATH override and
    CUDA_VISIBLE_DEVICES from; it is not the child process environment. `runner` replaces
    subprocess.run for the GNINA calls and is forwarded to gpu_backends for its nvidia-smi
    call, so a test has one injection point for both. `devices` overrides device detection.

    `on_error="raise"` (the default) turns a non-zero exit, or a GNINA that will not launch,
    into a RuntimeError carrying the command and the captured stderr. `on_error="skip"` leaves
    such a pose out of the mapping and appends a SkippedPose to the result's `.skipped`, which
    is the only place that loss is recorded: a caller that ignores it has silently lost a pose.

    A GNINA that exits 0 without printing a CNNaffinity yields a PoseScores whose cnn_affinity
    is None. That is reported, not substituted: callers ranking on the CNN must handle None.
    """
    if on_error not in ON_ERROR_MODES:
        raise ValueError(
            f"on_error must be one of {ON_ERROR_MODES}, got {on_error!r}")
    poses = list(pose_paths)
    if not poses:
        raise ValueError("No poses to rescore: pose_paths is empty.")

    if devices is None:
        device_indices = list(gpu_backends.visible_device_indices(env=env, runner=runner))
    else:
        # An explicit device list says which cards to use; it is not a way around the GPU-only
        # rule, so the hardware is still required. visible_device_indices() already calls
        # require_gpus() itself, which is why only this branch calls it.
        gpu_backends.require_gpus(runner=runner)
        device_indices = list(devices)
    if not device_indices:
        raise ValueError("devices is empty: there is no GPU index to run GNINA on.")

    binary = gpu_backends.require_binary("gnina", env=env)
    call = subprocess.run if runner is None else runner

    results = RescoreResults()
    for position, pose in enumerate(poses):
        # Round-robin rather than a pool: the calls are short and independent, so spreading
        # them over the cards keeps a second GPU from idling through the whole stage.
        device_index = device_indices[position % len(device_indices)]
        command = build_command(receptor, pose, binary=binary, device_index=device_index)
        try:
            completed = call(command, check=False, capture_output=True, text=True)
        except OSError as error:
            # The binary resolved moments ago, so this is it disappearing or not being
            # executable. Treated exactly like a non-zero exit: a pose with no score.
            failure = SkippedPose(pose, command, None, str(error))
        else:
            failure = None
            if completed.returncode != 0:
                failure = SkippedPose(pose, command, completed.returncode, completed.stderr)

        if failure is not None:
            if on_error == "raise":
                raise RuntimeError(failure.describe())
            results.skipped.append(failure)
            continue
        results[pose] = parse_scores(completed.stdout)
    return results
