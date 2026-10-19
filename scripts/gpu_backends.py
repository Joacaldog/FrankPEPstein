"""GPU discovery and preflight checks for the docking, rescoring and co-folding stages.

This path is GPU-only by design: there is no CPU fallback here. A missing GPU, a card too old
for the stage, or a missing Uni-Dock / GNINA binary raises GpuUnavailableError. That is the
point of the module. A run that silently degrades to CPU still finishes and still writes a TSV,
so its timings and its rankings end up compared against GPU runs as if they were the same
thing; failing at the start, with a message saying what to fix, is the cheaper outcome.

Two capability floors matter, and they are not the same number:

    MIN_COMPUTE_CAPABILITY     6.1   Uni-Dock and GNINA; Pascal (GTX 1080 / 1080 Ti) clears it
    MODERN_COMPUTE_CAPABILITY  8.0   bfloat16 and tensor cores, i.e. Boltz-2 and Chai-1

That split is why both exist: the lab machine's six GTX 1080 Ti cards can dock and rescore but
cannot run the co-folding models, while the desktop's RTX 5070 Ti (compute 12.0) can do both and
the GTX 1080 next to it cannot.

Binary locations can be overridden per machine, which is the normal case on a shared box where
neither tool is on PATH:

    FRANKPEP_UNIDOCK_PATH   Uni-Dock binary
    FRANKPEP_GNINA_PATH     GNINA binary

`describe_environment()` reports all of this for the run log and never raises; the require_*
functions are the ones that stop a run.
"""

import os
import shutil
import subprocess

# Floor for Uni-Dock and GNINA. Pascal sits exactly on it, which is the whole reason the GPU
# docking work could start on the hardware the lab already had.
MIN_COMPUTE_CAPABILITY = 6.1

# Floor for bfloat16 and tensor cores, i.e. the co-folding models. Ampere and newer clear it.
MODERN_COMPUTE_CAPABILITY = 8.0

_NVIDIA_SMI_QUERY = "index,name,compute_cap,memory.total"
_NVIDIA_SMI_FORMAT = "csv,noheader,nounits"


class GpuUnavailableError(RuntimeError):
    """Raised when a GPU or a GPU binary this run needs is not available.

    Nothing in the GPU path catches this: a stage that cannot reach the hardware it was asked
    for stops, it does not quietly run something else. RuntimeError is the base so a caller that
    already guards a stage with `except RuntimeError` keeps working.
    """


class Gpu:
    """One GPU as nvidia-smi reports it.

    `compute_cap` and `memory_mib` are None when nvidia-smi could not give a number for them
    (it prints "[N/A]" for both on some virtualised and very old-driver setups).
    """

    def __init__(self, index, name, compute_cap, memory_mib):
        self.index = index
        self.name = name
        self.compute_cap = compute_cap
        self.memory_mib = memory_mib

    @property
    def supports_modern_models(self):
        """True when this card has bfloat16 and tensor cores, i.e. can run Boltz-2 or Chai-1."""
        # An unknown capability counts as unsupported. Refusing a card we cannot identify costs
        # a clear error message; guessing costs a stage that dies inside CUDA an hour later.
        return self.compute_cap is not None and self.compute_cap >= MODERN_COMPUTE_CAPABILITY

    @property
    def supports_docking(self):
        """True when this card can run Uni-Dock and GNINA."""
        return self.compute_cap is not None and self.compute_cap >= MIN_COMPUTE_CAPABILITY

    def __repr__(self):
        cap = "unknown" if self.compute_cap is None else f"{self.compute_cap:.1f}"
        return (f"Gpu(index={self.index}, name={self.name!r}, compute_cap={cap}, "
                f"memory_mib={self.memory_mib})")


def _override_variable(tool):
    """Name of the environment variable that overrides where `tool` lives."""
    return f"FRANKPEP_{tool.strip().upper()}_PATH"


def _describe_gpu(gpu):
    """One line of fact about a card, shared by the run log and the error messages."""
    cap = "unknown" if gpu.compute_cap is None else f"{gpu.compute_cap:.1f}"
    memory = "unknown" if gpu.memory_mib is None else f"{gpu.memory_mib} MiB"
    return f"GPU {gpu.index}: {gpu.name}, compute {cap}, {memory}"


def _run_nvidia_smi(runner=None):
    """Return nvidia-smi's CSV output, or None when the tool is absent or fails."""
    if runner is None:
        # The PATH lookup belongs to the real path only: an injected runner *is* the stand-in
        # for nvidia-smi, so the tests must work on a machine where the executable is absent.
        if shutil.which("nvidia-smi") is None:
            return None
        runner = subprocess.run
    try:
        result = runner(
            ["nvidia-smi", f"--query-gpu={_NVIDIA_SMI_QUERY}", f"--format={_NVIDIA_SMI_FORMAT}"],
            check=False, capture_output=True, text=True)
    except OSError:
        return None
    if result.returncode != 0:
        return None
    return result.stdout


def parse_nvidia_smi(output):
    """Parse the CSV rows of `nvidia-smi --query-gpu` into Gpu objects.

    Separate from the call so the whole decision layer is testable without a GPU. Rows that do
    not parse are skipped rather than raising: one odd row should not hide the cards that did
    report properly, and a card that comes back with an unknown capability is refused later by
    its own capability properties anyway.
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
        # The numbers are read from the end and the name is whatever sits between: a model name
        # can itself contain a comma, and losing a real card to that would fail an otherwise
        # usable machine.
        try:
            compute_cap = float(fields[-2])
        except ValueError:
            compute_cap = None
        try:
            memory_mib = int(float(fields[-1]))
        except ValueError:
            memory_mib = None
        name = ", ".join(fields[1:-2])
        gpus.append(Gpu(index=index, name=name, compute_cap=compute_cap,
                        memory_mib=memory_mib))
    return gpus


def detect_gpus(runner=None):
    """List the GPUs nvidia-smi reports. Empty list when there are none.

    `runner` replaces subprocess.run, for tests. A missing or failing nvidia-smi is reported as
    "no GPU" here, not as an error: the require_* functions below own the decision to stop.
    """
    return parse_nvidia_smi(_run_nvidia_smi(runner=runner))


def require_gpus(runner=None):
    """The GPUs this run may dock on, or raise GpuUnavailableError.

    Only the cards that clear MIN_COMPUTE_CAPABILITY come back, so a caller never has to
    re-check: anything in the returned list can run Uni-Dock and GNINA.
    """
    gpus = detect_gpus(runner=runner)
    if not gpus:
        raise GpuUnavailableError(
            "No GPU is visible: nvidia-smi reported no device. This stage is GPU-only and has "
            "no CPU fallback. Check that `nvidia-smi --query-gpu="
            f"{_NVIDIA_SMI_QUERY} --format={_NVIDIA_SMI_FORMAT}` runs and lists your cards "
            "(an older driver that does not know the compute_cap field fails this query), that "
            "CUDA_VISIBLE_DEVICES is not hiding them, and otherwise run this stage on a machine "
            "with a GPU.")

    usable = [gpu for gpu in gpus if gpu.supports_docking]
    if not usable:
        found = "; ".join(_describe_gpu(gpu) for gpu in gpus)
        raise GpuUnavailableError(
            f"No visible GPU meets compute capability {MIN_COMPUTE_CAPABILITY:g}, the floor for "
            f"Uni-Dock and GNINA. Found: {found}. This stage is GPU-only and has no CPU "
            f"fallback: run it on a card of compute capability {MIN_COMPUTE_CAPABILITY:g} or "
            "newer (GTX 1080 Ti and up), or select such a card with CUDA_VISIBLE_DEVICES.")
    return usable


def require_binary(tool, env=None):
    """Path to `tool`'s executable ("unidock", "gnina"), or raise GpuUnavailableError.

    The override variable is consulted first and, when it is set, it is the only thing
    consulted: a stale FRANKPEP_*_PATH then fails loudly instead of resolving to some other
    build that happens to be on PATH, which would be the same silent substitution this module
    exists to prevent.
    """
    env = os.environ if env is None else env
    tool = tool.strip()
    variable = _override_variable(tool)
    override = (env.get(variable) or "").strip()

    if override:
        if os.path.exists(override):
            return override
        raise GpuUnavailableError(
            f"{variable} is set to {override!r}, which does not exist. Point {variable} at the "
            f"{tool} binary, or unset it to look {tool} up on PATH. This stage is GPU-only: "
            f"there is no CPU fallback for {tool}.")

    found = shutil.which(tool)
    if found:
        return found
    raise GpuUnavailableError(
        f"The {tool} binary was not found on PATH. Install {tool} or set {variable} to its full "
        f"path. This stage is GPU-only: there is no CPU fallback for {tool}.")


def visible_device_indices(env=None, runner=None):
    """Indices of the GPUs this process may use, or raise GpuUnavailableError.

    CUDA_VISIBLE_DEVICES lists *physical* indices, so it is intersected with what nvidia-smi
    reported and the result stays in that physical numbering: a child process launched with
    CUDA_VISIBLE_DEVICES=<one of these> sees that card as its own device 0. Entries that are not
    integers (UUID forms such as GPU-0f1e..., stray commas) are ignored, because there is
    nothing in nvidia-smi's index column to match them against.
    """
    env = os.environ if env is None else env
    usable = require_gpus(runner=runner)
    raw = env.get("CUDA_VISIBLE_DEVICES")

    if raw is None or not raw.strip():
        # CUDA itself reads an empty value as "no devices", but an empty value is also what a
        # shell leaves behind when an unset variable is expanded into it, so it is treated here
        # as "no restriction" rather than failing a run over a blank.
        return [gpu.index for gpu in usable]

    allowed = []
    for entry in raw.split(","):
        try:
            allowed.append(int(entry.strip()))
        except ValueError:
            continue

    indices = [gpu.index for gpu in usable if gpu.index in allowed]
    if not indices:
        found = ", ".join(str(gpu.index) for gpu in usable)
        raise GpuUnavailableError(
            f"CUDA_VISIBLE_DEVICES={raw!r} leaves no usable GPU. The cards that meet compute "
            f"capability {MIN_COMPUTE_CAPABILITY:g} are at indices {found}. Set "
            "CUDA_VISIBLE_DEVICES to a subset of those, or unset it to use all of them. This "
            "stage is GPU-only and has no CPU fallback.")
    return indices


def describe_environment(env=None, runner=None):
    """Lines describing what this machine can run. Printed at the start of a stage.

    A report, never a check: it does not raise, so the run log says what was seen before
    require_gpus() or require_binary() stops the run with the specific reason.
    """
    gpus = detect_gpus(runner=runner)
    if not gpus:
        return ["GPU: none visible (nvidia-smi reported no device). This path is GPU-only; "
                "there is no CPU fallback."]

    lines = []
    for gpu in gpus:
        if gpu.supports_modern_models:
            note = "can run co-folding models (bfloat16, tensor cores)"
        elif gpu.supports_docking:
            note = "docking and rescoring only, too old for co-folding models"
        else:
            note = f"unusable here, below compute capability {MIN_COMPUTE_CAPABILITY:g}"
        lines.append(f"{_describe_gpu(gpu)} ({note})")

    try:
        indices = visible_device_indices(env=env, runner=runner)
    except GpuUnavailableError as error:
        lines.append(f"Usable device indices: none. {error}")
    else:
        lines.append("Usable device indices: " + ", ".join(str(index) for index in indices))
    return lines
