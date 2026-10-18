# GPU backends

The GPU paths are opt-in. With no environment variable set, the pipeline runs exactly as before
on CPU, and every backend here falls back to that CPU path rather than failing a run.

Check what a machine would do before launching anything long:

```bash
python3 scripts/gpu_check.py
FRANKPEP_RESCORE_BACKEND=consensus python3 scripts/gpu_check.py
```

The second form is the useful one: it says whether the opt-in takes effect or silently falls
back, and why.

## Variables

| Variable | Values | Effect |
| --- | --- | --- |
| `FRANKPEP_RESCORE_BACKEND` | `vina` (default), `gnina`, `consensus` | Adds a `CNNaffinity_fullReceptor` column to the final TSV, from GNINA's CNN scoring the poses FrankVINA 2 already produced |
| `FRANKPEP_DOCK_BACKEND` | `vina` (default), `unidock` | Selects Uni-Dock for batch docking on GPU |
| `FRANKPEP_GNINA_PATH` | path | GNINA binary when it is not on `PATH` |
| `FRANKPEP_UNIDOCK_PATH` | path | Uni-Dock binary when it is not on `PATH` |

A backend is used only when all three hold: it was requested, `nvidia-smi` reports a GPU, and
its binary is found. Otherwise the run says why and continues on CPU.

## What each backend needs

| Backend | Compute capability | Runs on |
| --- | --- | --- |
| Uni-Dock, GNINA | 6.1 and up | GTX 1080, GTX 1080 Ti, and anything newer |
| Boltz-2, Chai-1 (not wired up here) | 8.0 and up, for bfloat16 and tensor cores | RTX 3080, RTX 5070 Ti, and newer |

That split is why the rescoring path came first: it needs no hardware the lab machine does not
already have.

## Reading the extra column

The Vina columns are binding energies in kcal/mol, so lower is better. `CNNaffinity_fullReceptor`
is a pKd-like number from a different scoring function, so **higher is better**. The column is
recorded alongside the Vina ranking and never used to drop a candidate: a scoring function the
pipeline has not been calibrated against informs the ranking, it does not decide it. Compare the
two orderings over a run you already trust before giving the CNN score any weight.

## Installing the tools

Neither tool is in `environment.yml`, because neither is needed for a CPU run and both pull a
CUDA stack. Install whichever you intend to use, into the activated `FrankPEPstein` environment:

- GNINA: prebuilt binary from <https://github.com/gnina/gnina/releases>, or build from source.
- Uni-Dock: <https://github.com/dptech-corp/Uni-Dock>.

Then point the pipeline at it, if it is not on `PATH`:

```bash
export FRANKPEP_GNINA_PATH=/opt/gnina/gnina
export FRANKPEP_RESCORE_BACKEND=consensus
```

## Status

`FRANKPEP_RESCORE_BACKEND` is wired into FrankVINA 2's final rescore. `FRANKPEP_DOCK_BACKEND`
resolves and builds the Uni-Dock command, but no stage calls it yet: the docking stages live in
files another change is rewriting, so wiring it is a separate change on top of this one.

## Tests

```bash
python3 -m unittest discover -s tests
```

They inject `nvidia-smi` output and binary paths, so they cover the fallback decisions on a
machine with no GPU and neither tool installed.
