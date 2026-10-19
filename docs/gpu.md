# GPU backends

The GPU paths are opt-in, and **GPU-only: there is no CPU fallback.** When a GPU path is asked
for and the card or the binary is missing, the stage stops with a message naming what to fix. A
run that quietly degrades to CPU still finishes and still writes a TSV, so its timings and its
rankings end up compared against GPU runs as if they were the same thing; failing at the start
is the cheaper outcome.

Not asking for a GPU path leaves the pipeline exactly as it was: `FRANKPEP_RESCORE_BACKEND`
unset means Vina on CPU, same columns, same numbers.

Check a machine before launching anything long:

```bash
python3 scripts/gpu_check.py
```

It reports every visible card with its compute capability, which of them this pipeline may use,
and never raises — it is a report, not a gate.

## Variables

| Variable | Values | Effect |
| --- | --- | --- |
| `FRANKPEP_RESCORE_BACKEND` | `vina` (default), `gnina` | `gnina` adds `CNNaffinity_fullReceptor` and `CNNscore_fullReceptor` to the final TSV, from GNINA scoring the poses FrankVINA 2 already produced. Requires a GPU. |
| `FRANKPEP_GNINA_PATH` | path | GNINA binary, when it is not on `PATH`. When set, it is the only thing consulted, so a stale value fails loudly instead of silently resolving a different build. |
| `FRANKPEP_UNIDOCK_PATH` | path | Uni-Dock binary, same rule. |
| `CUDA_VISIBLE_DEVICES` | e.g. `0,1` | Restricts which of the detected cards the pipeline may use. |

A misspelled backend name raises instead of defaulting to `vina`: a run launched expecting GNINA
columns that quietly produced Vina-only ones is the failure this path exists to avoid.

## Capability floors

| Stage | Compute capability | Cards that clear it |
| --- | --- | --- |
| Uni-Dock, GNINA | 6.1 | GTX 1080, GTX 1080 Ti, and newer |
| Co-folding models (Boltz-2, Chai-1), not wired up here | 8.0, for bfloat16 and tensor cores | RTX 3080, RTX 5070 Ti, and newer |

That split is why the rescoring work came first: it needs no hardware the lab machine does not
already have. A card whose capability `nvidia-smi` cannot report is refused rather than guessed.

## What is wired up

| Module | What it does | Called by |
| --- | --- | --- |
| `scripts/gpu_backends.py` | Discovery and preflight: which cards exist, which clear the floor, where the binaries are. Raises `GpuUnavailableError` rather than degrading. | both modules below, and FrankVINA 2 |
| `scripts/gnina_rescore.py` | Runs GNINA `--score_only` over existing poses, round-robin across the usable cards, and parses `CNNaffinity`, `CNNscore` and `Affinity`. | FrankVINA 2's final rescore |
| `scripts/unidock_runner.py` | Batches ligands into one Uni-Dock call per card, parses each output pose's affinity, reports every ligand that came back without a usable score. | nothing yet — see below |

FrankVINA 2 preflights the GPU at the start of the stage, not at the end: a missing card costs a
second rather than a whole run's MODELLER refinement.

`unidock_runner` is complete and tested but not yet called by any stage. The screening dock lives
in `scripts/frankVINA_1.py`, which another change is currently rewriting; wiring it is a separate
change on top of this one.

## Reading the extra columns

The Vina columns are binding energies in kcal/mol, so **lower is better**. `CNNaffinity` is a
pKd-like number, so **higher is better**. `CNNscore` is a pose-quality probability in [0, 1], not
an affinity at all.

An empty cell means GNINA reported no number for that pose. It never means zero: a 0.0
`CNNaffinity` is a real, very weak prediction, so filling a gap with it would read as a confident
answer.

The CNN columns are recorded alongside the Vina ranking and never used to drop a candidate. A
scoring function this pipeline has not been calibrated against informs the ranking; it does not
decide it. Compare the two orderings over a run you already trust before giving the CNN weight.

## Installing the tools

Neither tool is in `environment.yml`: neither is needed for a CPU run and both pull a CUDA stack.
Install whichever you intend to use, then point the pipeline at it.

- GNINA: prebuilt binary from <https://github.com/gnina/gnina/releases>, or build from source.
- Uni-Dock: <https://github.com/dptech-corp/Uni-Dock>.

```bash
export FRANKPEP_GNINA_PATH=/opt/gnina/gnina
export FRANKPEP_RESCORE_BACKEND=gnina
```

## A note on peptide-specific scoring

GNINA's CNN was trained on protein–small-molecule complexes, so for 5–15 residue peptides it is a
second opinion, not an authority. The scoring functions actually built for peptide–protein
interfaces are CPU-only:

- **Rosetta / PyRosetta** `FlexPepDocking -flexpep_score_only` with `ref2015`, ranked on
  `reweighted_sc`, which weights interface and peptide residues. Pip- or conda-installable, free
  for academic and non-profit use, paid for commercial use. Takes a PDB with the receptor chain
  first and the peptide second.
- **HADDOCK3** `[topoaa]` + `[emscoring]`, which scores an existing complex without re-docking.
  Apache-2.0, but its CNS dependency needs a free academic registration with a university email.
  Its score beat MM-PBSA at ranking on a protein–peptide benchmark.

Neither uses a GPU, so they parallelise across CPU cores instead. ADCP has no mode that scores an
externally produced pose, and InterPepRank's code could not be obtained.

## Tests

```bash
python3 -m unittest discover -s tests
```

They inject `nvidia-smi` output, binary paths and the subprocess runner, so the whole decision
layer is covered on a machine with no GPU and neither tool installed. What they cannot cover is a
real Uni-Dock or GNINA call; that has to be run on the hardware.
