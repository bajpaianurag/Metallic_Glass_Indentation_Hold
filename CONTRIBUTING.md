# Contributing

## Scope

This repository is primarily a research code archive (LAMMPS inputs + analysis utilities).
Contributions that improve reproducibility, documentation, portability (cluster/local), and
post-processing robustness are welcome.

## Suggested workflow

1. Create a new branch for your change.
2. Keep changes minimal and well-scoped.
3. If you modify analysis scripts, include a short note in the README (expected inputs/outputs).
4. Avoid committing large trajectories/outputs; use Git LFS or external archival (Zenodo) if needed.

## Style

- Python: keep scripts runnable as standalone CLIs; prefer explicit arguments over hard-coded paths.
- Shell: be conservative and portable (bash).
- LAMMPS: document any required packages/fixes used by an input.
