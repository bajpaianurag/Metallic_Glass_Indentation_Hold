# MD workflow for Cu–Zr-(Al-Ag-Ti) metallic glass: glass preparation, depth-hold nano-indentation, and OVITO-based post-processing

This repository contains LAMMPS input scripts and post-processing utilities used to (i) generate an amorphous Cu–Zr-(Al/Ag/Ti) metallic glass, (ii) perform constant-depth hold nano-indentation, and (iii) compute structural/kinetic descriptors such as Voronoi-based free volume and non-affine displacement (D2min) using OVITO’s Python API.

## Repository layout

- `MD codes/`
  - `1.glass/` — glass preparation (cooling/heating and anneals)
    - `lammps_script` — main LAMMPS input
    - `slurm_run.sh` — example SLURM submission script
    - `post-processing/` — Python utilities (OVITO) for RDF, Voronoi/free volume, bond analysis
  - `2.indentation/` — nano-indentation workflow
    - `lammps_script` — indentation input
    - `slurm_run.sh` — example SLURM submission script
    - `postprocess/` — Python utilities (OVITO + analysis) for D2min^2, STZ tagging, and Oliver–Pharr style post-processing
  - `ZrCu.meam`, `library.meam` — MEAM potential files used by the LAMMPS inputs
  - `move.sh`, `wait.sh` — helper shell scripts

- `matlab/` — MATLAB plotting/analysis scripts (D2min, free volume, and force metrics)

## Requirements

### Simulation

- **LAMMPS** built with the packages required by the inputs (commonly: `MEAM`, `MISC`, and standard MD fixes).
- A working MPI installation is recommended for production runs.

### Post-processing (Python)

The post-processing scripts in `post-processing/` and `postprocess/` rely on:

- `ovito` (OVITO Python module; tested with OVITO 3.x)
- `numpy`
- `pandas` (indentation curve analysis)
- `matplotlib` (plots)

Install via pip:

```bash
python -m pip install -r requirements.txt
```

If you work on clusters where `ovito` is more convenient via conda:

```bash
conda env create -f environment.yml
conda activate md-ovito
```

## Quick start

### 1) Glass preparation (LAMMPS)

```bash
cd "MD codes/1.glass"
# Run locally
lmp -in lammps_script

# Or submit on SLURM (edit modules, partitions, MPI launcher as needed)
sbatch slurm_run.sh
```

Outputs typically include `restart.*`, `log.lammps`, and trajectory/thermo files (names depend on the input).

### 2) Nano-indentation (LAMMPS)

```bash
cd "MD codes/2.indentation"
lmp -in lammps_script
# or
sbatch slurm_run.sh
```

The indentation input is written for LAMMPS `fix indent`-style workflows and generates averaged time series such as:

- `indentation_out.load`
- `indentation_out.unload`
- *(optional)* `indentation_out.hold`

### 3) Post-processing

#### 3A) D2min (OVITO)

From the indentation post-processing directory:

```bash
cd "MD codes/2.indentation/postprocess"
python d2min.py \
  --pre  ../dump.pre \
  --post ../dump.post \
  --cutoff 4.0 \
  --stz-thresh 1.5 \
  --hist-bins 100
```

Key outputs:

- `d2min.csv` — per-atom D2min with `STZ` flag
- `d2min_hist.csv` — histogram (optionally log-binned)
- `stz_ids.txt` — atom IDs tagged as STZ (D2min ≥ threshold)
- `stz_only.xyz` — STZ-only snapshot (can be visualized in OVITO)

#### 3B) Load–depth, hardness, and modulus (Oliver–Pharr style)

```bash
cd "MD codes/2.indentation/postprocess"
python load_depth_headness_modulus.py
```

Outputs:

- `indent_summary.csv`
- `indent_curve.csv`
- `indent_load_depth.png`
- `unloading_fit.png`

#### 3C) Voronoi/free-volume analysis (OVITO)

```bash
cd "MD codes/1.glass/post-processing"
python free.py final.cfg --out freeV --use-radii --radii "Cu=1.28,Zr=1.58" --export-per-atom
```

This computes Voronoi cell volume, occupied atomic volume, and free volume per atom (with optional clipping).

## Notes on data files in Git

LAMMPS outputs (e.g., large `dump.*`, `restart.*`, `log.lammps`, and `.cfg` snapshots) can be multiple GB and typically should **not** be committed to GitHub.

This repository includes a conservative `.gitignore` that excludes common heavy outputs by default. If you want to version large trajectories, use **Git LFS**.

## Citation

If you use this code in academic work, please cite the associated paper/preprint. A `CITATION.cff` template is included; update it with the final bibliographic details.

## License

This code is released under the MIT License (see `LICENSE`). If you require a different license for your project/paper, replace the license file accordingly.
