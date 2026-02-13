# OVITO 3.13 only: Save total + partial RDF
from ovito.io import import_file
from ovito.modifiers import CoordinationAnalysisModifier
import numpy as np
import pandas as pd
import os

# ===== User configuration =====
# Path to a structure file that OVITO can read (e.g., .cfg, .xyz, .lammpsdump)
input_file = "final.cfg"
# Radial cutoff in Å; histogram will span [0, r_max]
r_max = 10.0
# Bin width in Å; number of bins is derived from r_max/dr
dr    = 0.1
n_bins = max(1, int(round(r_max / dr)))

# Output directory for all CSVs
outdir = "RDF"
os.makedirs(outdir, exist_ok=True)


def sanitize(name: str) -> str:
    """Return a file/name-safe version of a label.
    Removes spaces and replaces slashes with dashes to avoid path issues.
    """
    return name.replace(" ", "").replace("/", "-").replace("\\", "-")


# ===== 1) Total RDF =====
# Build a data pipeline and compute a total (species-agnostic) g(r)
pipe_total = import_file(input_file)
pipe_total.modifiers.append(CoordinationAnalysisModifier(
    cutoff=r_max,
    number_of_bins=n_bins,
    partial=False  # total RDF only
))

data_total = pipe_total.compute()
if data_total.particles.count == 0:
    raise RuntimeError("No particles found in the input frame.")

# The RDF is exposed via the 'coordination-rdf' data table in OVITO
tbl_total = data_total.tables.get('coordination-rdf', None)
if tbl_total is None:
    raise RuntimeError(
        f"Could not find the total RDF table. Available tables: {list(data_total.tables.keys())}"
    )

# xy() returns a 2D array where the first column is r and the remaining columns are g(r)
xy_tot = tbl_total.xy()
if xy_tot.shape[1] < 2:
    raise RuntimeError("Unexpected total RDF format (number of columns < 2).")

r = xy_tot[:, 0]
g_total = xy_tot[:, 1]

# Save total RDF as a simple two-column CSV
pd.DataFrame({"r (Angstrom)": r, "g(r)": g_total}).to_csv(
    os.path.join(outdir, "total_rdf.csv"), index=False
)

# ===== 2) Partial RDF (per species pair) =====
# Recompute with partial=True to obtain all pairwise g_{αβ}(r)
pipe_part = import_file(input_file)
pipe_part.modifiers.append(CoordinationAnalysisModifier(
    cutoff=r_max,
    number_of_bins=n_bins,
    partial=True   # produce species pair-resolved RDFs
))

data_part = pipe_part.compute()
tbl_part = data_part.tables.get('coordination-rdf', None)
if tbl_part is None:
    raise RuntimeError(
        f"Could not find the partial RDF table. Available tables: {list(data_part.tables.keys())}"
    )

# First column is r; remaining columns are one series per species pair
xy_par = tbl_part.xy()
if xy_par.shape[1] < 2:
    raise RuntimeError("Partial RDF appears to be empty (number of columns < 2).")

r_p = xy_par[:, 0]
Y = xy_par[:, 1:]  # shape: (n_bins, n_series)

# Retrieve human-readable series labels for each pair (e.g., 'A-A', 'A-B', ...)
names = list(getattr(tbl_part.y, "component_names", []))
if not names or len(names) != Y.shape[1]:
    # Fallback to generic labels if names are missing or count mismatches
    names = [f"series_{i}" for i in range(Y.shape[1])]

# 2-1) Save all partial series into a single wide CSV (columns: r, pair1, pair2, ...)
df_all = pd.DataFrame({"r (Angstrom)": r_p})
for i, nm in enumerate(names):
    df_all[sanitize(nm)] = Y[:, i]
df_all.to_csv(os.path.join(outdir, "partial_rdf_all.csv"), index=False)

# 2-2) Also save each pair as a separate two-column CSV for convenient plotting
for i, nm in enumerate(names):
    safe = sanitize(nm)
    pd.DataFrame({"r (Angstrom)": r_p, "g(r)": Y[:, i]}).to_csv(
        os.path.join(outdir, f"{safe}_rdf.csv"), index=False
    )

print(
    f"Done: wrote {outdir}/total_rdf.csv, {outdir}/partial_rdf_all.csv, and {len(names)} per-pair CSV files"
)

