#!/usr/bin/env python3
# -*- coding: utf-8 -*-

import argparse, json, re, sys, csv
import numpy as np
from math import pi
from ovito.io import import_file, export_file
from ovito.modifiers import VoronoiAnalysisModifier, ComputePropertyModifier

def export_per_atom_csv(data, atomic_vol_key, out_csv_path):
    """
    Write a per-atom CSV with identifiers, types, positions, and volume-related properties.
    CSVs don't carry dtype metadata well, so we cast some values to float for portability.
    """
    P = data.particles

    # Safely fetch Particle ID / Type / Position
    if "Particle Identifier" in P:
        # Cast to float to avoid dtype issues in CSV consumers
        pid = np.asarray(P["Particle Identifier"], dtype=float)
    else:
        # If no explicit IDs, create 1..N consecutive IDs
        pid = np.arange(1, P.count + 1, dtype=float)

    ptype = np.asarray(P["Particle Type"], dtype=float) if "Particle Type" in P else np.full(P.count, -1.0)

    # Position is typically stored as an (N x 3) array named "Position" in OVITO
    if "Position" in P:
        pos = np.asarray(P["Position"])
        x, y, z = pos[:,0], pos[:,1], pos[:,2]
    else:
        # Rare case: positions are only present as component properties
        x = np.asarray(P["Position.X"])
        y = np.asarray(P["Position.Y"])
        z = np.asarray(P["Position.Z"])

    # Voronoi (atomic) volume, occupied volume, free volume, effective radius
    atomic_vol = np.asarray(P[atomic_vol_key], dtype=float)
    v_atom     = np.asarray(P["V_atom"], dtype=float)
    free_v     = np.asarray(P["free_volume"], dtype=float)
    rad_eff    = np.asarray(P["Radius_effective"], dtype=float)

    header = ["ParticleID","Type","X","Y","Z", atomic_vol_key, "V_atom","free_volume","Radius_effective"]

    with open(out_csv_path, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(header)
        for i in range(P.count):
            w.writerow([pid[i], ptype[i], x[i], y[i], z[i],
                        atomic_vol[i], v_atom[i], free_v[i], rad_eff[i]])

def parse_radii_map(radii_str):
    """
    Parse a mapping string like "Cu=1.28,Zr=1.58,Al=1.43" into {"Cu":1.28, "Zr":1.58, "Al":1.43}.
    Also accepts numeric ParticleType IDs, e.g., "1=1.35,2=1.25".
    """
    if not radii_str:
        return {}
    m = {}
    for token in re.split(r'[,\s]+', radii_str.strip()):
        if not token:
            continue
        if '=' not in token:
            raise ValueError(f"Malformed radius entry: {token} (e.g., use Cu=1.28)")
        k, v = token.split('=', 1)
        k = k.strip()
        v = float(v.strip())
        m[k] = v
    return m

import numpy as np

def build_radius_expression(data, radii_map, fallback=1.0):
    """
    Build an OVITO ComputeProperty expression for an effective radius per particle,
    based on the current frame's ParticleType (by name or numeric ID).

    - radii_map may contain entries like {"Cu":1.28, "Zr":1.58} or {"1":1.35,"2":1.25}.
    - Name matching uses the names of the particle types present in 'data'.
    - Any present types not listed in 'radii_map' receive the 'fallback' radius.
    """
    # 1) Determine which type IDs actually appear in the current data
    type_ids = []
    try:
        if "Particle Type" in data.particles:
            arr = np.asarray(data.particles["Particle Type"])
            type_ids = sorted({int(x) for x in arr.tolist()})
    except Exception:
        type_ids = []

    # If we cannot discover any types, just return the fallback scalar
    if not type_ids:
        return str(float(fallback))

    # 2) Build a robust name->id map by querying OVITO's type container
    name_to_id = {}
    try:
        types = data.particles.particle_types
        for tid in type_ids:
            try:
                t = types.type_by_id(int(tid))
                tname = getattr(t, "name", None) if t is not None else None
                if tname:
                    name_to_id[str(tname)] = int(tid)
            except Exception:
                # Skip any problematic id
                continue
    except Exception:
        # If the particle_types container is unavailable, numeric mapping still works
        name_to_id = {}

    # 3) Convert radii_map into a piecewise expression over ParticleType
    used_ids = set()
    terms = []

    # (a) Match by name first (case-insensitive fallback)
    folded = {k.casefold(): v for k, v in name_to_id.items()}
    for key, radius in radii_map.items():
        if key.isdigit():
            continue
        key_strip = key.strip()
        # Exact match first
        if key_strip in name_to_id:
            tid = name_to_id[key_strip]
            terms.append(f"(ParticleType=={tid})*({float(radius)})")
            used_ids.add(tid)
            continue
        # Case-insensitive match
        kf = key_strip.casefold()
        if kf in folded:
            tid = folded[kf]
            terms.append(f"(ParticleType=={tid})*({float(radius)})")
            used_ids.add(tid)

    # (b) Match by numeric id
    for key, radius in radii_map.items():
        if key.isdigit():
            tid = int(key)
            terms.append(f"(ParticleType=={tid})*({float(radius)})")
            used_ids.add(tid)

    # (c) Apply fallback to any present but unspecified types
    unmatched_ids = [tid for tid in type_ids if tid not in used_ids]
    if unmatched_ids:
        unmatched_sum = " + ".join([f"(ParticleType=={tid})" for tid in unmatched_ids])
        terms.append(f"({float(fallback)})*({unmatched_sum})")

    # If no terms materialized (should be rare), return fallback scalar
    if not terms:
        return str(float(fallback))

    return " + ".join(terms)

def main():
    ap = argparse.ArgumentParser(description="Compute per-atom Voronoi volume and free volume statistics from AtomEye CFG using OVITO.")
    ap.add_argument("cfg", help="Input CFG file (AtomEye format)")
    ap.add_argument("--out", default="freeV", help="Output file prefix (default: freeV)")
    ap.add_argument("--use-radii", action="store_true", help="Use radical (radius-weighted) Voronoi")
    ap.add_argument("--radii", default="", help='Atomic radii map, e.g. "Cu=1.28,Zr=1.58" or "1=1.35,2=1.25"')
    ap.add_argument("--fallback-radius", type=float, default=1.0, help="Fallback radius if a type is unspecified (default: 1.0)")
    ap.add_argument("--bins", type=int, default=150, help="Histogram bins (default: 150)")
    ap.add_argument("--hist-min", type=float, default=None, help="Histogram minimum (auto if omitted)")
    ap.add_argument("--hist-max", type=float, default=None, help="Histogram maximum (auto if omitted)")
    ap.add_argument("--clip-negative", action="store_true", help="Clamp negative free volume to 0")
    ap.add_argument("--export-per-atom", action="store_true", help="Export a per-atom CSV")
    args = ap.parse_args()

    # 1) Load input and determine dimensionality from the simulation cell
    pipeline = import_file(args.cfg)
    data = pipeline.compute()
    cell = data.cell
    is_3d = True
    if cell is None:
        is_3d = False
    elif hasattr(cell, "is2D"):
        is_3d = not cell.is2D
    else:
        try:
            is_3d = (getattr(cell, "matrix", None) is not None and len(cell.matrix) == 3)
        except Exception:
            is_3d = True

    if not is_3d:
        print("[WARNING] The periodic cell is missing or 2D. Voronoi results may be unreliable.", file=sys.stderr)

    # 2) Voronoi decomposition
    voro = VoronoiAnalysisModifier(
        use_radii=args.use_radii,   # Set True only for radical (radius-weighted) Voronoi
        # Optional knobs: compute_indices=False, edge_length_threshold=..., face_area_threshold=...
    )
    pipeline.modifiers.append(voro)

    # 3) Effective atomic radius & occupied atomic volume
    radii_map = parse_radii_map(args.radii)

    # Evaluate once to obtain ParticleType info present in this frame
    data_tmp = pipeline.compute()
    radius_expr = build_radius_expression(data_tmp, radii_map, fallback=args.fallback_radius)

    # Compute effective radius from the piecewise expression
    pipeline.modifiers.append(
        ComputePropertyModifier(
            output_property="Radius_effective",
            expressions=[radius_expr]
        )
    )
    # Compute occupied atomic volume assuming a sphere: (4/3)*pi*r^3
    pipeline.modifiers.append(
        ComputePropertyModifier(
            output_property="V_atom",
            expressions=[f"(4.0/3.0)*{pi}*(Radius_effective^3)"]   # volume of a sphere
        )
    )

    # 4) Free volume = Voronoi cell volume - occupied atomic volume
    free_expr = "AtomicVolume - V_atom"
    if args.clip_negative:
        free_expr = f"max({free_expr}, 0.0)"

    pipeline.modifiers.append(
        ComputePropertyModifier(
            output_property="free_volume",
            expressions=[free_expr]
        )
    )

    # 5) Compute & basic statistics
    data = pipeline.compute()
    free = np.asarray(data.particles["free_volume"])
    atomic_vol = np.asarray(data.particles["Atomic Volume"])
    v_atom = np.asarray(data.particles["V_atom"])

    N = len(free)
    mean_fv = float(free.mean())
    std_fv = float(free.std())
    min_fv = float(free.min())
    max_fv = float(free.max())

    print(f"# Atoms: {N}")
    print(f"# mean(free V) = {mean_fv:.6g}")
    print(f"# std (free V) = {std_fv:.6g}")
    print(f"# min / max    = {min_fv:.6g} / {max_fv:.6g}")

    # 6) Histogram over free volume (robust min/max by percentiles unless overridden)
    if args.hist_min is None:
        hist_min = float(np.percentile(free, 0.5))  # a bit more robust than raw min
    else:
        hist_min = args.hist_min
    if args.hist_max is None:
        hist_max = float(np.percentile(free, 99.5))
        if hist_max <= hist_min:
            hist_max = float(free.max())
    else:
        hist_max = args.hist_max

    hist, edges = np.histogram(free, bins=args.bins, range=(hist_min, hist_max), density=False)
    centers = 0.5 * (edges[:-1] + edges[1:])
    hist_out = np.column_stack([centers, hist])

    np.savetxt(f"{args.out}_freeV_hist.csv", hist_out, delimiter=",", header="center,count", comments="")
    with open(f"{args.out}_stats.json", "w") as f:
        json.dump({
            "N": int(N),
            "mean_free_volume": mean_fv,
            "std_free_volume": std_fv,
            "min_free_volume": min_fv,
            "max_free_volume": max_fv,
            "hist_bins": int(args.bins),
            "hist_min": hist_min,
            "hist_max": hist_max,
            "use_radii": bool(args.use_radii),
            "clip_negative": bool(args.clip_negative),
        }, f, indent=2)

    print(f"[Saved] Histogram: {args.out}_freeV_hist.csv")
    print(f"[Saved] Stats:     {args.out}_stats.json")

    # 7) Optional per-atom CSV dump
    if args.export_per_atom:
        keys = list(data.particles.keys())
        if "AtomicVolume" in keys:
            atomic_vol_key = "AtomicVolume"
        elif "Atomic Volume" in keys:
            atomic_vol_key = "Atomic Volume"
        else:
            raise RuntimeError(f"Could not find the Voronoi volume property. Available keys: {keys}")

        out_csv = f"{args.out}_per_atom.csv"
        export_per_atom_csv(data, atomic_vol_key, out_csv)
        print(f"[Saved] Per-atom CSV: {out_csv}")
        # Columns written: ID, Type, X,Y,Z, AtomicVolume, V_atom, free_volume, Radius_effective

if __name__ == "__main__":
    main()

