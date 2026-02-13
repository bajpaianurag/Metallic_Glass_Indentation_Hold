#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OVITO 3.13 — D2min + histogram + STZ tagging
+ STZ-only export to IMD/XYZ (native OVITO exporters)
+ Optional STZ-only export to AtomEye CFG via ASE (if installed)
"""
import argparse, csv, os, numpy as np
from ovito.io import import_file, export_file
from ovito.modifiers import AtomicStrainModifier, ExpressionSelectionModifier, DeleteSelectedModifier
from ovito.pipeline import FileSource
from ovito.data import DataCollection

def build_pipeline(post_path, pre_path, cutoff, mic, affine_mode, bottom_z):
    pipe = import_file(post_path)

    # (optional) remove bottom fixed layer by absolute Z
    if bottom_z is not None:
        sel = ExpressionSelectionModifier(expression=f"Position.Z < {bottom_z}")
        pipe.modifiers.append(sel)
        pipe.modifiers.append(DeleteSelectedModifier())

    mod = AtomicStrainModifier(
        cutoff=cutoff,
        minimum_image_convention=mic,
        output_nonaffine_squared_displacements=True,
        select_invalid_particles=True
    )

    amap = (affine_mode or "off").lower()
    if   amap == "off":
        mod.affine_mapping = AtomicStrainModifier.AffineMapping.Off
    elif amap in ("to_reference","to-reference","toref"):
        mod.affine_mapping = AtomicStrainModifier.AffineMapping.ToReference
    elif amap in ("to_current","to-current","tocur"):
        mod.affine_mapping = AtomicStrainModifier.AffineMapping.ToCurrent
    else:
        raise SystemExit(f"Unknown affine mode: {affine_mode}")

    # external reference snapshot
    mod.reference = FileSource(); mod.reference.load(pre_path)
    pipe.modifiers.append(mod)
    return pipe

def write_csv_with_stz(data: DataCollection, path_csv: str, stz_mask: np.ndarray):
    parts = data.particles
    has_id   = "Particle Identifier" in parts
    has_type = "Particle Type" in parts
    if "Nonaffine Squared Displacement" not in parts:
        raise SystemExit("D2min property not found.")
    N   = parts.count
    pos = parts.positions
    d2  = parts["Nonaffine Squared Displacement"].array
    ids = parts["Particle Identifier"].array if has_id else np.arange(N)
    ty  = parts["Particle Type"].array if has_type else None

    with open(path_csv, "w", newline="") as f:
        w = csv.writer(f)
        header = (["id"] + (["type"] if has_type else []) + ["x","y","z","D2min","STZ"])
        w.writerow(header)
        for i in range(N):
            row = [int(ids[i])]
            if has_type: row.append(int(ty[i]))
            row += [pos[i,0], pos[i,1], pos[i,2], float(d2[i]), int(stz_mask[i])]
            w.writerow(row)

def save_histogram(d2, out_csv, bins=80, log_bins=False, hi_percentile=99.9):
    d2 = np.asarray(d2, float)
    d2 = d2[np.isfinite(d2) & (d2 >= 0.0)]
    if d2.size == 0:
        raise SystemExit("No finite non-negative D2min values to histogram.")

    d2_max = np.percentile(d2, hi_percentile)
    if log_bins:
        d2_pos = d2[d2 > 0]
        vmin = d2_pos.min() if d2_pos.size else 1e-12
        vmin = max(vmin, d2_max*1e-12)
        edges = np.logspace(np.log10(vmin), np.log10(d2_max), bins+1)
    else:
        edges = np.linspace(0.0, d2_max, bins+1)

    counts, edges = np.histogram(d2, bins=edges)
    centers = 0.5*(edges[:-1] + edges[1:])
    width   = edges[1:] - edges[:-1]
    prob    = counts / counts.sum() if counts.sum() > 0 else np.zeros_like(counts)
    dens    = prob / width

    with open(out_csv, "w", newline="") as f:
        w = csv.writer(f)
        w.writerow(["bin_left","bin_right","bin_center","count","probability","density"])
        for i in range(len(counts)):
            w.writerow([edges[i], edges[i+1], centers[i], int(counts[i]), prob[i], dens[i]])

def export_stz_only(data: DataCollection, args):
    """Create STZ-only snapshot and export as requested formats."""
    # Clone to avoid touching the original collection
    stz_data = data.clone()
    # Keep only STZ==1 atoms
    stz_data.apply(ExpressionSelectionModifier(expression="STZ == 0"))
    stz_data.apply(DeleteSelectedModifier())

    # Export IMD
    if args.stz_imd is not None:
        export_file(
            stz_data, args.stz_imd, "imd",
            columns=[
                "Particle Identifier", "Particle Type",
                "Position.X", "Position.Y", "Position.Z",
                "Nonaffine Squared Displacement"
            ],
            multiple_frames=False
        )
    # Export XYZ
    if args.stz_xyz is not None:
        export_file(
            stz_data, args.stz_xyz, "xyz",
            columns=[
                "Particle Identifier", "Particle Type",
                "Position.X", "Position.Y", "Position.Z",
                "Nonaffine Squared Displacement"
            ],
            multiple_frames=False
        )
    # Optional: Export AtomEye CFG via ASE if available
    if args.stz_cfg is not None:
        try:
            from ovito.io.ase import ovito_to_ase
            from ase.io import write as ase_write
            atoms = ovito_to_ase(stz_data)
            # Attach useful arrays to ASE object (AtomEye can read auxiliary via extended CFG)
            d2 = stz_data.particles["Nonaffine Squared Displacement"].array
            atoms.new_array("D2min", np.asarray(d2))
            # Optionally pass IDs/types if present
            if "Particle Identifier" in stz_data.particles:
                atoms.new_array("id", np.asarray(stz_data.particles["Particle Identifier"].array, dtype=int))
            if "Particle Type" in stz_data.particles:
                atoms.new_array("type", np.asarray(stz_data.particles["Particle Type"].array, dtype=int))
            # Write CFG
            ase_write(args.stz_cfg, atoms, format="cfg")
        except ImportError as e:
            raise SystemExit("ASE not installed. Install with: pip install ase") from e

def main():
    ap = argparse.ArgumentParser(description="OVITO D2min + histogram + STZ tagging + STZ-only export")
    ap.add_argument("--pre",  default="../dump.pre")
    ap.add_argument("--post", default="../dump.post")
    ap.add_argument("--cutoff",  type=float, default=4.0)   # Å
    ap.add_argument("--no-mic",  action="store_true")
    ap.add_argument("--affine",  default="off")             # off|to_reference|to_current
    ap.add_argument("--bottom-z", type=float, default=None) # e.g., 10.0
    ap.add_argument("--csv",   default="d2min.csv")
    ap.add_argument("--dump",  default="d2min.post.dump")
    # histogram
    ap.add_argument("--hist",  default="d2min_hist.csv")
    ap.add_argument("--hist-bins", type=int, default=100)
    ap.add_argument("--hist-log",  action="store_true")
    ap.add_argument("--hist-hi",   type=float, default=99.9)
    # STZ tagging
    ap.add_argument("--stz-thresh", type=float, default=1.5, help="D2min threshold for STZ tagging (Å^2)")
    ap.add_argument("--stz-ids",   default="stz_ids.txt")
    # STZ-only export targets (any may be None)
    ap.add_argument("--stz-xyz", default="stz_only.xyz")
    ap.add_argument("--stz-imd", default=None)
    ap.add_argument("--stz-cfg", default=None)
    args = ap.parse_args()

    for p in (args.pre, args.post):
        if not os.path.exists(p):
            raise SystemExit(f"Missing file: {p}")

    pipe = build_pipeline(
        post_path=args.post,
        pre_path=args.pre,
        cutoff=args.cutoff,
        mic=(not args.no_mic),
        affine_mode=args.affine,
        bottom_z=args.bottom_z
    )

    # compute
    data = pipe.compute()
    parts = data.particles
    if "Nonaffine Squared Displacement" not in parts:
        raise SystemExit("D2min property missing. Check IDs & reference.")
    d2 = parts["Nonaffine Squared Displacement"].array
    N  = parts.count

    # STZ mask
    stz_mask = (d2 >= float(args.stz_thresh)).astype(np.int8)

    # Attach STZ property and export full dump (with STZ)
    data.particles_.create_property("STZ", data=stz_mask.astype(int))

    # CSV
    write_csv_with_stz(data, args.csv, stz_mask)

    # STZ id list
    has_id = "Particle Identifier" in parts
    ids = parts["Particle Identifier"].array if has_id else np.arange(N)
    stz_ids = [int(ids[i]) for i in range(N) if stz_mask[i] == 1]
    with open(args.stz_ids, "w") as f:
        for idv in stz_ids:
            f.write(f"{idv}\n")

    # histogram
    save_histogram(d2, args.hist, bins=args.hist_bins, log_bins=args.hist_log, hi_percentile=args.hist_hi)

    # Export full system dump (with STZ column)
    export_file(
        data, args.dump, "lammps/dump",
        columns=[
            "Particle Identifier", "Particle Type",
            "Position.X", "Position.Y", "Position.Z",
            "Nonaffine Squared Displacement", "STZ"
        ],
        multiple_frames=False
    )

    # Export STZ-only snapshots (IMD/XYZ/CFG)
    export_stz_only(data, args)

    print("Done.")
    print(f"  CSV         : {args.csv}")
    print(f"  LAMMPS dump : {args.dump}")
    print(f"  Hist        : {args.hist}")
    print(f"  STZ ids     : {args.stz_ids}")
    if args.stz_imd: print(f"  STZ-only IMD: {args.stz_imd}")
    if args.stz_xyz: print(f"  STZ-only XYZ: {args.stz_xyz}")
    if args.stz_cfg: print(f"  STZ-only CFG: {args.stz_cfg}")

if __name__ == "__main__":
    main()
