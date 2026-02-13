#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
OVITO 3.13 — Bond-length / bond-angle histograms from a .cfg structure.

- Uses CutoffNeighborFinder (PBC-aware) to find neighbors and explicitly compute
  distances and angles.
- Writes histogram data to CSV files; you can plot them later as you like.
"""

import numpy as np
import itertools
import argparse
import math

from ovito.io import import_file
from ovito.data import CutoffNeighborFinder


def compute_bond_lengths_and_angles(data, cutoff_len, cutoff_ang=None):
    """
    Compute bond lengths (pair distances) and bond angles (triplets) from an OVITO DataCollection.

    Parameters
    ----------
    data : ovito.data.DataCollection
        Result of pipeline.compute().
    cutoff_len : float
        Cutoff radius (Å) for collecting pairs used in bond-length statistics.
    cutoff_ang : float or None
        Cutoff radius (Å) for angle neighbors; if None, falls back to cutoff_len.

    Returns
    -------
    lengths : np.ndarray
        1D array of pair distances.
    angles_deg : np.ndarray
        1D array of angles in degrees for all neighbor-vector pairs around each center atom.
    """
    if cutoff_ang is None:
        cutoff_ang = cutoff_len

    positions = np.asarray(data.particles.positions, dtype=float)
    n = len(positions)

    # --- 1) Bond lengths (unique i<j pairs) ---
    finder_len = CutoffNeighborFinder(cutoff_len, data)
    lengths = []

    # For each atom i, iterate neighbors j and keep only i<j to avoid duplicates.
    for i in range(n):
        for nb in finder_len.find(i):
            j = int(nb.index)
            if i < j:
                # nb.distance is the minimum-image distance with PBC applied.
                lengths.append(float(nb.distance))

    lengths = np.array(lengths, dtype=float)

    # --- 2) Bond angles (angles between two neighbor vectors around the same center) ---
    finder_ang = CutoffNeighborFinder(cutoff_ang, data)
    angles = []

    for i in range(n):
        # Collect neighbor displacement vectors (minimum-image R_j - R_i).
        neigh_vectors = []
        for nb in finder_ang.find(i):
            v = np.asarray(nb.delta, dtype=float)  # minimum-image vector
            d = np.linalg.norm(v)
            if d > 1e-12:  # numerical guard against zero-length vectors
                neigh_vectors.append(v)

        # For each unordered pair of neighbor vectors (v_a, v_b), compute the angle.
        for va, vb in itertools.combinations(neigh_vectors, 2):
            na = np.linalg.norm(va)
            nb_ = np.linalg.norm(vb)
            if na < 1e-12 or nb_ < 1e-12:
                continue
            cosang = np.dot(va, vb) / (na * nb_)
            # Clamp due to floating-point noise before arccos.
            cosang = max(-1.0, min(1.0, float(cosang)))
            ang = math.degrees(math.acos(cosang))
            angles.append(ang)

    angles = np.array(angles, dtype=float)

    return lengths, angles


def main():
    ap = argparse.ArgumentParser(
        description="Bond length / angle histograms from cfg (OVITO 3.13)."
    )
    ap.add_argument("input", help="Input .cfg file")
    ap.add_argument(
        "--cutoff-length",
        type=float,
        default=3.2,
        help="Cutoff for bond length pairs (Å). Default: 3.2",
    )
    ap.add_argument(
        "--cutoff-angle",
        type=float,
        default=None,
        help="Cutoff for angle neighbors (Å). Default: same as --cutoff-length",
    )
    ap.add_argument(
        "--bins-length",
        type=int,
        default=100,
        help="Number of bins for the length histogram. Default: 100",
    )
    ap.add_argument(
        "--bins-angle",
        type=int,
        default=90,
        help="Number of bins for the angle histogram (0–180°). Default: 90",
    )
    args = ap.parse_args()

    # 0) Load the file and compute a single animation frame (if any).
    pipeline = import_file(args.input)
    data = pipeline.compute()

    # 1) Compute bond lengths and angles using the requested cutoffs.
    lengths, angles = compute_bond_lengths_and_angles(
        data, cutoff_len=args.cutoff_length, cutoff_ang=args.cutoff_angle
    )

    print(f"# pairs (lengths): {len(lengths)}")
    print(f"# triplets (angles): {len(angles)}")

    # 2) Build histograms and save as CSV.
    import matplotlib.pyplot as plt  # imported here so the script can run headless if needed

    # Length histogram.
    len_counts, len_edges = np.histogram(lengths, bins=args.bins_length)
    len_centers = 0.5 * (len_edges[:-1] + len_edges[1:])

    # Angle histogram (recommended 0–180°).
    ang_counts, ang_edges = np.histogram(angles, bins=args.bins_angle, range=(0.0, 180.0))
    ang_centers = 0.5 * (ang_edges[:-1] + ang_edges[1:])

    # Save CSVs: each row = (bin_center, count).
    np.savetxt(
        "hist_bond_length.csv",
        np.c_[len_centers, len_counts],
        delimiter=",",
        header="length_A,count",
        comments="",
    )
    np.savetxt(
        "hist_bond_angle.csv",
        np.c_[ang_centers, ang_counts],
        delimiter=",",
        header="angle_deg,count",
        comments="",
    )
    print("Saved: hist_bond_length.csv, hist_bond_angle.csv")


if __name__ == "__main__":
    main()

