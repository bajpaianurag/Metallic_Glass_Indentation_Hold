#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Post-processing for nano-indentation with LAMMPS fix indent (rigid indenter)
- Handles negative load after lift-off: slope fit uses only positive-load, high-load
  initial-unloading segment; plots clamp P<0 to 0.
- Optionally reads a 'hold' segment between load and unload if '../indentation_out.hold' exists.

Outputs:
  - indent_summary.csv
  - indent_curve.csv (includes Load_clamped for plotting)
  - indent_load_depth.png, unloading_fit.png
"""

from typing import Tuple, List
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import os

# ================== User Parameters ==================
R_ang = 30.0           # Å, spherical tip radius
nu_sample = 0.33       # Poisson's ratio of sample
epsilon = 0.75         # Oliver–Pharr geometry factor (spherical ~0.75)
beta = 1.0             # Shape factor (~1.0 for spherical)

# --- Unloading slope selection (robust to lift-off) ---
min_frac_of_Pmax = 0.20   # use only points with P >= (this)*Pmax (and P>0)
unload_head_frac  = 0.20  # among those, keep only the head (closest to Pmax) fraction
min_unload_points = 12    # require at least this many points for linear fit

# --- Smoothing (for plotting & Pmax detection; fit uses raw by default) ---
smooth_window = 9         # odd integer; set to 1 to disable
use_smoothed_for_Pmax = True
use_smoothed_for_display = True
# =====================================================

# -------- Unit conversions --------
EV_PER_ANG3_TO_GPA = 160.21766208  # 1 eV/Å^3 = 160.21766208 GPa
EV_PER_ANG_TO_N    = 1.602176634e-9
ANG_TO_M           = 1e-10

# -------- I/O paths --------
load_path = "../indentation_out.load"
hold_path = "../indentation_out.hold"      # <-- optional
unload_path = "../indentation_out.unload"
summary_path = "indent_summary.csv"
curve_path = "indent_curve.csv"
ph_path = "indent_load_depth.png"
un_path = "unloading_fit.png"

def read_ave_time_vector(path: str) -> pd.DataFrame:
    df = pd.read_csv(path, sep=r'\s+', header=None, comment="#", engine="python")
    names = ["StepOrTime","Load","Depth","Zmax","Lbox","Eindent","Temp","Fxsum","Fysum","Fzsum"]
    df.columns = names[:df.shape[1]]
    return df

def simple_smooth(y: np.ndarray, win: int) -> np.ndarray:
    if win is None or win <= 1 or win % 2 == 0:
        return y.copy()
    k = win // 2
    ypad = np.r_[np.full(k, y[0]), y, np.full(k, y[-1])]
    kern = np.ones(win) / win
    return np.convolve(ypad, kern, mode="valid")

def contact_area_spherical(R: float, hc: float) -> float:
    disc = 2.0*R*hc - hc*hc
    if disc <= 0.0: 
        return 0.0
    a = np.sqrt(disc)
    return np.pi * a * a

def select_unloading_segment(df_un: pd.DataFrame, Pmax: float) -> pd.DataFrame:
    """Return initial-unloading, positive-load, high-load segment for slope fit."""
    u = df_un.sort_values("Depth", ascending=False).reset_index(drop=True)  # start near max depth
    # Keep only positive loads and above fraction of Pmax
    mask_pos = (u["Load"] > 0.0) & (u["Load"] >= min_frac_of_Pmax * Pmax)
    u_pos = u.loc[mask_pos].copy()
    if len(u_pos) < min_unload_points:
        raise RuntimeError(
            f"Positive unloading points above {min_frac_of_Pmax:.0%}*Pmax are too few "
            f"({len(u_pos)} < {min_unload_points}). Consider lowering 'min_frac_of_Pmax' "
            f"or 'min_unload_points'."
        )
    # Take only the head fraction (closest to Pmax)
    m = max(min_unload_points, int(np.ceil(len(u_pos) * unload_head_frac)))
    return u_pos.iloc[:m].copy()

def concat_series_for_max(dfs: List[pd.DataFrame]) -> Tuple[np.ndarray, np.ndarray]:
    """Concatenate Depth/Load series in order for Pmax detection."""
    loads = []
    depths = []
    for df in dfs:
        if df is None or df.empty:
            continue
        loads.append(df["Load"].to_numpy())
        depths.append(df["Depth"].to_numpy())
    if not loads:
        return np.array([]), np.array([])
    return np.concatenate(loads), np.concatenate(depths)

def main():
    # ---------- Read ----------
    # require load & unload; hold is optional
    for p in (load_path, unload_path):
        if not os.path.exists(p):
            raise SystemExit(f"Missing '{p}'. Run the modified LAMMPS script first.")
    dfL = read_ave_time_vector(load_path)
    dfU = read_ave_time_vector(unload_path)
    dfH = read_ave_time_vector(hold_path) if os.path.exists(hold_path) else None

    # ---------- Optional smoothing for Pmax detection & display ----------
    if use_smoothed_for_Pmax and smooth_window > 1:
        L_all, D_all = concat_series_for_max(
            [dfL[["Load","Depth"]], (dfH[["Load","Depth"]] if dfH is not None else None), dfU[["Load","Depth"]]]
        )
        if L_all.size == 0:
            raise SystemExit("No data available to detect Pmax.")
        P_series_for_max = simple_smooth(L_all, smooth_window)
        D_series_for_max = simple_smooth(D_all, smooth_window)
        idx_pmax = int(np.argmax(P_series_for_max))
        Pmax = float(P_series_for_max[idx_pmax])
        hmax = float(D_series_for_max[idx_pmax])
    else:
        pieces = [dfL[["Depth","Load"]]]
        if dfH is not None:
            pieces.append(dfH[["Depth","Load"]])
        pieces.append(dfU[["Depth","Load"]])
        df_all = pd.concat(pieces, ignore_index=True)
        idx_pmax = int(df_all["Load"].idxmax())
        Pmax = float(df_all.loc[idx_pmax, "Load"])
        hmax = float(df_all.loc[idx_pmax, "Depth"])

    # ---------- Select unloading fit segment (robust to lift-off) ----------
    ufit = select_unloading_segment(dfU[["Depth","Load"]], Pmax)

    # ---------- Linear fit on raw data in that segment ----------
    x = ufit["Depth"].to_numpy()
    y = ufit["Load"].to_numpy()
    a, b = np.polyfit(x, y, 1)    # P ≈ a*h + b  => S0 = a
    S0_eV_per_A2 = float(a)

    # ---------- Oliver–Pharr ----------
    hc = hmax - epsilon * (Pmax / S0_eV_per_A2)
    hc = max(hc, 0.0)
    A = contact_area_spherical(R_ang, hc)

    H_eV_per_A3 = (Pmax / A) if A > 0 else np.nan
    H_GPa = H_eV_per_A3 * EV_PER_ANG3_TO_GPA if np.isfinite(H_eV_per_A3) else np.nan

    S0_SI = S0_eV_per_A2 * (EV_PER_ANG_TO_N / ANG_TO_M)  # N/m
    sqrtA_SI = np.sqrt(A) * ANG_TO_M                     # m
    Er_Pa = (np.sqrt(np.pi) / (2.0 * beta)) * (S0_SI / sqrtA_SI) if sqrtA_SI > 0 else np.nan
    Er_GPa = Er_Pa / 1e9
    Es_GPa = (1.0 - nu_sample**2) * Er_GPa if np.isfinite(Er_GPa) else np.nan

    # ---------- Export curve (with clamped load for plotting) ----------
    def add_clamped(df: pd.DataFrame) -> pd.DataFrame:
        out = df.copy()
        if use_smoothed_for_display and smooth_window > 1:
            out["Load_s"]  = simple_smooth(out["Load"].values,  smooth_window)
            out["Depth_s"] = simple_smooth(out["Depth"].values, smooth_window)
        else:
            out["Load_s"]  = out["Load"]
            out["Depth_s"] = out["Depth"]
        out["Load_clamped"] = np.maximum(out["Load_s"].values, 0.0)
        return out

    dL = add_clamped(dfL); dL["Phase"] = "load"
    dU = add_clamped(dfU); dU["Phase"] = "unload"
    parts = [dL]
    if dfH is not None:
        dH = add_clamped(dfH); dH["Phase"] = "hold"
        parts.append(dH)
    parts.append(dU)
    df_plot = pd.concat(parts, ignore_index=True)

    pd.DataFrame({
        "Depth_A": df_plot["Depth_s"],
        "Load_eV_per_A": df_plot["Load_s"],
        "Load_clamped_eV_per_A": df_plot["Load_clamped"],
        "Phase": df_plot["Phase"]
    }).to_csv(curve_path, index=False)

    # ---------- Save summary ----------
    """
    pd.DataFrame([{
        "R_ang": R_ang, "nu_sample": nu_sample, "epsilon": epsilon, "beta": beta,
        "Pmax_eV_per_A": Pmax, "hmax_A": hmax, "S0_eV_per_A2": S0_eV_per_A2,
        "hc_A": hc, "A_contact_A2": A,
        "Hardness_GPa": H_GPa, "Er_GPa": Er_GPa, "Es_GPa": Es_GPa,
        "unload_points_used": len(ufit),
        "min_frac_of_Pmax": min_frac_of_Pmax, "unload_head_frac": unload_head_frac
    }]).to_csv(summary_path, index=False)
    """

    # ---------- Plots ----------
    # 1) Load–Depth (clamped for visibility)
    plt.figure(figsize=(6,4))
    plt.plot(dL["Depth_s"], dL["Load_clamped"], label="Loading")
    if dfH is not None:
        plt.plot(dH["Depth_s"], dH["Load_clamped"], label="Hold")
    plt.plot(dU["Depth_s"], dU["Load_clamped"], label="Unloading")
    plt.scatter([hmax], [max(Pmax,0)], s=30, label="Pmax")
    plt.xlabel("Depth h (Å)"); plt.ylabel("Load P (eV/Å)")
    plt.title("Load–Depth Curve")
    plt.grid(True, alpha=0.3); plt.legend(); plt.tight_layout()
    plt.savefig(ph_path, dpi=200); plt.close()

    # 2) Unloading fit region (raw + fit line)
    plt.figure(figsize=(6,4))
    plt.plot(dfU["Depth"], dfU["Load"], label="Unloading (raw)")
    plt.scatter(ufit["Depth"], ufit["Load"], s=12, label="Fit segment (P>0, high-load)")
    plt.plot(ufit["Depth"], a*ufit["Depth"] + b, label="Linear fit")
    plt.xlabel("Depth h (Å)"); plt.ylabel("Load P (eV/Å)")
    plt.title("Unloading initial slope")
    plt.grid(True, alpha=0.3); plt.legend(); plt.tight_layout()
    plt.savefig(un_path, dpi=200); plt.close()

    # Console summary
    print(f"Pmax = {Pmax:.4f} eV/Å @ hmax = {hmax:.4f} Å")
    print(f"S0   = {S0_eV_per_A2:.6f} eV/Å²  (fit points: {len(ufit)})")
    print(f"hc   = {hc:.4f} Å,  A = {A:.4f} Å²")
    print(f"H    = {H_GPa:.3f} GPa")
    print(f"Er   = {Er_GPa:.3f} GPa")
    print(f"Es   = {Es_GPa:.3f} GPa  (rigid indenter)")
    present = "yes" if os.path.exists(hold_path) else "no"
    print(f"'hold' segment present: {present}")
    print(f"Saved: {summary_path}, {curve_path}, {ph_path}, {un_path}")

if __name__ == "__main__":
    main()

