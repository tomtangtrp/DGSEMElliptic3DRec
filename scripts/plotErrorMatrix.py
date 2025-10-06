#!/usr/bin/env python3
"""
Plot 3D DG convergence directly from saved error matrices.

Accepts either:
  1) METHOD + CASENAME  (will read: plots/<case>_3D_convergence_<METHOD>.h5)
  2) --file <path-to-summary.h5>

The summary file must contain datasets:
  Nel_1d, ks, L2, brokenL2, brokenH1
and (optionally) group /Meta with attrs 'case' and 'method'.

Outputs (same naming as plotConvergence.py):
  plots/<case>_3D_L2_<METHOD>.png
  plots/<case>_3D_brokenL2_<METHOD>.png
  plots/<case>_3D_brokenH1_<METHOD>.png
"""

import argparse
import os
import sys
import h5py
import numpy as np
import matplotlib.pyplot as plt

PLOT_DIR = "../plots"
os.makedirs(PLOT_DIR, exist_ok=True)

def plot_one(method, casename, xvals, M, title, ylabel, suffix, ks):
    # choose colors for however many ks we have
    cmap = plt.get_cmap("tab10")
    colors = [cmap(i % 10) for i in range(len(ks))]

    fig, ax = plt.subplots(figsize=(6.1, 4.2), dpi=130)
    drew = False
    for i, k in enumerate(ks):
        y = M[i, :]
        mask = np.isfinite(y) & (y > 0)
        if not np.any(mask):
            continue
        drew = True
        ax.loglog(np.array(xvals)[mask], y[mask], "-o",
                  color=colors[i], linewidth=1.8, markersize=5, label=f"k={int(k)}")
        # slope guide: match your plotConvergence.py (SIPG/IIPG: i+2; NIPG/NIPG0: i+1)
        if method in ("SIPG", "IIPG"):
            expo = i + 2
            scale = (2.5 ** (i + 2)) * y[mask][0]
        else:
            expo = i + 1
            scale = (2.0 ** (i + 2)) * y[mask][0]
        x = np.array(xvals, dtype=float)
        ax.loglog(x, scale * (1.0 / x) ** expo, "--", color=colors[i], linewidth=1.0,
                  label=fr"$O(h^{expo})$")

    ax.set_xlabel("Nel_1d")
    ax.set_ylabel(ylabel)
    ax.set_title(f"{method} · {casename} · {title}")
    if drew:
        ax.legend(frameon=True, loc="right", bbox_to_anchor=(1.28, 0.5))
    fig.tight_layout()

    out_png = os.path.join(PLOT_DIR, f"{casename}_3D_{suffix}_{method}.png")
    fig.savefig(out_png)
    plt.close(fig)
    print(f"[saved] {out_png}")

def main():
    ap = argparse.ArgumentParser(description="Plot 3D DG convergence from saved HDF5 matrices")
    ap.add_argument("method", nargs="?", help="SIPG|IIPG|NIPG|NIPG0 (omit if using --file)")
    ap.add_argument("casename", nargs="?", help="case name (omit if using --file)")
    ap.add_argument("--file", help="Direct path to summary HDF5 (plots/<case>_3D_convergence_<METHOD>.h5)")
    ap.add_argument("--show", action="store_true", help="Show figures after saving")
    args = ap.parse_args()

    # Resolve input file
    if args.file:
        summary_path = args.file
        if not os.path.isfile(summary_path):
            print(f"ERROR: file not found: {summary_path}")
            sys.exit(1)
        method = args.method
        casename = args.casename
    else:
        if not (args.method and args.casename):
            print("Usage: plotErrorMatrix.py METHOD CASENAME  OR  plotErrorMatrix.py --file <summary.h5>")
            sys.exit(1)
        method = args.method.upper()
        casename = args.casename
        summary_path = os.path.join(PLOT_DIR, f"{casename}_3D_convergence_{method}.h5")
        if not os.path.isfile(summary_path):
            print(f"ERROR: summary file not found: {summary_path}")
            print("Hint: run plotConvergence.py with --save summary to create it.")
            sys.exit(1)

    # Load matrices
    with h5py.File(summary_path, "r") as h:
        Nel_1d = np.array(h["Nel_1d"])
        ks     = np.array(h["ks"])
        M_L2   = np.array(h["L2"])
        M_bL2  = np.array(h["brokenL2"])
        M_bH1  = np.array(h["brokenH1"])
        if method is None or casename is None:
            # try to read from /Meta if not provided
            if "Meta" in h:
                meta = h["Meta"]
                if method is None:  method  = str(meta.attrs.get("method", "UNKNOWN"))
                if casename is None: casename = str(meta.attrs.get("case", "case"))
            else:
                # fallback generic names
                method = method or "METHOD"
                casename = casename or "case"

    # Plots (match plotConvergence.py)
    plot_one(method, casename, Nel_1d, M_L2,
             "convergence (relative vector L2)", "relative ||u-uh||_2", "L2", ks)
    plot_one(method, casename, Nel_1d, M_bL2,
             "convergence (relative broken L2)", "relative broken L2", "brokenL2", ks)
    plot_one(method, casename, Nel_1d, M_bH1,
             "convergence (relative broken H1)", "relative broken H1 (semi)", "brokenH1", ks)

    if args.show:
        plt.show()

if __name__ == "__main__":
    main()
