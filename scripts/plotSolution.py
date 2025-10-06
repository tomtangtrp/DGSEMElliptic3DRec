#!/usr/bin/env python3
"""
3D DG solution plotter for rectangular grids.

- Reads HDF5 written by testSolver_sp (groups: /Grid, /NumericalSolution, /ExactSolution)
- Reconstructs per-element nodal coordinates (uniform nodes, (k+1)^3)
- Plots 3x subplots (Exact, Numerical, Error) as 3D scatters
- CLI mirrors the 2D plotter (plus Nelz & point size options)

Similar in spirit to DG3DElliptic_RecElement_sparse.py's 3D visualization.
"""

import argparse
import h5py
import numpy as np
import matplotlib.pyplot as plt
import os
from typing import Tuple

# ---------- small helpers ----------
def read_attr(attrs, key, default=None):
    """Robustly read an HDF5 attribute."""
    if key in attrs:
        v = attrs[key]
        try:
            if isinstance(v, (bytes, bytearray)):
                try:
                    return v.decode()
                except Exception:
                    return str(v)
            if hasattr(v, "shape") and v.shape == ():
                return v[()]
            if hasattr(v, "item"):
                return v.item()
            return v
        except Exception:
            return v
    return default

def best_factor_triple(N: int, Lx: float, Ly: float, Lz: float) -> Tuple[int,int,int]:
    """
    Pick (Nelx, Nely, Nelz) so Nelx*Nely*Nelz == N and ratios ~ domain aspect.
    Brute-force over divisors (OK for moderate N).
    """
    if N <= 0:
        return (1, 1, 1)
    # normalize target aspect ratios
    Ax, Ay, Az = Lx, Ly, Lz
    if Ax <= 0: Ax = 1.0
    if Ay <= 0: Ay = 1.0
    if Az <= 0: Az = 1.0
    rx = Ax / Az
    ry = Ay / Az

    best = (N, 1, 1)
    best_err = float("inf")
    for nely in range(1, N + 1):
        if N % nely: 
            continue
        M = N // nely
        for nelx in range(1, M + 1):
            if M % nelx:
                continue
            nelz = M // nelx
            # compare (nelx/nelz) to (Lx/Lz) and (nely/nelz) to (Ly/Lz)
            ex = (nelx / max(1, nelz)) / rx
            ey = (nely / max(1, nelz)) / ry
            # closeness to 1 for both
            err = abs(ex - 1.0) + abs(ey - 1.0)
            if err < best_err:
                best_err = err
                best = (nelx, nely, nelz)
    return best

def element_xyz(ex, ey, ez, Nelx, Nely, Nelz, xLa, xLb, yLa, yLb, zLa, zLb, k):
    """
    Return element-local nodal grids (X,Y,Z) on the physical element [ax,bx]x[ay,by]x[az,bz]
    using uniform nodes (k+1 per axis). Shapes: (k+1, k+1, k+1) with indexing='xy' -> (Ny,Nx,Nz).
    """
    hx = (xLb - xLa) / Nelx
    hy = (yLb - yLa) / Nely
    hz = (zLb - zLa) / Nelz
    ax = xLa + ex * hx; bx = ax + hx
    ay = yLa + ey * hy; by = ay + hy
    az = zLa + ez * hz; bz = az + hz
    k1 = k + 1
    xloc = np.linspace(ax, bx, k1)
    yloc = np.linspace(ay, by, k1)
    zloc = np.linspace(az, bz, k1)
    # indexing='xy' -> shapes (len(y), len(x), len(z)), matching y-x-z flatten (z fastest)
    X, Y, Z = np.meshgrid(xloc, yloc, zloc, indexing="xy")
    return X, Y, Z

def eid_from_exeyez(ex, ey, ez, Nelx, Nely, Nelz) -> int:
    """Row-major, z-fastest element id: e = (ey*Nelx + ex)*Nelz + ez."""
    return (ey * Nelx + ex) * Nelz + ez

# ---------- main plotting ----------
def main():
    ap = argparse.ArgumentParser(description="Plot DG3D rectangular solution (exact, numerical, error) as 3D scatters.")
    ap.add_argument("h5file", help="Path to HDF5 written by the solver")
    ap.add_argument("--nelx", type=int, help="Override Nelx (optional)")
    ap.add_argument("--nely", type=int, help="Override Nely (optional)")
    ap.add_argument("--nelz", type=int, help="Override Nelz (optional)")
    ap.add_argument("--size", type=float, default=10.0, help="Marker size for scatter")
    ap.add_argument("--elev", type=float, default=22.0, help="View elev")
    ap.add_argument("--azim", type=float, default=-58.0, help="View azim")
    ap.add_argument("--save", help="Save figure to this path (PNG)")
    args = ap.parse_args()

    with h5py.File(args.h5file, "r") as hf:
        # --- Grid / domain
        G = hf["/Grid"]
        xLa = float(read_attr(G.attrs, "xLa", 0.0))
        xLb = float(read_attr(G.attrs, "xLb", 1.0))
        yLa = float(read_attr(G.attrs, "yLa", 0.0))
        yLb = float(read_attr(G.attrs, "yLb", 1.0))
        zLa = float(read_attr(G.attrs, "zLa", 0.0))
        zLb = float(read_attr(G.attrs, "zLb", 1.0))
        Lx, Ly, Lz = (xLb - xLa), (yLb - yLa), (zLb - zLa)

        # --- Numerical solution + meta
        NG = hf["/NumericalSolution"]
        Nel  = int(read_attr(NG.attrs, "Nel"))
        k    = int(read_attr(NG.attrs, "k"))
        k1   = k + 1
        locdim = k1 * k1 * k1
        method = read_attr(NG.attrs, "method", "")
        sigma0 = read_attr(NG.attrs, "sigma0", None)

        # retrieve u_h vector
        u_h = np.array(NG["u_h"])
        assert u_h.size == Nel * locdim, f"u_h length {u_h.size} != Nel*(k+1)^3={Nel*locdim}"

        # Determine Nelx/Nely/Nelz: CLI overrides -> attrs -> /Grid -> factor triple
        if args.nelx and args.nely and args.nelz:
            Nelx, Nely, Nelz = args.nelx, args.nely, args.nelz
        else:
            # Prefer attributes written by solver if present
            Nelx_attr = read_attr(NG.attrs, "Nelx")
            Nely_attr = read_attr(NG.attrs, "Nely")
            Nelz_attr = read_attr(NG.attrs, "Nelz")
            if Nelx_attr is not None and Nely_attr is not None and Nelz_attr is not None:
                Nelx, Nely, Nelz = int(Nelx_attr), int(Nely_attr), int(Nelz_attr)
            else:
                Nel_x_grid = read_attr(G.attrs, "Nel_x")
                Nel_y_grid = read_attr(G.attrs, "Nel_y")
                Nel_z_grid = read_attr(G.attrs, "Nel_z")
                if Nel_x_grid is not None and Nel_y_grid is not None and Nel_z_grid is not None:
                    Nelx, Nely, Nelz = int(Nel_x_grid), int(Nel_y_grid), int(Nel_z_grid)
                else:
                    Nelx, Nely, Nelz = best_factor_triple(Nel, Lx, Ly, Lz)

        # --- Exact (optional)
        EX = None
        if "/ExactSolution" in hf and "u_exact" in hf["/ExactSolution"]:
            EX = np.array(hf["/ExactSolution/u_exact"])
            if EX.size != u_h.size:
                print(f"[warn] exact length {EX.size} != numerical {u_h.size}; exact will be ignored.")
                EX = None

        # =========================
        # Assemble global nodal coords & values (per-element)
        # =========================
        pts = Nel * locdim
        Xall = np.empty(pts, dtype=np.float64)
        Yall = np.empty(pts, dtype=np.float64)
        Zall = np.empty(pts, dtype=np.float64)
        Uh   = np.asarray(u_h, dtype=np.float64)
        Uex  = np.asarray(EX,  dtype=np.float64) if EX is not None else None

        idx = 0
        for ey in range(Nely):
            for ex in range(Nelx):
                for ez in range(Nelz):
                    e = eid_from_exeyez(ex, ey, ez, Nelx, Nely, Nelz)
                    # local coords (Ny,Nx,Nz) with z fastest on ravel()
                    X, Y, Z = element_xyz(ex, ey, ez, Nelx, Nely, Nelz, xLa, xLb, yLa, yLb, zLa, zLb, k)
                    nloc = X.size  # k1^3
                    seg = slice(e*locdim, (e+1)*locdim)

                    Xall[idx:idx+nloc] = X.ravel(order="C")
                    Yall[idx:idx+nloc] = Y.ravel(order="C")
                    Zall[idx:idx+nloc] = Z.ravel(order="C")
                    # u_h already in the solver's local ordering (z fastest): use as-is
                    idx += nloc

        # shared color scale across exact & numerical
        if Uex is not None:
            vmin = float(min(Uex.min(), Uh.min()))
            vmax = float(max(Uex.max(), Uh.max()))
        else:
            vmin = float(Uh.min()); vmax = float(Uh.max())

        # =========================
        # Figure: exact, numerical, error scatters
        # =========================
        fig = plt.figure(figsize=(14, 4.6), dpi=140)

        # Exact (if present)
        if Uex is not None:
            ax1 = fig.add_subplot(1, 3, 1, projection="3d")
            ax1.view_init(elev=args.elev, azim=args.azim)
            p1 = ax1.scatter(Xall, Yall, Zall, c=Uex, s=args.size, marker="o", vmin=vmin, vmax=vmax)
            ax1.set_title("Exact")
            ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.set_zlabel("z")
            fig.colorbar(p1, ax=ax1, shrink=0.7, pad=0.05)
        else:
            ax1 = fig.add_subplot(1, 3, 1, projection="3d")
            ax1.view_init(elev=args.elev, azim=args.azim)
            p1 = ax1.scatter(Xall, Yall, Zall, c=Uh, s=args.size, marker="o")
            ax1.set_title("Numerical")
            ax1.set_xlabel("x"); ax1.set_ylabel("y"); ax1.set_zlabel("z")
            fig.colorbar(p1, ax=ax1, shrink=0.7, pad=0.05)

        # Numerical
        ax2 = fig.add_subplot(1, 3, 2, projection="3d")
        ax2.view_init(elev=args.elev, azim=args.azim)
        p2 = ax2.scatter(Xall, Yall, Zall, c=Uh, s=args.size, marker="o", vmin=vmin, vmax=vmax)
        ax2.set_title("Numerical")
        ax2.set_xlabel("x"); ax2.set_ylabel("y"); ax2.set_zlabel("z")
        fig.colorbar(p2, ax=ax2, shrink=0.7, pad=0.05)

        # Error
        ax3 = fig.add_subplot(1, 3, 3, projection="3d")
        ax3.view_init(elev=args.elev, azim=args.azim)
        if Uex is not None:
            err = Uh - Uex
            p3 = ax3.scatter(Xall, Yall, Zall, c=err, s=args.size, marker="o")
            ax3.set_title("Error (u_h - u)")
            fig.colorbar(p3, ax=ax3, shrink=0.7, pad=0.05)
        else:
            p3 = ax3.scatter(Xall, Yall, Zall, c=Uh, s=args.size, marker="o")
            ax3.set_title("Numerical (no exact)")
            fig.colorbar(p3, ax=ax3, shrink=0.7, pad=0.05)
        ax3.set_xlabel("x"); ax3.set_ylabel("y"); ax3.set_zlabel("z")

        # Title banner
        bits = [f"{method}", f"k={k}", f"Nel={Nelx}×{Nely}×{Nelz}",
                f"[{xLa},{xLb}]×[{yLa},{yLb}]×[{zLa},{zLb}]"]
        if sigma0 is not None:
            try:
                bits.insert(1, f"σ₀={float(sigma0)}")
            except Exception:
                bits.insert(1, f"sigma0={sigma0}")
        fig.suptitle("  ".join(bits), y=0.99, fontsize=11)

        plt.tight_layout(rect=[0, 0.00, 1, 0.97])

        # ---------- save/show ----------
        if args.save:
            fig.savefig(args.save, dpi=160, bbox_inches="tight")
        plt.show()

if __name__ == "__main__":
    main()
