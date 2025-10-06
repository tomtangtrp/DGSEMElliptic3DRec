#!/usr/bin/env python3
import os, sys, glob, argparse, re
from datetime import datetime
import numpy as np
import h5py
import matplotlib.pyplot as plt
import shutil

DATA_DIR = "./data"
PLOT_DIR = "./plots"
os.makedirs(PLOT_DIR, exist_ok=True)

K_LINES = [1,2,3,4,5]
COLORS  = ["#1f77b4","#ff7f0e","#2ca02c","#d62728","#9467bd"]

def _attr(attrs, key, default=None):
    if key in attrs:
        v = attrs[key]
        try:
            if hasattr(v, "shape") and v.shape == ():
                return v[()]
            if isinstance(v, (bytes, bytearray)):
                try: return v.decode()
                except Exception: return str(v)
            return v
        except Exception:
            return v
    return default

def _ds(grp, name, default=np.nan):
    return float(np.array(grp[name])) if (grp is not None and name in grp) else float(default)

def _parse_sigma_from_filename(path):
    m = re.search(r"_sigma0=([0-9eE+\-\.]+)\.h5$", os.path.basename(path))
    return float(m.group(1)) if m else np.nan

def read_record(path):
    import h5py
    with h5py.File(path, "r") as f:
        gnum = f["/NumericalSolution"]
        method = str(_attr(gnum.attrs, "method", "UNKNOWN")).upper()
        nel    = int(_attr(gnum.attrs, "Nel", 0))
        k      = int(_attr(gnum.attrs, "k",  -1))
        l2  = _ds(gnum, "l2_error")
        bL2 = _ds(gnum, "broken_L2_rel")
        bH1 = _ds(gnum, "broken_H1_rel")
        sigma0_attr = _attr(gnum.attrs, "sigma0", None)
        try:    sigma0 = float(sigma0_attr) if sigma0_attr is not None else _parse_sigma_from_filename(path)
        except: sigma0 = _parse_sigma_from_filename(path)
        G = f["/Grid"] if "/Grid" in f else None
        nelx = int(_attr(G.attrs, "Nel_x")) if G and "Nel_x" in G.attrs else None
        nely = int(_attr(G.attrs, "Nel_y")) if G and "Nel_y" in G.attrs else None
        nelz = int(_attr(G.attrs, "Nel_z")) if G and "Nel_z" in G.attrs else None
    return method, nel, nelx, nely, nelz, k, l2, bL2, bH1, sigma0, path

def infer_nel1d(nel, nelx, nely, nelz):
    if all(isinstance(v,int) for v in (nelx,nely,nelz)) and nelx==nely==nelz:
        return nelx
    return int(round(nel ** (1.0/3.0))) if nel>0 else None

def archive_data_dir():
    ts = datetime.now().strftime("runat_%y%m%d_%H%M")
    dst = os.path.join(DATA_DIR, ts)
    os.makedirs(dst, exist_ok=True)
    for name in os.listdir(DATA_DIR):
        src = os.path.join(DATA_DIR, name)
        if os.path.abspath(src) == os.path.abspath(dst):
            continue
        shutil.move(src, dst)
    print(f"Archived data files to {dst}")

def build_matrix(records, method, nel1d_values, metric):
    idx = {"l2":6, "bL2":7, "bH1":8}[metric]
    M = np.full((len(K_LINES), len(nel1d_values)), np.nan, dtype=float)
    for rec in records:
        m, nel, nx, ny, nz, kval, *_ = rec
        if m != method or kval not in K_LINES: continue
        n1d = infer_nel1d(nel, nx, ny, nz)
        if n1d is None or n1d not in nel1d_values: continue
        i = K_LINES.index(kval); j = nel1d_values.index(n1d)
        M[i,j] = float(rec[idx])
    return M

def slope_guides(ax, xvals, y0, method, kidx):
    if not np.isfinite(y0) or y0 <= 0: return
    if method in ("SIPG","IIPG"):
        expo = kidx + 2; scale = (2.5 ** (kidx+2)) * y0
    else:
        expo = kidx + 1; scale = (2.0 ** (kidx+2)) * y0
    x = np.array(xvals, dtype=float)
    ax.loglog(x, scale * (1.0/x)**expo, "--", color=COLORS[kidx], linewidth=1.0,
              label=fr"$O(h^{expo})$")

def plot_one(method, casename, xvals, M, title, ylabel, suffix):
    fig, ax = plt.subplots(figsize=(6.1,4.2), dpi=130)
    drew = False
    for i, k in enumerate(K_LINES):
        y = M[i,:]
        mask = np.isfinite(y) & (y>0)
        if not np.any(mask): continue
        drew = True
        ax.loglog(np.array(xvals)[mask], y[mask], "-o",
                  color=COLORS[i], linewidth=1.8, markersize=5, label=f"k={k}")
        slope_guides(ax, xvals, y[mask][0], method, i)
    ax.set_xlabel("Nel_1d"); ax.set_ylabel(ylabel)
    ax.set_title(f"{method} · {casename} · {title}")
    if drew:
        ax.legend(frameon=True, loc="right", bbox_to_anchor=(1.28,0.5))
    fig.tight_layout()
    out_png = os.path.join(PLOT_DIR, f"conv_{casename}_{suffix}_{method}.png")
    fig.savefig(out_png); plt.close(fig)
    print(f"[saved] {out_png}")

def save_summary_h5(casename, method, nel1d_vals, M_L2, M_bL2, M_bH1, records):
    out_path = os.path.join(PLOT_DIR, f"{casename}_3D_convergence_{method}.h5")
    with h5py.File(out_path, "w") as h:
        h.create_dataset("Nel_1d", data=np.array(nel1d_vals, dtype=float))
        h.create_dataset("ks", data=np.array(K_LINES, dtype=int))
        h.create_dataset("L2", data=np.array(M_L2,  order="C"))
        h.create_dataset("brokenL2", data=np.array(M_bL2, order="C"))
        h.create_dataset("brokenH1", data=np.array(M_bH1, order="C"))
        meta = h.create_group("Meta")
        meta.attrs["case"] = casename
        meta.attrs["method"] = method
        meta.attrs["n_files"] = len(records)
        meta.attrs["timestamp"] = datetime.now().isoformat(timespec="seconds")
        tbl = meta.create_dataset(
            "samples", shape=(len(records),),
            dtype=h5py.string_dtype(encoding="utf-8")
        )
        tbl[...] = [os.path.basename(r[10]) for r in records]
    print(f"[saved] {out_path}")

def save_per_sample_h5(casename, method, recs):
    outdir = os.path.join(PLOT_DIR, "samples")
    os.makedirs(outdir, exist_ok=True)
    for r in recs:
        method_r, nel, nx, ny, nz, k, l2, bL2, bH1, sigma0, src = r
        sigtxt = ("nan" if (sigma0 is None or not np.isfinite(sigma0)) else f"{sigma0:g}")
        fname = f"sample_{casename}_method={method}_Nel={nel}_k={k}_sigma={sigtxt}.h5"
        out_path = os.path.join(outdir, fname)
        with h5py.File(out_path, "w") as h:
            gE = h.create_group("Errors")
            gE.create_dataset("rel_vec_L2", data=np.array(l2, dtype=float))
            gE.create_dataset("rel_broken_L2", data=np.array(bL2, dtype=float))
            gE.create_dataset("rel_broken_H1", data=np.array(bH1, dtype=float))
            gM = h.create_group("Meta")
            gM.attrs["case"]     = casename
            gM.attrs["method"]   = method
            gM.attrs["Nel"]      = int(nel)
            if nx is not None: gM.attrs["Nel_x"] = int(nx)
            if ny is not None: gM.attrs["Nel_y"] = int(ny)
            if nz is not None: gM.attrs["Nel_z"] = int(nz)
            gM.attrs["k"]        = int(k)
            if np.isfinite(sigma0): gM.attrs["sigma0"] = float(sigma0)
            gM.attrs["source_h5"]= os.path.basename(src)
            gM.attrs["timestamp"]= datetime.now().isoformat(timespec="seconds")
        print(f"[saved] {out_path}")

def main():
    ap = argparse.ArgumentParser(description="3D DG convergence plots + optional HDF5 export")
    ap.add_argument("method", help="SIPG|IIPG|NIPG|NIPG0")
    ap.add_argument("casename", help="case name (without 'cases/' and '.json')")
    ap.add_argument("--save", choices=["summary","samples","both"], help="write HDF5 files")
    args = ap.parse_args()

    method = args.method.upper()
    casename = args.casename
    if method not in ("SIPG","IIPG","NIPG","NIPG0"):
        print("METHOD must be one of SIPG|IIPG|NIPG|NIPG0"); sys.exit(1)

    pattern = os.path.join(DATA_DIR, f"DG3DRec_{casename}_method={method}_*.h5")
    files = sorted(glob.glob(pattern))
    if not files:
        print(f"[no files] pattern: {pattern}")
        for p in sorted(glob.glob(os.path.join(DATA_DIR, '*.h5'))):
            print("  ", os.path.basename(p))
        sys.exit(0)

    # print("[files]")
    # for p in files: print("  ", os.path.basename(p))

    recs = []
    for p in files:
        try:    recs.append(read_record(p))
        except Exception as e: print(f"[skip] {os.path.basename(p)}: {e}")

    if not recs:
        print("[no usable records]"); sys.exit(0)

    nel1d_vals = sorted({infer_nel1d(r[1], r[2], r[3], r[4]) for r in recs
                         if infer_nel1d(r[1], r[2], r[3], r[4]) is not None})
    if not nel1d_vals:
        print("[no Nel_1d inferred]"); sys.exit(0)

    M_L2  = build_matrix(recs, method, nel1d_vals, "l2")
    M_bL2 = build_matrix(recs, method, nel1d_vals, "bL2")
    M_bH1 = build_matrix(recs, method, nel1d_vals, "bH1")

    plot_one(method, casename, nel1d_vals, M_L2,  "convergence (relative vector L2)", "relative ||u-uh||_2", "L2")
    plot_one(method, casename, nel1d_vals, M_bL2, "convergence (relative broken L2)", "relative broken L2", "brokenL2")
    plot_one(method, casename, nel1d_vals, M_bH1, "convergence (relative broken H1)", "relative broken H1 (semi)", "brokenH1")

    if args.save in ("summary","both"):
        save_summary_h5(casename, method, nel1d_vals, M_L2, M_bL2, M_bH1, recs)
    if args.save in ("samples","both"):
        save_per_sample_h5(casename, method, recs)
    archive_data_dir()

if __name__ == "__main__":
    main()
