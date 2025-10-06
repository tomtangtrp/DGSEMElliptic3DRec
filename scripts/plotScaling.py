#!/usr/bin/env python3
# Plot log-log OMP scaling from one or many CSVs produced by ompScaling.py (new schema).
# - Accepts CSV columns: threads, time, speedup, run1, run2, run3
# - Backward-compatible with old column 'mean_wall_time'
# - Filename pattern example:
#   testSolver_sp_omp_scaling_testcase0_method==SIPG_Nel=16x16x16_k=3.csv
#
# Title is built from detected key==val tokens (method, Nel, k, testcase if present).
# For multiple files, we draw one curve per file; if their metadata differ, we still plot them
# and show short labels derived from filename stems.
#
# Usage:
#   python plotScaling_v2.py file1.csv [file2.csv ...] [--out out.png] [--title-prefix "OMP scaling"]
#
import argparse
import csv
import re
from pathlib import Path
import matplotlib.pyplot as plt

def read_series(csv_path: Path):
    threads, times = [], []
    with open(csv_path, "r", newline="") as f:
        r = csv.DictReader(f)
        for row in r:
            try:
                t = int(row["threads"])
            except Exception:
                continue
            # prefer 'time' if present, else fall back
            val = None
            for key in ("time", "mean_wall_time", "wall_time", "walltime"):
                if key in row and row[key] not in (None, ""):
                    try:
                        val = float(row[key])
                        break
                    except Exception:
                        pass
            if val is None:
                continue
            if str(val).lower() == "nan":
                continue
            threads.append(t)
            times.append(val)
    z = sorted(zip(threads, times))
    return [t for t, _ in z], [m for _, m in z]

def parse_metadata_from_stem(stem: str):
    # find tokens after 'omp_scaling_' if present
    md = {}
    m = re.search(r'omp_scaling_(.*)$', stem)
    tail = m.group(1) if m else stem
    # split by underscore, look for key==value tokens
    for tok in tail.split('_'):
        if '==' in tok:
            k, v = tok.split('==', 1)
            md[k] = v
        elif '=' in tok:
            k, v = tok.split('=', 1)
            md[k] = v
        elif tok.startswith('testcase'):
            md.setdefault('testcase', tok.replace('testcase', ''))
    return md

def compact_title(md: dict):
    parts = []
    if 'testcase' in md and md['testcase'] != '':
        parts.append(f"testcase{md['testcase']}")
    if 'method' in md:
        parts.append(md['method'])
    if 'Nel' in md:
        parts.append(f"Nel={md['Nel']}")
    if 'k' in md:
        parts.append(f"k={md['k']}")
    return ", ".join(parts) if parts else ""

def main():
    ap = argparse.ArgumentParser(
        description="Plot log-log OMP scaling from CSV(s) with new schema")
    ap.add_argument("csvs", nargs="+",
                    help="CSV files like '..._omp_scaling_testcase0_method==SIPG_...csv'")
    ap.add_argument("--out", help="Output PNG filename (optional)")
    ap.add_argument("--title-prefix", default="OMP scaling",
                    help="Title prefix (optional)")
    args = ap.parse_args()

    files = [Path(x) for x in args.csvs]
    if len(files) < 1:
        raise SystemExit("Provide at least one CSV file.")

    plt.figure(figsize=(7, 5))
    all_threads = set()

    # Collect metadata to build a shared title when possible
    md_list = [parse_metadata_from_stem(p.stem) for p in files]
    base_title = compact_title(md_list[0]) if md_list else ""

    for p, md in zip(files, md_list):
        th, tm = read_series(p)
        if not th:
            print(f"[WARN] No valid data in {p}")
            continue
        all_threads.update(th)
        # Short label: prefer explicit 'label' key in md, else filename stem
        label = md.get('label') or p.stem
        plt.loglog(th, tm, "--o", label=label)

    xt = sorted(all_threads)
    ax = plt.gca()
    ax.set_xscale("log"); ax.set_yscale("log")
    if xt:
        ax.set_xticks(xt)
        ax.set_xticklabels([str(t) for t in xt])
    plt.xlabel("Threads")
    plt.ylabel("Wall time (seconds)")
    full_title = f"{args.title_prefix}: {base_title}" if base_title else args.title_prefix
    plt.title(full_title)
    plt.grid(True, which="both", alpha=0.3)
    if len(files) > 1:
        plt.legend()

    # Default output name
    if args.out:
        outpng = args.out
    else:
        tail = base_title.replace(", ", "_") if base_title else files[0].stem
        outpng = f"scaling_{tail}.png"

    plt.savefig(outpng, dpi=150, bbox_inches="tight")
    print(f"Wrote {outpng}")

if __name__ == "__main__":
    main()
