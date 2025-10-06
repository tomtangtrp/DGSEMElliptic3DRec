#!/usr/bin/env python3
# OpenMP scaling harness for DG3D testSolver_sp (3D solver)
# - Sweeps thread counts (physical-core aware with --no-ht)
# - Parses "DG Solver Wall time: <sec> seconds"
# - Tees original solver output by default (preserves single-run info)
# - Writes CSV + PNGs

import argparse
import os
import platform
import re
import subprocess
from pathlib import Path
from statistics import mean
from typing import Optional, List, Dict
import matplotlib.pyplot as plt
import csv

WALL_RE = re.compile(r'DG\s+Solver\s+Wall\s*time:\s*([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)\s*(?:seconds?|s)?', re.IGNORECASE)
ALT_WALL_RE = re.compile(r'Wall\s*time[^:]*:\s*([0-9]*\.?[0-9]+(?:[eE][-+]?\d+)?)', re.IGNORECASE)

def exe_tag(exe_path: str) -> str:
    b = os.path.basename(exe_path.rstrip("/"))
    return re.sub(r'[^A-Za-z0-9_.-]+', '_', b)

def detect_logical_cores() -> int:
    try:
        return os.cpu_count() or 1
    except Exception:
        return 1

def threads_per_core() -> int:
    try:
        import platform, subprocess
        if 'linux' in platform.system().lower():
            try:
                out = subprocess.check_output(['lscpu'], text=True, stderr=subprocess.DEVNULL)
                for line in out.splitlines():
                    if line.strip().startswith('Thread(s) per core:'):
                        return int(float(line.split(':',1)[1].strip()))
            except Exception:
                pass
        elif 'darwin' in platform.system().lower():
            # Apple Silicon: effectively 1
            return 1
    except Exception:
        pass
    return 1

def detect_physical_cores() -> int:
    """
    Physical core detection that works on Linux (e.g., NCI Gadi) and macOS.
    Linux: prefer 'lscpu' fields; fallback to unique (CORE,SOCKET) pairs.
    macOS: sysctl hw.physicalcpu.
    Final fallback: logical//2 (but never <1).
    """
    try:
        sysname = platform.system().lower()
        if 'linux' in sysname:
            try:
                out = subprocess.check_output(['lscpu'], text=True, stderr=subprocess.DEVNULL)
                sockets = None
                cores_per_socket = None
                for line in out.splitlines():
                    if line.startswith('Socket(s):'):
                        sockets = int(line.split(':',1)[1].strip())
                    elif line.startswith('Core(s) per socket:'):
                        # some distros show decimals; be defensive
                        cores_per_socket = int(float(line.split(':',1)[1].strip()))
                if sockets and cores_per_socket:
                    return sockets * cores_per_socket
                # CSV fallback
                out = subprocess.check_output(['lscpu', '-p=core,socket,online'], text=True, stderr=subprocess.DEVNULL)
                cores = set()
                for ln in out.splitlines():
                    if ln.startswith('#'):
                        continue
                    parts = ln.split(',')
                    if len(parts) >= 2:
                        # ONLINE may be 'Y'/'N'; if present and 'N', skip
                        if len(parts) >= 3 and parts[2] not in ('', 'Y'):
                            continue
                        cores.add((parts[0], parts[1]))
                if cores:
                    return len(cores)
            except Exception:
                pass
        elif 'darwin' in sysname:
            try:
                out = subprocess.check_output(['sysctl', '-n', 'hw.physicalcpu'], text=True)
                return max(1, int(out.strip()))
            except Exception:
                pass
        logical = detect_logical_cores()
        return max(1, logical // 2)
    except Exception:
        return 1

def default_threads(no_ht: bool) -> List[int]:
    logical = detect_logical_cores()
    physical = detect_physical_cores()
    max_threads = physical if no_ht else logical
    #candidates = [1,2,3,4,6,8,12,16,20,24,28,32,36,40,44,48,56,64,72,80,88,96,112,128]
    candidates = [1,2,4,8,12,16,20,24,28,32,36,40,44,48,56,64,72,80,88,96,112,128]
    threads = [t for t in candidates if t <= max_threads]
    if not threads or threads[-1] != max_threads:
        threads.append(max_threads)
    threads = sorted(set(t for t in threads if t > 0))
    return threads

def parse_threads(s: str) -> List[int]:
    arr: List[int] = []
    for part in s.split(','):
        part = part.strip()
        if not part:
            continue
        if '-' in part:
            a,b = part.split('-',1)
            a = int(a); b = int(b)
            step = 1 if b>=a else -1
            arr.extend(list(range(a,b+step,step)))
        else:
            arr.append(int(part))
    arr = [t for t in arr if t>0]
    return sorted(set(arr))

def build_solver_args(ns: argparse.Namespace) -> List[str]:
    if ns.solver_args:
        return ns.solver_args
    bc = ns.bc or 'l=dc,r=dc,b=dc,t=dc,back=dc,front=dc'
    return [
        'w' if ns.write else 'nw',
        ns.method.upper(),
        str(ns.k),
        f'{ns.sigma0:g}',
        str(ns.nelx),
        str(ns.nely),
        str(ns.nelz),
        f'{ns.xla:g}', f'{ns.xlb:g}',
        f'{ns.yla:g}', f'{ns.ylb:g}',
        f'{ns.zla:g}', f'{ns.zlb:g}',
        bc,
        ns.case
    ]

def run_once(exe: str, solver_args: List[str], extra_env: Dict[str,str], *, tee: bool = True, log_file: Optional[str] = None) -> float:
    env = os.environ.copy()
    env.update(extra_env or {})
    cmd = [exe] + solver_args
    p = subprocess.Popen(cmd, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True, env=env)
    out_lines = []
    log_fh = open(log_file, 'w', encoding='utf-8', errors='replace') if log_file else None
    try:
        for ln in p.stdout:  # type: ignore
            out_lines.append(ln)
            if tee:
                print(ln, end='')
            if log_fh:
                log_fh.write(ln)
        p.wait()
    finally:
        if log_fh:
            log_fh.flush(); log_fh.close()
    out_text = ''.join(out_lines)

    matches = list(WALL_RE.finditer(out_text))
    if not matches:
        matches = list(ALT_WALL_RE.finditer(out_text))
    if not matches:
        tail = ''.join(out_lines[-30:])
        raise RuntimeError(f"Could not parse wall time from solver output.\nLast lines:\n{tail}")
    val = matches[-1].group(1)
    try:
        return float(val)
    except Exception:
        tail = ''.join(out_lines[-30:])
        raise RuntimeError(f"Parsed wall time token '{val}' but could not convert to float.\nLast lines:\n{tail}")

def main():
    ap = argparse.ArgumentParser(description="OpenMP scaling harness for DG3D testSolver_sp (3D).")
    ap.add_argument('--exe', required=True, help='Path to solver executable (e.g., ./testSolver_sp)')
    ap.add_argument('--case', required=True, help='Path to JSON case file')
    ap.add_argument('--method', default='SIPG', choices=['SIPG','IIPG','NIPG','NIPG0'])
    ap.add_argument('-k', '--k', type=int, default=3)
    ap.add_argument('--sigma0', type=float, help='Penalty base (default: k*k+1 for SIPG/IIPG)')
    ap.add_argument('--nelx', type=int, default=8)
    ap.add_argument('--nely', type=int, default=8)
    ap.add_argument('--nelz', type=int, default=8)
    ap.add_argument('--xla', type=float, default=0.0)
    ap.add_argument('--xlb', type=float, default=1.0)
    ap.add_argument('--yla', type=float, default=0.0)
    ap.add_argument('--ylb', type=float, default=1.0)
    ap.add_argument('--zla', type=float, default=0.0)
    ap.add_argument('--zlb', type=float, default=1.0)
    ap.add_argument('--bc', default='l=dc,r=dc,b=dc,t=dc,back=dc,front=dc',
                    help='Comma/semicolon separated, e.g. "l=dc,r=nc,t=dc,back=dc,front=nc"')
    ap.add_argument('--write', action='store_true', help='Ask solver to write HDF5 output (default: false)')
    ap.add_argument('--threads', default='auto', help='Comma list or ranges (e.g. "1,2,4,8-16"); or "auto"')
    ap.add_argument('--reps', type=int, default=3, help='Repetitions per thread count')
    ap.add_argument('--mode', choices=['min','avg'], default='min', help='Aggregate over repetitions')
    ap.add_argument('--no-ht', action='store_true', help='Limit to physical cores only')
    ap.add_argument('--pin', action='store_true', help='Set a conservative OMP/MKL pinning env')
    ap.add_argument('--show', action='store_true', help='Show plots interactively')
    grp = ap.add_mutually_exclusive_group()
    grp.add_argument('--tee', dest='tee', action='store_true', help='Print solver stdout live (default)')
    grp.add_argument('--no-tee', dest='tee', action='store_false', help='Suppress solver stdout (quieter)')
    ap.set_defaults(tee=True)
    ap.add_argument('--log-dir', default=None, help="If set, save each run's solver stdout to this directory")
    ns, rest = ap.parse_known_args()
    ns.solver_args = []
    if rest and rest[0] == '--':
        ns.solver_args = rest[1:]

    if ns.sigma0 is None:
        ns.sigma0 = ns.k*ns.k + 1.0

    if ns.threads == 'auto':
        threads = default_threads(ns.no_ht)
    else:
        threads = parse_threads(ns.threads)
        if ns.no_ht:
            physical = detect_physical_cores()
            threads = [t for t in threads if t <= physical]
    if not threads:
        raise SystemExit("No threads to run.")

    solver_args = build_solver_args(ns)
    exe = ns.exe
    tag = exe_tag(exe)
    case_name = Path(ns.case).stem

    nel_tag = f"Nel={ns.nelx}x{ns.nely}x{ns.nelz}"
    png_dir = "./plots"
    png1 = os.path.join(png_dir, f"{tag}_walltime_vs_threads_{case_name}_method={ns.method}_{nel_tag}_k={ns.k}.png")
    png2 = os.path.join(png_dir, f"{tag}_speedup_vs_threads_{case_name}_method={ns.method}_{nel_tag}_k={ns.k}.png")
    csv_dir = "./scripts"
    csv_path = os.path.join(csv_dir,f"{tag}_omp_scaling_{case_name}_method=={ns.method}_{nel_tag}_k={ns.k}.csv")

    results = []
    for t in threads:
        per = []
        for r in range(ns.reps):
            env = {
                'OMP_NUM_THREADS': str(t),
                'MKL_NUM_THREADS': str(t),
            }
            if ns.pin:
                env.update({
                    'OMP_WAIT_POLICY': 'PASSIVE',
                    'OMP_PROC_BIND': 'true',
                    'OMP_PLACES': 'cores',
                    'KMP_BLOCKTIME': '0',
                })
            log_file = None
            if ns.log_dir:
                Path(ns.log_dir).mkdir(parents=True, exist_ok=True)
                log_stem = f"{tag}_{case_name}_threads={t}_run={r+1}.log"
                log_file = str(Path(ns.log_dir) / log_stem)
            try:
                sec = run_once(exe, solver_args, env, tee=ns.tee, log_file=log_file)
            except Exception as e:
                print(f"[threads={t}] run {r+1}/{ns.reps} failed: {e}")
                continue
            per.append(sec)
            print(f"[threads={t}] run {r+1}/{ns.reps}: {sec:.6f} s")
        if not per:
            continue
        agg = min(per) if ns.mode=='min' else mean(per)
        results.append({'threads': t, 'times': per, 'time': agg})

    if not results:
        raise SystemExit("No successful runs.")

    base = results[0]['time']
    for row in results:
        row['speedup'] = base / row['time'] if row['time']>0 else 0.0

    max_reps = max(len(r['times']) for r in results)
    with open(csv_path, 'w', newline='') as f:
        w = csv.writer(f)
        header = ['threads','time','speedup'] + [f'run{i+1}' for i in range(max_reps)]
        w.writerow(header)
        for r in results:
            row = [r['threads'], f"{r['time']:.9f}", f"{r['speedup']:.6f}"]
            row += [f"{x:.9f}" for x in r['times']]
            w.writerow(row)
    print(f"Wrote {csv_path}")

    xs = [r['threads'] for r in results]
    ys = [r['time'] for r in results]
    plt.figure()
    plt.plot(xs, ys, marker='o')
    plt.xlabel("Threads")
    plt.ylabel("Wall time (s)")
    plt.title(f"OMP scaling — {tag}\ncase={case_name}, method={ns.method}, k={ns.k}, {nel_tag}")
    plt.xticks(xs, [str(x) for x in xs])
    plt.grid(True, alpha=0.3)
    plt.savefig(png1, dpi=150, bbox_inches="tight")
    print(f"Wrote {png1}")

    su = [r['speedup'] for r in results]
    plt.figure()
    plt.plot(xs, su, marker='o')
    plt.xlabel("PhysicalCores")
    plt.ylabel("Speedup (vs first point)")
    plt.title(f"OMP speedup — {tag}\ncase={case_name}, method={ns.method}, k={ns.k}, {nel_tag}")
    plt.xticks(xs, [str(x) for x in xs])
    plt.grid(True, alpha=0.3)
    plt.savefig(png2, dpi=150, bbox_inches="tight")
    print(f"Wrote {png2}")

    if ns.show:
        plt.show()

if __name__ == '__main__':
    main()
