#!/usr/bin/env bash
set -euo pipefail

# Usage: ./runSolver.sh METHOD CASENAME [N1D] [K] [BC_SPEC]
#   METHOD   : SIPG|IIPG|NIPG|NIPG0
#   CASENAME : e.g. gaussian_iso  (the file is cases/gaussian_iso.json)
#   N1D      : optional, default 4  (Nelx=Nely=Nelz=N1D)
#   K        : optional, default 5
#   BC_SPEC  : optional, default "L=D,R=D,B=D,T=D,back=D,front=D"

if [[ $# -lt 2 || $# -gt 5 ]]; then
  echo "Usage: $0 METHOD CASENAME [N1D] [K] [BC_SPEC]"
  exit 1
fi

METHOD="$1"
CASENAME="$2"
N1D="${3:-3}"
K="${4:-3}"
BC_SPEC="${5:-L=D,R=D,B=D,T=D,back=D,front=D}"

case "$METHOD" in SIPG|IIPG|NIPG|NIPG0) ;; *) echo "Bad METHOD: $METHOD"; exit 1;; esac

BIN="./testSolver_sp"
WRITE="w"

# domain
XLA=0; XLB=1; YLA=0; YLB=1; ZLA=0; ZLB=1

# grid
NELX="$N1D"; NELY="$N1D"; NELZ="$N1D"

# penalty base (direct but safe default; solver may override for NIPG/NIPG0)
SIGMA0=$(( K*K + 1 ))

CASE_JSON="cases/${CASENAME}.json"
if [[ ! -f "$CASE_JSON" ]]; then
  echo "ERROR: missing case json: $CASE_JSON"
  exit 1
fi
OUTDIR="./data"

echo ">>> OMP_NUM_THREADS=1 $BIN $WRITE $METHOD $K $SIGMA0 $NELX $NELY $NELZ $XLA $XLB $YLA $YLB $ZLA $ZLB \"$BC_SPEC\" $CASE_JSON"
"$BIN" "$WRITE" "$METHOD" "$K" "$SIGMA0" \
  "$NELX" "$NELY" "$NELZ" \
  "$XLA" "$XLB" "$YLA" "$YLB" "$ZLA" "$ZLB" \
  "$BC_SPEC" "$CASE_JSON"

# ---- plot only if WRITE == "w" ----
if [[ "$WRITE" == "w" ]]; then
  case_name="$(basename "$CASE_JSON")"; case_name="${case_name%.*}"
  nel=$(( NELX * NELY * NELZ ))   # 3D total elements

  # Effective sigma0 for filename (if your solver overrides by method)
  SIGMA0_EFF="$SIGMA0"
  case "$METHOD" in
    NIPG)  SIGMA0_EFF="1" ;;
    NIPG0) SIGMA0_EFF="0" ;;
  esac

  fname="DG3DRec_${case_name}_method=${METHOD}_Nel=${nel}_k=${K}_sigma0=${SIGMA0_EFF}.h5"
  [[ -n "$OUTDIR" ]] && h5path="${OUTDIR%/}/$fname" || h5path="$fname"

  echo ">>> python3 /scripts/plotSolution.py \"$h5path\""
  python3 ./scripts/plotSolution.py "$h5path"
fi
