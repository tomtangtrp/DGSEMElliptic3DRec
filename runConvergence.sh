#!/usr/bin/env bash
set -euo pipefail

# Usage: ./runConvergence.sh METHOD CASENAME [K_LIST] [N1D_LIST] [BC_SPEC] [SAVE_MODE]
#   METHOD    : SIPG|IIPG|NIPG|NIPG0
#   CASENAME  : e.g. gaussian_iso  (uses cases/<CASENAME>.json)
#   K_LIST    : optional space-separated (default: "1 2 3 4 5")
#   N1D_LIST  : optional space-separated (default: "2 3 4 5 6")
#   BC_SPEC   : optional (default: "L=D,R=D,B=D,T=D,back=D,front=D")
#   SAVE_MODE : optional {summary|samples|both} → forwarded to plotConvergence.py as --save

if [[ $# -lt 2 || $# -gt 6 ]]; then
  echo "Usage: $0 METHOD CASENAME [\"K_LIST\"] [\"N1D_LIST\"] [BC_SPEC] [SAVE_MODE]"
  exit 1
fi

METHOD="$1"
CASENAME="$2"
K_LIST=${3:-"1 2 3 4 5"}
N1D_LIST=${4:-"2 3 4 5 6"}
BC_SPEC=${5:-"L=D,R=D,B=D,T=D,back=D,front=D"}
#SAVE_MODE=${6:-}   # "", or one of: summary | samples | both
SAVE_MODE=${6:-"summary"}    # "", or one of: summary | samples | both
case "$METHOD" in SIPG|IIPG|NIPG|NIPG0) ;; *)
  echo "ERROR: METHOD must be one of SIPG|IIPG|NIPG|NIPG0"
  exit 1
  ;;
esac

if [[ -n "$SAVE_MODE" && "$SAVE_MODE" != "summary" && "$SAVE_MODE" != "samples" && "$SAVE_MODE" != "both" ]]; then
  echo "WARNING: SAVE_MODE must be {summary|samples|both}. Ignoring: '$SAVE_MODE'"
  SAVE_MODE=""
fi

BIN="./testSolver_sp"
WRITE="w"
#SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
#DATA_DIR="$SCRIPT_DIR/data"
#PLOTS_DIR="$SCRIPT_DIR/plots"

SCRIPT_DIR="./scripts"
DATA_DIR="./data"
PLOTS_DIR="./plots"
#mkdir -p "$DATA_DIR" "$PLOTS_DIR"

# Fixed plot script path (no auto-detect)
PLOT_SCRIPT="$SCRIPT_DIR/plotConvergence.py"

# Domain
XLA=0; XLB=1; YLA=0; YLB=1; ZLA=0; ZLB=1

# Case file
#CASE_JSON="$SCRIPT_DIR/cases/${CASENAME}.json"
CASE_JSON="./cases/${CASENAME}.json"
if [[ ! -f "$CASE_JSON" ]]; then
  echo "ERROR: missing case file: $CASE_JSON"
  exit 1
fi

echo "=== 3D convergence sweep: METHOD=${METHOD}, CASE=${CASENAME} ==="
for K in $K_LIST; do
  SIGMA0=$(( K*K + 1 ))   # base penalty; solver may override for NIPG/NIPG0
  for N1D in $N1D_LIST; do
    NELX=$N1D; NELY=$N1D; NELZ=$N1D
    echo ">>> k=${K}, Nelx=Nely=Nelz=${N1D} (Nel=$((N1D*N1D*N1D)))"
    "$BIN" "$WRITE" "$METHOD" "$K" "$SIGMA0" \
      "$NELX" "$NELY" "$NELZ" \
      "$XLA" "$XLB" "$YLA" "$YLB" "$ZLA" "$ZLB" \
      "$BC_SPEC" "$CASE_JSON"
  done
done

# Plot (only if the script exists here)
if [[ -f "$PLOT_SCRIPT" ]]; then
  echo "=== Plotting with: $PLOT_SCRIPT  (METHOD=${METHOD}, CASE=${CASENAME}) ==="
  if ! command -v python3 >/dev/null 2>&1; then
    echo "WARNING: python3 not found; skipping plotting."
    exit 0
  fi

  SAVE_ARGS=()
  if [[ -n "$SAVE_MODE" ]]; then
    SAVE_ARGS+=(--save "$SAVE_MODE")
  fi

  python3 "$PLOT_SCRIPT" "$METHOD" "$CASENAME" "${SAVE_ARGS[@]}" || true

  echo "=== Plots generated (if any): ==="
  ls -1 "$PLOTS_DIR/${CASENAME}_3D_"*"_${METHOD}.png" 2>/dev/null || {
    echo "No plots found. Expected HDF5 files like:"
    echo "  $DATA_DIR/DG3DRec_${CASENAME}_method=${METHOD}_*.h5"
  }
else
  echo "NOTE: $PLOT_SCRIPT not found; skipping plotting."
fi

echo "=== Done. ==="
