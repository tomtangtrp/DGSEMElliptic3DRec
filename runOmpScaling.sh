#!/usr/bin/env bash
set -euo pipefail

# Wrapper for DG3D/ompScaling.py targeting the 3D solver testSolver_sp
# Usage: override env vars as needed: EXE, CASE_JSON, THREADS, REPS, MODE, NO_HT, PIN, TEE, LOG_DIR

EXE="${EXE:-./testSolver_sp}"
CASE_JSON="${CASE_JSON:-./cases/testcase0.json}"

WRITE="${WRITE:-nw}"                 # w | nw
METHOD="${METHOD:-SIPG}"             # SIPG | IIPG | NIPG | NIPG0
K="${K:-3}"
SIGMA0="${SIGMA0:-$((K*K+1))}"
NELX="${NELX:-4}"
NELY="${NELY:-4}"
NELZ="${NELZ:-4}"
XLA="${XLA:-0}"; XLB="${XLB:-1}"
YLA="${YLA:-0}"; YLB="${YLB:-1}"
ZLA="${ZLA:-0}"; ZLB="${ZLB:-1}"
BC_SPEC="${BC_SPEC:-l=dc,r=dc,b=dc,t=dc,back=dc,front=dc}"

THREADS="${THREADS:-auto}"
REPS="${REPS:-3}"
MODE="${MODE:-min}"                  # min | avg
NO_HT="${NO_HT:-1}"                  # 1 = physical cores only
PIN="${PIN:-0}"                      # 1 = set OMP pinning env
SHOW="${SHOW:-0}"
TEE="${TEE:-1}"
LOG_DIR="${LOG_DIR:-}"
PYTHON="${PYTHON:-python3}"
export DG_OMP_DEFAULT_NO_HT=1

CMD=( "${PYTHON}" "$(dirname "$0")/scripts/ompScaling.py"
  --exe "${EXE}"
  --case "${CASE_JSON}"
  --method "${METHOD}"
  -k "${K}"
  --sigma0 "${SIGMA0}"
  --nelx "${NELX}" --nely "${NELY}" --nelz "${NELZ}"
  --xla "${XLA}" --xlb "${XLB}" --yla "${YLA}" --ylb "${YLB}" --zla "${ZLA}" --zlb "${ZLB}"
  --bc "${BC_SPEC}"
  --threads "${THREADS}"
  --reps "${REPS}"
  --mode "${MODE}"
)
if [[ "${WRITE}" == "w" ]]; then
  CMD+=( --write )
fi
if [[ "${NO_HT}" == "1" ]]; then
  CMD+=( --no-ht )
fi
if [[ "${PIN}" == "1" ]]; then
  CMD+=( --pin )
fi
if [[ "${SHOW}" == "1" ]]; then
  CMD+=( --show )
fi
if [[ "${TEE}" == "1" ]]; then
  CMD+=( --tee )
else
  CMD+=( --no-tee )
fi
if [[ -n "${LOG_DIR}" ]]; then
  CMD+=( --log-dir "${LOG_DIR}" )
fi

echo ">>> ${CMD[*]}"
exec "${CMD[@]}"
