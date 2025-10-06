#!/usr/bin/env bash
# perf_auto_diag_safe.sh
# Per-core performance + automatic bound verdict for Intel 14900K (Raptor Lake)
# Uses only architectural events so it works on any stock kernel

# Allow perf to access all events (run once per boot if needed):
#   sudo sysctl -w kernel.perf_event_paranoid=-1
#   sudo modprobe msr

OUT=$(sudo perf stat -a --per-core \
    -e cycles \
    -e instructions \
    -e cache-references \
    -e cache-misses \
    -e bus-cycles \
    ./testSolver_sp n NIPG 3 8 8 8 2>&1)

echo "$OUT"

# Extract totals across all cores
cycles=$(echo "$OUT" | awk '/^\s*[0-9,.]+\s+cycles/ {gsub(",","",$1); sum+=$1} END{print sum}')
instr=$(echo "$OUT" | awk '/^\s*[0-9,.]+\s+instructions/ {gsub(",","",$1); sum+=$1} END{print sum}')
cref=$(echo "$OUT" | awk '/^\s*[0-9,.]+\s+cache-references/ {gsub(",","",$1); sum+=$1} END{print sum}')
cmiss=$(echo "$OUT" | awk '/^\s*[0-9,.]+\s+cache-misses/ {gsub(",","",$1); sum+=$1} END{print sum}')

# Compute metrics
if [[ -n "$cycles" && -n "$instr" && "$cycles" -gt 0 ]]; then
    ipc=$(awk -v i="$instr" -v c="$cycles" 'BEGIN{printf "%.2f", i/c}')
else
    ipc="N/A"
fi

if [[ -n "$cref" && "$cref" -gt 0 ]]; then
    miss_rate=$(awk -v m="$cmiss" -v r="$cref" 'BEGIN{printf "%.2f", (m/r)*100}')
else
    miss_rate="N/A"
fi

echo
echo "=== Derived Metrics ==="
echo "IPC (Instructions per Cycle): $ipc"
echo "Cache Miss Rate: $miss_rate %"

# Simple verdict logic
if [[ "$ipc" != "N/A" && "$miss_rate" != "N/A" ]]; then
    if (( $(echo "$ipc < 1.0" | bc -l) )) && (( $(echo "$miss_rate > 30.0" | bc -l) )); then
        echo "Verdict: Likely MEMORY-BOUND (low IPC + high miss rate)"
    elif (( $(echo "$ipc > 1.5" | bc -l) )) && (( $(echo "$miss_rate < 20.0" | bc -l) )); then
        echo "Verdict: Likely COMPUTE-BOUND (high IPC + low miss rate)"
    else
        echo "Verdict: Mixed or workload-dependent"
    fi
else
    echo "Verdict: Unable to determine (missing data)"
fi

