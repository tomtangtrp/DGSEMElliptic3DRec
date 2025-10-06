python3 ./scripts/scanSigma.py \
  --exe ./testSolver_sp \
  --method SIPG -k 3 \
  --nelx 2 --nely 2 --nelz 2 \
  --xla 0 --xlb 1 --yla 0 --ylb 1 --zla 0 --zlb 1 \
  --bc "L=D,R=D,B=D,T=D,back=D,front=D" \
  --case ./cases/testcase0.json \
  --samples 50 --sigma_max_factor 3.0 \
  --write nw --outdir ./plots --logy

