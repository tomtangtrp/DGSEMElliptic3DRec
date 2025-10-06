# DGSEMElliptic2DRec
practice project on DGSEM 3D rectangular grid for Elliptic PDEs, using Eigen's intel Pardiso sparse direct solver with OPENMP

# machine.mk
For different hardware platform (example included are ubuntu(x86)+gcc, gadi+gcc/intel, macos(arm64 apple silicon)+clang, just select different `BUILD_ENV` in machine.mk


# Dependencies
### install highfive
```bash
cd /usr/local/include
git clone --recursive https://github.com/BlueBrain/HighFive.git
cd HighFive
cmake -DCMAKE_INSTALL_PREFIX=build/install -DHIGHFIVE_USE_BOOST=Off -B build .
cmake --build build --parallel
cmake --install build
```

### ubuntu
```bash
sudo apt install libeigen3-dev
sudo apt install libhdf5-dev
```

### gadi
```bash
module load eigen/3.3.7
module load hdf5/1.10.7
module load intel-mkl/2025.2.0
odule load intel-compiler-llvm/2025.2.0
```

# Build
```bash
make testSolver_sp
```

# Run
### simple run
```bash
./testSolver_sp n SIPG 3 10 8 8 8 0 1 0 1 0 1 "L=D,R=D,B=D,T=D,back=N,front=N" cases/testCase0.json
```

### Convergence test
```bash
./runConvergence.sh SIPG testcase0
```

### Sigma penalty scan
```bash
python3 ./scripts/scanSigma.py \
  --exe ./testSolver_sp \
  --method SIPG -k 3 \
  --nelx 2 --nely 2 --nelz 2 \
  --xla 0 --xlb 1 --yla 0 --ylb 1 --zla 0 --zlb 1 \
  --bc "L=D,R=D,B=D,T=D,back=D,front=D" \
  --case ./cases/testcase0.json \
  --samples 50 --sigma_max_factor 3.0 \
  --write nw --outdir ./plots --logy
```

### omp scaling test
```bash
./runOmpScaling.sh
```

# Demonstration
### Solution: testcase0
<!-- plain mark down fig size is too large: ![solution testcase0 DC](./plots/solution_exponential_DC_SIPG.png) -->
<img src="./plots/solution_testcase0_DC_SIPG.png" alt="solution testcase0 DC" width="600"/>

### Convergence: testcase0
<img src="./plots/conv_testcase0_brokenL2_SIPG.png" alt="conv testcase0 DC" width="400"/>
<img src="./plots/conv_testcase0_brokenH1_SIPG.png" alt="conv testcase0 DC" width="400"/>

### DG(SIPG) penalty parameter scan: testcase0
<img src="./plots/scanSigma_testcase0_k3_SIPG.png" alt="sigma testcase0 k3" width="400"/>

### Pardiso OMP scaling on Gadi single node intel 8268 processor:testcase0
<img src="./plots/testcase0_SIPG_k2_Nel4096_walltime_vs_cores.png" alt="omp scaling" width="400"/>





