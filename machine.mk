# (1)g++ compiler linked to generic BLAS and LAPACK ( do not add comment after variable definition)
# following sytanx is better than BUILD_ENV?=gcc_ubuntu_blas, also !!! make sure there is no trailing space !!!
#BUILD_ENV ?= gcc_ubuntu_blas
# (2) g++ compiler linked to mkl BLAS and LAPACK, also use <Eigen/PardisoSupport> for intel's Pardiso sparse library
BUILD_ENV ?= gcc_ubuntu_intel
# (3)
#BUILD_ENV ?= gcc_gadiy08_intel
#BUILD_ENV ?= intel_gadiy08_intel
# (4) g++ compiler linked to Apple Accelerate
#BUILD_ENV ?= gcc_macos_accelerate

OPT_FLAGS = -Ofast
MATH_FLAGS = -ffast-math -fno-math-errno -fassociative-math -fno-signed-zeros
MP_FLAGS = -fopenmp

ifeq ($(BUILD_ENV),gcc_ubuntu_blas)
CXX = g++
CPPFLAGS = -I/usr/include/eigen3 -I/usr/include/mkl -I/usr/include/hdf5/serial -I/usr/local/include/HighFive/include -I./include
#ARCH_FLAGS = -march=alderlake -mtune=alderlake
ARCH_FLAGS = -march=alderlake -mavx2 -mtune=alderlake
# LAPACKE is c API for LAPACK
EIGEN_FLAGS = -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
CXXFLAGS = -std=c++17 -funroll-loops -mfma -DNDEBUG $(MP_FLAGS) $(OPT_FLAGS) $(MATH_FLAGS) $(ARCH_FLAGS) $(EIGEN_FLAGS)
LDFLAGS = -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LDLIBS = -lblas -llapacke -lhdf5 -lhdf5_hl -lpthread -lm
endif

ifeq ($(BUILD_ENV),gcc_ubuntu_intel)
CXX = g++
CPPFLAGS = -I/usr/include/eigen3 -I/usr/include/mkl -I/usr/include/hdf5/serial -I/usr/local/include/HighFive/include -I./include
#ARCH_FLAGS = -march=alderlake -mtune=alderlake
ARCH_FLAGS = -march=alderlake -mavx2 -mtune=alderlake
MKL_FLAGS = -m64
# LAPACKE is c API for LAPACK
# EIGEN_FLAGS = -DEIGEN_USE_MKL_ALL
# DEBUG_FLAGS = -g -Wall
DEBUG_FLAGS = -DNDEBUG -w
CXXFLAGS = -std=c++17 -funroll-loops -mfma $(DEBUG_FLAGS) $(MKL_FLAGS) $(MP_FLAGS) $(OPT_FLAGS) $(MATH_FLAGS) $(ARCH_FLAGS) $(EIGEN_FLAGS)
LDFLAGS = -L/usr/lib/x86_64-linux-gnu -L/usr/lib/x86_64-linux-gnu/hdf5/serial
LDLIBS = -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lhdf5 -lhdf5_hl -lpthread -lm
endif


ifeq ($(BUILD_ENV),gcc_gadiy08_intel)
CXX = g++
CPPFLAGS = -I${EIGEN_BASE}/eigen3 -I${MKLROOT}/include -I$(HDF5_BASE)/include -I../minclude/HighFive/include -I./include
#ARCH_FLAGS = -march=skylake -mtune=skylake
ARCH_FLAGS = -march=skylake-avx512 -mtune=skylake
MKL_FLAGS = -m64
# LAPACKE is c API for LAPACK
EIGEN_FLAGS = -DEIGEN_USE_MKL_ALL
CXXFLAGS = -std=c++17 -funroll-loops -mfma -DNDEBUG $(MKL_FLAGS) $(MP_FLAGS) $(OPT_FLAGS) $(MATH_FLAGS) $(ARCH_FLAGS) $(EIGEN_FLAGS)
LDFLAGS = -L${MKLROOT}/lib/intel64 -L${HDF5_BASE}/lib
LDLIBS = -lmkl_intel_lp64 -lmkl_gnu_thread -lmkl_core -lhdf5 -lhdf5_hl -lpthread -lm
endif

ifeq ($(BUILD_ENV),intel_gadiy08_intel)
CXX = icpx
CPPFLAGS = -I${EIGEN_BASE}/eigen3 -I${MKLROOT}/include -I$(HDF5_BASE)/include -I../minclude/HighFive/include -I./include
OPT_FLAGS = -O3
#ARCH_FLAGS = -march=skylake-avx512 -mtune=skylake
ARCH_FLAGS = -march=cascadelake -mavx512f -mtune=cascadelake
MKL_FLAGS = -m64
MP_FLAGS = -qopenmp
# LAPACKE is c API for LAPACK
EIGEN_FLAGS = -DEIGEN_USE_MKL_ALL
CXXFLAGS = -std=c++17 -funroll-loops -mfma -DNDEBUG $(MKL_FLAGS) $(MP_FLAGS) $(OPT_FLAGS) $(MATH_FLAGS) $(ARCH_FLAGS) $(EIGEN_FLAGS)
LDFLAGS = -L${MKLROOT}/lib/intel64 -L${HDF5_BASE}/lib
LDLIBS = -lmkl_intel_lp64 -lmkl_intel_thread -lmkl_core -liomp5 -lhdf5 -lhdf5_hl -lpthread -lm -lstdc++fs
endif

ifeq ($(BUILD_ENV),gcc_macos_accelerate)
#CXX = /usr/local/bin/g++
#CXX = /Library/Developer/CommandLineTools/usr/bin/g++
SDK := /Library/Developer/CommandLineTools/SDKs/MacOSX.sdk
CXX = /Library/Developer/CommandLineTools/usr/bin/clang++ -isysroot $(SDK)
CPPFLAGS = -I/usr/local/include/eigen3 -I/opt/homebrew/opt/hdf5/include -I/usr/local/include/HighFive/include -I./in
clude
#ARCH_FLAGS = -mcpu=apple-m3 -mtune=apple-m3
ARCH_FLAGS = -mcpu=apple-m4 -mtune=apple-m4
# LAPACKE is c API for LAPACK
#EIGEN_FLAGS = -DEIGEN_USE_BLAS -DEIGEN_USE_LAPACKE
EIGEN_FLAGS = -DEIGEN_NO_DEBUG -DEIGEN_VECTORIZE=1 -DEIGEN_VECTORIZE_NEON=1 -DEIGEN_UNROLLING_LIMIT=256
# To use <Eigen/AccelerateSupprt>, disabled EIGEN_FLAGS
#DEBUG_FLAGS = -w -DNDEBUG -DEIGEN_NO_DEBUG
# DEBUG_FLAGS = -g -Wall
CXXFLAGS = -std=c++17 -funroll-loops -flto=thin -ffp-contract=fast -fstrict-aliasing $(DEBUG_FLAGS) $(OPT_FLAGS) $(E
IGEN_FLAGS) $(MATH_FLAGS) $(ARCH_FLAGS)
LDFLAGS += -flto=thin -framework Accelerate -L/opt/homebrew/opt/hdf5/lib
LDLIBS = -lhdf5 -lhdf5_hl -lpthread -lm
endif


