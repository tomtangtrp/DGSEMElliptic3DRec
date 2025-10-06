#pragma once
#include <Eigen/Sparse>
#include <stdexcept>
#include <vector>

#if defined(__APPLE__)
  #include <Accelerate/Accelerate.h>
#endif

namespace accelwrap {

// Solve A x = b via Apple Accelerate Sparse LU
// A: Eigen::SparseMatrix<double, ColMajor> (square, CSC, 0-based)
// b: VectorXd (size n)
// returns x
inline Eigen::VectorXd solve_with_accelerate_lu(
    const Eigen::SparseMatrix<double, Eigen::ColMajor>& A,
    const Eigen::VectorXd& b)
{
#if !defined(__APPLE__)
    throw std::runtime_error("Accelerate LU unavailable: not on Apple platform.");
#else
    using SpMat = Eigen::SparseMatrix<double, Eigen::ColMajor>;
    if (A.rows() != A.cols())  throw std::runtime_error("A must be square.");
    if (b.size() != A.rows())  throw std::runtime_error("b.size() must equal A.rows().");

    // Ensure CSC invariants (compressed)
    if (!A.isCompressed()) const_cast<SpMat&>(A).makeCompressed();

    const int n   = static_cast<int>(A.rows());
    const int nnz = static_cast<int>(A.nonZeros());

    // Eigen CSC pointers
    const auto* outer = A.outerIndexPtr();   // length n+1 (col starts)
    const auto* inner = A.innerIndexPtr();   // length nnz  (row indices)
    const double* vals = A.valuePtr();       // length nnz

    // Apple types: columnStarts = long*, rowIndices = int*
    std::vector<long> colStarts(n + 1);
    for (int j = 0; j <= n; ++j) colStarts[j] = static_cast<long>(outer[j]);

    // If Eigen::StorageIndex is already int, we can alias; else, copy to int32
    const int* rowIdx_ptr = nullptr;
    std::vector<int> rowIdx_copy;
    if constexpr (std::is_same_v<typename SpMat::StorageIndex,int>) {
        rowIdx_ptr = inner;
    } else {
        rowIdx_copy.resize(nnz);
        for (int k = 0; k < nnz; ++k) rowIdx_copy[k] = static_cast<int>(inner[k]);
        rowIdx_ptr = rowIdx_copy.data();
    }

    // Build SparseMatrix_Double (CSC)
    SparseAttributes_t matAttrs{};            // Ordinary general matrix
    matAttrs.kind = SparseOrdinary;
    matAttrs.transpose = false;

    SparseMatrixStructure S{};
    S.rowCount     = n;
    S.columnCount  = n;
    S.columnStarts = colStarts.data();        // long*
    S.rowIndices   = const_cast<int*>(rowIdx_ptr); // int*
    S.attributes   = matAttrs;
    S.blockSize    = 1;

    SparseMatrix_Double A_accel{};
    A_accel.structure = S;
    A_accel.data      = const_cast<double*>(vals);

    // Factor (LU)
    SparseOpaqueFactorization_Double F = SparseFactor(SparseFactorizationLU, A_accel);
    if (F.status != SparseStatusOK) {
        SparseCleanup(F); // safe even if partially constructed
        throw std::runtime_error("SparseFactor(LU) failed, status = " + std::to_string(F.status));
    }

    // Wrap RHS/solution as DenseMatrix_Double (1 column)
    Eigen::VectorXd x = Eigen::VectorXd::Zero(n);
    Eigen::VectorXd bcopy = b; // keep b const; could also alias if you don't need b after

    SparseAttributes_t dmAttrs{}; // defaults
    DenseMatrix_Double B{};
    B.attributes   = dmAttrs;
    B.rowCount     = n;
    B.columnCount  = 1;
    B.columnStride = n;          // column-major, contiguous
    B.data         = const_cast<double*>(bcopy.data());

    DenseMatrix_Double X{};
    X.attributes   = dmAttrs;
    X.rowCount     = n;
    X.columnCount  = 1;
    X.columnStride = n;
    X.data         = x.data();

    // Direct solve (matrix RHS form, no workspace)
    SparseSolve(F, B, X);  // This is the 3-arg direct-solve overload Eigen also uses.  [oai_citation:5‡Eigen](https://www.eigen.tuxfamily.org/dox/AccelerateSupport_8h_source.html)

    SparseCleanup(F);
    return x;
#endif
}

} // namespace accelwrap