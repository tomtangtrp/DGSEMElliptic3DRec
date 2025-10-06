#include <iostream>
#include <Eigen/Dense>
#include "NumUtilities.hpp"

using namespace Eigen;

VectorXd NumUtilities::MapPhysical(VectorXd& xi_1d, double& a, double& b)
{
    int size = xi_1d.size();
    VectorXd x(size);
    for (int i=0; i<xi_1d.size(); i++)
    {
        x[i] = (b-a)*0.5*xi_1d[i] + (b+a)*0.5;
    } 
    return x;
}

// ---- symmetry helpers ----

double NumUtilities::SkewRelFrobenius(const SparseMatrix<double>& A)
{
    if (A.rows() != A.cols()) return std::numeric_limits<double>::infinity();

    // Work with compressed copies for stable ops
    SparseMatrix<double> Ac = A;
    Ac.makeCompressed();

    // Skew part: A - A^T
    SparseMatrix<double> Skew = Ac - SparseMatrix<double>(Ac.transpose());

    const double nA   = Ac.norm();
    const double nSk  = Skew.norm();
    const double denom = (nA > 0.0 ? nA : 1.0);

    return nSk / denom;
}

bool NumUtilities::IsSymmetric(const SparseMatrix<double>& A, double rel_tol)
{
    if (A.rows() != A.cols()) return false;
    const double rel = SkewRelFrobenius(A);
    return std::isfinite(rel) && (rel <= rel_tol);
}