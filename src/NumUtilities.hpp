#ifndef NUMUTILITIESHEADERDEF
#define NUMUTILITIESHEADERDEF

#include <iostream>
#include <Eigen/Dense>
#include <Eigen/Sparse>   // (new)

using namespace Eigen;

class NumUtilities
{
  public:
    VectorXd MapPhysical(VectorXd& xi_1d, double& a, double& b);

    // --- Simple sparse symmetry checks ---
    // Relative test: returns true if ||A - A^T||_F <= rel_tol * max(1, ||A||_F)
    static bool IsSymmetric(const Eigen::SparseMatrix<double>& A, double rel_tol = 1e-12);

    // Returns the relative skew: ||A - A^T||_F / max(1, ||A||_F).
    // (Use this to print a quick diagnostic value.)
    static double SkewRelFrobenius(const Eigen::SparseMatrix<double>& A);
};

#endif
