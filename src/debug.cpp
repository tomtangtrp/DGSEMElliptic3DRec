#include <Eigen/Dense>
#include <iostream>

// Self-defined Kronecker product
template<typename DerivedA, typename DerivedB>
Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, Eigen::Dynamic>
kron(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedB>& B)
{
    using Scalar = typename DerivedA::Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> K(
        A.rows() * B.rows(), A.cols() * B.cols());

    for (Eigen::Index i = 0; i < A.rows(); ++i)
        for (Eigen::Index j = 0; j < A.cols(); ++j)
            K.block(i * B.rows(), j * B.cols(), B.rows(), B.cols()) = A(i, j) * B;

    return K;
}

int main() {
    using namespace Eigen;

    RowVector3d dphi_n;
    dphi_n << 0.5, -2, 1.5;

    // W_1 = diag([1/3, 4/3, 1/3])
    Matrix3d W_1 = Matrix3d::Zero();
    W_1(0, 0) = 1.0 / 3.0;
    W_1(1, 1) = 4.0 / 3.0;
    W_1(2, 2) = 1.0 / 3.0;

    // gR is a 27-vector in row-major ordering (z-fastest)
    VectorXd gR(9);
    gR << 2.71828183, 12.18249396, 54.59815003,
    7.3890561 , 33.11545196, 148.4131591,
    20.08553692, 90.0171313 , 403.42879349;

    // Step 1: fastest axis pair first: axis₁ (dphi_n) ⊗ axis₂ (W1)
    MatrixXd kron12 = kron(dphi_n, W_1);   // (3×9)

    // Step 2: prepend slowest axis₀ (W1)
    MatrixXd kron_final = kron(W_1, kron12); // (9×27)

    // Step 3: multiply
    RowVectorXd result = gR.transpose() * kron_final;
    VectorXd result_vec = result.transpose();

    std::cout << std::fixed << result_vec << "\n";

    return 0;
}
