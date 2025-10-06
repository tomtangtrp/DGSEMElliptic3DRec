#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;

int main(int argc, char* argv[])
{
    // ------------------ Vector----------------------------------------------
    // (1) special vector
    // N requires to be constexpr
    constexpr int N = 5;
    Vector<double, N> v0 = Vector<double, N>::Zero();
    std::cout << v0 << std::endl;
    Vector<double, N> v1 = Vector<double, N>::Ones();
    std::cout << v1 << std::endl;

    // Known size
    Vector3d v3(3);
    v3 << 1.5, 2.5, 3.5;
    std::cout << "fixed size use VectorNd and VectorXi = \n" << v3 << "\n";
    VectorXi v6(6);
    v6 << 0.5, 1.5, 2.5, 3.5, 4.5, 5.5;
    std::cout << "fixed size use 'VectorNd' and VectorXi = \n" << v6 << "\n";
    std::cout << "use v6.transpose to print horizontally (only for visual)= " << v6.transpose() << "\n";


    //Dynamic size use 
    VectorXd v;
    
    // outer product
    Eigen::VectorXd e0(3);
    Eigen::VectorXd dphi_0(3);
    e0 << 1.0, 0.0, 0.0;
    dphi_0 << -1.5, 2.0, -0.5;
    Eigen::MatrixXd outer = e0*dphi_0.transpose();
    std::cout << "outer(e0, dphi_0) = \n" << outer << "\n";



    std::cout << std::endl;

    return 0;
}