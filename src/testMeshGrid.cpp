#include <iostream>
#include <Eigen/Dense>
#include <type_traits>
#include "MeshGrid3D.hpp"
#include "MeshGrid2D.hpp"
#include "gllQuadrature.hpp"

using namespace std;
using namespace Eigen;

Eigen::VectorXd exact_sol(
    const Eigen::Tensor<double,3,Eigen::RowMajor>& X,
    const Eigen::Tensor<double,3,Eigen::RowMajor>& Y,
    const Eigen::Tensor<double,3,Eigen::RowMajor>& Z);

Eigen::VectorXd g_back(
        const Eigen::ArrayXXd& fX,
        const Eigen::ArrayXXd& fY,
        double zLa);

VectorXd MapPhysical(VectorXd& xi_1d, double& a, double& b);

int main(int argc, char* argv[])
{   
    int k = 2;
    GLL mGLL(k);
    VectorXd xi_1d = mGLL.getGLLNodes();
    VectorXd eta_1d = xi_1d;
    VectorXd zeta_1d = xi_1d;
    double ax = 0.0;
    double bx = 0.5;
    double ay = 0.0;
    double by = 0.5;
    double az = 0.0;
    double bz = 0.5;
    VectorXd x_1d = MapPhysical(xi_1d, ax, bx);
    VectorXd y_1d = MapPhysical(eta_1d, ay, by);
    VectorXd z_1d = MapPhysical(eta_1d, az, bz);


    MeshGrid3D<double> mesh(x_1d, y_1d, z_1d);
    // grids are built when you first request them
    // auto& X = mesh.X();
    // auto& Y = mesh.Y();
    const Eigen::Tensor<double, 3, Eigen::RowMajor>& X = mesh.X();
    const Eigen::Tensor<double, 3, Eigen::RowMajor>& Y = mesh.Y();
    const Eigen::Tensor<double, 3, Eigen::RowMajor>& Z = mesh.Z();

    
    std::cout << "X = \n" << X << std::endl;
    std::cout << "Y = \n" << Y << std::endl;
    std::cout << "Z = \n" << Y << std::endl;

    
    VectorXd exact_ary = exact_sol(X, Y, Z);

    MeshGrid2D<double> facemesh(x_1d, y_1d);
    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fX = facemesh.X();
    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fY = facemesh.Y();
    VectorXd gBa = g_back(fX, fY, az);

    std::cout << "exact_sol = \n" << exact_ary << std::endl;
    std::cout << "gBa = \n" << gBa << std::endl;
    
    return 0;
}

Eigen::VectorXd exact_sol(
    const Eigen::Tensor<double,3,Eigen::RowMajor>& X,
    const Eigen::Tensor<double,3,Eigen::RowMajor>& Y,
    const Eigen::Tensor<double,3,Eigen::RowMajor>& Z)
{
    const Eigen::Index Ny = X.dimension(0);
    const Eigen::Index Nx = X.dimension(1);
    const Eigen::Index Nz = X.dimension(2);
    const Eigen::Index N  = Ny * Nx * Nz;

    // Vectorized element-wise compute using flat Maps (respect RowMajor layout).
    Eigen::Map<const Eigen::ArrayXd> x(X.data(), N);
    Eigen::Map<const Eigen::ArrayXd> y(Y.data(), N);
    Eigen::Map<const Eigen::ArrayXd> z(Z.data(), N);

    Eigen::VectorXd out(N);
    out = (x + 2.0*y + 3.0*z).exp().matrix();  // element-wise, then copy to VectorXd
    return out;
}

Eigen::VectorXd g_back(
    const Eigen::ArrayXXd& fX,
    const Eigen::ArrayXXd& fY,
    double zLa)
{
    // Element-wise face evaluation
    const Eigen::ArrayXXd F = (fX + 2.0 * fY + 3.0 * zLa).exp(); // shape: Ny x Nx

    // Force row-major flatten (row-by-row) without loops:
    // transpose -> column-major flatten -> equals row-major of original
    const Eigen::ArrayXXd Ft = F.transpose().eval();
    Eigen::VectorXd out = Eigen::Map<const Eigen::VectorXd>(Ft.data(), Ft.size());
    return out;
}



VectorXd MapPhysical(VectorXd& xi_1d, double& a, double& b)
{   
    int size = xi_1d.size();
    VectorXd x(size);
    for (int i=0; i<xi_1d.size(); i++)
    {
        x[i] = (b-a)*0.5*xi_1d[i] + (b+a)*0.5;
    } 
    return x;
}
