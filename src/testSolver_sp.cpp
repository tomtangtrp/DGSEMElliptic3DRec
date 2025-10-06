// C++ stdlib
#include <iostream>
#include <string>
#include <string_view>
#include <vector>
#include <chrono>
#include <cmath>
#include <sstream>
#include <cctype>

// Eigen
#include <Eigen/PardisoSupport>
#include <Eigen/Dense>
#include <Eigen/Sparse>
#include <Eigen/SparseLU>
#include <Eigen/SparseCholesky>

// HighFive
#include <highfive/H5File.hpp>
#include <highfive/eigen.hpp>

// Project headers
#include "gllQuadrature.hpp"
#include "MeshGrid3D.hpp"
#include "MeshGrid2D.hpp"
#include "RecGrid.hpp"
#include "ConnectionTable.hpp"
#include "RecElement.hpp"
#include "JsonCase.hpp"
#include "DGError.hpp"
#include "NumUtilities.hpp"

using namespace Eigen;


// ----- small helpers (unchanged patterns) -----
template <typename SparseMat> std::size_t sparseMemoryBytes(const SparseMat& M) {
    using Scalar = typename SparseMat::Scalar;
    using StorageIndex = typename SparseMat::StorageIndex;

    std::size_t bytes = 0;
    bytes += sizeof(Scalar)      * M.nonZeros();       // values
    bytes += sizeof(StorageIndex)* M.nonZeros();       // inner indices
    bytes += sizeof(StorageIndex)* (M.outerSize() + 1);// outer ptrs
    return bytes;
}

template<typename Derived>
Eigen::VectorXd flatten_transpose_rowmajor(const Eigen::MatrixBase<Derived>& m) {
#if EIGEN_VERSION_AT_LEAST(3,4,0)
    return m.transpose().reshaped();
#else
    Eigen::Matrix<typename Derived::Scalar, Eigen::Dynamic, Eigen::Dynamic> tmp = m.transpose();
    return Eigen::Map<const Eigen::VectorXd>(tmp.data(), tmp.size());
#endif
}

template<typename DerivedA, typename DerivedB>
Eigen::Matrix<typename DerivedA::Scalar, Eigen::Dynamic, Eigen::Dynamic>
kron(const Eigen::MatrixBase<DerivedA>& A, const Eigen::MatrixBase<DerivedB>& B) {
    using Scalar = typename DerivedA::Scalar;
    Eigen::Matrix<Scalar, Eigen::Dynamic, Eigen::Dynamic> K(
        A.rows()*B.rows(), A.cols()*B.cols());
    for (Eigen::Index i = 0; i < A.rows(); ++i)
        for (Eigen::Index j = 0; j < A.cols(); ++j)
            K.block(i*B.rows(), j*B.cols(), B.rows(), B.cols()) = A(i,j)*B;
    return K;
}

// ----- DG method enum & parser (same names as 2D) -----
enum class Method { SIPG, IIPG, NIPG, NIPG0 };
static bool parseMethod(std::string_view s, Method& out) {
    if (s=="SIPG")  { out=Method::SIPG;  return true; }
    if (s=="IIPG")  { out=Method::IIPG;  return true; }
    if (s=="NIPG")  { out=Method::NIPG;  return true; }
    if (s=="NIPG0") { out=Method::NIPG0; return true; }
    return false;
}

// ----- 3D boundary spec (compact-string) -----
enum class BC { DC, NC }; // Dirichlet, Neumann
struct BCspec3D {
    BC left=BC::DC, right=BC::DC, bottom=BC::DC, top=BC::DC, back=BC::DC, front=BC::DC;
};
static bool parseBCToken3D(std::string tok, BC& out) {
    for (auto& ch: tok) ch = std::tolower(ch);
    if (tok=="d" || tok=="dc" || tok=="dirichlet") { out=BC::DC; return true; }
    if (tok=="n" || tok=="nc" || tok=="neumann")   { out=BC::NC; return true; }
    return false;
}
static BCspec3D parseBCspec3D(const std::string& s) {
    BCspec3D spec;
    auto trim = [](std::string& z){ z.erase(0, z.find_first_not_of(" \t")); z.erase(z.find_last_not_of(" \t")+1); };
    auto apply = [&](std::string kv){
        auto p = kv.find('=');
        if (p==std::string::npos) return;
        std::string k = kv.substr(0,p), v = kv.substr(p+1);
        trim(k); trim(v);
        for (auto& ch: k) ch = std::tolower(ch);
        BC tmp;
        if (!parseBCToken3D(v, tmp)) return;
        if (k=="l" || k=="left")        spec.left  = tmp;
        else if (k=="r" || k=="right")  spec.right = tmp;
        else if (k=="b" || k=="bottom") spec.bottom= tmp;
        else if (k=="t" || k=="top")    spec.top   = tmp;
        else if (k=="back" || k=="bk")  spec.back  = tmp;
        else if (k=="front"|| k=="f")   spec.front = tmp;
    };
    std::string acc;
    for (char c: s) { if (c==',' || c==';') { apply(acc); acc.clear(); } else acc.push_back(c); }
    if (!acc.empty()) apply(acc);
    return spec;
}

int main(int argc, char* argv[])
{
    // Unified 3D CLI (same pattern as 2D, with z added)
    // {w|nw} METHOD k sigma0 Nel_x Nel_y Nel_z xLa xLb yLa yLb zLa zLb BC_SPEC CASE_JSON
    if (argc != 16) {
        std::cerr << "Usage: " << argv[0]
                  << " {w|nw} METHOD k sigma0 Nel_x Nel_y Nel_z xLa xLb yLa yLb zLa zLb BC_SPEC CASE_JSON\n";
        return 1;
    }

    const std::string write_h5 = argv[1];    // 'w' or 'nw'
    Method method;
    if (!parseMethod(argv[2], method)) {
        std::cerr << "Unknown METHOD: " << argv[2] << " (expected SIPG|IIPG|NIPG|NIPG0)\n";
        return 1;
    }
    const int    k      = std::stoi(argv[3]);
    double       sigma0 = std::stod(argv[4]);           // NIPG/NIPG0 override later
    const int    Nel_x  = std::stoi(argv[5]);
    const int    Nel_y  = std::stoi(argv[6]);
    const int    Nel_z  = std::stoi(argv[7]);
    const double xLa    = std::stod(argv[8]);
    const double xLb    = std::stod(argv[9]);
    const double yLa    = std::stod(argv[10]);
    const double yLb    = std::stod(argv[11]);
    const double zLa    = std::stod(argv[12]);
    const double zLb    = std::stod(argv[13]);
    const BCspec3D mBC   = parseBCspec3D(argv[14]);
    const std::string case_path = argv[15];

    if (!(xLb > xLa && yLb > yLa && zLb > zLa)) {
        std::cerr << "ERROR: require xLb>xLa, yLb>yLa, and zLb>zLa.\n";
        return 1;
    }

    // JSON case (alpha, f, u, faces)
    DG::JsonCase mCase(case_path);
    const double alpha = mCase.alpha();

    // One-time method parameterization (same semantics as your 2D code)
    int eps = 0; double sigma0b = 0.0; double beta0 = 1.0;
    switch (method) {
        case Method::SIPG:
            eps = -1; sigma0b = 2.0 * sigma0; beta0 = 1.0; break;
        case Method::IIPG:
            eps =  0; sigma0b = 2.0 * sigma0; beta0 = 1.0; break;
        case Method::NIPG:
            eps =  1; sigma0  = 1.0; sigma0b = 1.0; beta0 = 1.0; break;
        case Method::NIPG0:
            eps =  1; sigma0  = 0.0; sigma0b = 0.0; beta0 = 0.0; break;
    }

    auto start = std::chrono::high_resolution_clock::now();
    int locdim_1d = k+1;
    int locdim = locdim_1d*locdim_1d*locdim_1d;
    RecGrid mRecGrid(Nel_x, Nel_y, Nel_z);
    ConnectionTable mRecConnectionTable = mRecGrid.getRecConnectionTable();
    // mRecConnectionTable.printConnectionTable();

    int Nel = Nel_x*Nel_y*Nel_z;
    // std::cout << "Nel = " << Nel << "\n";
    int dim = Nel*locdim;

    double h_n_x = (xLb-xLa)/double(Nel_x);
    double h_n_y = (yLb-yLa)/double(Nel_y);
    double h_n_z = (zLb-zLa)/double(Nel_z);
    double J1d_x = h_n_x/2.0;
    double J1d_y = h_n_y/2.0;
    double J1d_z = h_n_z/2.0;
    double Jyx = J1d_y*J1d_x;
    double Jyz = J1d_y*J1d_z;
    double Jxz = J1d_x*J1d_z;
    double dfac_x = 2.0/h_n_x;
    double dfac_y = 2.0/h_n_y;
    double dfac_z = 2.0/h_n_z;

    GLL mGLL(k);
    // 1D operators
    VectorXd xi_1d = mGLL.getGLLNodes();
    VectorXd eta_1d = xi_1d;
    VectorXd zeta_1d = xi_1d;
    VectorXd w_1d = mGLL.getGLLWeights(); 
    // be careful about off-diagonal blow up
    MatrixXd W_1 = mGLL.get1dMass(); // Diagonal, important to initialize W1.setZero() or w_1d.asDiagonal(), 
    VectorXd bw_1d = mGLL.getBaryWeights(); 
    // std::cout << "W_1=\n" << W_1 << "\n";
    MatrixXd D_1 = mGLL.get1dDiff();
    MatrixXd M_x = J1d_x*W_1;
    // std::cout << "M_x=\n" << M_x << "\n";
    MatrixXd M_y= J1d_y*W_1;
    MatrixXd M_z = J1d_z*W_1;
    MatrixXd D_x = dfac_x*D_1;
    MatrixXd D_y= dfac_y*D_1;
    MatrixXd D_z = dfac_z*D_1;
    MatrixXd I = Eigen::MatrixXd::Identity(locdim_1d, locdim_1d);
    MatrixXd S_x = D_x.transpose()*M_x*D_x;
    MatrixXd S_y = D_y.transpose()*M_y*D_y;
    MatrixXd S_z = D_z.transpose()*M_z*D_z;
    VectorXd e0 = VectorXd::Zero(locdim_1d);
    e0[0] = 1.0;
    VectorXd en = VectorXd::Zero(locdim_1d);
    en[k] = 1.0;
    VectorXd dphi_0 = D_1.row(0);
    VectorXd dphi_n = D_1.row(k);
    // std::cout << "S_x" << S_x << "\n";

    // projection operators
    MatrixXd proj_0 = e0 * e0.transpose();
    MatrixXd proj_n = en * en.transpose();
    MatrixXd proj_0n = e0*en.transpose();
    MatrixXd proj_n0 = en*e0.transpose();
    MatrixXd d0e0T = dphi_0*e0.transpose();
    MatrixXd e0d0T = e0*dphi_0.transpose();
    MatrixXd dnenT = dphi_n*en.transpose();
    MatrixXd endnT = en*dphi_n.transpose();
    MatrixXd d0enT = dphi_0*en.transpose();
    // std::cout << "np.outer(dphi_0, en)" << d0enT << "\n";
    MatrixXd end0T = en*dphi_0.transpose();
    MatrixXd e0dnT = e0*dphi_n.transpose();
    MatrixXd dne0T = dphi_n*e0.transpose();
    // MatrixXd Diag_d0 = dphi_0.asDiagonal();
    // MatrixXd Diag_dn = dphi_n.asDiagonal();
    // 2D operators:
    MatrixXd M_yz = kron(M_y, M_z);
    MatrixXd M_xz = kron(M_x, M_z);
    MatrixXd M_yx = kron(M_y, M_x);

    // 3D operator: Mass and differentiation matrix on reference element
    MatrixXd M = kron(kron(M_y, M_x), M_z);
    // std::cout << "M" << M << "\n";
    MatrixXd K_x =  kron(kron(M_y, S_x), M_z);
    MatrixXd K_y =  kron(kron(S_y, M_x), M_z);
    MatrixXd K_z =  kron(kron(M_y, M_x), S_z);
    MatrixXd K = K_x + K_y + K_z;
    // std::cout << "K" << K << "\n";

    // Manufactured exact solution
    VectorXd u_exact(dim);

    // === Sparse global matrix and assembly ===
    using Trip = Eigen::Triplet<double>;
    std::vector<Trip> triplets;
    triplets.reserve(std::size_t( (size_t)Nel * (size_t)locdim * (size_t)locdim * 6 / 4 )); // rough guess

    auto addDenseBlock = [&](int r0, int c0, const Eigen::Ref<const MatrixXd>& B){
        const int Br = (int)B.rows();
        const int Bc = (int)B.cols();
        for (int r = 0; r < Br; ++r) {
            for (int c = 0; c < Bc; ++c) {
                double v = B(r,c);
                if (v != 0.0) triplets.emplace_back(r0 + r, c0 + c, v);
            }
        }
    };

    VectorXd b = VectorXd::Zero(dim);

    MatrixXd C11_FB = - 0.5*Jyx*dfac_z*kron(kron(W_1, W_1), endnT)
                      + 0.5*double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1), dnenT)
                      + Jyx*(sigma0/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1), proj_n);
    MatrixXd C22_FB = + 0.5*Jyx*dfac_z*kron(kron(W_1, W_1), e0d0T)
                      - 0.5*double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1), d0e0T)
                      + Jyx*(sigma0/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1), proj_0);                   
    MatrixXd C12_FB = - 0.5*Jyx*dfac_z*kron(kron(W_1, W_1), end0T)
                      - 0.5*double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1), dne0T)
                      - Jyx*(sigma0/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1), proj_n0);
    MatrixXd C21_FB = + 0.5*Jyx*dfac_z*kron(kron(W_1, W_1), e0dnT)
                      + 0.5*double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1), d0enT)
                      - Jyx*(sigma0/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1), proj_0n);

    MatrixXd C11_RL = - 0.5*Jyz*dfac_x*kron(kron(W_1, endnT), W_1)
                      + 0.5*double(eps)*Jyz*dfac_x*kron(kron(W_1, dnenT), W_1)
                      + Jyz*(sigma0/(std::pow(h_n_x, beta0)))*kron(kron(W_1, proj_n), W_1);
    MatrixXd C22_RL = + 0.5*Jyz*dfac_x*kron(kron(W_1, e0d0T), W_1)
                      - 0.5*double(eps)*Jyz*dfac_x*kron(kron(W_1, d0e0T), W_1)
                      + Jyz*(sigma0/(std::pow(h_n_x, beta0)))*kron(kron(W_1, proj_0), W_1);                   
    MatrixXd C12_RL = - 0.5*Jyz*dfac_x*kron(kron(W_1, end0T), W_1)
                      - 0.5*double(eps)*Jyz*dfac_x*kron(kron(W_1, dne0T), W_1)
                      - Jyz*(sigma0/(std::pow(h_n_x, beta0)))*kron(kron(W_1, proj_n0), W_1);
    MatrixXd C21_RL = + 0.5*Jyz*dfac_x*kron(kron(W_1, e0dnT), W_1)
                      + 0.5*double(eps)*Jyz*dfac_x*kron(kron(W_1, d0enT), W_1)
                      - Jyz*(sigma0/(std::pow(h_n_x, beta0)))*kron(kron(W_1, proj_0n), W_1);

    MatrixXd C11_TB = - 0.5*Jxz*dfac_y*kron(endnT, kron(W_1, W_1))
                      + 0.5*double(eps)*Jxz*dfac_y*kron(dnenT, kron(W_1, W_1))
                      + Jxz*(sigma0/(std::pow(h_n_y, beta0)))*kron(proj_n, kron(W_1, W_1));
    MatrixXd C22_TB = + 0.5*Jxz*dfac_y*kron(e0d0T, kron(W_1, W_1))
                      - 0.5*double(eps)*Jxz*dfac_y*kron(d0e0T, kron(W_1, W_1))
                      + Jxz*(sigma0/(std::pow(h_n_y, beta0)))*kron(proj_0, kron(W_1, W_1));                   
    MatrixXd C12_TB = - 0.5*Jxz*dfac_y*kron(end0T, kron(W_1, W_1))
                      - 0.5*double(eps)*Jxz*dfac_y*kron(dne0T, kron(W_1, W_1))
                      - Jxz*(sigma0/(std::pow(h_n_y, beta0)))*kron(proj_n0, kron(W_1, W_1));
    MatrixXd C21_TB = + 0.5*Jxz*dfac_y*kron(e0dnT, kron(W_1, W_1))
                      + 0.5*double(eps)*Jxz*dfac_y*kron(d0enT, kron(W_1, W_1))
                      - Jxz*(sigma0/(std::pow(h_n_y, beta0)))*kron(proj_0n, kron(W_1, W_1));                      
    
    // int nsign1 = 1;
    // int nsign2 = -1;

    // Corresponds to C11_FR
    MatrixXd Bfront = - Jyx*dfac_z*kron(kron(W_1, W_1), endnT)
                      + double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1), dnenT)
                      + Jyx*(sigma0b/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1), proj_n);
    // Corresponds to C22_FR
    MatrixXd Bback  = + Jyx*dfac_z*kron(kron(W_1, W_1), e0d0T)
                      - double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1), d0e0T)
                      + Jyx*(sigma0b/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1), proj_0);       

    MatrixXd Bright = - Jyz*dfac_x*kron(kron(W_1, endnT), W_1)
                      + double(eps)*Jyz*dfac_x*kron(kron(W_1, dnenT), W_1)
                      + Jyz*(sigma0b/(std::pow(h_n_x, beta0)))*kron(kron(W_1, proj_n), W_1);

    MatrixXd Bleft = + Jyz*dfac_x*kron(kron(W_1, e0d0T), W_1)
                     - double(eps)*Jyz*dfac_x*kron(kron(W_1, d0e0T), W_1)
                     + Jyz*(sigma0b/(std::pow(h_n_x, beta0)))*kron(kron(W_1, proj_0), W_1);    

    MatrixXd Btop =  - Jxz*dfac_y*kron(endnT, kron(W_1, W_1))
                     + double(eps)*Jxz*dfac_y*kron(dnenT, kron(W_1, W_1))
                     + Jxz*(sigma0b/(std::pow(h_n_y, beta0)))*kron(proj_n, kron(W_1, W_1));

    MatrixXd Bbottom = + Jxz*dfac_y*kron(e0d0T, kron(W_1, W_1))
                       - double(eps)*Jxz*dfac_y*kron(d0e0T, kron(W_1, W_1))
                       + Jxz*(sigma0b/(std::pow(h_n_y, beta0)))*kron(proj_0, kron(W_1, W_1));           


    // std::cout << "alpha=" << alpha << ", k=" << k << ", Nel=" << Nel << " (" << Nel_x << "x" << Nel_y << "x" << Nel_z << ")\n";
    // std::cout << "dim = " << dim << "\n";
    // std::cout << "sigma0 = " << sigma0 << ", sigma0b = " << sigma0b << ", beta0 = " << beta0 << ", eps = " << eps << "\n";
    // std::cout << "Assembling the global matrix...\n";

    // Volume terms
    for (int i=0; i<Nel; i++)
    {   
        int m_z = i%Nel_z;
        int m_x = (i/Nel_z)%Nel_x;
        int m_y = i/(Nel_x*Nel_z);
        // std::cout << "m_x = " << m_x << ", m_y = "<< m_y << ", m_z = " << m_z << "\n";
        double a_x = (double (m_x))*h_n_x;
        double b_x = a_x + h_n_x;
        double a_y = (double (m_y))*h_n_y;
        double b_y = a_y + h_n_y;
        double a_z = (double (m_z))*h_n_z;
        double b_z = a_z + h_n_z;
        VectorXd x_1d = (b_x-a_x)/2.0*xi_1d.array() + (b_x+a_x)/2.0;
        VectorXd y_1d = (b_y-a_y)/2.0*eta_1d.array() + (b_y+a_y)/2.0;
        VectorXd z_1d = (b_z-a_z)/2.0*zeta_1d.array() + (b_z+a_z)/2.0;

        MeshGrid3D<double> mesh(x_1d, y_1d, z_1d);
        const Eigen::Tensor<double, 3, Eigen::RowMajor>& X = mesh.X();
        const Eigen::Tensor<double, 3, Eigen::RowMajor>& Y = mesh.Y();
        const Eigen::Tensor<double, 3, Eigen::RowMajor>& Z = mesh.Z();

        // Assemble A with volume terms.
        addDenseBlock(i*locdim, i*locdim, K + alpha*M);

        // Source contribution to b
        // MatrixXd bf_mat = M*(eval_3D(X, Y, Z, source1).matrix());
        // VectorXd bf = M*source(X, Y, Z);
        // VectorXd bf = M*mCase.f(X, Y, Z);
        // VectorXd bf_local = M*flatten_transpose_rowmajor(mCase.f(X, Y, Z).matrix());
        VectorXd bf_local = M*(mCase.f(X, Y, Z));
        b.segment(i*locdim, locdim) = bf_local;
        // std::cout << "bf= \n" << bf << "\n";

        // Compute u_exact
        // VectorXd u_exact_local= exact_sol(X, Y, Z);
        // VectorXd u_exact_local= mCase.u(X, Y, Z);
        // VectorXd u_exact_local = flatten_transpose_rowmajor(mCase.u(X, Y, Z).matrix());
        VectorXd u_exact_local =mCase.u(X, Y, Z);
        u_exact.segment(i*locdim, locdim) = u_exact_local;

        RecElement elm_local(a_x, b_x, a_y, b_y, a_z, b_z);
        for (int ie=0; ie<6; ie++)
        {
            Face face = elm_local.getFaces()[ie];
            face.check_bdry(xLa, xLb, yLa, yLb, zLa, zLb);
            if (bool bool_bdry = face.get_is_bdry()) {
                int face_lid = face.get_face_lid();
                if (face_lid == 0) {
                    // Back face
                    MeshGrid2D<double> facemesh(x_1d, y_1d);
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fX = facemesh.X();
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fY = facemesh.Y(); 
                    if (mBC.back == BC::DC && alpha>=0.0) {
                    // std::cout << "Boundary Back" << "\n";
                    addDenseBlock(i*locdim, i*locdim, Bback);
                    VectorXd gBa =  mCase.g_back(fX, fY, zLa);
                    // std::cout << "gBa = \n" << gBa << "\n"; 
                    MatrixXd bB_mat = -double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1)*gBa, dphi_0) + Jyx*(sigma0b/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1)*gBa, e0);
                    // VectorXd bB = bB_mat.transpose().reshaped();
                    VectorXd bB = flatten_transpose_rowmajor(bB_mat);
                    // std::cout << "bB = \n" << bB << "\n"; 
                    b.segment(i*locdim, locdim) += bB;
                    }
                    else if (mBC.back == BC::NC) {
                        VectorXd gNBa = mCase.gN_back(fX, fY, zLa);
                        // Natural term: +∫ g_N v on the back face  -> (Jyx * (W⊗W) * gN) ⊗ e0
                        MatrixXd bNB_mat = Jyx * kron(kron(W_1, W_1)*gNBa, e0);
                        VectorXd bNB = flatten_transpose_rowmajor(bNB_mat);
                        b.segment(i*locdim, locdim) += bNB;
                    }
                }
                if (face_lid == 1) { 
                    // Front edge
                    MeshGrid2D<double> facemesh(x_1d, y_1d);
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fX = facemesh.X();
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fY = facemesh.Y(); 
                    if (mBC.front == BC::DC && alpha>=0.0) {
                    // std::cout << "Boundary Front" << "\n";
                    addDenseBlock(i*locdim, i*locdim, Bfront);
                    VectorXd gF =  mCase.g_front(fX, fY, zLb);
                    // std::cout << "gF = \n" << gF << "\n"; 
                    MatrixXd bF_mat = double(eps)*Jyx*dfac_z*kron(kron(W_1, W_1)*gF, dphi_n) + Jyx*(sigma0b/(std::pow(h_n_z, beta0)))*kron(kron(W_1, W_1)*gF, en);
                    VectorXd bF = flatten_transpose_rowmajor(bF_mat);
                    // VectorXd bF = bF_mat.transpose().reshaped();
                    // std::cout << "bF = \n" << bF << "\n"; 
                    b.segment(i*locdim, locdim) += bF;
                    }
                    else if (mBC.front == BC::NC) {
                        VectorXd gNF = mCase.gN_front(fX, fY, zLb);
                        // Natural term: +∫ g_N v on the front face  -> (Jyx * (W⊗W) * gN) ⊗ en
                        MatrixXd bNF_mat = Jyx * kron(kron(W_1, W_1)*gNF, en);
                        VectorXd bNF = flatten_transpose_rowmajor(bNF_mat);
                        b.segment(i*locdim, locdim) += bNF;
                    }
                }
                if (face_lid == 2) { 
                    // Left edge
                    MeshGrid2D<double> facemesh(z_1d, y_1d);
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fZ = facemesh.X();
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fY = facemesh.Y(); 
                    if (mBC.left == BC::DC && alpha>=0.0) {
                        // std::cout << "Boundary Left" << "\n";
                        addDenseBlock(i*locdim, i*locdim, Bleft);
                        VectorXd gL =  mCase.g_left(fY, fZ, xLa);
                        // std::cout << "gL = \n" << gL << "\n"; 
                        // MatrixXd bL_mat = -double(eps)*Jyz*dfac_x*(gL.transpose()*(kron(W_1, kron(dphi_0, W_1)))) + Jyz*(sigma0b/(std::pow(h_n_x, beta0)))*(gL.transpose()*kron(W_1, kron(e0, W_1)));
                        // VectorXd bL = bL_mat.transpose().reshaped();
                        // VectorXd bL = flatten_transpose_rowmajor(bL_mat);
                        RowVectorXd dphi_0_row = dphi_0.transpose();
                        RowVectorXd e0_row = e0.transpose();
                        RowVectorXd bL_row = -double(eps)*Jyz*dfac_x*(gL.transpose()*(kron(W_1, kron(dphi_0_row, W_1)))) + Jyz*(sigma0b/(std::pow(h_n_x, beta0)))*(gL.transpose()*kron(W_1, kron(e0_row, W_1)));
                        VectorXd bL = bL_row.transpose();
                        // std::cout << "bL = \n" << bL << "\n"; 
                        b.segment(i*locdim, locdim) += bL;
                    }
                    else if (mBC.left == BC::NC) {
                        VectorXd gNL = mCase.gN_left(fY, fZ, xLa);
                        RowVectorXd e0_row = e0.transpose();
                        RowVectorXd bNL_row = Jyz * (gNL.transpose() * kron(W_1, kron(e0_row, W_1)));
                        b.segment(i*locdim, locdim) += bNL_row.transpose();
                    }

                }
                if (face_lid == 3) { 
                    // Right edge
                    MeshGrid2D<double> facemesh(z_1d, y_1d);
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fZ = facemesh.X();
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fY = facemesh.Y(); 
                    if (mBC.right == BC::DC && alpha>=0.0) {
                    // std::cout << "Boundary Right" << "\n";
                    addDenseBlock(i*locdim, i*locdim, Bright);
                    VectorXd gR =  mCase.g_right(fY, fZ, xLb);
                    // std::cout << "gR = \n" << gR << "\n"; 
                    // MatrixXd bR_mat = double(eps)*Jyz*dfac_x*(gR.transpose()*(kron(W_1, kron(dphi_n, W_1)))) + Jyz*(sigma0b/(std::pow(h_n_x, beta0)))*(gR.transpose()*kron(W_1, kron(en, W_1)));
                    // VectorXd bR = bR_mat.transpose().reshaped();
                    // VectorXd bR = flatten_transpose_rowmajor(bR_mat);
                    RowVectorXd dphi_n_row = dphi_n.transpose();
                    RowVectorXd en_row = en.transpose();
                    RowVectorXd bR_row = double(eps)*Jyz*dfac_x*(gR.transpose()*(kron(W_1, kron(dphi_n_row, W_1)))) + Jyz*(sigma0b/(std::pow(h_n_x, beta0)))*(gR.transpose()*kron(W_1, kron(en_row, W_1)));
                    VectorXd bR = bR_row.transpose();
                    // std::cout << "bR = \n" << bR << "\n"; 
                    b.segment(i*locdim, locdim) += bR;
                }
                    else if (mBC.right == BC::NC) {
                        VectorXd gNR = mCase.gN_right(fY, fZ, xLb);
                        RowVectorXd en_row = en.transpose();
                        RowVectorXd bNR_row = Jyz * (gNR.transpose() * kron(W_1, kron(en_row, W_1)));
                        b.segment(i*locdim, locdim) += bNR_row.transpose();
                    }
                }   
                if (face_lid == 4) { 
                    // Bottom face
                    MeshGrid2D<double> facemesh(z_1d, x_1d);
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fZ = facemesh.X();
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fX = facemesh.Y(); 
                    if (mBC.bottom == BC::DC && alpha>=0.0) {
                    // std::cout << "Boundary Bottom" << "\n";
                    addDenseBlock(i*locdim, i*locdim, Bbottom);
                    VectorXd gBo =  mCase.g_bottom(fX, fZ, yLa);
                    // std::cout << "gBo = \n" << gBo << "\n"; 
                    MatrixXd bBo_mat = -double(eps)*Jxz*dfac_y*kron(dphi_0, kron(W_1, W_1)*gBo) + Jxz*(sigma0b/(std::pow(h_n_y, beta0)))*kron(e0, kron(W_1, W_1)*gBo);
                    // VectorXd bBo = bBo_mat.transpose().reshaped();
                    VectorXd bBo = flatten_transpose_rowmajor(bBo_mat);
                    // std::cout << "bBo = \n" << bBo << "\n"; 
                    b.segment(i*locdim, locdim) += bBo;
                    }
                    else if (mBC.bottom == BC::NC) {
                        VectorXd gNBo = mCase.gN_bottom(fX, fZ, yLa);
                        MatrixXd bNBo_mat = Jxz * kron(e0, kron(W_1, W_1)*gNBo);
                        VectorXd bNBo = flatten_transpose_rowmajor(bNBo_mat);
                        b.segment(i*locdim, locdim) += bNBo;
                    }
                }
                if (face_lid == 5) { 
                    // Top edge
                    MeshGrid2D<double> facemesh(z_1d, x_1d);
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fZ = facemesh.X();
                    const Eigen::Array<double, Eigen::Dynamic, Eigen::Dynamic>& fX = facemesh.Y();  
                    // std::cout << "Boundary top" << "\n";
                    // std::cout << "a_y = " << a_y << "\n"; 
                    // std::cout << "b_y = " << b_y << "\n"; 
                    if (mBC.top == BC::DC && alpha>=0.0) {
                    addDenseBlock(i*locdim, i*locdim, Btop);
                    VectorXd gT =  mCase.g_top(fX, fZ, yLb);
                    // std::cout << "g_top = \n" << gT << "\n"; 
                    MatrixXd bT_mat = double(eps)*Jxz*dfac_y*kron(dphi_n, kron(W_1, W_1)*gT) + Jxz*(sigma0b/(std::pow(h_n_y, beta0)))*kron(en, kron(W_1, W_1)*gT);
                    // VectorXd bT = bT_mat.transpose().reshaped();
                    VectorXd bT = flatten_transpose_rowmajor(bT_mat);
                    // std::cout << "bT = \n" << bT << "\n"; 
                    b.segment(i*locdim, locdim) += bT;
                    }
                    else if (mBC.top == BC::NC) {
                        VectorXd gNT = mCase.gN_top(fX, fZ, yLb);
                        MatrixXd bNT_mat = Jxz * kron(en, kron(W_1, W_1)*gNT);
                        VectorXd bNT = flatten_transpose_rowmajor(bNT_mat);
                        b.segment(i*locdim, locdim) += bNT;
                    } 
                }
            }
        }
    }

    // Contributions of A from Numerical flux from connection table
    for (int i=0; i<mRecConnectionTable.getSize(); i++)
    {
        Connection connect = mRecConnectionTable[i];
        int N_E1 = std::get<0>(connect.ElmConnect); 
        int N_E2 = std::get<1>(connect.ElmConnect);

        int face_minus_id = std::get<0>(connect.FaceConnect);
        if (face_minus_id == 1 && alpha >= 0.0){ // Front-Back connect
            addDenseBlock(N_E1*locdim, N_E1*locdim, C11_FB);
            addDenseBlock(N_E2*locdim, N_E2*locdim, C22_FB);
            addDenseBlock(N_E1*locdim, N_E2*locdim, C12_FB);
            addDenseBlock(N_E2*locdim, N_E1*locdim, C21_FB);      
        }
        else if (face_minus_id == 3){ // Right-Left connect
            addDenseBlock(N_E1*locdim, N_E1*locdim, C11_RL);
            addDenseBlock(N_E2*locdim, N_E2*locdim, C22_RL);
            addDenseBlock(N_E1*locdim, N_E2*locdim, C12_RL);
            addDenseBlock(N_E2*locdim, N_E1*locdim, C21_RL);      
        }
        else if (face_minus_id == 5 && alpha >=0.0){ // Top-Bottom connect
            addDenseBlock(N_E1*locdim, N_E1*locdim, C11_TB);
            addDenseBlock(N_E2*locdim, N_E2*locdim, C22_TB);
            addDenseBlock(N_E1*locdim, N_E2*locdim, C12_TB);
            addDenseBlock(N_E2*locdim, N_E1*locdim, C21_TB);      
        }
    }


    // Build sparse matrix (sum duplicate triplets)
    SparseMatrix<double> A(dim, dim);
#if EIGEN_VERSION_AT_LEAST(3,3,0)
    A.setFromTriplets(triplets.begin(), triplets.end(), std::plus<double>());
#else
    A.setFromTriplets(triplets.begin(), triplets.end());
#endif
    A.makeCompressed();

    // // init numerical solution and solve with sparse solvers
    // VectorXd u_h;
    // if (method == Method::SIPG){
    //     // SPD -> Cholesky
    //     SimplicialLDLT<SparseMatrix<double>> solver;
    //     solver.compute(A);
    //     if (solver.info() != Success) {
    //         std::cerr << "SimplicialLDLT factorization failed, falling back to SparseLU" << std::endl;
    //         SparseLU<SparseMatrix<double>> lu;
    //         lu.analyzePattern(A);
    //         lu.factorize(A);
    //         u_h = lu.solve(b);
    //     } else {
    //         u_h = solver.solve(b);
    //     }
    // }
    // else{
    //     // Possibly nonsymmetric -> SparseLU
    //     SparseLU<SparseMatrix<double>> solver;
    //     solver.analyzePattern(A);
    //     solver.factorize(A);
    //     u_h = solver.solve(b);
    // }

    // init numerical solution and solve with MKL PARDISO
    VectorXd u_h;

    // PardisoLDLT is for SPD; PardisoLU is for general
    // alpha<0 is Helmholtz, even SIPG is indefinite
    if (method == Method::SIPG){
        if (alpha < 0.0) {
        // symmetric but indefinite -> LDLT cholesky
        Eigen::PardisoLDLT<SparseMatrix<double>> solver;
        std::cout << "SIPG and alpha<0, Symmetric+Indefinite, trying Eigen's MKL Pardiso's LDLT first"<< "\n";
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            std::cerr << "PardisoLDLT factorization failed, falling back to PardisoLU" << std::endl;
            Eigen::PardisoLU<SparseMatrix<double>> lu;
            lu.compute(A);
            if (lu.info() != Eigen::Success) {
                throw std::runtime_error("PARDISO factorization failed");
            }
            u_h = lu.solve(b);
        } else {
            u_h = solver.solve(b); 
        }
        }
        else{ 
        // SPD -> LLT cholesky
        Eigen::PardisoLLT<SparseMatrix<double>> solver;
        std::cout << "SIPG and alpha>0, SPD, using Eigen's MKL Pardiso's LLT"<< "\n";
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            std::cerr << "PardisoLDLT factorization failed, falling back to PardisoLU" << std::endl;
            Eigen::PardisoLU<SparseMatrix<double>> lu;
            lu.compute(A);
            if (lu.info() != Eigen::Success) {
                throw std::runtime_error("PARDISO factorization failed");
            }
            u_h = lu.solve(b);
        } else {
            u_h = solver.solve(b);
        }
    }
    } else {
        std::cout << "IIPG or NIPG or NIPG0, Non-Symmetric, using Eigen's MKL Pardiso's LU"<< "\n";
        // Possibly nonsymmetric -> MKL LU
        Eigen::PardisoLU<SparseMatrix<double>> solver;
        solver.compute(A);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("PARDISO factorization failed");
        }
        u_h = solver.solve(b);
    }
    
// (Optional) Check solve status
if ((A * u_h).isApprox(b) == false) {
    std::cerr << "Warning: solution residual is large\n";
}

    // std::cout << "A = \n" << Eigen::MatrixXd(A) << "\n";
    // std::cout << "b = \n" << b << "\n";
    // std::cout << "u_h = \n" << u_h << "\n";
    // std::cout << "u_exact = \n" << u_exact << "\n";


    std::string method_str = argv[2];
    std::cout << "3D Rectangle grid [ " << xLa << ", " << xLb << " ] x [ " << yLa << ", " << yLb << " ] x [ " << zLa << ", " << zLb << " ]"\
     " with Nel_x=" << Nel_x << ", Nel_y=" << Nel_y << ", Nel_z=" << Nel_z << ", and total Nel=" << Nel << " elements:" << "\n";
    std::cout << "DG method=" << method_str << ", order k=" << k << ", sigma0=" << sigma0 << "\n";
    std::cout << "Solver method=" <<  "Sparse Direct" << ", package=" << "Eigen intel MKL Pardiso LLT+LDLT+PardisoLU" << "\n";


    auto end = std::chrono::high_resolution_clock::now();
    auto duration = end - start;
    double seconds = std::chrono::duration<double>(duration).count();
    std::cout << "DG Solver Wall time: " << seconds << " seconds" << std::endl;


    // -----------------------------------------------------------------------
    //                     Compute 3D errors (relative)
    // -----------------------------------------------------------------------
    // Check the preliminary vector l2 error
    // VectorXd error = u_exact - u_h;
    // double l2_error = std::sqrt(error.dot(error))/std::sqrt(u_exact.dot(u_exact));
    // std::cout << "relative vector l2 error = " << l2_error << std::endl;
    DG::DGError3D err(k, W_1, D_1);
    err.reset();
    err.add_from_flat(u_h, u_exact, Nel_x, Nel_y, Nel_z, h_n_x, h_n_y, h_n_z);
    // Vector L2 (flat)
    // const double E_vec   = DG::DGError::vector_L2(u_h, u_exact);
    const double E_vec_rel = DG::DGError3D::vector_L2_rel(u_h, u_exact);

    // Broken norms
    // const double E_L2   = err.broken_L2();
    const double E_L2_rel = err.broken_L2_rel();
    //const double E_H1   = err.broken_H1();
    const double E_H1_rel = err.broken_H1_rel();

    std::cout << "Errors:\n"
            << "  relative ||u-uh||_2        = " << E_vec_rel << "\n"
            << "  relative broken L2         = " << E_L2_rel << "\n"
            << "  relative broken H1 (semi)  = " << E_H1_rel << "\n";


    // Check sparse matrix memory consumption
    std::cout << "Sparse global A approx memory consumption = " << sparseMemoryBytes(A)/1024.0/1024.0 << " MB" << std::endl;
    // -----------------------------------------------------------------------
    //                            HDF5 output
    // -----------------------------------------------------------------------
    if (write_h5=="w" || write_h5=="W") {
        // derive case name from JSON filename
        std::string case_name = case_path;
        {
            auto p = case_name.find_last_of("/\\");
            if (p != std::string::npos) case_name = case_name.substr(p+1);
            auto dot = case_name.find_last_of('.');
            if (dot != std::string::npos) case_name = case_name.substr(0, dot);
        }

        std::ostringstream oss;
        oss << "DG3DRec"
            << "_"  << case_name
            << "_method=" << method_str
            << "_Nel=" << Nel
            << "_k=" << k
            << "_sigma0=" << sigma0
            << ".h5";
        const std::string h5path = std::string("data/") + oss.str();
        std::cout << "Writing solution to file: \"" << h5path << "\"\n";

        HighFive::File file(h5path, HighFive::File::Overwrite);

        // Grid group (3D attrs)
        auto grid = file.createGroup("Grid");
        grid.createAttribute("xLa", xLa); grid.createAttribute("xLb", xLb);
        grid.createAttribute("yLa", yLa); grid.createAttribute("yLb", yLb);
        grid.createAttribute("zLa", zLa); grid.createAttribute("zLb", zLb);
        grid.createAttribute("Nel_x", Nel_x);
        grid.createAttribute("Nel_y", Nel_y);
        grid.createAttribute("Nel_z", Nel_z);

        // Exact solution (optional but useful for plotting)
        auto gex = file.createGroup("ExactSolution");
        gex.createDataSet("u_exact", u_exact);

        // Numerical solution + errors
        auto gnum = file.createGroup("NumericalSolution");
        gnum.createAttribute("method", method_str);
        gnum.createAttribute("Nel", Nel);
        gnum.createAttribute("k", k);
        gnum.createAttribute("sigma0", sigma0);
        gnum.createDataSet("u_h", u_h);
        gnum.createDataSet("l2_error",       E_vec_rel);
        gnum.createDataSet("broken_L2_rel",  E_L2_rel);
        gnum.createDataSet("broken_H1_rel",  E_H1_rel);
    }

    return 0;
}
