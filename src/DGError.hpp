// DGError3D.hpp
#pragma once
#include <Eigen/Dense>
#include <cmath>
#include <limits>

namespace DG {

// Error metrics for DG on a rectangular, affine 3D mesh.
// Assumptions:
//  - Per-element nodal layout is (k+1)^3 with x-fastest local ordering:
//        idx_local = (kz*(k+1) + jy)*(k+1) + ix, where ix, jy, kz in [0..k]
//  - Global vectors u_h, u_exact are concatenations of element tiles in element-major order:
//        e = ((ez * Nel_y) + ey) * Nel_x + ex     (x fastest across elements)
//  - Quadrature: GLL nodes/weights; M1_ref is the (diagonal) 1D mass on reference [-1,1],
//    D1_ref is the 1D diff on reference. Physical derivatives use 2/h scaling.
//  - Mesh is axis-aligned; each element’s sizes are (hx, hy, hz).
//
// Provided metrics (absolute and relative):
//   - broken L2
//   - broken H1 (seminorm)
//   - vector L2 on flat vectors
class DGError3D {
public:
  DGError3D(int k, const Eigen::MatrixXd& M1_ref, const Eigen::MatrixXd& D1_ref)
  : k_(k), k1_(k+1), W1_(M1_ref.diagonal()), D1_(D1_ref)
  {
    reset();
    // No need to materialize a full 3D weight tensor; we multiply W1_ factors on the fly.
  }

  void reset() {
    sumL2_ = sumL2_exact_ = 0.0L;
    sumH1_ = sumH1_exact_ = 0.0L;
  }

  // Accumulate from one element cube.
  // Ue_flat: numerical; Ux_flat: exact, both length (k+1)^3 with x-fastest ordering.
  // hx, hy, hz are THIS element's physical sizes.
  void add_element(const Eigen::Ref<const Eigen::VectorXd>& Ue_flat,
                   const Eigen::Ref<const Eigen::VectorXd>& Ux_flat,
                   double hx, double hy, double hz)
  {
    const int nloc = k1_ * k1_ * k1_;
    // ---- L2 part ----
    {
      long double l2_num = 0.0L, l2_den = 0.0L;
      const long double J = 0.125L * hx * hy * hz; // |∂(ξ,η,ζ)/∂(x,y,z)|^{-1}
      for (int kz = 0; kz < k1_; ++kz) {
        const double wz = W1_(kz);
        for (int jy = 0; jy < k1_; ++jy) {
          const double wy = W1_(jy);
          const long double wyz = static_cast<long double>(wy * wz);
          for (int ix = 0; ix < k1_; ++ix) {
            const double wx = W1_(ix);
            const long double w = static_cast<long double>(wx) * wyz;
            const int idx = (kz*k1_ + jy)*k1_ + ix;
            const double d  = Ux_flat(idx) - Ue_flat(idx);
            l2_num += static_cast<long double>(d*d) * w;
            l2_den += static_cast<long double>(Ux_flat(idx)*Ux_flat(idx)) * w;
          }
        }
      }
      sumL2_       += l2_num * J;
      sumL2_exact_ += l2_den * J;
    }

    // ---- H1(semi) part ----
    {
      const double sx = 2.0 / hx;
      const double sy = 2.0 / hy;
      const double sz = 2.0 / hz;
      const long double J = 0.125L * hx * hy * hz;

      long double h1_num = 0.0L, h1_den = 0.0L;

      // We compute directional derivatives by lines and D1_:
      // x-lines (vary ix for fixed (jy,kz)), y-lines, z-lines.
      Eigen::VectorXd u(k1_), v(k1_), du(k1_), dv(k1_);

      // d/dx lines: for each (jy,kz), contiguous ix-line
      for (int kz = 0; kz < k1_; ++kz) {
        const double wz = W1_(kz);
        for (int jy = 0; jy < k1_; ++jy) {
          const double wy = W1_(jy);
          const long double wyz = static_cast<long double>(wy * wz);
          // load x-line
          for (int ix = 0; ix < k1_; ++ix) {
            const int idx = (kz*k1_ + jy)*k1_ + ix;
            u(ix) = Ue_flat(idx);
            v(ix) = Ux_flat(idx);
          }
          du = sx * (D1_ * u);
          dv = sx * (D1_ * v);
          for (int ix = 0; ix < k1_; ++ix) {
            const double wx = W1_(ix);
            const long double w = static_cast<long double>(wx) * wyz;
            const double de = dv(ix) - du(ix);
            h1_num += static_cast<long double>(de*de) * w;
            h1_den += static_cast<long double>(dv(ix)*dv(ix)) * w;
          }
        }
      }

      // d/dy lines: for each (ix,kz), j-line (stride k1_)
      for (int kz = 0; kz < k1_; ++kz) {
        const double wz = W1_(kz);
        for (int ix = 0; ix < k1_; ++ix) {
          const double wx = W1_(ix);
          for (int jy = 0; jy < k1_; ++jy) {
            const int idx = (kz*k1_ + jy)*k1_ + ix;
            u(jy) = Ue_flat(idx);
            v(jy) = Ux_flat(idx);
          }
          du = sy * (D1_ * u);
          dv = sy * (D1_ * v);
          for (int jy = 0; jy < k1_; ++jy) {
            const double wy = W1_(jy);
            const long double w = static_cast<long double>(wx * wy * wz);
            const double de = dv(jy) - du(jy);
            h1_num += static_cast<long double>(de*de) * w;
            h1_den += static_cast<long double>(dv(jy)*dv(jy)) * w;
          }
        }
      }

      // d/dz lines: for each (ix,jy), k-line (hop k1_*k1_)
      for (int jy = 0; jy < k1_; ++jy) {
        const double wy = W1_(jy);
        for (int ix = 0; ix < k1_; ++ix) {
          const double wx = W1_(ix);
          for (int kz = 0; kz < k1_; ++kz) {
            const int idx = (kz*k1_ + jy)*k1_ + ix;
            u(kz) = Ue_flat(idx);
            v(kz) = Ux_flat(idx);
          }
          du = sz * (D1_ * u);
          dv = sz * (D1_ * v);
          for (int kz = 0; kz < k1_; ++kz) {
            const double wz = W1_(kz);
            const long double w = static_cast<long double>(wx * wy * wz);
            const double de = dv(kz) - du(kz);
            h1_num += static_cast<long double>(de*de) * w;
            h1_den += static_cast<long double>(dv(kz)*dv(kz)) * w;
          }
        }
      }

      sumH1_       += h1_num * J;
      sumH1_exact_ += h1_den * J;
    }
  }

  // Accumulate for ALL elements from flat global vectors (x-fastest local tiles).
  // u_h, u_exact: length Nel*(k+1)^3. Nel_x, Nel_y, Nel_z element counts; hx, hy, hz sizes.
  void add_from_flat(const Eigen::VectorXd& u_h,
                     const Eigen::VectorXd& u_exact,
                     int Nel_x, int Nel_y, int Nel_z,
                     double hx, double hy, double hz)
  {
    const int Nel = Nel_x * Nel_y * Nel_z;
    const int locdim = k1_ * k1_ * k1_;
    for (int e = 0; e < Nel; ++e) {
      Eigen::Map<const Eigen::VectorXd>
          Uh(u_h.data()     + static_cast<Eigen::Index>(e) * locdim, locdim),
          Ux(u_exact.data() + static_cast<Eigen::Index>(e) * locdim, locdim);
      add_element(Uh, Ux, hx, hy, hz);
    }
  }

  // Absolute broken norms
  double broken_L2() const { return std::sqrt(static_cast<double>(std::max(0.0L, sumL2_))); }
  double broken_H1() const { return std::sqrt(static_cast<double>(std::max(0.0L, sumH1_))); }

  // Relative broken norms (NaN if exact integral is zero)
  double broken_L2_rel() const {
    return (sumL2_exact_ > 0.0L)
      ? std::sqrt(static_cast<double>(sumL2_ / sumL2_exact_))
      : std::numeric_limits<double>::quiet_NaN();
  }
  double broken_H1_rel() const {
    return (sumH1_exact_ > 0.0L)
      ? std::sqrt(static_cast<double>(sumH1_ / sumH1_exact_))
      : std::numeric_limits<double>::quiet_NaN();
  }

  // Vector L2 (Euclidean) on flat stacked vectors
  static double vector_L2(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    return (a-b).norm();
  }
  static double vector_L2_rel(const Eigen::VectorXd& a, const Eigen::VectorXd& b) {
    const double denom = b.norm();
    return (denom > 0.0) ? (a-b).norm() / denom : std::numeric_limits<double>::quiet_NaN();
  }

private:
  int k_, k1_;
  Eigen::VectorXd W1_;     // (k+1) 1D weights on reference (diag(M1_ref))
  Eigen::MatrixXd D1_;     // (k+1)x(k+1) 1D diff on reference

  long double sumL2_ = 0.0L, sumL2_exact_ = 0.0L;
  long double sumH1_ = 0.0L, sumH1_exact_ = 0.0L;
};

} // namespace DG
