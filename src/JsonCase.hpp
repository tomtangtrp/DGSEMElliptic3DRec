#pragma once
#include <memory>
#include <string>
#include <Eigen/Core>
#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>  // Tensor used in public API

namespace DG {

class JsonCase {
public:
  explicit JsonCase(const std::string& path);
  ~JsonCase();
  JsonCase(JsonCase&&) noexcept;
  JsonCase& operator=(JsonCase&&) noexcept;
  JsonCase(const JsonCase&) = delete;
  JsonCase& operator=(const JsonCase&) = delete;

  double alpha() const;

  // ===== Volume fields on 3D nodal grids (Z fastest) =====
  Eigen::VectorXd f(const Eigen::Tensor<double,3,Eigen::RowMajor>& X,
                    const Eigen::Tensor<double,3,Eigen::RowMajor>& Y,
                    const Eigen::Tensor<double,3,Eigen::RowMajor>& Z) const;

  Eigen::VectorXd u(const Eigen::Tensor<double,3,Eigen::RowMajor>& X,
                    const Eigen::Tensor<double,3,Eigen::RowMajor>& Y,
                    const Eigen::Tensor<double,3,Eigen::RowMajor>& Z) const;

  bool has_exact() const;

  // ===== Dirichlet data on faces (flattened row-major vectors) =====
  // x = const faces (YZ grids)
  Eigen::VectorXd g_left (const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLa) const;
  Eigen::VectorXd g_right(const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLb) const;

  // y = const faces (XZ grids)
  Eigen::VectorXd g_bottom(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLa) const;
  Eigen::VectorXd g_top   (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLb) const;

  // z = const faces (XY grids)
  Eigen::VectorXd g_back (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLa) const;
  Eigen::VectorXd g_front(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLb) const;

  // ===== Neumann data on faces (same shapes/ordering) =====
  Eigen::VectorXd gN_left (const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLa) const;
  Eigen::VectorXd gN_right(const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLb) const;

  Eigen::VectorXd gN_bottom(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLa) const;
  Eigen::VectorXd gN_top   (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLb) const;

  Eigen::VectorXd gN_back (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLa) const;
  Eigen::VectorXd gN_front(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLb) const;

private:
  struct Impl;                    // heavy internals hidden
  std::unique_ptr<Impl> p_;
};

} // namespace DG
