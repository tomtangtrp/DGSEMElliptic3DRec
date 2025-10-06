#include "JsonCase.hpp"
#include <fstream>
#include <unordered_map>
#include <optional>
#include <stdexcept>
#include <sstream>

// heavy includes live here
#include "json.hpp"      // nlohmann::json
#include "exprtk.hpp"    // ExprTk
#include <unsupported/Eigen/CXX11/Tensor>

using json = nlohmann::json;

namespace DG {

// ========== Impl with 3D expressions ==========
struct JsonCase::Impl {
  struct Expr3D {
    exprtk::symbol_table<double> sym;
    exprtk::expression<double>   ex;
    exprtk::parser<double>       parser;
    mutable double x=0.0, y=0.0, z=0.0;
    bool compiled=false;

    Expr3D(){
      sym.add_variable("x",x);
      sym.add_variable("y",y);
      sym.add_variable("z",z);
      ex.register_symbol_table(sym);
    }
    void bind_consts(const std::unordered_map<std::string,double>& P){
      for (const auto& kv: P) sym.add_constant(kv.first, kv.second);
    }
    void compile_or_throw(const std::string& code, const char* tag){
      if (code.empty()) throw std::runtime_error(std::string("empty expr: ")+tag);
      compiled = parser.compile(code, ex);
      if (!compiled) {
        std::ostringstream oss; oss << "bad expr: " << tag;
        if (parser.error_count()) {
          oss << " (";
          for (std::size_t i=0;i<parser.error_count();++i){
            auto er = parser.get_error(i);
            if (i) oss << "; ";
            oss << "#" << i << ": " << er.diagnostic
                << " [token='" << er.token.value << "' pos=" << er.token.position << "]";
          }
          oss << ")";
        }
        throw std::runtime_error(oss.str());
      }
    }
    inline double eval(double xx, double yy, double zz) const {
      x=xx; y=yy; z=zz;
      return ex.value();
    }

    // Face evaluators that return 2D arrays (row-major)
    Eigen::ArrayXXd evalYZ(const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xfix){
      Eigen::ArrayXXd out(Y.rows(), Y.cols());
      for (int i=0;i<Y.rows();++i)
        for (int j=0;j<Y.cols();++j)
          out(i,j) = eval(xfix, Y(i,j), Z(i,j));
      return out;
    }
    Eigen::ArrayXXd evalXZ(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yfix){
      Eigen::ArrayXXd out(X.rows(), X.cols());
      for (int i=0;i<X.rows();++i)
        for (int j=0;j<X.cols();++j)
          out(i,j) = eval(X(i,j), yfix, Z(i,j));
      return out;
    }
    Eigen::ArrayXXd evalXY(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zfix){
      Eigen::ArrayXXd out(X.rows(), X.cols());
      for (int i=0;i<X.rows();++i)
        for (int j=0;j<X.cols();++j)
          out(i,j) = eval(X(i,j), Y(i,j), zfix);
      return out;
    }
  };

  std::unordered_map<std::string,double> P;
  double alpha = 1.0;

  // Volume expressions
  Expr3D f;                    // required
  std::optional<Expr3D> u;     // optional exact solution

  // Dirichlet faces
  std::optional<Expr3D> dL,dR,dB,dT,dBack,dFront;

  // Neumann faces
  std::optional<Expr3D> nL,nR,nB,nT,nBack,nFront;

  // ----- utils -----
  static std::string slurp(const std::string& path){
    std::ifstream in(path);
    if (!in) throw std::runtime_error("Cannot open JSON case file: "+path);
    return std::string((std::istreambuf_iterator<char>(in)), std::istreambuf_iterator<char>());
  }
  static std::optional<std::string> get_opt(const json& j, const char* key){
    if (!j.contains(key) || j.at(key).is_null()) return std::nullopt;
    return j.at(key).get<std::string>();
  }

  static void maybe_compile(std::optional<Expr3D>& dst,
                            const std::optional<std::string>& code,
                            const std::unordered_map<std::string,double>& P,
                            const char* tag){
    if (!code) return;
    dst.emplace();
    dst->bind_consts(P);
    dst->compile_or_throw(*code, tag);
  }

  explicit Impl(const std::string& path){
    auto txt = slurp(path);
    json j = json::parse(txt, nullptr, true, /*ignore_comments=*/true);

    if (j.contains("parameters")) {
      for (auto it=j["parameters"].begin(); it!=j["parameters"].end(); ++it)
        P[it.key()] = it.value().get<double>();
    }
    if (j.contains("alpha")) alpha = j.at("alpha").get<double>();

    f.bind_consts(P);
    f.compile_or_throw(j.at("source").get<std::string>(), "source");

    if (j.contains("exact") && j["exact"].contains("u")){
      u.emplace(); u->bind_consts(P);
      u->compile_or_throw(j["exact"]["u"].get<std::string>(), "exact.u");
    }

    if (j.contains("dirichlet")){
      const auto& d=j["dirichlet"];
      maybe_compile(dL,     get_opt(d,"left"),   P, "dirichlet.left");
      maybe_compile(dR,     get_opt(d,"right"),  P, "dirichlet.right");
      maybe_compile(dB,     get_opt(d,"bottom"), P, "dirichlet.bottom");
      maybe_compile(dT,     get_opt(d,"top"),    P, "dirichlet.top");
      maybe_compile(dBack,  get_opt(d,"back"),   P, "dirichlet.back");
      maybe_compile(dFront, get_opt(d,"front"),  P, "dirichlet.front");
    }
    if (j.contains("neumann")){
      const auto& n=j["neumann"];
      maybe_compile(nL,     get_opt(n,"left"),   P, "neumann.left");
      maybe_compile(nR,     get_opt(n,"right"),  P, "neumann.right");
      maybe_compile(nB,     get_opt(n,"bottom"), P, "neumann.bottom");
      maybe_compile(nT,     get_opt(n,"top"),    P, "neumann.top");
      maybe_compile(nBack,  get_opt(n,"back"),   P, "neumann.back");
      maybe_compile(nFront, get_opt(n,"front"),  P, "neumann.front");
    }
  }
};

// ----- local helpers -----
static Eigen::VectorXd flatten_rowmajor(const Eigen::ArrayXXd& A){
  const int R = A.rows(), C = A.cols();
  Eigen::VectorXd v(static_cast<Eigen::Index>(R)*C);
  Eigen::Index p = 0;
  for (int i=0;i<R;++i)
    for (int j=0;j<C;++j)
      v(p++) = A(i,j);
  return v;
}

// z-fastest flatten: p = (iy*nx + ix)*nz + iz
template <typename Expr>
static Eigen::VectorXd eval_on_tensor_rowmajor(const Expr& ex,
    const Eigen::Tensor<double,3,Eigen::RowMajor>& X,
    const Eigen::Tensor<double,3,Eigen::RowMajor>& Y,
    const Eigen::Tensor<double,3,Eigen::RowMajor>& Z)
{
  const int ny = static_cast<int>(X.dimension(0));
  const int nx = static_cast<int>(X.dimension(1));
  const int nz = static_cast<int>(X.dimension(2));
  Eigen::VectorXd out(static_cast<Eigen::Index>(ny)*nx*nz);
  for (int iy=0; iy<ny; ++iy)
    for (int ix=0; ix<nx; ++ix)
      for (int iz=0; iz<nz; ++iz){
        const Eigen::Index p = static_cast<Eigen::Index>((iy*nx + ix)*nz + iz);
        out(p) = ex.eval(X(iy,ix,iz), Y(iy,ix,iz), Z(iy,ix,iz));
      }
  return out;
}

// ===== JsonCase public API =====

JsonCase::JsonCase(const std::string& path) : p_(std::make_unique<Impl>(path)) {}
JsonCase::~JsonCase() = default;
JsonCase::JsonCase(JsonCase&&) noexcept = default;
JsonCase& JsonCase::operator=(JsonCase&&) noexcept = default;

double JsonCase::alpha() const { return p_->alpha; }

Eigen::VectorXd JsonCase::f(const Eigen::Tensor<double,3,Eigen::RowMajor>& X,
                            const Eigen::Tensor<double,3,Eigen::RowMajor>& Y,
                            const Eigen::Tensor<double,3,Eigen::RowMajor>& Z) const {
  return eval_on_tensor_rowmajor(p_->f, X, Y, Z);
}

Eigen::VectorXd JsonCase::u(const Eigen::Tensor<double,3,Eigen::RowMajor>& X,
                            const Eigen::Tensor<double,3,Eigen::RowMajor>& Y,
                            const Eigen::Tensor<double,3,Eigen::RowMajor>& Z) const {
  if (!p_->u) {
    const int ny = static_cast<int>(X.dimension(0));
    const int nx = static_cast<int>(X.dimension(1));
    const int nz = static_cast<int>(X.dimension(2));
    return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(ny)*nx*nz);
  }
  return eval_on_tensor_rowmajor(*p_->u, X, Y, Z);
}

bool JsonCase::has_exact() const { return static_cast<bool>(p_->u); }

// ----- Dirichlet faces -----
Eigen::VectorXd JsonCase::g_left (const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLa) const {
  if (p_->dL) return flatten_rowmajor(p_->dL->evalYZ(Y,Z,xLa));
  if (p_->u)  return flatten_rowmajor(p_->u ->evalYZ(Y,Z,xLa));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(Y.rows()*Y.cols()));
}
Eigen::VectorXd JsonCase::g_right(const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLb) const {
  if (p_->dR) return flatten_rowmajor(p_->dR->evalYZ(Y,Z,xLb));
  if (p_->u)  return flatten_rowmajor(p_->u ->evalYZ(Y,Z,xLb));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(Y.rows()*Y.cols()));
}

Eigen::VectorXd JsonCase::g_bottom(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLa) const {
  if (p_->dB) return flatten_rowmajor(p_->dB->evalXZ(X,Z,yLa));
  if (p_->u)  return flatten_rowmajor(p_->u ->evalXZ(X,Z,yLa));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}
Eigen::VectorXd JsonCase::g_top   (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLb) const {
  if (p_->dT) return flatten_rowmajor(p_->dT->evalXZ(X,Z,yLb));
  if (p_->u)  return flatten_rowmajor(p_->u ->evalXZ(X,Z,yLb));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}

Eigen::VectorXd JsonCase::g_back (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLa) const {
  if (p_->dBack) return flatten_rowmajor(p_->dBack->evalXY(X,Y,zLa));
  if (p_->u)     return flatten_rowmajor(p_->u    ->evalXY(X,Y,zLa));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}
Eigen::VectorXd JsonCase::g_front(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLb) const {
  if (p_->dFront) return flatten_rowmajor(p_->dFront->evalXY(X,Y,zLb));
  if (p_->u)      return flatten_rowmajor(p_->u    ->evalXY(X,Y,zLb));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}

// ----- Neumann faces -----
Eigen::VectorXd JsonCase::gN_left (const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLa) const {
  if (p_->nL) return flatten_rowmajor(p_->nL->evalYZ(Y,Z,xLa));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(Y.rows()*Y.cols()));
}
Eigen::VectorXd JsonCase::gN_right(const Eigen::ArrayXXd& Y, const Eigen::ArrayXXd& Z, double xLb) const {
  if (p_->nR) return flatten_rowmajor(p_->nR->evalYZ(Y,Z,xLb));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(Y.rows()*Y.cols()));
}

Eigen::VectorXd JsonCase::gN_bottom(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLa) const {
  if (p_->nB) return flatten_rowmajor(p_->nB->evalXZ(X,Z,yLa));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}
Eigen::VectorXd JsonCase::gN_top   (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Z, double yLb) const {
  if (p_->nT) return flatten_rowmajor(p_->nT->evalXZ(X,Z,yLb));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}

Eigen::VectorXd JsonCase::gN_back (const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLa) const {
  if (p_->nBack) return flatten_rowmajor(p_->nBack->evalXY(X,Y,zLa));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}
Eigen::VectorXd JsonCase::gN_front(const Eigen::ArrayXXd& X, const Eigen::ArrayXXd& Y, double zLb) const {
  if (p_->nFront) return flatten_rowmajor(p_->nFront->evalXY(X,Y,zLb));
  return Eigen::VectorXd::Zero(static_cast<Eigen::Index>(X.rows()*X.cols()));
}

} // namespace DG
