#ifndef MESHGRID3DHEADERDEF
#define MESHGRID3DHEADERDEF

#include <Eigen/Dense>
#include <unsupported/Eigen/CXX11/Tensor>

template<typename Scalar = double>
class MeshGrid3D {
public:
    using Array1 = Eigen::Array<Scalar, Eigen::Dynamic, 1>;
    using Array3 = Eigen::Tensor<Scalar, 3, Eigen::RowMajor>;

    // x size Nx, y size Ny, z size Nz
    MeshGrid3D(const Array1& x, const Array1& y, const Array1& z)
        : x_(x), y_(y), z_(z), computed_(false) {}

    const Array3& X() { compute(); return XX_; }
    const Array3& Y() { compute(); return YY_; }
    const Array3& Z() { compute(); return ZZ_; }

private:
    Array1 x_, y_, z_;
    Array3 XX_, YY_, ZZ_;
    bool computed_;

    void compute() {
        if (computed_) return;

        const int Nx = x_.size();
        const int Ny = y_.size();
        const int Nz = z_.size();

        XX_.resize(Ny, Nx, Nz); // dimensions ordered: y, x, z
        YY_.resize(Ny, Nx, Nz);
        ZZ_.resize(Ny, Nx, Nz);

        for (int iy = 0; iy < Ny; ++iy) {
            for (int ix = 0; ix < Nx; ++ix) {
                for (int iz = 0; iz < Nz; ++iz) {
                    XX_(iy, ix, iz) = x_(ix);
                    YY_(iy, ix, iz) = y_(iy);
                    ZZ_(iy, ix, iz) = z_(iz);
                }
            }
        }
        computed_ = true;
    }
};

#endif
