
#include <ranges>
#include <algorithm>
#include <vector>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef std::vector<double> value_t;
typedef std::vector<std::vector<value_t>> vvalue_t;
typedef std::vector<value_t> v1Dvalue_t;

typedef std::array<double, 2> point2D_t;
typedef std::vector<std::vector<point2D_t>> vpoint2D_t;

std::vector<double> operator+(const std::vector<double> &v1, const std::vector<double> &v2)
{
    assert(v1.size() == v2.size());
    std::vector<double> res(v1.size());
    std::ranges::transform(v1, v2, res.begin(), std::plus<double>{});
    return res;
}

std::vector<double> operator-(const std::vector<double> &v1, const std::vector<double> &v2)
{
    assert(v1.size() == v2.size());
    std::vector<double> res(v1.size());
    std::ranges::transform(v1, v2, res.begin(), std::minus<double>{});
    return res;
}

std::vector<double> operator/(const std::vector<double> &v1, const std::vector<double> &v2) {
    assert(v1.size() == v2.size());
    std::vector<double> res(v1.size());
    std::ranges::transform(v1, v2, res.begin(), std::divides<double>{});
    return res;
}

std::vector<double> operator*(const std::vector<double> &v1, double alpha)
{
    std::vector<double> res(v1.size());
    std::ranges::transform(v1, res.begin(), [alpha](double el)
                           { return alpha * el; });
    return res;
}

std::vector<double> operator*(double alpha, const std::vector<double> &v)
{
    return v * alpha;
}



namespace utils {
    void resize_vec2D(auto& vec, size_t nx, size_t ny) {
        vec.resize(nx);
        for(auto& row: vec) {
            row.resize(ny);
        }
    }

    void transform_2d(const auto& source, auto& des, auto f) {
        const size_t nx = source.size();
        const size_t ny = source[0].size();
        for (int i = 0; i < nx; i++)
        {
            for (int j = 0; j < ny; j++)
            {
                des[i][j] = f(source[i][j]);
            }
        }
    }

    void transform2d_without_ghostcells(const auto& source, auto& des, const size_t nx_g, const size_t ny_g, auto f) {
        const size_t nx_tot = source.size();
        const size_t ny_tot = source[0].size();
        const size_t nx = nx_tot - nx_g;
        const size_t ny = ny_tot - ny_g;

        assert(des.size() == nx);
        assert(des[0].size() == ny);

        for(int i = 0; i < nx; ++i) {
            for(int j = 0; j < ny; ++j) {
                des[i][j] = f(source[i + nx_g][j + ny_g]);
            }
        }
    }

    class ColumnView
    {
    public:
        ColumnView(vvalue_t &matrix, size_t col)
            : matrix(matrix), col(col) {}

        auto as_range()
        {
            return std::ranges::views::transform(matrix, [&](auto &row)
                                                { return row[col]; });
        }

    private:
        vvalue_t &matrix;
        size_t col;
    };
}

// Initial conditions are always given in primitive variables
namespace initialconditions {

    value_t test0_2d(point2D_t x)
    {
        return x[0] <= 0.25 ? value_t({1.0, 1.0, 0.0, 1.0}) : value_t({0.1, 1.0, 0.0, 1.});
    }

    value_t test3_2d(point2D_t x)
    {
        double x_0 = 0.5;
        double y_0 = 0.5;
        double R = 0.4;
        if ((x[0] - x_0) * (x[0] - x_0) + (x[1] - y_0) * (x[1] - y_0) < R * R)
        {
            return {1, 0, 0, 1};
        }
        else
        {
            return {0.125, 0, 0, 0.1};
        }
    }

}

namespace boundaryconditions {

    void transmissive(vvalue_t& solution, size_t nx_g, size_t ny_g) {
        const size_t nx_total = solution.size();
        const size_t nx = nx_total - 2 * nx_g;

        const size_t ny_total = solution[0].size();
        const size_t ny = ny_total - 2 * ny_g;


        for (int j = 0; j < ny + 2 * ny_g; ++j)
        {
            for(int i = 0; i < nx_g; ++i) {
                solution[i][j] = solution[nx_g][j];
            }

            for(int i = nx + nx_g; i < nx + 2 * nx_g; ++i) {
                solution[i][j] = solution[nx + nx_g - 1][j];
            }
        }

        for (int i = 0; i < nx + 2 * nx_g; ++i)
        {
            for(int j = 0; j < ny_g; ++j) {
                solution[i][j] = solution[i][ny_g];
            }
            for(int j = ny + ny_g; j < ny + 2 * ny_g; ++j) {
                solution[i][j] = solution[i][ny + ny_g - 1];
            }
        }
    }
}

namespace FluxDirection {
    struct F{};
    struct G{};
}

namespace NumericalFlux {
    struct FORCE{};
    struct BACKWARD{};
    struct FORWARD{};
    struct CENTRAL{};
    struct RI{};
    struct LF{};
    struct SLIC{};
}

namespace FluxFunction {
    struct Euler{};
};

namespace EOS {
    struct IdealGas{};
};

namespace SOURCE {
    struct NullSource{};
    struct CylindricalGeometric{};
}

template <typename NumericalFluxT, typename FluxFunctionT, typename EOST, typename SourceT>
class ConservationProblem {
public:
    ConservationProblem(NumericalFluxT numericalFluxT, FluxFunctionT fluxFunctionT, EOST eosT, SourceT sourceT, double gamma, double dx, double dy, double C = 0.8):
        _numericalFluxT(numericalFluxT), _fluxFunctionT(fluxFunctionT), _eosT(eosT), _sourceT(sourceT), _gamma(gamma), _dx(dx), _dy(dy), _C(C) {};

    double getMaxDt(std::ranges::input_range auto &&solution) {
        return _getMaxDT(this->_fluxFunctionT, solution);
    };

    double getSpeedOfSound(const value_t& u) {
        return _getSpeedOfSound(this->_fluxFunctionT, u);
    }
    
    v1Dvalue_t getFlux(std::ranges::input_range auto const& solution, size_t nk_g, double dt, auto fdt) {
        // Compile Time optimised
        if constexpr (std::is_same_v<NumericalFluxT, NumericalFlux::SLIC>) {
            return getSlopeLimitedFlux(solution, nk_g, dt, fdt);
        } else {
            return getFirstOrderFlux(solution, nk_g, dt, fdt);
        }
    }

    v1Dvalue_t getFirstOrderFlux(std::ranges::input_range auto const& solution, size_t nk_g, double dt, auto fdt) {
        const size_t nk_tot = solution.size();
        const size_t nk = nk_tot - 2 * nk_g;

        // This shouldn't be too expensive if nx/ny are around the same dimensionality
        this->_flux.resize(nk + 1);
        this->_fluxFunctionAtCells.resize(nk_tot);

        for (size_t i = 0; i < nk_tot; ++i)
        {
            _fluxFunctionAtCells[i] = f(fdt, solution[i]);
        }

        for (size_t i = 0; i <= nk; ++i)
        {
            _flux[i] = _getNumericalFlux(this->_numericalFluxT, fdt, solution, i + nk_g, dt);
        }

        return _flux;
    }


    v1Dvalue_t getSlopeLimitedFlux(std::ranges::input_range auto const& solution, size_t nk_g, double dt, auto fdt) {
        const size_t nk_tot = solution.size();
        const size_t nk = nk_tot - 2 * nk_g;

        // This shouldn't be too expensive if nx/ny are around the same dimensionality
        this->_solution_buffer.resize(nk_tot);
        this->_vSlopeLimiterR.resize(nk + 2);
        this->_vSlopeLimiterPhi.resize(nk + 2);

        this->_flux.resize(nk + 1);
        // this->_fluxFunctionAtCells.resize(nk_tot);
        v1Dvalue_t delta(nk + 2);

        // Compute delta
        for(size_t i = 0; i < nk + 2; ++i) {
            delta[i] = 0.5 * (solution[nk_g + i] - solution[nk_g + i - 2]);
        }

        // Compute r
        for(size_t i = 0; i < nk + 2; ++i) {
            _vSlopeLimiterR[i] = (solution[nk_g + i - 1] - solution[nk_g + i - 2]) / (solution[nk_g + i] - solution[nk_g + i - 1]) ;
        }
    
        // Compute phi
        std::transform(_vSlopeLimiterR.begin(), _vSlopeLimiterR.end(), _vSlopeLimiterPhi.begin(), [&](value_t r){
            // limiting only on Energy
            double rE = r[3];
            
            // MinBee
            if(rE <= 0) return 0.;
            else if(0 < rE && rE <= 1) return rE;
            return std::min<double>(1., 2. / (1. + rE));
        });

        const double _dtdx = dt / _dx;
        for (size_t i = 0; i <= nk; ++i)
        {
            value_t uL = solution[nk_g + i - 1];
            double phiL = _vSlopeLimiterPhi[i];
            value_t deltaL = delta[i];
    
            value_t uR = solution[nk_g + i];
            double phiR = _vSlopeLimiterPhi[i + 1];
            value_t deltaR = delta[i + 1];

            // data reconstruction
            uL = uL + 0.5 * phiL * deltaL;
            uR = uR + 0.5 * phiR * deltaR;

            // Half timestep update
            value_t update = 0.5 * _dtdx * (f(fdt, uR) - f(fdt, uL));
            uL = uL - update;
            uR = uR - update;

            value_t fL = f(fdt, uL);
            value_t fR = f(fdt, uR);

            _flux[i] = _getNumericalFlux(NumericalFlux::FORCE{}, fdt, uL, uR, fL, fR, dt);
        }

        return _flux;
    }

    // Source term
    value_t getSource(const value_t& u, const point2D_t& cell_centroid)
    {
        return _getSource(this->_sourceT, u, cell_centroid);
    };

    // Flux Fonctions
    value_t f(auto fdt, const value_t& u) {
        return _f(this->_fluxFunctionT, fdt, u);
    }
    value_t _f(FluxFunction::Euler, FluxDirection::F, const value_t& u);
    value_t _f(FluxFunction::Euler, FluxDirection::G, const value_t& u);

    // EOS
    vvalue_t conservativeToPrimitive(const vvalue_t& solution);
    value_t  conservativeToPrimitive(const value_t&  u);
    vvalue_t primitiveToConservative(const vvalue_t& primitive_solution);
    value_t  primitiveToConservative(const value_t&  q);

private:
    // Specialisation
    NumericalFluxT _numericalFluxT;
    FluxFunctionT _fluxFunctionT;
    EOST _eosT;
    SourceT _sourceT;

    v1Dvalue_t _flux;
    v1Dvalue_t _fluxFunctionAtCells;

    v1Dvalue_t _solution_buffer;
    v1Dvalue_t _vSlopeLimiterR;
    std::vector<double> _vSlopeLimiterPhi;


    double _dx;
    double _dy;
    double _gamma;
    double _C;
    
    // CFL conditions
    double _getMaxDT(FluxFunction::Euler, std::ranges::input_range auto &&u);

    // Speed of sound
    double _getSpeedOfSound(FluxFunction::Euler, const value_t& u);

    // Numerical Fluxes
    value_t _getNumericalFlux(NumericalFlux::BACKWARD, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::FORWARD,  auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::FORCE,    auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::CENTRAL,  auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::RI,       auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::LF,       auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);

    value_t _getNumericalFlux(NumericalFlux::LF,       auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt);
    value_t _getNumericalFlux(NumericalFlux::RI,       auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt);
    value_t _getNumericalFlux(NumericalFlux::FORCE,    auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt);
    

    // EOS
    value_t  _conservativeToPrimitive(EOS::IdealGas, const value_t& u);
    value_t  _primitiveToConservative(EOS::IdealGas, const value_t& u);

    // Source terms
    value_t _getSource(SOURCE::NullSource,           const value_t& u, const point2D_t& cell_centroid);
    value_t _getSource(SOURCE::CylindricalGeometric, const value_t& u, const point2D_t& cell_centroid);
};

// --------------------- Speed of sound ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::_getSpeedOfSound(FluxFunction::Euler, const value_t& u) {
    return sqrt(this->_gamma * (this->_gamma - 1.) * (u[3] - 0.5 * u[1] * u[1] / u[3]) / u[0]);
}

// --------------------- CFL Condition ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::_getMaxDT(FluxFunction::Euler, std::ranges::input_range auto &&solution) {
    double max_a = 0;
    for(const auto& row: solution) {
        for (const auto &el : row)
        {
            const double a = sqrt((el[1] * el[1] + el[2] * el[2]) / (el[0]*el[0])) + this->getSpeedOfSound(el);
            max_a = std::max<double>(max_a, a);
        }
    }
    return _C * std::min(_dx, _dy) / max_a;
}

// --------------------- Numerical Fluxes ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::BACKWARD, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    return _fluxFunctionAtCells[face_idx - 1];
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::FORWARD, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    return _fluxFunctionAtCells[face_idx];
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::CENTRAL, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    return 0.5 * (_fluxFunctionAtCells[face_idx - 1] + _fluxFunctionAtCells[face_idx]);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::LF, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    return 0.5 * (_dx / dt) * (u[face_idx - 1] - u[face_idx]) + 0.5 * (_fluxFunctionAtCells[face_idx - 1] + _fluxFunctionAtCells[face_idx]);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::RI, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    const value_t u_n = 0.5 * (u[face_idx - 1] + u[face_idx] + dt / _dx * (_fluxFunctionAtCells[face_idx - 1] - _fluxFunctionAtCells[face_idx]));
    return f(fdt, u_n);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::FORCE, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    return 0.5 * (this->_getNumericalFlux(NumericalFlux::LF{}, fdt, u, face_idx, dt) + this->_getNumericalFlux(NumericalFlux::RI{}, fdt, u, face_idx, dt));
}


// Need to implement this interface for 2nd order schemes
template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::LF, auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt) {
    return 0.5 * (_dx / dt) * (uL - uR) + 0.5 * (fL + fR);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::RI, auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt) {
    const value_t u_n = 0.5 * (uL + uR + dt / _dx * (fL - fR));
    return f(fdt, u_n);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::FORCE, auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt) {
    return 0.5 * (this->_getNumericalFlux(NumericalFlux::LF{}, fdt, uL, uR, fL, fR, dt) + this->_getNumericalFlux(NumericalFlux::RI{}, fdt, uL, uR, fL, fR, dt));
}

// --------------------- Flux Fonctions ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_f(FluxFunction::Euler, FluxDirection::F, const value_t& u)
{
    auto q = this->conservativeToPrimitive(u);
    const auto &rho = q[0];
    const auto &vx = q[1];
    const auto &vy = q[2];
    const auto &p = q[3];
    const auto &E = u[3];

    return {
        rho * vx,
        rho * vx * vx + p,
        rho * vx * vy,
        (E + p) * vx,
    };
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_f(FluxFunction::Euler, FluxDirection::G, const value_t& u) {
    auto q = this->conservativeToPrimitive(u);
    const auto &rho = q[0];
    const auto &vx = q[1];
    const auto &vy = q[2];
    const auto &p = q[3];
    const auto &E = u[3];

    return {
        rho * vy,
        rho * vx * vy,
        rho * vy * vy + p,
        (E + p) * vy,
    };
};


// --------------------- EOS ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
vvalue_t ConservationProblem<NF, FF, EOST, SourceT>::conservativeToPrimitive(const vvalue_t& solution)
{
        const size_t nx_tot = solution.size();
        const size_t ny_tot = solution[0].size();
        vvalue_t res(nx_tot);
        for (size_t i = 0; i < nx_tot; ++i)
        {
            res[i].resize(ny_tot);
            for (size_t j = 0; j < ny_tot; ++j)
            {
                res[i][j] = conservativeToPrimitive(solution[i][j]);
            }
        }
        return res;
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::conservativeToPrimitive(const value_t& u) {
    return _conservativeToPrimitive(this->_eosT, u);
};


template <typename NF, typename FF, typename EOST, typename SourceT>
vvalue_t ConservationProblem<NF, FF, EOST, SourceT>::primitiveToConservative(const vvalue_t& primitive_solution)
{
        const size_t nx_tot = primitive_solution.size();
        const size_t ny_tot = primitive_solution[0].size();
        vvalue_t res(nx_tot);
        for (size_t i = 0; i < nx_tot; ++i)
        {
            res[i].resize(ny_tot);
            for (size_t j = 0; j < ny_tot; ++j)
            {
                res[i][j] = primitiveToConservative(primitive_solution[i][j]);
            }
        }
        return res;
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::primitiveToConservative(const value_t& q) {
    return _primitiveToConservative(this->_eosT, q);
};


// ---- EOS EULER IDEAL GAS ----
template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_conservativeToPrimitive(EOS::IdealGas, const value_t& q) {
    return {
            q[0],
            q[1] / q[0],
            q[2] / q[0],
            (this->_gamma - 1.) * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0]),
    };
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_primitiveToConservative(EOS::IdealGas, const value_t& q) {
return {
        q[0],
        q[0] * q[1],
        q[0] * q[2],
        q[3] / (this->_gamma - 1.) + 0.5 * q[0] * (q[1] * q[1] + q[2] * q[2]),
    };
};

// ---- END EOS EULER IDEAL GAS ----
// ---------------------    END EOS   ---------------------

// --------------------- Source Terms ---------------------
template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getSource(SOURCE::NullSource, const value_t& u, const point2D_t& cell_centroid) {
    value_t res(u.size());
    std::ranges::fill(res, 0.);
    return res;
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getSource(SOURCE::CylindricalGeometric, const value_t& u, const point2D_t& cell_centroid) {
    const auto &r = cell_centroid[0];
    const auto &z = cell_centroid[1];
    const auto &q = this->conservativeToPrimitive(u);

    return {-u[1] / r, -u[1] * u[1] / (u[0] * r), -u[1] * u[2] / (u[0] * r), -(u[3] + q[3]) * u[1] / (u[0] * r)};
};

// --------------------- END Source Term ---------------------


class Mesh
{
public:
    // n is the number of cells
    Mesh(int nx, int ny, int nx_ghostcell, int ny_ghostcell, double min_x, double max_x, double min_y, double max_y) :
        _nx(nx), _ny(ny), _min_x(min_x), _max_x(max_x), _min_y(min_y), _max_y(max_y),
        _dx((max_x - min_x) / (double)nx), _dy((max_y - min_y) / (double)ny),
        _nx_ghostcells(nx_ghostcell), _ny_ghostcells(ny_ghostcell)
    {
        const size_t nPointsx = nx + 1;
        const size_t nPointsy = ny + 1;

        this->_lx = (_max_x - _min_x);
        this->_ly = (_max_y - _min_y);

        this->_points.resize(nPointsx);
        for (int i = 0; i < nPointsx; ++i)
        {
            this->_points[i].resize(nPointsy);
            for (int j = 0; j < nPointsy; j++)
            {
                this->_points[i][j] = {this->_lx * i / (double)nPointsx, this->_ly * j / (double)nPointsy};
            }
        }

        this->_cellPoints.resize(nx);
        for (int i = 0; i < nx; ++i)
        {
            this->_cellPoints[i].resize(ny);
            for (int j = 0; j < ny; j++)
            {
                const point2D_t &p = this->_points[i][j];
                this->_cellPoints[i][j] = {p[0] + _dx * 0.5, p[1] + _dy * 0.5};
            }
        }

        utils::resize_vec2D(this->_fullCellPoints, nx + nx_ghostcell, ny + ny_ghostcell);

        for (size_t i = 0; i < nx; ++i)
        {
            for (int j = 0; j < ny; j++)
            {
                const point2D_t &p = this->_points[i][j];
                this->_fullCellPoints[i+nx_ghostcell][j+ny_ghostcell] = {p[0] + _dx * 0.5, p[1] + _dy * 0.5};
            }
        }
    }

    const vpoint2D_t &getPoints()
    {
        return _points;
    }

    const vpoint2D_t &getCellPoints()
    {
        return _cellPoints;
    }

    const vpoint2D_t &getFullCellPoints()
    {
        return _fullCellPoints;
    }

    const double &getDx()
    {
        return _dx;
    }

    const double &getDy()
    {
        return _dy;
    }

    const int _nx;
    const int _ny;
    const int _nx_ghostcells;
    const int _ny_ghostcells;

private:
    const double _dx;
    const double _dy;
    const double _min_x;
    const double _max_x;
    const double _min_y;
    const double _max_y;
    double _lx{};
    double _ly{};

    vpoint2D_t _points;
    vpoint2D_t _cellPoints;
    vpoint2D_t _fullCellPoints; // Include ghost cells
};

template <typename NF, typename FF, typename EOST, typename SourceT>
class Simulation {
public:
    Simulation(double final_time, ConservationProblem<NF, FF, EOST, SourceT> problem, Mesh mesh, std::function<value_t(point2D_t)> ic, std::function<void(vvalue_t&, size_t, size_t)> bc) : 
        _final_time(final_time), _problem(problem), _mesh(mesh), _ic(ic), _bc(bc), _nx(mesh._nx), _ny(mesh._ny), _nx_g(mesh._nx_ghostcells), _ny_g(mesh._ny_ghostcells)
        {
            const vpoint2D_t &mesh_points = mesh.getPoints();
            const vpoint2D_t &mesh_fullCellPoints = mesh.getFullCellPoints();
            
            utils::resize_vec2D(solution, _nx + 2 * _nx_g, _ny +  2 *_ny_g);
            utils::resize_vec2D(_solution_buffer, _nx + 2 * _nx_g, _ny +  2 * _ny_g);

            utils::transform_2d(mesh_fullCellPoints, solution, [&](point2D_t x){
                return this->_problem.primitiveToConservative(ic(x));
            });
            
            this->_bc(solution, _nx_g, _ny_g);
        };
    
    void F_sweep(double dt) {
        const double _dtdx = dt / this->_mesh.getDx();

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            const v1Dvalue_t &flux = this->_problem.getFlux(solution[i], this->_ny_g, dt, FluxDirection::F{});
            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                this->_solution_buffer[i][j] = solution[i][j] - _dtdx * (flux[j - _ny_g + 1] - flux[j - _ny_g]);
            }
        }

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                this->solution[i][j] = this->_solution_buffer[i][j];
            }
        }
    }

    void G_sweep(double dt) {
        const double _dtdy = dt / this->_mesh.getDy();

        for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
        {
            utils::ColumnView col(solution, j);
            const v1Dvalue_t &flux = this->_problem.getFlux(col.as_range(), this->_nx_g, dt, FluxDirection::G{});
            for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
            {
                this->_solution_buffer[i][j] = solution[i][j] - _dtdy * (flux[i - _nx_g + 1] - flux[i - _nx_g]);
            }
        }

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                this->solution[i][j] = this->_solution_buffer[i][j];
            }
        }

    }

    bool evolve(int i) {
        const double rest_dt = (_final_time - _ct);
        double min_dt = std::min(rest_dt, this->_problem.getMaxDt(solution));

        if(i%2 == 0) {
            this->F_sweep(min_dt);
            this->_bc(solution, _nx_g, _ny_g);
            this->G_sweep(min_dt);
            this->_bc(solution, _nx_g, _ny_g);
        } else {
            this->G_sweep(min_dt);
            this->_bc(solution, _nx_g, _ny_g);
            this->F_sweep(min_dt);
            this->_bc(solution, _nx_g, _ny_g);
        }

        _ct += min_dt;
        return _ct < _final_time;
    }

    void plot() {
        _plot(FF{}, SourceT{});
    }

    vvalue_t solution;

private:
    std::function<value_t(point2D_t)> _ic;
    std::function<void(vvalue_t&, size_t, size_t)> _bc;

    Mesh _mesh;
    ConservationProblem<NF, FF, EOST, SourceT> _problem;

    vvalue_t _solution_buffer;

    const size_t _nx;
    const size_t _ny;
    const size_t _nx_g;
    const size_t _ny_g;

    double _ct = 0.;
    const double _final_time;

    // plotters
    void _plot(FluxFunction::Euler, SOURCE::NullSource);
};

// Plotters
template <typename NF, typename FF, typename EOST, typename SourceT>
void Simulation<NF, FF, EOST, SourceT>::_plot(FluxFunction::Euler, SOURCE::NullSource) {
    vvalue_t primitive_solution = this->_problem.conservativeToPrimitive(solution);
    std::vector<float> density(_nx * _ny);
    std::vector<float> velocity_x(_nx * _ny);
    std::vector<float> velocity_y(_nx * _ny);
    std::vector<float> pressure(_nx * _ny);

    for(int i = 0; i < _nx; ++i) {
        for(int j = 0; j < _ny; ++j) {
            density[i + _nx * j]    = primitive_solution[i + _nx_g][j + _ny_g][0];
            velocity_x[i + _nx * j] = primitive_solution[i + _nx_g][j + _ny_g][1];
            velocity_y[i + _nx * j] = primitive_solution[i + _nx_g][j + _ny_g][2];
            pressure[i + _nx * j]   = primitive_solution[i + _nx_g][j + _ny_g][3];
        }
    }

    plt::clf();
    plt::suptitle(boost::lexical_cast<std::string>(this->_ct));
    const float *pVx = &(velocity_x[0]);
    const float *pVy = &(velocity_y[0]);
    const float *pD = &(density[0]);
    const float *pP = &(pressure[0]);

    plt::subplot(2, 2, 1);
    plt::title("Pressure");
    plt::imshow(pP, _ny, _nx, 1);

    plt::subplot(2, 2, 2);
    plt::title("Density");
    plt::imshow(pD, _ny, _nx, 1);

    plt::subplot(2, 2, 3);
    plt::title("Velocity X");
    plt::imshow(pVx, _ny, _nx, 1);

    plt::subplot(2, 2, 4);
    plt::title("Velocity Y");
    plt::imshow(pVy, _ny, _nx, 1);

    plt::pause(0.01);
    // plt::show();
}

// Full case description
namespace cases  {
    struct EulerForceIdealGAS {
        NumericalFlux::FORCE numerical_flux;
        FluxFunction::Euler flux_function;
        EOS::IdealGas eos;
        SOURCE::NullSource source;
        double gamma = 1.4;
        std::function<value_t(point2D_t)> ic = initialconditions::test3_2d;
        std::function<void(vvalue_t&, size_t, size_t)> bc = boundaryconditions::transmissive;
    };

    struct EulerSlicIdealGAS {
        NumericalFlux::SLIC numerical_flux;
        FluxFunction::Euler flux_function;
        EOS::IdealGas eos;
        SOURCE::NullSource source;
        double gamma = 1.4;
        std::function<value_t(point2D_t)> ic = initialconditions::test3_2d;
        std::function<void(vvalue_t&, size_t, size_t)> bc = boundaryconditions::transmissive;
    };
}

int main(int argc, char *argv[])
{
    const int nx = 200;
    const int ny = 200;
    const int nx_ghostcell = 2;
    const int ny_ghostcell = 2;
    const double min_x = 0.;
    const double max_x = 1.;
    const double min_y = 0.;
    const double max_y = 1.;
    const double dx = (max_x - min_x) / nx;
    const double dy = (max_y - min_y) / ny;


    Mesh mesh {
        nx, ny, nx_ghostcell, ny_ghostcell, min_x, max_x, min_y, max_y
    };

    const double gamma = 1.4;
    const double C = 0.7;


    // cases::EulerForceIdealGAS test_case;
    cases::EulerSlicIdealGAS test_case;


    ConservationProblem problem{
        test_case.numerical_flux, 
        test_case.flux_function, 
        test_case.eos, 
        test_case.source, 
        test_case.gamma,
        dx,
        dy,
        C
    };

    const double final_time = 10.0;

    Simulation sim(final_time, problem, mesh, test_case.ic, test_case.bc);
    
    // plt::figure_size(1000, 600);
    sim.plot();
    int i = 0;
    while (sim.evolve(i))
    {
        std::cout << "step #" << i << std::endl;
        if (i % 5 == 0)
            sim.plot();
        ++i;
    }

    return 0;
}
