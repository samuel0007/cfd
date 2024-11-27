
#include <ranges>
#include <algorithm>
#include <vector>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include <optional>
#include <set>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef std::vector<double> value_t;
typedef std::vector<value_t> v1Dvalue_t;
typedef std::vector<v1Dvalue_t> vvalue_t;

typedef std::vector<double> v1Dscalar;
typedef std::vector<std::vector<double>> vscalar;

typedef std::array<double, 2> point2D_t;
typedef std::vector<point2D_t> v1Dpoint2D_t;
typedef std::vector<v1Dpoint2D_t> vpoint2D_t;

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
    inline double normsq(double a, double b, double c) {
        return a * a + b * b + c * c;
    }

    inline double norm(double a, double b, double c) {
        return sqrt(normsq(a, b, c));
    }

    inline double dot(double a1, double b1, double c1, double a2, double b2, double c2) {
        return a1 * a2 + b1 * b2 + c1 * c2;
    }

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
    namespace euler {
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
        };

        value_t RP_test1_1D(point2D_t x) {
            return x[1] <= 0.5 ? value_t({1.0, 0.0, 0.0, 1.0}) : value_t({0.125, 0.0, 0.0, 0.1});
        }

        value_t RP_test2_1D(point2D_t x) {
            return x[1] <= 0.5 ? value_t({1.0, 0.0, -2.0, 0.4}) : value_t({1.0, 0.0, 2.0, 0.4});
        }

        value_t RP_test3_1D(point2D_t x) {
            return x[1] <= 0.5 ? value_t({1.0, 0.0, 0.0, 1000.0}) : value_t({1.0, 0.0, 0.0, 0.01});
        }

        value_t RP_test4_1D(point2D_t x) {
            return x[1] <= 0.5 ? value_t({1.0, 0.0, 0.0, 0.01}) : value_t({1.0, 0.0, 0.0, 100.0});
        }

        value_t RP_test5_1D(point2D_t x) {
            return x[1] <= 0.5 ? value_t({5.99924, 0., 19.5975, 460.894}) : value_t({5.99242, 0., -6.19633, 46.095});
        }
    }

    namespace mhd {
        // Primitive variables:
        // 0: rho, 1: vx, 2: vy, 3: vz, 4: p, 5: Bx, 6: By, 7: Bz
        // Domain: [0, 1]
        value_t toro_test_1d(point2D_t x) {
            return x[1] <= 0.5 ? value_t({1.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0}) : value_t({0.125, 0.0, 0.0, 0.0, 0.1, 0.0, 0.0, 0.0});
        }
        // Domain: [0, 800]
        value_t brio_wu_test_1d(point2D_t x) {
            return x[1] <= 400 ? value_t({1.0, 0.0, 0.0, 0.0, 1.0, 0.75, 1.0, 0.0}) : value_t({0.125, 0.0, 0.0, 0.0, 0.1, 0.75, -1.0, 0.0});
        }
    }

    namespace levelset {
        // Domain: [0, 1]
        double mid_interface(point2D_t x)
        {
            return (x[1] - 0.5);
        }

        double helium_slab_interface(point2D_t x) {
            return 0.1 - std::abs(x[1] - 0.5);
        }

        double null_interface(point2D_t x) {
            return -1;
        }
    }

    namespace multimaterial {
        value_t stationary_contact_discontinuity_L(point2D_t x) {
            return {1.0, 0., 0., 1.0};
        }
        value_t stationary_contact_discontinuity_R(point2D_t x) {
            return {0.5, 0., 0., 1.0};
        }

        value_t moving_contact_discontinuity_L(point2D_t x) {
            return {1.0, 0., 0.5, 1.0};
        }
        value_t moving_contact_discontinuity_R(point2D_t x) {
            return {0.5, 0., 0.5, 1.0};
        }

        value_t RP_test1_1D_L(point2D_t x) {
            return {1.0, 0.0, 0.0, 1.0};
        }
        value_t RP_test1_1D_R(point2D_t x) {
            return {0.125, 0.0, 0.0, 0.1};
        }

        value_t RP_test2_1D_L(point2D_t x) {
            return {1.0, 0.0, -2.0, 0.4};
        }
        value_t RP_test2_1D_R(point2D_t x) {
            return {1.0, 0.0, 2.0, 0.4};
        }

        value_t RP_test3_1D_L(point2D_t x) {
            return {1.0, 0.0, 0.0, 1000.0};
        }
        value_t RP_test3_1D_R(point2D_t x) {
            return {1.0, 0.0, 0.0, 0.01};
        }

        value_t RP_test4_1D_L(point2D_t x) {
            return {1.0, 0.0, 0.0, 0.01};
        }
        value_t RP_test4_1D_R(point2D_t x) {
            return {1.0, 0.0, 0.0, 100.0};
        }

        value_t RP_test5_1D_L(point2D_t x) {
            return {5.99924, 0., 19.5975, 460.894};
        }
        value_t RP_test5_1D_R(point2D_t x) {
            return {5.99242, 0., -6.19633, 46.095};
        }

        
        value_t Fedkiws_testA_1D_L(point2D_t x) { // gamma = 1.4
            return {1, 0., 0, 1e5};
        }
        value_t Fedkiws_testA_1D_R(point2D_t x) { // gamma = 1.2
            return {0.125, 0., 0, 1e4};
        }

        value_t Fedkiws_testB_1D_L(point2D_t x) { // gamma = 1.4
            return x[1] < 0.05 ? value_t({1.333, 0., 0.3535 * sqrt(1e5), 1.5 * 1e5}) : value_t({1, 0., 0, 1e5});
        }
        value_t Fedkiws_testB_1D_R(point2D_t x) { // gamma = 1.67
            return {0.1379, 0., 0, 1e5};
        }

        value_t Helium_slab_L(point2D_t x) { // gamma = 1.4
            return x[1] < 0.25 ? value_t({1.3765, 0., 0.3948, 1.57}) : value_t({1., 0., 0., 1.});
        }
        value_t Helium_slab_R(point2D_t x) { // gamma = 1.67
            return {0.138, 0., 0, 1};
        }
    }

}

namespace boundaryconditions {

    template <typename T>
    void transmissive(T& solution, size_t nx_g, size_t ny_g) {
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
    struct GODEXACT{};
    struct GODHLLC{};
    struct UPWIND{};
}

namespace FluxFunction {
    struct Euler{};
    struct MHD1D{};
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
    ConservationProblem(NumericalFluxT numericalFluxT, FluxFunctionT fluxFunctionT, EOST eosT, SourceT sourceT, double gamma, double dx, double dy, double C = 0.8, double eps = 1e-10, int max_it = 100):
        _numericalFluxT(numericalFluxT), _fluxFunctionT(fluxFunctionT), _eosT(eosT), _sourceT(sourceT), _dx(dx), _dy(dy), _C(C), _eps(eps), _max_it(max_it),_eos(gamma, fluxFunctionT, eosT) {};

    double getMaxDt(std::ranges::input_range auto &&solution) {
        return _getMaxDT(this->_fluxFunctionT, solution);
    };

    // expose some interface for the eos, this probably show that the IC/BC should be owned by the problem!
    vvalue_t conservativeToPrimitive(const vvalue_t& solution) {
        return this->_eos.conservativeToPrimitive(solution);
    }
    value_t  conservativeToPrimitive(const value_t&  u) {
        return this->_eos.conservativeToPrimitive(u);
    }
    vvalue_t primitiveToConservative(const vvalue_t& primitive_solution) {
        return this->_eos.primitiveToConservative(primitive_solution);
    }
    value_t  primitiveToConservative(const value_t&  q) {
        return this->_eos.primitiveToConservative(q);
    }
    double getGamma() const {
        return this->_eos.getGamma();
    }

    class eos {
    public:
        eos(double gamma, FluxFunctionT fluxFunctionT, EOST eosT): _gamma(gamma) {};
    
        double getSpeedOfSound(const value_t& u) {
            const auto q = this->conservativeToPrimitive(u);
            return _getSpeedOfSound(this->_fluxFunctionT, q);
        }

        double getSpeedOfSoundFromPrimitive(const value_t& q) {
            return _getSpeedOfSound(this->_fluxFunctionT, q);
        }

        vvalue_t conservativeToPrimitive(const vvalue_t& solution);
        value_t  conservativeToPrimitive(const value_t&  u);
        vvalue_t primitiveToConservative(const vvalue_t& primitive_solution);
        value_t  primitiveToConservative(const value_t&  q);

        double getGamma() const {
            return _gamma;
        }

    private:
        double _gamma;
        FluxFunctionT _fluxFunctionT;
        EOST _eosT;

        // Speed of sound
        double _getSpeedOfSound(FluxFunction::Euler, const value_t& q);
        double _getSpeedOfSound(FluxFunction::MHD1D, const value_t& q);

        // EOS
        value_t  _conservativeToPrimitive(EOS::IdealGas, FluxFunction::Euler, const value_t& u);
        value_t  _primitiveToConservative(EOS::IdealGas, FluxFunction::Euler, const value_t& u);
        value_t  _conservativeToPrimitive(EOS::IdealGas, FluxFunction::MHD1D, const value_t& u);
        value_t  _primitiveToConservative(EOS::IdealGas, FluxFunction::MHD1D, const value_t& u);
    };


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

    constexpr double getDelta(auto fdt) {
        if constexpr (std::is_same_v<decltype(fdt), FluxDirection::F>) {
            return _dx;
        } else {
            return _dy;
        }
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

        const double d = getDelta(fdt);
        const double _dtdx = dt / d;
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
    // Flux Fonctions
    value_t f_fromPrimitive(auto fdt, const value_t& q) {
        const value_t u = this->primitiveToConservative(q);
        return _f(this->_fluxFunctionT, fdt, u);
    }

    // Riemann Solvers
    value_t ExactRiemannSolver(auto fdt, const value_t& uL, const value_t& uR) {
        return _ExactRiemannSolver(this->_fluxFunctionT, fdt, uL, uR);
    }
    value_t HLLCFluxSolver(auto fdt, const value_t& uL, const value_t& uR) {
        return _HLLCFluxSolver(this->_fluxFunctionT, fdt, uL, uR);
    }

   
private:
    // Specialisation
    NumericalFluxT _numericalFluxT;
    FluxFunctionT _fluxFunctionT;
    EOST _eosT;
    SourceT _sourceT;

    eos _eos;

    v1Dvalue_t _flux;
    v1Dvalue_t _fluxFunctionAtCells;

    v1Dvalue_t _solution_buffer;
    v1Dvalue_t _vSlopeLimiterR;
    std::vector<double> _vSlopeLimiterPhi;

    double _dx;
    double _dy;
    double _C;

    const double _eps;
    const int _max_it;
    
    // Flux functions
    value_t _f(FluxFunction::Euler, FluxDirection::F, const value_t& u);
    value_t _f(FluxFunction::Euler, FluxDirection::G, const value_t& u);

    value_t _f(FluxFunction::MHD1D, FluxDirection::G, const value_t& u);

    // CFL conditions
    double _getMaxDT(FluxFunction::Euler, std::ranges::input_range auto &&u);
    double _getMaxDT(FluxFunction::MHD1D, std::ranges::input_range auto &&u);

    
    // Numerical Fluxes
    value_t _getNumericalFlux(NumericalFlux::BACKWARD, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::FORWARD,  auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::FORCE,    auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::CENTRAL,  auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::RI,       auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::LF,       auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::GODEXACT, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);
    value_t _getNumericalFlux(NumericalFlux::GODHLLC,  auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt);

    value_t _getNumericalFlux(NumericalFlux::LF,       auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt);
    value_t _getNumericalFlux(NumericalFlux::RI,       auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt);
    value_t _getNumericalFlux(NumericalFlux::FORCE,    auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt);

    // Source terms
    value_t _getSource(SOURCE::NullSource,           const value_t& u, const point2D_t& cell_centroid);
    value_t _getSource(SOURCE::CylindricalGeometric, const value_t& u, const point2D_t& cell_centroid);

    // Riemann Solvers

    // Exact Euler
    value_t _ExactRiemannSolver(FluxFunction::Euler, auto fdt, const value_t& uL, const value_t& uR);
    std::pair<value_t, double> _ExactRiemannSolver1D(FluxFunction::Euler, const value_t& uL1D, const value_t& uR1D);
    double newtons_raphson_pressure(double rhoL, double vL, double pL, double rhoR, double vR, double pR);
    double pressure_function(double p_star, double rhoK, double vK, double pK, double A_K, double B_K, double cs_K);
    double pressure_function_derivative(double p_star, double rhoK, double pK, double A_K, double B_K, double cs_K);
    double compute_S_K(double p_star, double pK, double vK, double cs_K_hat);
    value_t rarefaction_fan_K(double rhoK, double vK, double pK, double cs_K_hat);
    double compute_rho_star_K(double p_star, double pK, double rhoK);

    // HLLC Euler Solver
    value_t _HLLCFluxSolver(FluxFunction::Euler, auto fdt, const value_t& uL, const value_t& uR);
};

// --------------------- Speed of sound ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::eos::_getSpeedOfSound(FluxFunction::Euler, const value_t& q) {
    // cs = sqrt(gamma rho / p)
    return sqrt(this->_gamma * q[0] / q[3]);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::eos::_getSpeedOfSound(FluxFunction::MHD1D, const value_t& q) {
    return sqrt(this->_gamma * q[0] / q[4]);
}


// --------------------- CFL Condition ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::_getMaxDT(FluxFunction::Euler, std::ranges::input_range auto &&solution) {
    double max_a = 0;
    for(const auto& row: solution) {
        for (const auto &el : row)
        {
            const double a = sqrt((el[1] * el[1] + el[2] * el[2]) / (el[0]*el[0])) + this->_eos.getSpeedOfSound(el);
            max_a = std::max(max_a, a);
        }
    }
    return _C * std::min(_dx, _dy) / max_a;
}

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::_getMaxDT(FluxFunction::MHD1D, std::ranges::input_range auto &&solution) {
    double max_a = 0;
    for(const auto& row: solution) {
        for (const auto &el : row)
        {
            // max_a = v + fast wave
            const double& rho = el[0];
            const double& Bx  = el[5];
            const double& By  = el[6];
            const double& Bz  = el[7];
            const double cs = this->getSpeedOfSound(el);
            const double cs2 = cs * cs;
            const double B2 = utils::normsq(Bx, By, Bz);
            const double cf = sqrt(0.5 * sqrt(cs2 + B2 / rho + sqrt((cs2 + B2 / rho)*(cs2 + B2 / rho) - 4 * (cs2 * Bx * Bx) / rho)));
            const double a = sqrt(utils::normsq(el[1], el[2], el[3]) / (el[0]*el[0])) + cf;
            max_a = std::max(max_a, a);
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
    double d = getDelta(fdt);
    return 0.5 * (d / dt) * (u[face_idx - 1] - u[face_idx]) + 0.5 * (_fluxFunctionAtCells[face_idx - 1] + _fluxFunctionAtCells[face_idx]);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::RI, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    double d = getDelta(fdt);
    const value_t u_n = 0.5 * (u[face_idx - 1] + u[face_idx] + dt / d * (_fluxFunctionAtCells[face_idx - 1] - _fluxFunctionAtCells[face_idx]));
    return f(fdt, u_n);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::FORCE, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    return 0.5 * (this->_getNumericalFlux(NumericalFlux::LF{}, fdt, u, face_idx, dt) + this->_getNumericalFlux(NumericalFlux::RI{}, fdt, u, face_idx, dt));
}


template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::GODEXACT, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    value_t res = this->ExactRiemannSolver(fdt, u[face_idx - 1], u[face_idx]);
    return f(fdt, res);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::GODHLLC, auto fdt, std::ranges::input_range auto &&u, size_t face_idx, double dt) {
    return this->HLLCFluxSolver(fdt, u[face_idx - 1], u[face_idx]);;
}


// Need to implement this interface for 2nd order schemes
template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::LF, auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt) {
    double d = getDelta(fdt);
    return 0.5 * (d / dt) * (uL - uR) + 0.5 * (fL + fR);
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_getNumericalFlux(NumericalFlux::RI, auto fdt, const value_t& uL, const value_t& uR, const value_t& fL, const value_t& fR, double dt) {
    double d = getDelta(fdt);
    const value_t u_n = 0.5 * (uL + uR + dt / d * (fL - fR));
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
    const auto &vx  = q[1];
    const auto &vy  = q[2];
    const auto &p   = q[3];
    const auto &E   = u[3];

    return {
        rho * vx,
        rho * vx * vx + p,
        rho * vx * vy,
        (E + p) * vx,
    };
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_f(FluxFunction::Euler, FluxDirection::G, const value_t& u) {
    auto q = this->_eos.conservativeToPrimitive(u);
    const auto &rho = q[0];
    const auto &vx  = q[1];
    const auto &vy  = q[2];
    const auto &p   = q[3];
    const auto &E   = u[3];

    return {
        rho * vy,
        rho * vx * vy,
        rho * vy * vy + p,
        (E + p) * vy,
    };
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_f(FluxFunction::MHD1D, FluxDirection::G, const value_t& u) {
    auto q = this->conservativeToPrimitive(u);

    const double& rho = q[0];
    const double& vx  = q[1];
    const double& vy  = q[2];
    const double& vz  = q[3];
    const double& p   = q[4];
    const double& Bx  = q[5];
    const double& By  = q[6];
    const double& Bz  = q[7];
    const double& U   = u[4];

    const double B2 = utils::normsq(Bx, By, Bz);

    return {
        rho * vx,
        rho * vx * vx + p + 0.5 * B2 - Bx * Bx,
        rho * vx * vy - Bx * By,
        rho * vx * vz - Bx * Bz,
        (U + p + 0.5 * B2) * vx - utils::dot(vx, vy, vz, Bx, By, Bz) * Bx,
        0, // DBx/Dt = 0
        By * vx - Bx * vy,
        Bz * vx - Bx * vz,
    };
}

// --------------------- EOS ---------------------

template <typename NF, typename FF, typename EOST, typename SourceT>
vvalue_t ConservationProblem<NF, FF, EOST, SourceT>::eos::conservativeToPrimitive(const vvalue_t& solution)
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
value_t ConservationProblem<NF, FF, EOST, SourceT>::eos::conservativeToPrimitive(const value_t& u) {
    return _conservativeToPrimitive(this->_eosT, this->_fluxFunctionT, u);
};


template <typename NF, typename FF, typename EOST, typename SourceT>
vvalue_t ConservationProblem<NF, FF, EOST, SourceT>::eos::primitiveToConservative(const vvalue_t& primitive_solution)
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
value_t ConservationProblem<NF, FF, EOST, SourceT>::eos::primitiveToConservative(const value_t& q) {
    return _primitiveToConservative(this->_eosT, this->_fluxFunctionT, q);
};


// ---- EOS EULER IDEAL GAS ----
template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::eos::_conservativeToPrimitive(EOS::IdealGas, FluxFunction::Euler, const value_t& q) {
    return {
        q[0],
        q[1] / q[0],
        q[2] / q[0],
        (this->_gamma - 1.) * (q[3] - 0.5 * (q[1] * q[1] + q[2] * q[2]) / q[0]),
    };
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::eos::_primitiveToConservative(EOS::IdealGas, FluxFunction::Euler, const value_t& u) {
return {
        u[0],
        u[0] * u[1],
        u[0] * u[2],
        u[3] / (this->_gamma - 1.) + 0.5 * u[0] * (u[1] * u[1] + u[2] * u[2]),
    };
};

// ---- EOS MHD1D IDEAL GAS ----
template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::eos::_conservativeToPrimitive(EOS::IdealGas, FluxFunction::MHD1D, const value_t& u) {
    const double& rho    = u[0];
    const double& rhovx  = u[1];
    const double& rhovy  = u[2];
    const double& rhovz  = u[3];
    const double& U      = u[4];
    const double& Bx     = u[5];
    const double& By     = u[6];
    const double& Bz     = u[7];

    const double vx = rhovx / rho;
    const double vy = rhovy / rho;
    const double vz = rhovz / rho;

    return {
        rho,
        vx,
        vy,
        vz,
        (this->_gamma - 1.) * (U - 0.5 * rho * utils::normsq(vx, vy, vz) - 0.5 * utils::normsq(Bx, By, Bz)),
        Bx,
        By,
        Bz
    };
};

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::eos::_primitiveToConservative(EOS::IdealGas, FluxFunction::MHD1D, const value_t& q) {
    const double& rho = q[0];
    const double& vx  = q[1];
    const double& vy  = q[2];
    const double& vz  = q[3];
    const double& p   = q[4];
    const double& Bx  = q[5];
    const double& By  = q[6];
    const double& Bz  = q[7];

    return {
        rho,
        rho * vx,
        rho * vy,
        rho * vz,
        p / (this->_gamma - 1.) + 0.5 * rho * utils::normsq(vx, vy, vz) + 0.5 * utils::normsq(Bx, By, Bz),
        Bx,
        By,
        Bz,
    };
};


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

// --------------------- Riemann Solvers ---------------------

// ---------------- Exact Euler Solver ----------------
template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::compute_rho_star_K(double p_star, double pK, double rhoK)
{
    double gamma = this->_eos.getGamma();
    if (p_star > pK)
    { // Shock
        double r = p_star / pK;
        double gr = (gamma - 1) / (gamma + 1);
        return rhoK * (r + gr) / (gr * r + 1);
    }
    else
    { // Rarefaction
        return rhoK * std::pow(p_star / pK, (1. / gamma));
    }
}

// PRE: cs_k_hat = cs_L or -cs_R
template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::compute_S_K(double p_star, double pK, double vK, double cs_K_hat)
{
    double gamma = this->_eos.getGamma();
    return vK - cs_K_hat * sqrt((gamma + 1) * p_star / (2 * gamma * pK) + (gamma - 1) / (2 * gamma));
}

// PRE: cs_k_hat = cs_L or -cs_R
template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::rarefaction_fan_K(double rhoK, double vK, double pK, double cs_K_hat)
{
    double gamma = this->_eos.getGamma();
    const double rho = rhoK * std::pow(2. / (gamma + 1.) + vK * (gamma - 1.) / ((gamma + 1.) * cs_K_hat), 2 / (gamma - 1));
    const double v = (2 / (gamma + 1)) * (cs_K_hat + vK * (gamma - 1) / 2.);
    const double p = pK * std::pow(2. / (gamma + 1) + vK * (gamma - 1.) / ((gamma + 1) * cs_K_hat), (2. * gamma) / (gamma - 1));
    return {rho, v, p};
}

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::pressure_function(double p_star, double rhoK, double vK, double pK, double A_K, double B_K, double cs_K)
{
    double gamma = this->_eos.getGamma();
    if (p_star > pK)
    {
        // Shock
        return (p_star - pK) * sqrt(A_K / (p_star + B_K));
    }
    else
    {
        // Rarefaction
        return (2. * cs_K / (gamma - 1.)) * (std::pow(p_star / pK, (gamma - 1) / (2 * gamma)) - 1);
    }
}

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::pressure_function_derivative(double p_star, double rhoK, double pK, double A_K, double B_K, double cs_K)
    {
        double gamma = this->_eos.getGamma();
        if (p_star > pK)
        {
            // Shock
            return sqrt(A_K / (B_K + p_star)) * (1. - (p_star - pK) / (2 * (B_K + p_star)));
        }
        else
        {
            // Rarefaction
            return (1. / (rhoK * cs_K)) * std::pow(p_star / pK, -(gamma + 1) / (2 * gamma));
        }
    }

template <typename NF, typename FF, typename EOST, typename SourceT>
double ConservationProblem<NF, FF, EOST, SourceT>::newtons_raphson_pressure(double rhoL, double vL, double pL, double rhoR, double vR, double pR) {
    double gamma = this->_eos.getGamma();

    const double A_L = 2. / ((gamma + 1) * rhoL);
    const double A_R = 2. / ((gamma + 1) * rhoR);
    const double B_L = pL * (gamma - 1) / (gamma + 1);
    const double B_R = pR * (gamma - 1) / (gamma + 1);

    const double cs_L = sqrt(gamma * pL / rhoL);
    const double cs_R = sqrt(gamma * pR / rhoR);

    double p_star_old;
    // double p_star = 0.5 * (pL + pR);
    double p_star = std::max(_eps, 0.5 * (pL + pR) - 0.125 * (vR - vL) * (rhoL + rhoR) * (cs_L + cs_R));

    int i = 0;
    do
    {
        p_star_old = p_star;
        const double f_L = pressure_function(p_star_old, rhoL, vL, pL, A_L, B_L, cs_L);
        const double f_R = pressure_function(p_star_old, rhoR, vR, pR, A_R, B_R, cs_R);
        const double f = f_L + f_R + (vR - vL);
        const double f_d_L = pressure_function_derivative(p_star_old, rhoL, pL, A_L, B_L, cs_L);
        const double f_d_R = pressure_function_derivative(p_star_old, rhoR, pR, A_R, B_R, cs_R);
        const double f_d = f_d_L + f_d_R;
        p_star = p_star_old - f / f_d;
        ++i;
    } while ((fabs(p_star - p_star_old) / p_star_old) > _eps && i < this->_max_it);

    return p_star;
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_ExactRiemannSolver(FluxFunction::Euler, auto fdt, const value_t& uL, const value_t& uR) {
    int velocity_id;
    int trans_velocity_id;
    if constexpr (std::is_same_v<decltype(fdt), FluxDirection::F>) {
        velocity_id = 1;
        trans_velocity_id = 2;
    } else {
        trans_velocity_id = 1;
        velocity_id = 2;
    }

    const value_t qL = this->_eos.conservativeToPrimitive(uL);
    const value_t qR = this->_eos.conservativeToPrimitive(uR);

    const value_t qL1D{qL[0], qL[velocity_id], qL[3]};
    const value_t qR1D{qR[0], qR[velocity_id], qR[3]};
    
    // value_t q1D;
    auto [q1D, S_star] = _ExactRiemannSolver1D(FluxFunction::Euler{}, qL1D, qR1D);
    // this could be done for any advected quantity by the solver
    double trans_velocity = S_star < 0 ? qL[trans_velocity_id] : qR[trans_velocity_id];
    
    value_t primitive_res;
    if constexpr (std::is_same_v<decltype(fdt), FluxDirection::F>) {
        primitive_res = {q1D[0], q1D[1], trans_velocity, q1D[2]};
    } else {
        primitive_res = {q1D[0], trans_velocity, q1D[1], q1D[2]};
    }
    return this->_eos.primitiveToConservative(primitive_res);
};

template <typename NF, typename FF, typename EOST, typename SourceT>
std::pair<value_t, double> ConservationProblem<NF, FF, EOST, SourceT>::_ExactRiemannSolver1D(FluxFunction::Euler, const value_t& qL1D, const value_t& qR1D) {
    const double gamma = this->_eos.getGamma();
    const double &rhoL = qL1D[0];
    const double &vL   = qL1D[1];
    const double &pL   = qL1D[2];

    const double &rhoR = qR1D[0];
    const double &vR   = qR1D[1];
    const double &pR   = qR1D[2];

    const double p_star = newtons_raphson_pressure(rhoL, vL, pL, rhoR, vR, pR);

    const double A_L = 2. / ((gamma + 1) * rhoL);
    const double A_R = 2. / ((gamma + 1) * rhoR);
    const double B_L = pL * (gamma - 1) / (gamma + 1);
    const double B_R = pR * (gamma - 1) / (gamma + 1);

    const double cs_L = sqrt(gamma * pL / rhoL);
    const double cs_R = sqrt(gamma * pR / rhoR);

    const double v_star = 0.5 * (vL + vR) + 0.5 * (pressure_function(p_star, rhoR, vR, pR, A_R, B_R, cs_R) - pressure_function(p_star, rhoL, vL, pL, A_L, B_L, cs_L));

    const double rho_star_L = compute_rho_star_K(p_star, pL, rhoL);
    const double rho_star_R = compute_rho_star_K(p_star, pR, rhoR);

    const double cs_star_L = sqrt(gamma * p_star / rho_star_L);
    const double cs_star_R = sqrt(gamma * p_star / rho_star_R);

    value_t res;
    // Left Wave is right of center -> uL
    if (p_star > pL)
    { // If shock, check if right of x = 0
        double S_L = compute_S_K(p_star, pL, vL, cs_L);
        if (S_L > 0) {
            res = {rhoL, vL, pL};
            return std::make_pair(res, v_star);
        }
    }
    else
    { // If rarefaction, check if tail right of x = 0
        if ((vL - cs_L) > 0) {
            res = {rhoL, vL, pL};
            return std::make_pair(res, v_star);
        }
    }

    // Right wave is left of center -> uR
    if (p_star > pR)
    { // If shock, check if left of x = 0
        double S_R = compute_S_K(p_star, pR, vR, -cs_R);
        if (S_R < 0) {
            res = {rhoR, vR, rhoR};
            return std::make_pair(res, v_star);
        }
    }
    else
    { // If rarefaction, check if tail left of x = 0
        if ((vR + cs_R) < 0) {
            res = {rhoR, vR, pR};
            return std::make_pair(res, v_star);
        }
    }

    // Left wave is a rarefaction fan that spans over the center
    if ((p_star < pL) && ((vL - cs_L) < 0) && ((v_star - cs_star_L) > 0))
    {
        res = rarefaction_fan_K(rhoL, vL, pL, cs_L);
        return std::make_pair(res, v_star);
    }

    // Right wave is a rarefaction fan that spans over the center
    if ((p_star < pR) && ((vR + cs_R) > 0) && ((v_star + cs_star_R) < 0))
    {
        res = rarefaction_fan_K(rhoR, vR, pR, -cs_R);
        return std::make_pair(res, v_star);
    }

    // Left wave is fully on the left, right wave is fully the right (rarefaction tails included)
    // if v_star < 0 -> interface is at the right of the contact discontinuity
    // else if v_star > 0 -> interface is at the left of the contact discontinuity
    if (v_star < 0)
    {
        res = {rho_star_R, v_star, p_star};
        return std::make_pair(res, v_star);
    }
    return std::make_pair(value_t{rho_star_L, v_star, p_star}, v_star);
    // return std::make_tuple(value_t{rho_star_L, v_star, p_star}, v_star);
}

// ---------------- HLLC Euler Solver ----------------
value_t EulerMixedRPHLLCFluxSolver(auto fdt, const value_t& uL, const value_t& uR, const value_t fL, const value_t fR, auto&& eosL, auto&& eosR) {
    int velocity_id;
    int trans_velocity_id;
    if constexpr (std::is_same_v<decltype(fdt), FluxDirection::F>) {
        velocity_id = 1;
        trans_velocity_id = 2;
    } else {
        trans_velocity_id = 1;
        velocity_id = 2;
    }

    const value_t qL = eosL.conservativeToPrimitive(uL);
    const value_t qR = eosR.conservativeToPrimitive(uR);

    const double &rhoL = qL[0];
    const double &vL   = qL[velocity_id];
    const double &pL   = qL[3];

    const double &rhoR = qR[0];
    const double &vR   = qR[velocity_id];
    const double &pR   = qR[3];

    const double cs_L = eosL.getSpeedOfSoundFromPrimitive(qL);
    const double cs_R = eosR.getSpeedOfSoundFromPrimitive(qR);

    // 1. Estimate pressure
    const double p_avg = 0.5 * (pL + pR);
    const double cs_avg = 0.5 * (cs_L + cs_R);
    const double p_star = std::max(0., p_avg - 0.5 *(vR - vL) * p_avg * cs_avg);

    // 2. Estimate wave speeds
    double gammaL = eosL.getGamma();
    double gammaR = eosR.getGamma();

    auto q_k = [&](double pK, double gammaK) -> double {
        return p_star <= pK ? 1. : sqrt(1 + ((gammaK + 1) / (2 * gammaK)) * (p_star / pK - 1));
    };

    const double SL = vL - cs_L * q_k(pL, gammaL);
    // const double SL = vL - cs_L;

    const double SR = vR + cs_R * q_k(pR, gammaR);
    // const double SR = vR + cs_R;

    const double S_star = (pR - pL + rhoL * vL * (SL - vL) - rhoR * vR * (SR - vR)) / (rhoL * (SL - vL) - rhoR * (SR - vR));
    
    // 3. Compute HLLC Fluxes
    if (SL >= 0) { // All the waves are on the right of the interface, solution is left state
        return fL;
    }
    if(0 >= SR) {  // All the waves are on the left of the interface, solution is right state
        return fR;
    }

    value_t interface_state(4);

    auto prefactor_k = [&](double rhoK, double vK, double SK, double S_star) -> double {
        return rhoK * ((SK - vK) / (SK - S_star));
    };
    auto rho_star_k = [&](double rhoK, double vK, double SK, double S_star) -> double {
        return prefactor_k(rhoK, vK, SK, S_star);
    };
    auto v_star_k = [&](double rhoK, double vK, double SK, double S_star) -> double {
        return prefactor_k(rhoK, vK, SK, S_star) * S_star;
    };
    auto E_star_k = [&](double rhoK, double pK, double vK, double EK, double SK, double S_star) -> double {
        return prefactor_k(rhoK, vK, SK, S_star) * (EK / rhoK + (S_star - vK)*(S_star + pK / (rhoK * (SK - vK))));
    };
    
    if(SL <= 0 && 0 <= S_star) { // left star state
        const double EL = uL[3];
        interface_state[0] = rho_star_k(rhoL, vL, SL, S_star);
        interface_state[velocity_id] = v_star_k(rhoL, vL, SL, S_star);
        interface_state[trans_velocity_id] = prefactor_k(rhoL, vL, SL, S_star) * qL[trans_velocity_id];
        interface_state[3] = E_star_k(rhoL, pL, vL, EL, SL, S_star);
        return fL + SL * (interface_state - uL);
    } else { // right star state
        const double ER = uR[3];
        interface_state[0] = rho_star_k(rhoR, vR, SR, S_star);
        interface_state[velocity_id] = v_star_k(rhoR, vR, SR, S_star);
        interface_state[trans_velocity_id] = prefactor_k(rhoR, vR, SR, S_star) * qR[trans_velocity_id];
        interface_state[3] = E_star_k(rhoR, pR, vR, ER, SR, S_star);
        return fR + SR * (interface_state - uR);
    }
}

template <typename NF, typename FF, typename EOST, typename SourceT>
value_t ConservationProblem<NF, FF, EOST, SourceT>::_HLLCFluxSolver(FluxFunction::Euler, auto fdt, const value_t& uL, const value_t& uR) {
    value_t fL = this->f(fdt, uL);
    value_t fR = this->f(fdt, uR);
    return EulerMixedRPHLLCFluxSolver(fdt, uL, uR, fL, fR, this->_eos, this->_eos);
};


template <typename NumericalFluxT, typename FluxFunctionT1, typename FluxFunctionT2>
class TracerProblem {
public:
    TracerProblem(NumericalFluxT numericalFlux, FluxFunctionT1 fluxFunctionT1, FluxFunctionT2 fluxFunctionT2, double dx, double dy):
        _numericalFluxT(numericalFlux), _fluxFunctionT1(fluxFunctionT1), _fluxFunctionT2(fluxFunctionT2), _dx(dx), _dy(dy) {};

    // This so far doesn't consider the direction? only upwind sceme
    v1Dscalar getDDT(std::ranges::input_range auto const& tracer, std::ranges::input_range auto const& solution1, std::ranges::input_range auto const& solution2, size_t nk_g, double dt, auto fdt) {
        return this->_getDDT(this->_numericalFluxT, this->_fluxFunctionT1, this->_fluxFunctionT2, tracer, solution1, solution2, nk_g, dt, fdt);
    }

    constexpr double getDelta(auto fdt) {
        if constexpr (std::is_same_v<decltype(fdt), FluxDirection::F>) {
            return _dx;
        } else {
            return _dy;
        }
    }

private:
    NumericalFluxT _numericalFluxT;
    FluxFunctionT1 _fluxFunctionT1;
    FluxFunctionT2 _fluxFunctionT2;

    double _dx;
    double _dy;

    v1Dscalar _ddt;

    // Probably would be better to implement a get_velocity(FluxFunction) utils function
    const v1Dscalar& _getDDT(NumericalFlux::UPWIND, FluxFunction::Euler, FluxFunction::Euler, std::ranges::input_range auto const& tracer, std::ranges::input_range auto const& solution1, std::ranges::input_range auto const& solution2, size_t nk_g, double dt, auto fdt) {
        assert(tracer.size() == solution1.size() && solution1.size() == solution2.size());
        // Only implemented in 1D for y direction so far
        static_assert(std::is_same_v<decltype(fdt), FluxDirection::G>);
        const size_t nk_total = solution1.size();
        const size_t nk = nk_total - 2 * nk_g;
        _ddt.resize(nk);

        const double d = getDelta(fdt);

        for(size_t i = 0; i < nk; ++i) {
            const double& phi = tracer[i + nk_g];
            // if phi < 0, then the actual material in this cell is the solution1. Otherwise, solution2.
            const value_t& u = phi < 0 ? solution1[i + nk_g] : solution2[i + nk_g];
            const double& vy = u[2] / u[0];
            const double dphi = vy > 0 ?  phi - tracer[i + nk_g - 1] : tracer[i + nk_g + 1] - phi;
            _ddt[i] = - dphi / d * vy;
        }

        return _ddt;
    }
};

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

    // obvious TODO
    std::vector<double> getCellPoints1D()
    {
        std::vector<double> res(_ny);
        for(size_t i = 0; i < _ny; ++i) {
            res[i] = _cellPoints[0][i][1];
        }
        return res;
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


namespace GhostFluidScheme {
    struct original{};
    struct RP{};
}

// Probably the Problem should own the boundary conditions and initial conditions
template <typename NF1, typename FF1, typename EOST1, typename SourceT1,
          typename NF2, typename FF2, typename EOST2, typename SourceT2,
          typename NFL, typename GhostFluidT>
class Simulation {
public:
    Simulation(
        double final_time,
        ConservationProblem<NF1, FF1, EOST1, SourceT1> problem_1, 
        ConservationProblem<NF2, FF2, EOST2, SourceT2> problem_2, 
        Mesh mesh,
        std::function<value_t(point2D_t)> ic1,
        std::function<value_t(point2D_t)> ic2,
        std::function<void(vvalue_t&, size_t, size_t)> bc1,
        std::function<void(vvalue_t&, size_t, size_t)> bc2,
        TracerProblem<NFL, FF1, FF2> levelSet_problem,
        std::function<double(point2D_t)> levelSet_ic,
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc,
        GhostFluidT ghostFluidT
    ) : 
        _final_time(final_time), _problem_1(problem_1), _problem_2(problem_2), _mesh(mesh), _ic1(ic1), _bc1(bc1), _ic2(ic2), _bc2(bc2), _levelSet_bc(levelSet_bc), _nx(mesh._nx), _ny(mesh._ny), _nx_g(mesh._nx_ghostcells), _ny_g(mesh._ny_ghostcells), _levelSet_problem(levelSet_problem),
        _ghostFluidT(ghostFluidT)
        // _final_time(final_time), _problem(problem), _mesh(mesh), _ic(ic), _bc(bc), _nx(mesh._nx), _ny(mesh._ny), _nx_g(mesh._nx_ghostcells), _ny_g(mesh._ny_ghostcells), _levelSet_problem(levelSet_problem)
        {
            const vpoint2D_t &mesh_points = mesh.getPoints();
            const vpoint2D_t &mesh_fullCellPoints = mesh.getFullCellPoints();

            utils::resize_vec2D(solution1,         _nx + 2 * _nx_g, _ny +  2 * _ny_g);
            utils::resize_vec2D(solution2,         _nx + 2 * _nx_g, _ny +  2 * _ny_g);
            utils::resize_vec2D(_solution_buffer1, _nx + 2 * _nx_g, _ny +  2 * _ny_g);
            utils::resize_vec2D(_solution_buffer2, _nx + 2 * _nx_g, _ny +  2 * _ny_g);
            utils::resize_vec2D(levelSet,          _nx + 2 * _nx_g, _ny +  2 * _ny_g);
            utils::resize_vec2D(_levelSet_buffer,  _nx + 2 * _nx_g, _ny +  2 * _ny_g);

            utils::transform_2d(mesh_fullCellPoints, solution1, [&](point2D_t x){
                return this->_problem_1.primitiveToConservative(ic1(x));
            });
            utils::transform_2d(mesh_fullCellPoints, solution2, [&](point2D_t x){
                return this->_problem_2.primitiveToConservative(ic2(x));
            });

            utils::transform_2d(mesh_fullCellPoints, levelSet, [&](point2D_t x) {
                return levelSet_ic(x);
            });
            
            this->_levelSet_bc(levelSet, _nx_g, _ny_g);
            
            this->ghostFluidBC();
            this->_bc1(solution1, _nx_g, _ny_g);
            this->_bc2(solution2, _nx_g, _ny_g);
        };
    
    void G_sweep(double dt) {
        const double dtdy = dt / this->_mesh.getDy();

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            const v1Dvalue_t &flux1 = this->_problem_1.getFlux(solution1[i], this->_ny_g, dt, FluxDirection::G{});
            const v1Dvalue_t &flux2 = this->_problem_2.getFlux(solution2[i], this->_ny_g, dt, FluxDirection::G{});

            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                if(this->levelSet[i][j] < 0) {
                    this->_solution_buffer1[i][j] = solution1[i][j] - dtdy * (flux1[j - _ny_g + 1] - flux1[j - _ny_g]);
                } else {
                    this->_solution_buffer2[i][j] = solution2[i][j] - dtdy * (flux2[j - _ny_g + 1] - flux2[j - _ny_g]);
                }
            }
        }

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                if(this->levelSet[i][j] < 0) {
                    this->solution1[i][j] = this->_solution_buffer1[i][j];
                } else {
                    this->solution2[i][j] = this->_solution_buffer2[i][j];
                }
            }
        }
    }

    void F_sweep(double dt) {
        const double dtdx = dt / this->_mesh.getDx();

        for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
        {
            utils::ColumnView col1(solution1, j);
            utils::ColumnView col2(solution2, j);
            const v1Dvalue_t &flux1 = this->_problem_1.getFlux(col1.as_range(), this->_nx_g, dt, FluxDirection::F{});
            const v1Dvalue_t &flux2 = this->_problem_2.getFlux(col2.as_range(), this->_nx_g, dt, FluxDirection::F{});
            for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
            {
                this->_solution_buffer1[i][j] = solution1[i][j] - dtdx * (flux1[i - _nx_g + 1] - flux1[i - _nx_g]);
                this->_solution_buffer2[i][j] = solution2[i][j] - dtdx * (flux2[i - _nx_g + 1] - flux2[i - _nx_g]);
            }
        }

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                this->solution1[i][j] = this->_solution_buffer1[i][j];
                this->solution2[i][j] = this->_solution_buffer2[i][j];
            }
        }
    }

    void evolveLevelSet1D(double dt) {
        const double dtdy = dt / this->_mesh.getDy();

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            const v1Dscalar &levelSet_ddt = this->_levelSet_problem.getDDT(levelSet[i], solution1[i], solution2[i], this->_ny_g, dt, FluxDirection::G{});
            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                this->_levelSet_buffer[i][j] = levelSet[i][j] + dt * levelSet_ddt[j - _ny_g];
            }
        }

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            for (size_t j = _ny_g; j < this->_ny + _ny_g; ++j)
            {
                this->levelSet[i][j] = this->_levelSet_buffer[i][j];
            }
        }
    }

    void reinitLevelSet1D() {
        // Find interface
        const double& dy = this->_mesh.getDy();

        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            size_t left_idx = _ny_g;
            double left_value = this->levelSet[i][left_idx];
            size_t right_idx;
            double right_value;

            for (size_t j = _ny_g + 1; j < this->_ny + _ny_g; ++j)
            {
                right_idx = j;
                right_value = this->levelSet[i][right_idx];
                if(left_value * right_value < 0) {
                    break;
                } // at the interface, both values have different signs
                left_value = right_value;
                left_idx = right_idx;
            }
            // const double sLdy = left_value - this->levelSet[i][left_idx - 1];
            // const double sRdy = this->levelSet[i][right_idx + 1] - right_value;

            for(size_t j = _ny_g; j < _ny + _ny_g; ++j) {
                if(j < left_idx) {
                    levelSet[i][j] = (left_idx - j)  * (-dy);
                } else if(j > right_idx) {
                    levelSet[i][j] = (j - right_idx) * dy;
                }
            }

        }
    }
    
    void ghostFluidBC() {
        if constexpr(std::is_same_v<decltype(this->_ghostFluidT), GhostFluidScheme::original>) {
            return originalGhostFluidBC();
        } else if(std::is_same_v<decltype(this->_ghostFluidT), GhostFluidScheme::RP>) {
            return RPGhostFluidBC();
        } else {
            assert(false && "missing implementation for selected ghost fluid scheme");
        }
    }
    // this only work for a euler-euler system
    void originalGhostFluidBC() {
        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            std::map<size_t, value_t> interfaces1_idx;
            std::map<size_t, value_t> interfaces2_idx;
            // Copy velocity pressure and find intervals
            for(size_t j = _ny_g + 1; j < _ny + _ny_g; ++j) {
                const double& left_phi = this->levelSet[i][j - 1];
                const double& right_phi = this->levelSet[i][j];
                if(left_phi * right_phi < 0) {
                    // we are at an interface
                    if(left_phi < 0) {
                        interfaces1_idx[j - 1] = this->_problem_1.conservativeToPrimitive(solution1[i][j - 1]);
                        interfaces2_idx[j]     = this->_problem_2.conservativeToPrimitive(solution2[i][j]);
                    } else {
                        interfaces1_idx[j]     = this->_problem_1.conservativeToPrimitive(solution1[i][j]);
                        interfaces2_idx[j - 1] = this->_problem_2.conservativeToPrimitive(solution2[i][j - 1]);
                    }
                }
            }

            // If they are no interfaces, we don't have to do any ghost fluid method
            if(interfaces1_idx.empty()) return;

            auto get_closest_interface_primitive_value = [&](size_t idx, std::map<size_t, value_t> interfaces) -> value_t {
                auto closest_it = std::ranges::min_element(interfaces, 
                [idx](const auto& a, const auto& b) {
                    return std::abs(static_cast<long>(a.first) - static_cast<long>(idx)) <
                        std::abs(static_cast<long>(b.first) - static_cast<long>(idx));
                });

                return closest_it->second;
            };

            // Constant Entropy Extrapolation:
            // If this is a ghost region, find closest interface and extrapolate
            // rho_g = ( P_g / P_i ) ^ (1/gamma_i)  * rho_i
            for(size_t j = _ny_g; j < _ny + _ny_g; ++j) {
                if(this->levelSet[i][j] < 0) {
                    // Real: 1, Ghost: 2
                    // Only velocity and pressure are correct in q_g
                    value_t q_g = this->_problem_1.conservativeToPrimitive(solution1[i][j]);

                    // Find q_i
                    value_t q_i = get_closest_interface_primitive_value(j, interfaces2_idx);

                    // Modify density
                    double gamma = this->_problem_2.getGamma();
                    q_g[0] = std::pow(q_g[3] / q_i[3], 1. / gamma) * q_i[0];
                    solution2[i][j] = this->_problem_2.primitiveToConservative(q_g);
                } else {
                    // Real: 2, Ghost: 1
                    value_t q_g = this->_problem_2.conservativeToPrimitive(solution2[i][j]);
                    value_t q_i = get_closest_interface_primitive_value(j, interfaces1_idx);

                    double gamma = this->_problem_1.getGamma();
                    q_g[0] = std::pow(q_g[3] / q_i[3], 1. / gamma) * q_i[0];
                    solution1[i][j] = this->_problem_1.primitiveToConservative(q_g);
                }
            }
        }
    }

    void RPGhostFluidBC() {
        for (size_t i = _nx_g; i < this->_nx + _nx_g; ++i)
        {
            std::map<size_t, value_t> interfaces1_idx;
            std::map<size_t, value_t> interfaces2_idx;
            // Copy velocity pressure and find intervals
            for(size_t j = _ny_g + 1; j < _ny + _ny_g; ++j) {
                const double& left_phi = this->levelSet[i][j - 1];
                const double& right_phi = this->levelSet[i][j];
                if(left_phi * right_phi < 0) {
                    // we are at an interface
                    if(left_phi < 0) {
                        interfaces1_idx[j - 1] = this->_problem_1.conservativeToPrimitive(solution1[i][j - 1]);
                        interfaces2_idx[j]     = this->_problem_2.conservativeToPrimitive(solution2[i][j]);
                    } else {
                        interfaces1_idx[j]     = this->_problem_1.conservativeToPrimitive(solution1[i][j]);
                        interfaces2_idx[j - 1] = this->_problem_2.conservativeToPrimitive(solution2[i][j - 1]);
                    }
                }
            }

            // If they are no interfaces, we don't have to do any ghost fluid method
            if(interfaces1_idx.empty()) return;

            auto get_closest_interface_primitive_value = [&](size_t idx, std::map<size_t, value_t> interfaces) -> value_t {
                auto closest_it = std::ranges::min_element(interfaces, 
                [idx](const auto& a, const auto& b) {
                    return std::abs(static_cast<long>(a.first) - static_cast<long>(idx)) <
                        std::abs(static_cast<long>(b.first) - static_cast<long>(idx));
                });

                return closest_it->second;
            };

            // Constant Entropy Extrapolation:
            // If this is a ghost region, find closest interface and extrapolate
            // rho_g = ( P_g / P_i ) ^ (1/gamma_i)  * rho_i
            for(size_t j = _ny_g; j < _ny + _ny_g; ++j) {
                if(this->levelSet[i][j] < 0) {
                    // Real: 1, Ghost: 2
                    // Only velocity and pressure are correct in q_g
                    value_t q_g = this->_problem_1.conservativeToPrimitive(solution1[i][j]);

                    // Find q_i
                    value_t q_i = get_closest_interface_primitive_value(j, interfaces2_idx);

                    // Modify density
                    double gamma = this->_problem_2.getGamma();
                    q_g[0] = std::pow(q_g[3] / q_i[3], 1. / gamma) * q_i[0];
                    solution2[i][j] = this->_problem_2.primitiveToConservative(q_g);
                } else {
                    // Real: 2, Ghost: 1
                    value_t q_g = this->_problem_2.conservativeToPrimitive(solution2[i][j]);
                    value_t q_i = get_closest_interface_primitive_value(j, interfaces1_idx);

                    double gamma = this->_problem_1.getGamma();
                    q_g[0] = std::pow(q_g[3] / q_i[3], 1. / gamma) * q_i[0];
                    solution1[i][j] = this->_problem_1.primitiveToConservative(q_g);
                }
            }

        }

    }

    bool evolve(int i) {
        const double rest_dt = (_final_time - _ct);
        double min_dt = std::min(rest_dt, this->_problem_1.getMaxDt(solution1));
        // min_dt        = std::min(min_dt,  this->_problem_2.getMaxDt(solution2));

        if(i%2 == 0) {
            this->G_sweep(min_dt);
            this->_bc1(solution1, _nx_g, _ny_g);
            this->F_sweep(min_dt);
            this->_bc1(solution1, _nx_g, _ny_g);
        } else {
            this->F_sweep(min_dt);
            this->_bc1(solution1, _nx_g, _ny_g);
            this->G_sweep(min_dt);
            this->_bc1(solution1, _nx_g, _ny_g);
        }

        _ct += min_dt;
        return _ct < _final_time;
    }

    bool evolve1D(int i) {
        const double rest_dt = (_final_time - _ct);
        double min_dt = std::min(rest_dt, this->_problem_1.getMaxDt(solution1));
        min_dt        = std::min(min_dt,  this->_problem_2.getMaxDt(solution2));
        // min_dt        = std::min(1e-4, min_dt);

        std::cout << min_dt << "\n";

        // Level set update
        this->evolveLevelSet1D(min_dt);
        this->_levelSet_bc(levelSet, _nx_g, _ny_g);
        // this->reinitLevelSet1D();

        this->G_sweep(min_dt);

        // Ghost Fluid Boundaries
        this->ghostFluidBC();
        // Physical Boundaries
        this->_bc1(solution1, _nx_g, _ny_g);
        this->_bc2(solution2, _nx_g, _ny_g);
            
        _ct += min_dt;
        return _ct < _final_time;
    }

    // TODO plot for two materials
    void plot() {
        _plot(FF1{}, SourceT1{});
    }

    // TODO
    void plot1D() {
        _plot1D(FF1{}, SourceT1{});
    }

    vvalue_t solution1;
    vvalue_t solution2;

    vscalar levelSet;

private:
    const GhostFluidT _ghostFluidT;
    std::function<value_t(point2D_t)> _ic1;
    std::function<value_t(point2D_t)> _ic2;
    std::function<void(vvalue_t&, size_t, size_t)> _bc1;
    std::function<void(vvalue_t&, size_t, size_t)> _bc2;
    std::function<void(vscalar&, size_t, size_t)> _levelSet_bc;

    Mesh _mesh;
    ConservationProblem<NF1, FF1, EOST1, SourceT1> _problem_1;
    ConservationProblem<NF2, FF2, EOST2, SourceT2> _problem_2;

    TracerProblem<NFL, FF1, FF2> _levelSet_problem;

    vvalue_t _solution_buffer1;
    vvalue_t _solution_buffer2;

    vscalar _levelSet_buffer;

    const size_t _nx;
    const size_t _ny;
    const size_t _nx_g;
    const size_t _ny_g;

    double _ct = 0.;
    const double _final_time;

    // plotters
    void _plot(FluxFunction::Euler, SOURCE::NullSource);
    void _plot1D(FluxFunction::Euler, SOURCE::NullSource);
    void _plot1D(FluxFunction::MHD1D, SOURCE::NullSource);
};

// Plotters
template <typename NF1, typename FF1, typename EOST1, typename SourceT1,
          typename NF2, typename FF2, typename EOST2, typename SourceT2,
          typename NFL, typename GFS>
void Simulation<NF1, FF1, EOST1, SourceT1, NF2, FF2, EOST2, SourceT2, NFL, GFS>::_plot(FluxFunction::Euler, SOURCE::NullSource) {
    vvalue_t primitive_solution = this->_problem_1.conservativeToPrimitive(solution1);
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

template <typename NF1, typename FF1, typename EOST1, typename SourceT1,
          typename NF2, typename FF2, typename EOST2, typename SourceT2,
          typename NFL, typename GFS>
void Simulation<NF1, FF1, EOST1, SourceT1, NF2, FF2, EOST2, SourceT2, NFL, GFS>::_plot1D(FluxFunction::Euler, SOURCE::NullSource) {
    vvalue_t primitive_solution1 = this->_problem_1.conservativeToPrimitive(solution1);
    vvalue_t primitive_solution2 = this->_problem_2.conservativeToPrimitive(solution2);

    assert(_nx_g == 0 && _nx == 1);
    std::vector<double> density1(_ny);
    std::vector<double> velocity1(_ny);
    std::vector<double> pressure1(_ny);

    std::vector<double> density2(_ny);
    std::vector<double> velocity2(_ny);
    std::vector<double> pressure2(_ny);

    std::vector<double> gamma(_ny);

    for(int j = 0; j < _ny; ++j) {
        gamma[j] = levelSet[0][_ny_g + j];
        
        density1[j]  = gamma[j] < 0 ? primitive_solution1[0][_ny_g + j][0]: INFINITY;
        velocity1[j] = gamma[j] < 0 ? primitive_solution1[0][_ny_g + j][2]: INFINITY;
        pressure1[j] = gamma[j] < 0 ? primitive_solution1[0][_ny_g + j][3]: INFINITY;

        density2[j]  = gamma[j] >= 0 ? primitive_solution2[0][_ny_g + j][0]: INFINITY;
        velocity2[j] = gamma[j] >= 0 ? primitive_solution2[0][_ny_g + j][2]: INFINITY;
        pressure2[j] = gamma[j] >= 0 ? primitive_solution2[0][_ny_g + j][3]: INFINITY;
    }

    auto cell_points = _mesh.getCellPoints1D();
    
    // int ny_plots = this->_has_levelset ? 3 : 2;
    int ny_plots = 3;
    plt::clf();
    plt::suptitle(boost::lexical_cast<std::string>(this->_ct));
    plt::subplot(2, ny_plots, 1);
    plt::title("Pressure");
    // plt::ylim(0, 2);
    // plt::plot(cell_points, pressure1, "r*");
    plt::plot(cell_points, pressure1, "r");
    plt::plot(cell_points, pressure2, "b");

    // plt::imshow(pP, _ny, _nx, 1);

    plt::subplot(2, ny_plots, 2);
    plt::title("Density");
    // plt::ylim(0, 2);
    // plt::plot(cell_points, density1, "r*");
    plt::plot(cell_points, density1, "r");
    plt::plot(cell_points, density2, "b");

    plt::subplot(2, ny_plots, 3);
    plt::title("Velocity");
    // plt::ylim(-2, 2);

    // plt::plot(cell_points, velocity1, "r*");
    plt::plot(cell_points, velocity1, "r");
    plt::plot(cell_points, velocity2, "b");

    // if(this->_has_levelset) {

    plt::subplot(2, ny_plots, 4);
    plt::title("Level Set");
    // plt::plot(cell_points, gamma, "r*");
    plt::ylim(-1, 1);
    plt::plot(cell_points, gamma);
    // }

    plt::pause(0.01);
    // plt::show();

}

template <typename NF1, typename FF1, typename EOST1, typename SourceT1,
          typename NF2, typename FF2, typename EOST2, typename SourceT2,
          typename NFL, typename GFS>
void Simulation<NF1, FF1, EOST1, SourceT1, NF2, FF2, EOST2, SourceT2, NFL, GFS>::_plot1D(FluxFunction::MHD1D, SOURCE::NullSource) {
    vvalue_t primitive_solution1 = this->_problem_1.conservativeToPrimitive(solution1);

    assert(_nx_g == 0 && _nx == 1);
    std::vector<double> density(_ny);
    std::vector<double> velocity_x(_ny);
    std::vector<double> velocity_y(_ny);
    std::vector<double> velocity_z(_ny);
    std::vector<double> v(_ny);
    std::vector<double> pressure(_ny);
    std::vector<double> By(_ny);
    std::vector<double> Bz(_ny);
    std::vector<double> E(_ny);


    for(int j = 0; j < _ny; ++j) {
        density[j]    = primitive_solution1[0][_ny_g + j][0];
        velocity_x[j] = primitive_solution1[0][_ny_g + j][1];
        velocity_y[j] = primitive_solution1[0][_ny_g + j][2];
        velocity_z[j] = primitive_solution1[0][_ny_g + j][3];
        pressure[j]   = primitive_solution1[0][_ny_g + j][4];
        By[j]         = primitive_solution1[0][_ny_g + j][6];
        Bz[j]         = primitive_solution1[0][_ny_g + j][7];

        v[j] = utils::norm(velocity_x[j], velocity_y[j], velocity_z[j]);
        E[j] = solution1[0][_ny_g + j][4];
    }

    auto cell_points = _mesh.getCellPoints1D();
  
    plt::clf();
    plt::suptitle(boost::lexical_cast<std::string>(this->_ct));
    plt::subplot(3, 3, 1);
    plt::title("Pressure");
    // plt::plot(cell_points, pressure, "r*");
    plt::plot(cell_points, pressure);

    // plt::imshow(pP, _ny, _nx, 1);

    plt::subplot(3, 3, 2);
    plt::title("Density");
    // plt::plot(cell_points, density, "r*");
    plt::plot(cell_points, density);

    plt::subplot(3, 3, 3);
    plt::title("Energy");
    // plt::plot(cell_points, E, "r*");
    plt::plot(cell_points, E);

    plt::subplot(3, 3, 4);
    plt::title("Velocity X");
    // plt::plot(cell_points, velocity_x, "r*");
    plt::plot(cell_points, velocity_x);

    plt::subplot(3, 3, 5);
    plt::title("Velocity Y");
    // plt::plot(cell_points, velocity_y, "r*");
    plt::plot(cell_points, velocity_y);

    plt::subplot(3, 3, 6);
    plt::title("Velocity Z");
    // plt::plot(cell_points, velocity_z, "r*");
    plt::plot(cell_points, velocity_z);

    plt::subplot(3, 3, 7);
    plt::title("V");
    // plt::plot(cell_points, v, "r*");
    plt::plot(cell_points, v);

    plt::subplot(3, 3, 8);
    plt::title("By");
    // plt::plot(cell_points, By, "r*");
    plt::plot(cell_points, By);

    plt::subplot(3, 3, 9);
    plt::title("Bz");
    // plt::plot(cell_points, Bz, "r*");
    plt::plot(cell_points, Bz);

    plt::pause(0.01);
}

// Full case description
namespace cases  {
    struct EulerForceIdealGAS {
        NumericalFlux::FORCE numerical_flux;
        FluxFunction::Euler flux_function;
        EOS::IdealGas eos;
        SOURCE::NullSource source;
        double gamma = 1.4;
        std::function<value_t(point2D_t)> ic = initialconditions::euler::RP_test1_1D;
        std::function<void(vvalue_t&, size_t, size_t)> bc = boundaryconditions::transmissive<vvalue_t>;
    };

    struct EulerSlicIdealGAS {
        NumericalFlux::SLIC numerical_flux;
        FluxFunction::Euler flux_function;
        EOS::IdealGas eos;
        SOURCE::NullSource source;
        double gamma = 1.4;
        std::function<value_t(point2D_t)> ic = initialconditions::euler::RP_test1_1D;
        std::function<void(vvalue_t&, size_t, size_t)> bc = boundaryconditions::transmissive<vvalue_t>;
    };

    struct EulerGODIdealGAS {
        NumericalFlux::FORCE numerical_flux;
        FluxFunction::Euler flux_function;
        EOS::IdealGas eos;
        SOURCE::NullSource source;
        double gamma = 1.4;
        std::function<value_t(point2D_t)> ic = initialconditions::euler::test3_2d;
        std::function<void(vvalue_t&, size_t, size_t)> bc = boundaryconditions::transmissive<vvalue_t>;
    };

    struct MHD1DForceIdealGas {
        NumericalFlux::FORCE numerical_flux;
        FluxFunction::MHD1D flux_function;
        EOS::IdealGas eos;
        SOURCE::NullSource source;
        double gamma = 2;
        std::function<value_t(point2D_t)> ic = initialconditions::mhd::brio_wu_test_1d;
        std::function<void(vvalue_t&, size_t, size_t)> bc = boundaryconditions::transmissive<vvalue_t>;
    };

    struct Euler1DGODIdealGasLevelSet {
        NumericalFlux::GODEXACT numerical_flux1;
        NumericalFlux::GODEXACT numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.4;
        std::function<value_t(point2D_t)> ic1 = initialconditions::euler::RP_test1_1D;
        std::function<value_t(point2D_t)> ic2 = initialconditions::euler::RP_test1_1D;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::mid_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;
    };

    struct StationaryContactDiscontinuityMM {
        NumericalFlux::GODEXACT numerical_flux1;
        NumericalFlux::GODEXACT numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.67;
        std::function<value_t(point2D_t)> ic1 = initialconditions::multimaterial::stationary_contact_discontinuity_L;
        std::function<value_t(point2D_t)> ic2 = initialconditions::multimaterial::stationary_contact_discontinuity_R;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::mid_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;
    };

    struct StationaryContactDiscontinuityMMHLLC {
        NumericalFlux::GODHLLC numerical_flux1;
        NumericalFlux::GODHLLC numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.67;
        std::function<value_t(point2D_t)> ic1 = initialconditions::multimaterial::stationary_contact_discontinuity_L;
        std::function<value_t(point2D_t)> ic2 = initialconditions::multimaterial::stationary_contact_discontinuity_R;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::mid_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;
        GhostFluidScheme::original ghostfluidscheme;
    };

    struct MovingContactDiscontinuityMM {
        NumericalFlux::GODEXACT numerical_flux1;
        NumericalFlux::GODEXACT numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.67;
        std::function<value_t(point2D_t)> ic1 = initialconditions::multimaterial::moving_contact_discontinuity_L;
        std::function<value_t(point2D_t)> ic2 = initialconditions::multimaterial::moving_contact_discontinuity_R;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::mid_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;
        GhostFluidScheme::original ghostfluidscheme;
    };

    struct RPGODMM {
        NumericalFlux::GODEXACT numerical_flux1;
        NumericalFlux::GODEXACT numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.67;
        std::function<value_t(point2D_t)> ic1 = initialconditions::multimaterial::Fedkiws_testB_1D_L;
        std::function<value_t(point2D_t)> ic2 = initialconditions::multimaterial::Fedkiws_testB_1D_R;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::mid_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;
        GhostFluidScheme::original ghostfluidscheme;
    };

    struct RPHLLCMM {
        NumericalFlux::GODHLLC numerical_flux1;
        NumericalFlux::GODHLLC numerical_flux2;
        // NumericalFlux::GODEXACT numerical_flux1;
        // NumericalFlux::GODEXACT numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.4;
        std::function<value_t(point2D_t)> ic1 = initialconditions::euler::RP_test1_1D;
        std::function<value_t(point2D_t)> ic2 = initialconditions::euler::RP_test1_1D;;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::null_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;

        GhostFluidScheme::original ghostfluidscheme;
    };

    struct HeliumSlabExact {
        NumericalFlux::GODEXACT numerical_flux1;
        NumericalFlux::GODEXACT numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.67;
        std::function<value_t(point2D_t)> ic1 = initialconditions::multimaterial::Helium_slab_L;
        std::function<value_t(point2D_t)> ic2 = initialconditions::multimaterial::Helium_slab_R;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::helium_slab_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;
        GhostFluidScheme::original ghostfluidscheme;
    };

    struct HeliumSlabHLLC {
        NumericalFlux::GODHLLC numerical_flux1;
        NumericalFlux::GODHLLC numerical_flux2;
        FluxFunction::Euler flux_function1;
        FluxFunction::Euler flux_function2;
        EOS::IdealGas eos1;
        EOS::IdealGas eos2;

        SOURCE::NullSource source1;
        SOURCE::NullSource source2;
        double gamma1 = 1.4;
        double gamma2 = 1.67;
        std::function<value_t(point2D_t)> ic1 = initialconditions::multimaterial::Helium_slab_L;
        std::function<value_t(point2D_t)> ic2 = initialconditions::multimaterial::Helium_slab_R;

        std::function<void(vvalue_t&, size_t, size_t)> bc1 = boundaryconditions::transmissive<vvalue_t>;
        std::function<void(vvalue_t&, size_t, size_t)> bc2 = boundaryconditions::transmissive<vvalue_t>;

        NumericalFlux::UPWIND levelSet_numericalFlux;
        std::function<double(point2D_t)> levelSet_ic = initialconditions::levelset::helium_slab_interface;
        std::function<void(vscalar&, size_t, size_t)> levelSet_bc = boundaryconditions::transmissive<vscalar>;
        GhostFluidScheme::original ghostfluidscheme;
    };
}

int main(int argc, char *argv[])
{
    const double C = 0.1;
    const double final_time = 100.0;

    const int nx = 1;
    const int nx_ghostcell = 0;
    const int ny = 1000;
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


    // cases::EulerForceIdealGAS test_case;
    // cases::EulerSlicIdealGAS test_case;
    // cases::EulerGODIdealGAS test_case;
    // cases::MHD1DForceIdealGas test_case;
    // cases::Euler1DGODIdealGasLevelSet test_case;
    // cases::StationaryContactDiscontinuityMM test_case;
    // cases::MovingContactDiscontinuityMM test_case;
    // cases::RPGODMM test_case;
    // cases::HeliumSlabExact test_case;
    // cases::StationaryContactDiscontinuityMMHLLC test_case;
    cases::RPHLLCMM test_case;
    // cases::HeliumSlabHLLC test_case;

    ConservationProblem problem_1{
        test_case.numerical_flux1, 
        test_case.flux_function1, 
        test_case.eos1, 
        test_case.source1, 
        test_case.gamma1,
        dx,
        dy,
    };

    ConservationProblem problem_2{
        test_case.numerical_flux2, 
        test_case.flux_function2, 
        test_case.eos2, 
        test_case.source2, 
        test_case.gamma2,
        dx,
        dy,
    };

    TracerProblem levelSet_problem{
        test_case.levelSet_numericalFlux,
        test_case.flux_function1,
        test_case.flux_function2,
        dx,
        dy
    };


    Simulation sim(
        final_time,
        problem_1,
        problem_2,
        mesh,
        test_case.ic1,
        test_case.ic2,
        test_case.bc1,
        test_case.bc2,
        levelSet_problem,
        test_case.levelSet_ic,
        test_case.levelSet_bc,
        test_case.ghostfluidscheme
    );

    
    // plt::figure_size(1000, 600);
    sim.plot1D();
    // sim.plot();

    int i = 0;
    // while (sim.evolve(i))
    while (sim.evolve1D(i))
    {
        std::cout << "step #" << i << std::endl;
        if (i % 30 == 0) {
            sim.plot1D();
            // sim.plot();
        }
        ++i;
    }

    return 0;
}
