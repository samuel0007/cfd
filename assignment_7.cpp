//
// Created by russo on 14.10.2024.
//

#include <iostream>
#include <utility>
#include <vector>
#include <ranges>
#include <algorithm>
#include <boost/lexical_cast.hpp>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef std::vector<double> value_t;
typedef std::vector<std::vector<value_t>> vvalue_t;
typedef std::vector<value_t> v1Dvalue_t;


typedef std::array<double, 2> point2D_t;
typedef std::vector<std::vector<point2D_t>> vpoint2D_t;


class ColumnView {
public:
    ColumnView(vvalue_t& matrix, size_t col)
        : matrix(matrix), col(col) {}

    auto as_range() {
        return std::ranges::views::transform(matrix, [this](auto& row) { return row[col]; });
    }

private:
    vvalue_t& matrix;
    size_t col;
};

std::vector<double> operator+(const std::vector<double>& v1, const std::vector<double>& v2)
{
    assert(v1.size() == v2.size());
    std::vector<double> res(v1.size());
    std::ranges::transform(v1, v2, res.begin(), std::plus<double>{});
    return res;
}

std::vector<double> operator-(const std::vector<double>& v1, const std::vector<double>& v2)
{
    assert(v1.size() == v2.size());
    std::vector<double> res(v1.size());
    std::ranges::transform(v1, v2, res.begin(), std::minus<double>{});
    return res;
}

std::vector<double> operator*(const std::vector<double>& v1, double alpha) {
    std::vector<double> res(v1.size());
    std::ranges::transform(v1, res.begin(), [alpha](double el){return alpha*el;});
    return res;
}

std::vector<double> operator*(double alpha, const std::vector<double>& v) {
    return v * alpha;
}


namespace utils
{

    value_t primitiveToConservative2D(double gamma, value_t u){
        return {
            u[0],
            u[0]*u[1],
            u[0]*u[2],
            u[3] / (gamma - 1.) + 0.5 * u[0] * (u[1] * u[1] + u[2] * u[2]),
        };
    }

    value_t conservativeToPrimitive(double gamma, value_t u) {
        return {
            u[0],
            u[1] / u[0],
            u[2] / u[0],
            (gamma - 1.) * (u[3] - 0.5 * (u[1] * u[1] + u[2] * u[2]) / u[0]),
        };
    }

    vvalue_t conservativeToPrimitive(double gamma, vvalue_t solution) {
        const size_t nx = solution.size();
        const size_t ny = solution[0].size();
        vvalue_t res(nx);
        for(size_t i = 0; i < nx; ++i) {
            res[i].resize(ny);
            for(size_t j = 0; j < ny; ++j) {
                res[i][j] = conservativeToPrimitive(gamma, solution[i][j]);
            }
        }
        return res;
    }

    value_t test0_2d(point2D_t x) {
        return x[0] <= 0.25 ? primitiveToConservative2D(1.4, {1.0, 1.0, 0.0, 1.0}) : primitiveToConservative2D(1.4, {0.1, 1.0, 0.0, 1.});
    }

    value_t test2_2d(point2D_t x) {
        return x[0] <= 0.5 ? primitiveToConservative2D(1.4, {1.0, -2.0, 0.5, 0.4}) : primitiveToConservative2D(1.4, {1.0, 2.0, 0.5, 0.4});
    }

    value_t test3_2d(point2D_t x) {
        double x_0 = 1.;
        double y_0 = 1.;
        double R = 0.4;
        if((x[0] - x_0) * (x[0] - x_0) + (x[1] - y_0) * (x[1] - y_0) < R*R) {
            return {1, 0, 0, 1};
        } else {
            return {0.125, 0, 0, 0.1};
        }
    }

    value_t test0(point2D_t x) {
        return x[0] <= 0.25 ? primitiveToConservative2D(1.4, {1.0, 1.0, 1.0}) : primitiveToConservative2D(1.4, {0.1, 1.0, 1.});
    }

    value_t test1(point2D_t x) {
        return x[0] <= 0.5 ? primitiveToConservative2D(1.4, {1.0, 0.0, 1.0}) : primitiveToConservative2D(1.4, {0.125, 0., 0.1});
    }

    value_t test2(double x) {
        return x <= 0.5 ? primitiveToConservative2D(1.4, {1.0, -2.0, 0.4}) : primitiveToConservative2D(1.4, {1.0, 2.0, 0.4});
    }

    value_t test3(double x) {
        return x <= 0.5 ? primitiveToConservative2D(1.4, {1.0, 0.0, 1000.}) : primitiveToConservative2D(1.4, {1.0, 0., 0.01});
    }

    value_t test4(double x) {
        return x <= 0.5 ? primitiveToConservative2D(1.4, {1.0, 0.0, 0.01}) : primitiveToConservative2D(1.4, {1.0, 0., 100.});
    }

    value_t test5(double x) {
        return x <= 0.5 ? primitiveToConservative2D(1.4, {5.99924, 19.5975, 460.894}) : primitiveToConservative2D(1.4, {5.99242, -6.19633, 46.095});
    }
}


// enum NumericalFlux
// {
//     Central,
//     Backward,
//     Forward,
//     LF,
//     RI,
//     FORCE,
//     GOD,
// };




template <typename ProblemT, typename NumericalFluxT>
class ConservationEquation1D
{   
public:
    ConservationEquation1D(size_t n, double dx, ProblemT& problem, NumericalFluxT& fluxPolicy, double C = 0.9) : _n(n), _nCells(n + 2), _dx(dx), _C(C), _problem(problem), _fluxPolicy(fluxPolicy)
    {
        this->_flux.resize(_n + 1);
        this->_fluxFunctionAtCells.resize(_n + 2);
    }
    double getMaxDt(std::ranges::input_range auto&& u) {
        return _problem.getMaxDt(u, this->_C, this->_dx);
    };

    const v1Dvalue_t &getFlux(std::ranges::input_range auto const& u, double dt)
    {
        // Precompute f(u) for every cell node
        for (size_t i = 0; i < this->_n + 2; ++i)
        {
            _fluxFunctionAtCells[i] = _problem.f(u[i]);
        }

        for (size_t i = 0; i < this->_n + 1; ++i)
        {
            _flux[i] = _fluxPolicy.getNumericalFlux(_fluxFunctionAtCells, u, i, dt, _dx);
        }

        return _flux;
    }

protected:
    size_t _nCells;
    size_t _n;
    v1Dvalue_t _flux;
    v1Dvalue_t _fluxFunctionAtCells;
    double _dx;
    double _dxdt;
    double _dtdx;
    double _C;
    ProblemT& _problem;
    NumericalFluxT& _fluxPolicy;
};

namespace NumericalFlux {
    template <typename F, typename ProblemT>
    class FluxPolicy {
    public:
        ProblemT& problem;
        FluxPolicy(ProblemT& problem_): problem(problem_) {

        };

        value_t getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return static_cast<F&>(*this)._getNumericalFlux(fluxFunctionAtCells, u, face_idx, dt, dx);
        }

        // Default is CD
        value_t _getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return 0.5 * (fluxFunctionAtCells[face_idx] + fluxFunctionAtCells[face_idx + 1]);
        };

        value_t _getNumericalFluxLF(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return 0.5 * dx / dt * (u[face_idx] - u[face_idx + 1]) + 0.5 * (fluxFunctionAtCells[face_idx] + fluxFunctionAtCells[face_idx + 1]);
        };

        value_t _getNumericalFluxRI(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx)
        {
            const value_t u_n = 0.5 * (u[face_idx] + u[face_idx + 1] + dt / dx * (fluxFunctionAtCells[face_idx] - fluxFunctionAtCells[face_idx + 1]));
            return problem.f(u_n);
        }
    };

    template <typename ProblemT>
    class Forward: public FluxPolicy<Forward<ProblemT>, ProblemT> {
    public:
        value_t _getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return fluxFunctionAtCells[face_idx];
        }

    };

    template <typename ProblemT>
    class Central: public FluxPolicy<Central<ProblemT>, ProblemT>  {
    public:
        value_t _getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return 0.5 * (fluxFunctionAtCells[face_idx] + fluxFunctionAtCells[face_idx + 1]);
        }
    };

    template <typename ProblemT>
    class LF: public FluxPolicy<LF<ProblemT>, ProblemT>  {
    public:
        value_t _getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return this->_getNumericalFluxLF(fluxFunctionAtCells, u, face_idx, dt, dx);;
        }
    };

    template <typename ProblemT>
    class RI: public FluxPolicy<RI<ProblemT>, ProblemT>  {
    public:
        value_t _getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return this->_getNumericalFluxRI(fluxFunctionAtCells, u, face_idx, dt, dx);;
        }
    };

    template <typename ProblemT>
    class GOD: public FluxPolicy<GOD<ProblemT>, ProblemT>  {
    public:
        value_t _getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return this->problem.f(this->problem.RiemannSolver(u[face_idx], u[face_idx + 1]));
        }
    };

    template <typename ProblemT>
    class FORCE: public FluxPolicy<FORCE<ProblemT>, ProblemT>  {
    public:
        value_t _getNumericalFlux(std::ranges::input_range auto&& fluxFunctionAtCells, std::ranges::input_range auto&& u, size_t face_idx, double dt, double dx) {
            return 0.5 * (this->_getNumericalFluxRI(fluxFunctionAtCells, u, face_idx, dt, dx) + this->_getNumericalFluxLF(fluxFunctionAtCells, u, face_idx, dt, dx));
        }
    };
}


template <typename DIRECTION>
class EulerProblem1D {
public:
    explicit EulerProblem1D(double gamma, double eps) : _gamma(gamma), _eps(eps) {
        this->_max_it = 100;
    };

    value_t f(value_t u)
    {   
        return static_cast<DIRECTION&>(*this)._f(u);
    }

    double pressure_function(double p_star, double rhoK, double vK, double pK, double A_K, double B_K, double cs_K) {
        if(p_star > pK) {
            // Shock
            return (p_star - pK) * sqrt(A_K / (p_star + B_K));
        } else {
            // Rarefaction
            return (2. * cs_K / (this->_gamma - 1.)) * (std::pow(p_star / pK, (this->_gamma - 1) / (2*this->_gamma)) - 1);
        }
    }

    double pressure_function_derivative(double p_star, double rhoK, double pK, double A_K, double B_K, double cs_K) {
        if(p_star > pK) {
            // Shock
            return sqrt(A_K / (B_K + p_star)) * (1. - (p_star - pK) / (2*(B_K + p_star)));
        } else {
            // Rarefaction
            return (1. / (rhoK*cs_K)) * std::pow(p_star / pK, -(_gamma + 1)/(2*_gamma)); 
        }
    }
    
    double newtons_raphson_pressure(const value_t& uL, const value_t& uR) {
        const double& rhoL = uL[0];
        const double& vL = uL[1];
        const double& pL = uL[2];

        const double& rhoR = uR[0];
        const double& vR = uR[1];
        const double& pR = uR[2];

        const double A_L = 2. / ((this->_gamma + 1) * rhoL);
        const double A_R = 2. / ((this->_gamma + 1) * rhoR);
        const double B_L  = pL * (this->_gamma - 1) / (this->_gamma + 1);
        const double B_R  = pR * (this->_gamma - 1) / (this->_gamma + 1);

        const double cs_L = sqrt(this->_gamma * pL / rhoL);
        const double cs_R = sqrt(this->_gamma * pR / rhoR);

        double p_star_old;
        // double p_star = 0.5 * (pL + pR);
        double p_star = std::max(_eps, 0.5 * (pL + pR) - 0.125 * (vR - vL) * (rhoL + rhoR)*(cs_L + cs_R));

        int i = 0;
        do {
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

        // if(uL != uR)
        //     std::cout << i << " " << pL << " " << pR << " " <<  p_star << "\n";
    
        return p_star;
    }

    double compute_rho_star_K(double p_star, double pK, double rhoK) {
        if(p_star > pK) { // Shock
            return rhoK * std::pow(p_star / pK, (1. / _gamma));
        } else {
            double r = p_star / pK;
            double gr = (_gamma - 1) / (_gamma + 1);
            return rhoK * (r + gr) / (gr * r + 1);
        }
    }

    // PRE: cs_k_hat = cs_L or -cs_R
    double compute_S_K(double p_star, double pK, double vK, double cs_K_hat){
        return vK - cs_K_hat * sqrt((_gamma + 1) * p_star / (2*_gamma * pK) + (_gamma - 1) / (2*_gamma));
    }

    // PRE: cs_k_hat = cs_L or -cs_R
    value_t rarefaction_fan_K(double rhoK, double vK, double pK, double cs_K_hat) {
        const double rho = rhoK * std::pow(2. / (_gamma + 1.) + vK * (_gamma - 1.) / ((_gamma + 1.)*cs_K_hat), 2 / (_gamma - 1));
        const double v = (2 / (_gamma + 1)) * (cs_K_hat + vK * (_gamma - 1) / 2.);
        const double p = pK * std::pow(2. / (_gamma + 1) + vK *(_gamma - 1.) / ((_gamma + 1) * cs_K_hat), (2. * _gamma) / (_gamma - 1));
        return {rho, v, p};
    }

    value_t RiemannSolver(const value_t& u0, const value_t& u1) {
        // get primitive variables
        const value_t uL = utils::conservativeToPrimitive(this->_gamma, u0);
        const value_t uR = utils::conservativeToPrimitive(this->_gamma, u1);

        return utils::primitiveToConservative2D(this->_gamma, RiemannSolverPrimitive(uL, uR));
    }

    double getMaxDt(std::ranges::input_range auto&& u, double C, double dx)
    {   
        double max_a = 0;
        for(const auto& el: u) {
            const double cs = sqrt(_gamma* (_gamma - 1.) * (el[3] - 0.5*el[1]*el[1]/el[3]) / el[0]);
            const double a = abs(el[1] / el[0]) + cs;
            max_a = std::max<double>(max_a, a);
        }
        return C * dx / max_a;
    }

protected:
    value_t RiemannSolverPrimitive(const value_t& uL, const value_t& uR)
    {
        const double& rhoL = uL[0];
        const double& vL = uL[1];
        const double& pL = uL[2];

        const double& rhoR = uR[0];
        const double& vR = uR[1];
        const double& pR = uR[2];

        const double p_star = newtons_raphson_pressure(uL, uR);

        const double A_L = 2. / ((this->_gamma + 1) * rhoL);
        const double A_R = 2. / ((this->_gamma + 1) * rhoR);
        const double B_L  = pL * (this->_gamma - 1) / (this->_gamma + 1);
        const double B_R  = pR * (this->_gamma - 1) / (this->_gamma + 1);

        const double cs_L = sqrt(this->_gamma * pL / rhoL);
        const double cs_R = sqrt(this->_gamma * pR / rhoR);

        const double v_star = 0.5 * (vL + vR) + 0.5 * (pressure_function(p_star, rhoR, vR, pR, A_R, B_R, cs_R) - pressure_function(p_star, rhoL, vL, pL, A_L, B_L, cs_L));

        const double rho_star_L = compute_rho_star_K(p_star, pL, rhoL);
        const double rho_star_R = compute_rho_star_K(p_star, pR, rhoR);

        const double cs_star_L = sqrt(this->_gamma * p_star / rho_star_L);
        const double cs_star_R = sqrt(this->_gamma * p_star / rho_star_R);


        // Left Wave is right of center -> uL
        if(p_star > pL) { // If shock, check if right of x = 0 
            double S_L = compute_S_K(p_star, pL, vL, cs_L);
            if(S_L > 0) return {rhoL, vL, pL};
        } else { // If rarefaction, check if tail right of x = 0
            if((vL - cs_L) > 0) return {rhoL, vL, pL}; 
        }

        // Right wave is left of center -> uR
        if(p_star > pR) { // If shock, check if left of x = 0 
            double S_R = compute_S_K(p_star, pR, vR, -cs_R);
            if(S_R < 0) return {rhoR, vR, rhoR};
        } else { // If rarefaction, check if tail left of x = 0
            if((vR + cs_R) < 0) return {rhoR, vR, pR}; 
        }

        // Left wave is a rarefaction fan that spans over the center
        if((p_star < pL) && ((vL - cs_L) < 0) && ((v_star - cs_star_L) > 0)) {
            return rarefaction_fan_K(rhoL, vL, pL, cs_L);
        }

        // Right wave is a rarefaction fan that spans over the center
        if((p_star < pR) && ((vR + cs_R) > 0) && ((v_star + cs_star_R) < 0)) {
            return rarefaction_fan_K(rhoR, vR, pR, -cs_R);
        }

        // Left wave is fully on the left, right wave is fully the right (rarefaction tails included)
        // if v_star < 0 -> interface is at the right of the contact discontinuity
        // else if v_star > 0 -> interface is at the left of the contact discontinuity
        if(v_star < 0) {
            return {rho_star_R, v_star, p_star};
        }

        return {rho_star_L, v_star, p_star};
    }

    double _gamma;
    double _eps;
    int _max_it;
};

class F: public EulerProblem1D<F> {
public:
    value_t _f(value_t u)
    {   
        auto q = utils::conservativeToPrimitive(this->_gamma, u);
        const auto& rho = q[0];
        const auto& vx = q[1];
        const auto& vy = q[2];
        const auto& p = q[3];
        const auto& E = u[3];

        return {
            rho * vx,
            rho * vx * vx + p,
            rho * vx * vy,
            (E + p)*vx,
        };
    }
};

class G: public EulerProblem1D<G> {
public:
    value_t _f(value_t u)
    {   
        auto q = utils::conservativeToPrimitive(this->_gamma, u);
        const auto& rho = q[0];
        const auto& vx = q[1];
        const auto& vy = q[2];
        const auto& p = q[3];
        const auto& E = u[3];

        return {
            rho * vy,
            rho * vx * vy,
            rho * vy * vy + p,
            (E + p)*vy,
        };
    }
};

class Mesh
{

public:
    // n is the number of cells
    explicit Mesh(int nx, int ny, double min_x, double max_x, double min_y, double max_y) : _nx(nx), _ny(ny), _min_x(min_x), _max_x(max_x), _min_y(min_y), _max_y(max_y), _dx((max_x - min_x) / (double)nx), _dy((max_y - min_y) / (double)ny)
    {
        const size_t nPointsx = nx + 1;
        const size_t nPointsy = ny + 1;

        this->_lx = (_max_x - _min_x);
        this->_ly = (_max_y - _min_y);

        this->_points.resize(nPointsx);
        for (int i = 0; i < nPointsx; ++i) {
            this->_points[i].resize(nPointsy);
            for(int j = 0; j < nPointsy; j++) {
                this->_points[i][j] = {this->_lx * i / (double)nPointsx, this->_ly * j / (double)nPointsy};
            }
        }

        this->_cellPoints.resize(nx);
        for (int i = 0; i < nx; ++i) {
            this->_cellPoints[i].resize(ny);
            for(int j = 0; j < ny; j++) {
                const point2D_t& p = this->_points[i][j];
                this->_cellPoints[i][j] = {p[0] + _dx * 0.5, p[1] + _dy * 0.5};
            }
        }

        this->_fullCellPoints.resize(nx + 2);
        for(size_t i = 0; i < nx + 2; ++i) {
            this->_fullCellPoints[i].resize(ny + 2);
        }

        // Copy closest cell point on ghost cells
        // 4 corners
        this->_fullCellPoints[0][0] = _cellPoints[0][0];
        this->_fullCellPoints[nx + 1][ny + 1] = _cellPoints[nx - 1][ny - 1];
        this->_fullCellPoints[nx + 1][0] = _cellPoints[nx - 1][0];
        this->_fullCellPoints[0][ny + 1] = _cellPoints[0][ny - 1];


        for(int j = 0; j < ny; ++j) {
            this->_fullCellPoints[0][j+1] = _cellPoints[0][j];
            this->_fullCellPoints[1][j+1] = _cellPoints[0][j];

            this->_fullCellPoints[nx][j+1] = _cellPoints[nx - 1][j];
            this->_fullCellPoints[nx + 1][j+1] = _cellPoints[nx - 1][j];
        }

        for(int i = 0; i < nx; ++i) {
            this->_fullCellPoints[i+1][0] = _cellPoints[i][0];
            this->_fullCellPoints[i+1][1] = _cellPoints[i][0];

            this->_fullCellPoints[i+1][ny] = _cellPoints[i][ny - 1];
            this->_fullCellPoints[i+1][ny + 1] = _cellPoints[i][ny - 1];
        }

        for (int i = 1; i < nx + 1; ++i) {
            for(int j = 1; j < ny + 1; ++j) {
                this->_fullCellPoints[i][j] = _cellPoints[i - 1][j - 1];
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

private:
    const double _dx;
    const double _dy;
    const double _min_x;
    const double _max_x;
    const double _min_y;
    const double _max_y;
    double _lx{};
    double _ly{};
    const int _nx;
    const int _ny;
    vpoint2D_t _points;
    vpoint2D_t _cellPoints;
    vpoint2D_t _fullCellPoints; // Include ghost cells
};


template <typename F1, typename F2, typename P1, typename P2>
class Simulation
{
public:
    Simulation(double final_time, ConservationEquation1D<P1, F1> &eq_f, ConservationEquation1D<P2, F2> &eq_g, Mesh mesh, std::function<value_t(point2D_t)> ic) : _final_time(
                                                                                                                  final_time),
                                                                                                              _eq_f(eq_f), _eq_g(eq_g),_mesh(mesh)
    {
        const vpoint2D_t &mesh_points = mesh.getPoints();
        const vpoint2D_t &mesh_fullCellPoints = mesh.getFullCellPoints();
        this->nx = mesh.getCellPoints().size();
        this->ny = mesh.getCellPoints()[0].size();

        this->solution.resize(this->nx + 2);
        this->_solution_buffer.resize(this->nx + 2);
        for (size_t i = 0; i < this->nx + 2; ++i) {
            this->solution[i].resize(this->ny+2);
            this->_solution_buffer[i].resize(this->ny+2);
        }

        for (int i = 0; i < this->nx + 2; i++)
        {
            for(int j = 0; j < this->ny + 2; j++) {
                this->solution[i][j] = ic(mesh_fullCellPoints[i][j]);
            }
        }

        // Transmissive BCs
        this->solution[0][0] = this->solution[1][1];
        this->solution[nx+1][ny+1] = this->solution[nx][ny];
        this->solution[nx + 1][0] = this->solution[nx][0];
        this->solution[0][ny + 1] = this->solution[0][ny];

        for(int j = 1; j < ny + 1; ++j) {
            this->solution[0][j] = this->solution[1][j];
            this->solution[nx + 1][j] = this->solution[nx][j];
        }

        for(int i = 1; i < nx + 1; ++i) {
            this->solution[i][0] = this->solution[i][1];
            this->solution[i][ny + 1] = this->solution[i][ny];
        }
    };

    bool evolve()
    {
        double rest_dt = (_final_time - _ct);
        double min_dt = std::min(rest_dt, this->_eq_g.getMaxDt(solution[0]));

        // Y dt sweep
        for(size_t i = 1; i < this->nx + 1; ++i) {
            min_dt = std::min(min_dt, this->_eq_g.getMaxDt(solution[i]));
        }

        for(size_t j = 1; j < this->ny + 1; ++j) {
            ColumnView col(solution, j);
            min_dt = std::min(min_dt, this->_eq_f.getMaxDt(col.as_range()));
        }

        const double _dtdx = min_dt / this->_mesh.getDx();

        // X sweep
        for (size_t i = 1; i < this->nx + 1; ++i) {
            const v1Dvalue_t &flux = this->_eq_f.getFlux(solution[i], min_dt);
            for(size_t j = 1; j < this->ny + 1; ++j)
            {
                this->_solution_buffer[i][j] = solution[i][j] - _dtdx * (flux[j] - flux[j - 1]);
            }
        }

        for(size_t i = 1; i < this->nx + 1; ++i) {
            for(size_t j = 1; j < this->ny + 1; ++j) {
                this->solution[i][j] = this->_solution_buffer[i][j];
            }
        }

        // Transmissive BCs
        this->solution[0][0] = this->solution[1][1];
        this->solution[nx+1][ny+1] = this->solution[nx][ny];
        this->solution[nx + 1][0] = this->solution[nx][0];
        this->solution[0][ny + 1] = this->solution[0][ny];

        for(int j = 1; j < ny + 1; ++j) {
            this->solution[0][j] = this->solution[1][j];
            this->solution[nx + 1][j] = this->solution[nx][j];
        }

        for(int i = 1; i < nx + 1; ++i) {
            this->solution[i][0] = this->solution[i][1];
            this->solution[i][ny + 1] = this->solution[i][ny];
        }

        const double _dtdy = min_dt / this->_mesh.getDy();

        // Y Sweep
        for(size_t j = 1; j < this->ny + 1; ++j) {
            ColumnView col(solution, j);
            const v1Dvalue_t &flux = this->_eq_g.getFlux(col.as_range(), min_dt);
            for (size_t i = 1; i < this->nx + 1; ++i)
            {
                this->_solution_buffer[i][j] = solution[i][j] - _dtdy * (flux[i] - flux[i - 1]);
            }
        }

        for(size_t i = 1; i < this->nx + 1; ++i) {
            for(size_t j = 1; j < this->ny + 1; ++j) {
                this->solution[i][j] = this->_solution_buffer[i][j];
            }
        }

        // Transmissive BCs
        this->solution[0][0] = this->solution[1][1];
        this->solution[nx+1][ny+1] = this->solution[nx][ny];
        this->solution[nx + 1][0] = this->solution[nx][0];
        this->solution[0][ny + 1] = this->solution[0][ny];

        for(int j = 1; j < ny + 1; ++j) {
            this->solution[0][j] = this->solution[1][j];
            this->solution[nx + 1][j] = this->solution[nx][j];
        }

        for(int i = 1; i < nx + 1; ++i) {
            this->solution[i][0] = this->solution[i][1];
            this->solution[i][ny + 1] = this->solution[i][ny];
        }

        // std::cout << min_dt << "\t" << _ct << "\r";
        _ct += min_dt;

        return _ct < _final_time;
    };


    void plot2d() {
        plt::clf();

        const size_t nx = solution.size();
        const size_t ny = solution[0].size();

        vvalue_t primitive_values = utils::conservativeToPrimitive(1.4, solution);
        std::vector<float> density(nx*ny);
        std::vector<float> velocity_x(nx*ny);
        std::vector<float> velocity_y(nx*ny);
        std::vector<float> pressure(nx*ny);

        for(int i = 0; i < nx; ++i) {
            for(int j = 0; j < ny; ++j) {
                density[j + i*ny] = primitive_values[i][j][0];
                velocity_x[j + i*ny] = primitive_values[i][j][1];
                velocity_y[j + i*ny] = primitive_values[i][j][2];
                pressure[j + i*ny] = primitive_values[i][j][3];
            }
        }

        std::cout << std::ranges::max(velocity_x) << std::endl;
        std::cout << std::ranges::max(velocity_y) << std::endl;

        plt::suptitle(boost::lexical_cast<std::string>(this->_ct));
        const float* pVx = &(velocity_x[0]);
        const float* pVy = &(velocity_y[0]);
        const float* pD = &(density[0]);
        const float* pP = &(pressure[0]);
        
        plt::subplot(2, 2, 1);
        plt::imshow(pP, nx, ny, 1);

        plt::subplot(2, 2, 2);
        plt::imshow(pD, nx, ny, 1);

        plt::subplot(2, 2, 3);
        plt::imshow(pVx, nx, ny, 1);

        plt::subplot(2, 2, 4);
        plt::imshow(pVy, nx, ny, 1);

        plt::pause(0.01);
        // plt::show();


    }
    void plot(double blim = -1.0, double tlim = 1.0)
    {
        plt::clf();
        
        std::vector<value_t> primitive_values(solution.size());
        std::vector<double> density(solution.size());
        std::vector<double> velocity_x(solution.size());
        std::vector<double> velocity_y(solution.size());
        std::vector<double> pressure(solution.size());
        // plt::ylim(std::ranges::min(density), std::ranges::max(density));
        
        // TODO
        std::ranges::transform(solution, primitive_values.begin(), [&](value_t u){return utils::conservativeToPrimitive(1.4, u);});
        std::ranges::transform(primitive_values, density.begin(), [&](value_t u){return u[0];});
        std::ranges::transform(primitive_values, velocity_x.begin(), [&](value_t u){return u[1];});
        std::ranges::transform(primitive_values, velocity_y.begin(), [&](value_t u){return u[2];});
        std::ranges::transform(primitive_values, pressure.begin(), [&](value_t u){return u[3];});

        
        // plt::ylim(std::ranges::min(density), std::ranges::max(density));

        // for(const auto& el: density) {
        //     std::cout << el << " ";
        // }
        // std::cout << "\n";
        plt::suptitle(boost::lexical_cast<std::string>(this->_ct));
        plt::subplot(2, 2, 1);
        plt::plot(_mesh.getFullCellPoints(), density, "r*");
        plt::subplot(2, 2, 2);
        plt::plot(_mesh.getFullCellPoints(), pressure, "r*");
        plt::subplot(2, 2, 3);
        plt::plot(_mesh.getFullCellPoints(), velocity_x, "r*");
        plt::subplot(2, 2, 3);
        plt::plot(_mesh.getFullCellPoints(), velocity_x, "r*");
        // plt::show();
        plt::pause(0.01);
    }

    vvalue_t solution;

private:
    size_t nx;
    size_t ny;

    double _dt;
    double _ct = 0.;
    double _final_time;
    ConservationEquation1D<P1, F1> &_eq_f;
    ConservationEquation1D<P2, F2> &_eq_g;
    Mesh _mesh;
    vvalue_t _solution_buffer;
};


int main(int argc, char *argv[])
{
    int nx = 200;
    int ny = 200;

    double min_x = 0;
    double max_x = 2;
    double min_y = 0.;
    double max_y = 2.;
    double final_time = 10.0;
    double C = 0.8;

    Mesh mesh(nx, ny, min_x, max_x, min_y, max_y);

    // value_t a = 1.0;
    // F
    EulerProblem1D<F> problem_F(1.4, 1e-10);
    EulerProblem1D<G> problem_G(1.4, 1e-10);

    NumericalFlux::FluxPolicy<NumericalFlux::FORCE<EulerProblem1D<F>>, EulerProblem1D<F>> fluxPolicy_F(problem_F);
    NumericalFlux::FluxPolicy<NumericalFlux::FORCE<EulerProblem1D<G>>, EulerProblem1D<G>> fluxPolicy_G(problem_G);

    ConservationEquation1D equation_f(nx, mesh.getDx(), problem_F, fluxPolicy_F, C);
    // G
    ConservationEquation1D equation_g(ny, mesh.getDy(), problem_G, fluxPolicy_G, C);


    Simulation sim(final_time, equation_f, equation_g, mesh, utils::test3_2d);
    // sim.evolve();
    // sim.plot();
    // sim.evolve();
    plt::figure_size(1000, 600);

    sim.plot2d();
    int i = 0;
    while (sim.evolve())
    {
        if(i%5 == 0)
            sim.plot2d();
        ++i;
    }

    return 0;
}
