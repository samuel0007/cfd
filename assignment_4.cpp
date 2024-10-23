//
// Created by russo on 14.10.2024.
//

#include <iostream>
#include <utility>
#include <vector>
#include <ranges>
#include <algorithm>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef std::vector<double> value_t;

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

    value_t primitiveToConservative(double gamma, value_t u){
        return {
            u[0],
            u[0]*u[1],
            u[2] / (gamma - 1.) + 0.5 * u[0] * u[1] * u[1],
        };
    }

    value_t conservativeToPrimitive(double gamma, value_t u) {
        return {
            u[0],
            u[1] / u[0],
            (gamma - 1.) * (u[2] - 0.5 * u[1] * u[1] / u[0]),
        };
    }

    value_t test0(double x) {
        return x <= 0.25 ? primitiveToConservative(1.4, {1.0, 1.0, 1.0}) : primitiveToConservative(1.4, {0.1, 1.0, 1.});
    }

    value_t test1(double x) {
        return x <= 0.5 ? primitiveToConservative(1.4, {1.0, 0.0, 1.0}) : primitiveToConservative(1.4, {0.125, 0., 0.1});
    }
}

enum NumericalFlux
{
    Central,
    Backward,
    Forward,
    LF,
    RI,
    FORCE,
    GOD,
};

template <NumericalFlux F>
class ConservationEquation
{
public:
    ConservationEquation(size_t n, double dx, double C = 0.9) : _n(n), _nCells(n + 2), _dx(dx), _C(C)
    {
        this->_flux.resize(_n + 1);
        this->_fluxFunctionAtCells.resize(_n + 2);
    }

    virtual value_t f(value_t u) = 0;
    virtual value_t RiemannSolver(const value_t& u0, const value_t& u1) = 0;
    virtual double getMaxDt(const std::vector<value_t> &u) = 0;

    // u.size() = n + 2
    const std::vector<value_t> &getFlux(const std::vector<value_t> &u, double dt)
    {
        this->_dtdx = dt / this->_dx;
        this->_dxdt = this->_dx / dt;

        // Precompute f(u) for every cell node
        for (size_t i = 0; i < this->_n + 2; ++i)
        {
            _fluxFunctionAtCells[i] = f(u[i]);
        }

        for (size_t i = 0; i < this->_n + 1; ++i)
        {
            _flux[i] = this->_getNumericalFlux(u, i);
        }

        return _flux;
    }

protected:
    // Default is CD
    value_t _getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
    {
        return 0.5 * (this->_fluxFunctionAtCells[face_idx] + this->_fluxFunctionAtCells[face_idx + 1]);
    };

    value_t _getNumericalFluxLF(const std::vector<value_t> &u, size_t face_idx)
    {
        return 0.5 * this->_dxdt * (u[face_idx] - u[face_idx + 1]) + 0.5 * (_fluxFunctionAtCells[face_idx] + _fluxFunctionAtCells[face_idx + 1]);
    }

    value_t _getNumericalFluxRI(const std::vector<value_t> &u, size_t face_idx)

    {
        const value_t u_n = 0.5 * (u[face_idx] + u[face_idx + 1] + _dtdx * (_fluxFunctionAtCells[face_idx] - _fluxFunctionAtCells[face_idx + 1]));
        return f(u_n);
    }

    size_t _nCells;
    size_t _n;
    std::vector<value_t> _flux;
    std::vector<value_t> _fluxFunctionAtCells;
    double _dx;
    double _dxdt;
    double _dtdx;
    double _C;
};

template <>
value_t ConservationEquation<NumericalFlux::Forward>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return this->_fluxFunctionAtCells[face_idx];
}

template <>
value_t ConservationEquation<NumericalFlux::Central>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return 0.5 * (this->_fluxFunctionAtCells[face_idx] + this->_fluxFunctionAtCells[face_idx + 1]);
}

template <>
value_t ConservationEquation<NumericalFlux::Backward>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return this->_fluxFunctionAtCells[face_idx + 1];
}

template <>
value_t ConservationEquation<NumericalFlux::LF>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return _getNumericalFluxLF(u, face_idx);
}

template <>
value_t ConservationEquation<NumericalFlux::RI>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return _getNumericalFluxRI(u, face_idx);
}

template <>
value_t ConservationEquation<NumericalFlux::FORCE>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return 0.5 * (_getNumericalFluxRI(u, face_idx) + _getNumericalFluxLF(u, face_idx));
}

template <>
value_t ConservationEquation<NumericalFlux::GOD>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return f(RiemannSolver(u[face_idx], u[face_idx + 1]));
}


template <NumericalFlux T>
class EulerEquation : public ConservationEquation<T>
{
public:
    explicit EulerEquation(double gamma, double eps, auto... args) : _gamma(gamma),_eps(eps), ConservationEquation<T>(args...) {
        this->_max_it = 100;
    };

    value_t f(value_t u)
    {   
        const auto& u1 = u[0];
        const auto& u2 = u[1];
        const auto& u3 = u[2];

        return {
            u2,
            0.5*(3-this->_gamma) * u2 * u2 / u1 + (this->_gamma - 1) * u3,
            this->_gamma * u2 * u3 / u1 - 0.5 * (this->_gamma - 1) * u2 * u2 * u2 / (u1 * u1)
        };
    }

    double pressure_function(double p_star, double rhoK, double vK, double pK, double A_K, double B_K, double cs_K) {
        if(p_star > pK) {
            // Shock
            return (p_star - pK) * sqrt(A_K / (p_star + B_K));
        } else {
            // Rarefaction
            return 2. * cs_K / (this->_gamma - 1.) * (std::pow(p_star / pK, (this->_gamma - 1) / 2*this->_gamma) - 1);
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
        const double A_R = 2. / ((this->_gamma + 1) * rhoL);
        const double B_K  = (this->_gamma - 1) / (this->_gamma + 1);

        const double cs_L = sqrt(this->_gamma * pL / rhoL);
        const double cs_R = sqrt(this->_gamma * pR / rhoR);

        double p_star_old;
        double p_star = 0.5 * (pL + pR);

        int i = 0;
        do {
            p_star_old = p_star;
            const double f_L = pressure_function(p_star_old, rhoL, vL, pL, A_L, B_K, cs_L);
            const double f_R = pressure_function(p_star_old, rhoR, vR, pR, A_R, B_K, cs_R);
            const double f = f_R + f_R + (vR - vL);
            const double f_d_L = pressure_function_derivative(p_star_old, rhoL, pL, A_L, B_K, cs_L);
            const double f_d_R = pressure_function_derivative(p_star_old, rhoR, pR, A_R, B_K, cs_R);
            const double f_d = f_d_L + f_d_R;
            p_star = p_star_old - f / f_d;
            ++i;
        } while ((fabs(p_star - p_star_old) / p_star_old) > _eps && i < this->_max_it);

        // if(i != 1)
        std::cout << i << " " << pL<< " " << pR << " " << p_star << std::endl;

        return p_star;
    }

    value_t RiemannSolver(const value_t& u0, const value_t& u1)
    {
        // get primitive variables
        const value_t uL = utils::conservativeToPrimitive(this->_gamma, u0);
        const value_t uR = utils::conservativeToPrimitive(this->_gamma, u1);
        // Find p*
        const double p_star = newtons_raphson_pressure(uL, uR);

        // Find v*
        // Find rho_K*
        // 

       return 0.5*(u0 + u1);
    }



    double getMaxDt(const std::vector<value_t> &u)
    {   
        double max_a = 0;
        for(const auto& el: u) {
            const double cs = sqrt(_gamma* (_gamma - 1.) * (el[2] - 0.5*el[1]*el[1]/el[2]) / el[0]);
            const double a = abs(el[1] / el[0]) + cs;
            max_a = std::max<double>(max_a, a);
        }
        return this->_C * this->_dx / max_a;
    }

private:
    double _gamma;
    double _eps;
    int _max_it;
};

class Mesh
{
public:
    // n is the number of cells
    explicit Mesh(int n, double min_x, double max_x) : _n(n), _min_x(min_x), _max_x(max_x), _dx((max_x - min_x) / (double)n)
    {
        const size_t nPoints = n + 1;
        this->_l = (_max_x - _min_x);

        this->_points.resize(nPoints);
        for (int i = 0; i < nPoints; ++i)
        {
            this->_points[i] = this->_l * i / (double)nPoints;
        }

        this->_cellPoints.resize(n);
        for (int i = 0; i < n; ++i)
        {
            this->_cellPoints[i] = this->_points[i] + _dx * 0.5;
        }

        this->_fullCellPoints.resize(n + 2);
        this->_fullCellPoints[0] = this->_points[0] - _dx * 0.5; // avoid mess up of the plotting
        for (int i = 1; i < n + 1; ++i)
        {
            this->_fullCellPoints[i] = this->_cellPoints[i - 1];
        }
        this->_fullCellPoints[n + 1] = this->_points[n] + _dx * 0.5;
    }

    const std::vector<double> &getPoints()
    {
        return _points;
    }

    const std::vector<double> &getCellPoints()
    {
        return _cellPoints;
    }

    const std::vector<double> &getFullCellPoints()
    {
        return _fullCellPoints;
    }

    const double &getDx()
    {
        return _dx;
    }

private:
    double _dx;
    double _min_x;
    double _max_x;
    double _l{};
    int _n;
    std::vector<double> _points;
    std::vector<double> _cellPoints;
    std::vector<double> _fullCellPoints; // Include ghost cells
};

template <NumericalFlux T>
class Simulation
{
public:
    Simulation(double final_time, ConservationEquation<T> &eq, Mesh mesh, std::function<value_t(double)> ic) : _final_time(
                                                                                                                  final_time),
                                                                                                              _eq(eq), _mesh(mesh)
    {
        const std::vector<double> &mesh_points = mesh.getPoints();
        const std::vector<double> &mesh_fullCellPoints = mesh.getFullCellPoints();
        this->n = mesh.getCellPoints().size();

        this->solution.resize(n + 2);
        this->_solution_buffer.resize(n + 2);

        for (int i = 1; i < n + 1; ++i)
        {
            solution[i] = ic(mesh_fullCellPoints[i]);
        }
        // Transmissive BCs
        solution[0] = solution[1];
        solution[n + 1] = solution[n];
    };

    bool evolve()
    {
        double rest_dt = (_final_time - _ct);
        double min_dt = std::min(rest_dt, this->_eq.getMaxDt(solution));

        const std::vector<value_t> &flux = this->_eq.getFlux(solution, min_dt);
        const double _dtdx = min_dt / this->_mesh.getDx();
        for (size_t i = 1; i < n + 1; ++i)
        {
            this->_solution_buffer[i] = solution[i] - _dtdx * (flux[i] - flux[i - 1]);
        }
        this->solution = _solution_buffer;

        // Transmissive BCs
        this->solution[0] = this->solution[1];
        this->solution[n + 1] = this->solution[n];

        std::cout << min_dt << "\t" << _ct << "\r";
        _ct += min_dt;

        return _ct < _final_time;
    };

    void plot(double blim = -1.0, double tlim = 1.0)
    {
        plt::clf();
        // plt::ylim(blim, tlim);
        std::vector<double> density;
        density.resize(solution.size());
        std::ranges::transform(solution, density.begin(), [&](value_t u){return u[0];});
        plt::plot(_mesh.getFullCellPoints(), density, "r*");
        plt::show();
        // plt::pause(0.01);
    }

    std::vector<value_t> solution;

private:
    size_t n;
    double _dt;
    double _ct = 0.;
    double _final_time;
    ConservationEquation<T> &_eq;
    Mesh _mesh;
    std::vector<value_t> _solution_buffer;
};




int main(int argc, char *argv[])
{
    int n = 200;
    double min_x = 0;
    double max_x = 1;
    double final_time = 10.0;
    double C = 0.8;

    Mesh mesh(n, min_x, max_x);

    // value_t a = 1.0;
    // AdvectionEquation<NumericalFlux::LF> equation(a, n, mesh.getDx(), C);
    EulerEquation<NumericalFlux::GOD> equation(1.4, 1e-6, n, mesh.getDx(), C);

    Simulation sim(final_time, equation, mesh, utils::test1);
    sim.evolve();

    sim.plot();
    // while (sim.evolve())
    // {
    //     sim.plot();
    // }

    return 0;
}
