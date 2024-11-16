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
        return x <= 0.25 ? primitiveToConservative2D(1.4, {1.0, 1.0, 1.0}) : primitiveToConservative2D(1.4, {0.1, 1.0, 1.});
    }

    value_t test1(double x) {
        return x <= 0.5 ? primitiveToConservative2D(1.4, {1.0, 0.0, 1.0}) : primitiveToConservative2D(1.4, {0.125, 0., 0.1});
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
class ConservationEquation1D
{
public:
    ConservationEquation1D(size_t n, double dx, double C = 0.9) : _n(n), _nCells(n + 2), _dx(dx), _C(C)
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
value_t ConservationEquation1D<NumericalFlux::Forward>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return this->_fluxFunctionAtCells[face_idx];
}

template <>
value_t ConservationEquation1D<NumericalFlux::Central>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return 0.5 * (this->_fluxFunctionAtCells[face_idx] + this->_fluxFunctionAtCells[face_idx + 1]);
}

template <>
value_t ConservationEquation1D<NumericalFlux::Backward>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return this->_fluxFunctionAtCells[face_idx + 1];
}

template <>
value_t ConservationEquation1D<NumericalFlux::LF>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return _getNumericalFluxLF(u, face_idx);
}

template <>
value_t ConservationEquation1D<NumericalFlux::RI>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return _getNumericalFluxRI(u, face_idx);
}

template <>
value_t ConservationEquation1D<NumericalFlux::FORCE>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return 0.5 * (_getNumericalFluxRI(u, face_idx) + _getNumericalFluxLF(u, face_idx));
}

template <>
value_t ConservationEquation1D<NumericalFlux::GOD>::_getNumericalFlux(const std::vector<value_t> &u, size_t face_idx)
{
    return f(RiemannSolver(u[face_idx], u[face_idx + 1]));
}


template <NumericalFlux T>
class EulerEquation1D : public ConservationEquation1D<T>
{
public:
    explicit EulerEquation1D(double gamma, double eps, auto... args) : _gamma(gamma),_eps(eps), ConservationEquation1D<T>(args...) {
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
    Simulation(double final_time, ConservationEquation1D<T> &eq, Mesh mesh, std::function<value_t(double)> ic) : _final_time(
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

        // std::cout << min_dt << "\t" << _ct << "\r";
        _ct += min_dt;

        return _ct < _final_time;
    };

    void plot(double blim = -1.0, double tlim = 1.0)
    {
        plt::clf();
        
        std::vector<value_t> primitive_values(solution.size());
        std::vector<double> density(solution.size());
        std::vector<double> velocity(solution.size());
        std::vector<double> pressure(solution.size());
        // plt::ylim(std::ranges::min(density), std::ranges::max(density));
        
        // TODO
        std::ranges::transform(solution, primitive_values.begin(), [&](value_t u){return utils::conservativeToPrimitive(1.4, u);});
        std::ranges::transform(primitive_values, density.begin(), [&](value_t u){return u[0];});
        std::ranges::transform(primitive_values, velocity.begin(), [&](value_t u){return u[1];});
        std::ranges::transform(primitive_values, pressure.begin(), [&](value_t u){return u[2];});

        // plt::ylim(std::ranges::min(density), std::ranges::max(density));

        // for(const auto& el: density) {
        //     std::cout << el << " ";
        // }
        // std::cout << "\n";
        plt::suptitle(boost::lexical_cast<std::string>(this->_ct));
        plt::subplot(1, 3, 1);
        plt::plot(_mesh.getFullCellPoints(), density, "r*");
        plt::subplot(1, 3, 2);
        plt::plot(_mesh.getFullCellPoints(), velocity, "r*");
        plt::subplot(1, 3, 3);
        plt::plot(_mesh.getFullCellPoints(), pressure, "r*");
        // plt::show();
        plt::pause(0.01);
    }

    std::vector<value_t> solution;

private:
    size_t n;
    double _dt;
    double _ct = 0.;
    double _final_time;
    ConservationEquation1D<T> &_eq;
    Mesh _mesh;
    std::vector<value_t> _solution_buffer;
};




int main(int argc, char *argv[])
{
    int n = 1000;
    double min_x = 0;
    double max_x = 1;
    double final_time = 1.0;
    double C = 0.8;

    Mesh mesh(n, min_x, max_x);

    // value_t a = 1.0;
    // AdvectionEquation<NumericalFlux::LF> equation(a, n, mesh.getDx(), C);
    EulerEquation1D<NumericalFlux::GOD> equation(1.4, 1e-10, n, mesh.getDx(), C);

    Simulation sim(final_time, equation, mesh, utils::test1);
    // sim.evolve();
    // sim.plot();
    // sim.evolve();
    plt::figure_size(1000, 600);

    sim.plot();
    int i = 0;
    while (sim.evolve())
    {
        if(i%10 == 0)
            sim.plot();
        ++i;
    }

    return 0;
}
