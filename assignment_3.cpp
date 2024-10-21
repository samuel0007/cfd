//
// Created by russo on 14.10.2024.
//

#include <iostream>
#include <utility>
#include <vector>
#include <ranges>
#include <algorithm>
#include "matplotlibcpp.h"
#include <Eigen/Dense>

namespace plt = matplotlibcpp;

typedef double value_t;

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

template <NumericalFlux T>
class ConservationEquation
{
public:
    ConservationEquation(size_t n, double dx, double C = 0.9) : _n(n), _nCells(n + 2), _dx(dx), _C(C)
    {
        this->_flux.resize(_n + 1);
        this->_fluxFunctionAtCells.resize(_n + 2);
    }

    virtual value_t f(value_t u) = 0;
    virtual value_t RiemannSolver(value_t u0, value_t u1) = 0;
    virtual value_t getMaxDt(const std::vector<value_t> &u) = 0;

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
class AdvectionEquation : public ConservationEquation<T>
{
public:
    explicit AdvectionEquation(value_t a, auto... args) : _a(a), ConservationEquation<T>(args...) {};

    value_t f(value_t u)
    {
        return this->_a * u;
    }

    value_t RiemannSolver(value_t u0, value_t u1)
    {
        return this->_a > 0 ? u0 : u1;
    }

    double getMaxDt(const std::vector<value_t> &u)
    {
        return this->_C * this->_dx / abs(this->_a);
    }

private:
    value_t _a;
};

template <NumericalFlux T>
class BurgersEquation : public ConservationEquation<T>
{
public:
    explicit BurgersEquation(auto... args) : ConservationEquation<T>(args...) {};

    value_t f(value_t u)
    {
        return 0.5 * u * u;
    }

    value_t RiemannSolver(value_t u0, value_t u1)
    {
        const value_t S = 0.5 * (u0 + u1);
        if (u0 > u1)
        {
            return S > 0 ? u0 : u1;
        }
        else
        {
            if (u0 > 0)
            {
                return u0;
            }
            else if (u1 < 0)
            {
                return u1;
            }
            else
            {
                return 0.;
            }
        }
    }

    double getMaxDt(const std::vector<value_t> &u)
    {
        return this->_C * this->_dx / std::ranges::max(std::views::transform(u, [](value_t el)
                                                                             { return fabs(el); }));
    }
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
    Simulation(double final_time, ConservationEquation<T> &eq, Mesh mesh, std::function<double(double)> ic) : _final_time(
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
        plt::ylim(blim, tlim);
        plt::plot(_mesh.getFullCellPoints(), solution, "r*");
        //        plt::show();
        plt::pause(0.01);
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

namespace IC
{
    value_t shock_wave_RP(value_t x)
    {
        return x < 0.5 ? 2 : 1;
    }

    value_t rarefaction_RP(value_t x)
    {
        return x < 0.5 ? 2 : 1;
    }

    value_t cosine_IVP(value_t x)
    {
        return std::cos(2 * M_PI * x);
    }

}

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
    BurgersEquation<NumericalFlux::GOD> equation(n, mesh.getDx(), C);

    Simulation sim(final_time, equation, mesh, IC::rarefaction_RP);
    sim.plot();
    while (sim.evolve())
    {
        sim.plot(0, 3);
    }

    return 0;
}
