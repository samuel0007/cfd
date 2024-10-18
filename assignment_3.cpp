//
// Created by russo on 14.10.2024.
//

#include <iostream>
#include <utility>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef double value_t;


enum NumericalFlux {
    Central, Backward, Forward
};

template <NumericalFlux T>
class ConservationEquation {
public:
    virtual value_t f(value_t u) = 0;

    const std::vector<value_t>& getFlux(const std::vector<value_t>& u) {
        // Evaluate f(u) of the conservation form
        _flux[0] = this->_getNumericalFlux(u[this->_nCells-1], u[0]);
        for(size_t i = 1; i < this->_nCells; ++i) {
            _flux[i] = this->_getNumericalFlux(u[i-1], u[i]);
        }
        _flux[this->_nCells] = _flux[0];
        return _flux;
    }

    size_t _nCells;
    std::vector<value_t> _flux;

private:
    // Default is CD
    value_t _getNumericalFlux(value_t u0, value_t u1){
        return 0.5*(f(u0)+f(u1));
    };
};

template <>
value_t ConservationEquation<NumericalFlux::Forward>::_getNumericalFlux(value_t u0, value_t u1) {
    return f(u0);
}

template <>
value_t ConservationEquation<NumericalFlux::Central>::_getNumericalFlux(value_t u0, value_t u1) {
    return 0.5*(f(u0)+f(u1));
}

template <NumericalFlux T>
class AdvectionEquation: public ConservationEquation<T> {
public:
    explicit AdvectionEquation(value_t a, size_t nCells) : _a(a) {
        this->_nCells = nCells;
        this->_flux.resize(nCells+1);
    };
    
    value_t f(value_t u) {
        return this->_a*u;
    }

private:
    value_t _a;
};

class Mesh {
public:
    // n is the number of cells 
    explicit Mesh(int n, double min_x, double max_x) : _n(n), _min_x(min_x), _max_x(max_x), _dx((max_x - min_x) / (double) n) {
        const size_t nPoints = n + 1;
        this->_l = (_max_x - _min_x);

        this->_points.resize(nPoints);
        for (int i = 0; i < nPoints; ++i) {
            this->_points[i] = this->_l * i / (double) nPoints;
        }

        this->_cellPoints.resize(n);
        for (int i = 0; i < n; ++i) {
            this->_cellPoints[i] = this->_points[i] + _dx * 0.5;
        }

    }

    const std::vector<value_t>& getPoints() {
        return _points;
    }

    const std::vector<value_t>& getCellPoints() {
        return _cellPoints;
    }

    const double& getDx() {
        return _dx;
    }

private:
    double _dx;
    double _min_x;
    double _max_x;
    double _l{};
    int _n;
    std::vector<value_t> _points;
    std::vector<value_t> _cellPoints;

};

template <NumericalFlux T>
class Simulation {
public:
    Simulation(double final_time, double dt, ConservationEquation<T>& eq, Mesh mesh, std::function<double(double)> ic) : _final_time(
            final_time), _dt(dt), _eq(eq), _mesh(mesh) {
        const std::vector<value_t>& mesh_points = mesh.getPoints();
        const std::vector<value_t>& mesh_cellPoints = mesh.getCellPoints();

        this->solution.resize(mesh_cellPoints.size());
        this->_solution_buffer.resize(mesh_cellPoints.size());
        std::ranges::transform(mesh_cellPoints, solution.begin(), [&ic](double x) { return ic(x); });
    };


    bool evolve() {
        double rest_dt = (_final_time - _ct);
        double min_dt = std::min(rest_dt, _dt);

        const std::vector<double>& flux = this->_eq.getFlux(solution);
        size_t nCells = this->_mesh.getCellPoints().size();
        const auto& dx = this->_mesh.getDx();
        for(size_t i = 0; i < nCells; ++i) {
            this->_solution_buffer[i] = solution[i] - (_dt / dx) * (flux[i+1] - flux[i]);
        }
        this->solution = _solution_buffer;
        _ct += min_dt;

        return _ct < _final_time;
    };

    void plot() {
        plt::clf();
        plt::plot(_mesh.getCellPoints(), solution);
//        plt::show();
        plt::pause(0.001);
    }

    std::vector<value_t> solution;

private:
    double _dt;
    double _ct = 0.;
    double _final_time;
    ConservationEquation<T>& _eq;
    Mesh _mesh;
    std::vector<value_t> _solution_buffer;
};

int main(int argc, char *argv[]) {
    value_t a = 0.2;
    int n = 100;
    double min_x = 0;
    double max_x = 2 * M_PI;
    double final_time = 10.0;
    double dt = 0.05;

    AdvectionEquation<NumericalFlux::Forward> equation(a, n);
    Mesh mesh(n, min_x, max_x);
    Simulation sim(final_time, dt, equation, mesh, [](double x) { return std::sin(x); });
    sim.plot();
    while (sim.evolve()) {
        sim.plot();
    }

    return 0;
}
