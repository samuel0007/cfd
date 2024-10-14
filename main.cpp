//
// Created by russo on 14.10.2024.
//

#include <iostream>
#include <utility>
#include <vector>
#include "matplotlibcpp.h"

namespace plt = matplotlibcpp;

typedef double value_t;


enum SpaceDiscretisation {
    Central, Backward, Forward
};

template<SpaceDiscretisation T>
class AdvectionEquation {
public:
    explicit AdvectionEquation(value_t a, int n, double min_x, double max_x) : _a(a), _n(n), _min_x(min_x),
                                                                               _max_x(max_x),
                                                                               _dx((max_x - min_x) / (double) n) {
        // Init mesh
        this->_mesh.resize(n);
        this->_l = (_max_x - _min_x);
        std::cout << this->_l << std::endl;
        for (int i = 0; i < n; ++i) {
            this->_mesh[i] = this->_l * i / (double) n;
        }

    };

    inline std::vector<value_t> dudt(const std::vector<value_t> &u) {
        std::vector<value_t> res(this->_n);

        if constexpr(T == SpaceDiscretisation::Central) {
            res[0] = -this->_a * (u[1] - u[this->_n - 1]) / (2 * this->_dx);
            for (int i = 1; i < this->_n - 1; ++i) {
                res[i] = -this->_a * (u[(i + 1)] - u[(i - 1)]) / (2 * this->_dx);
            }
            res[this->_n - 1] = -this->_a * (u[0] - u[this->_n - 2]) / (2 * this->_dx);
        } else if constexpr(T == SpaceDiscretisation::Backward) {
            for (int i = 1; i < this->_n - 1; ++i) {
                res[i] = -this->_a * (u[i + 1] - u[i - 1]) / this->_dx;
            }
        } else if constexpr(T == SpaceDiscretisation::Forward) {
            for (int i = 1; i < this->_n - 1; ++i) {
                res[i] = -this->_a * (u[i + 1] - u[i - 1]) / this->_dx;
            }
        }
        return res;
    }

    const std::vector<value_t> &get_mesh() {
        return _mesh;
    }


private:
    value_t _a;
    double _dx;
    double _min_x;
    double _max_x;
    double _l{};
    int _n;
    std::vector<value_t> _mesh;

};

template<SpaceDiscretisation T>
class Simulation {
public:
    Simulation(double final_time, double dt, AdvectionEquation<T> eq, std::function<double(double)> ic) : _final_time(
            final_time), _dt(dt), _eq(std::move(eq)) {
        std::vector<double> mesh = _eq.get_mesh();
        solution.resize(mesh.size());
        std::transform(mesh.begin(), mesh.end(), solution.begin(), [&ic](double x) { return ic(x); });
    };


    bool evolve() {
        double rest_dt = (_final_time - _ct);
        double min_dt = std::min(rest_dt, _dt);

        std::vector<double> dudt = _eq.dudt(solution);
        std::transform(solution.begin(), solution.end(), dudt.begin(), solution.begin(), [&](double u, double dudt) {
            return u + dudt * min_dt;
        });
        _ct += min_dt;

        return _ct < _final_time;
    };

    void plot() {
        plt::clf();
        plt::plot(_eq.get_mesh(), solution);
//        plt::show();
        plt::pause(0.001);
    }

    std::vector<value_t> solution;

private:
    double _dt;
    double _ct = 0.;
    double _final_time;
    AdvectionEquation<T> _eq;
};

int main(int argc, char *argv[]) {
    value_t a = 0.2;
    int n = 100;
    double min_x = 0;
    double max_x = 2 * M_PI;
    double final_time = 10.0;
    double dt = 0.05;

    AdvectionEquation<SpaceDiscretisation::Central> equation(a, n, min_x, max_x);
    Simulation sim(final_time, dt, equation, [](double x) { return std::sin(x); });
    sim.plot();
    while (sim.evolve()) {
        sim.plot();
    }

    return 0;
}
