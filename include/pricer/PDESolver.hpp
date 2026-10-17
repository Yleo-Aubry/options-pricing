#pragma once

#include <cstddef>
#include <vector>

#include "pricer/BlackScholes.hpp"
#include "pricer/Payoff.hpp"

namespace pricer {

// Discretisation parameters.
struct GridParams {
    std::size_t spaceNodes = 401;    // number of nodes in x = ln S, odd so that the spot is a node
    std::size_t timeSteps = 400;     // number of time steps between maturity and today
    double theta = 0.5;              // 0 explicit, 0.5 Crank-Nicolson, 1 fully implicit
    double stdDevs = 6.0;            // half-width of the grid in standard deviations sigma sqrt(T)
    std::size_t rannacherSteps = 2;  // initial steps replaced by two implicit half-steps each
};

// Finite-difference solver of the Black-Scholes PDE in log-spot.
//
// With x = ln S and tau = T - t, the price V(x, tau) solves the constant-coefficient equation
//     dV/dtau = sigma^2 / 2 d2V/dx2 + (r - sigma^2 / 2) dV/dx - r V,   V(x, 0) = payoff(e^x).
// The grid is uniform in x and centred on ln S0, so the price and the Greeks are read at a
// node without interpolation.
class PDESolver {
public:
    PDESolver(double maturity, double rate, double vol, GridParams grid = {});

    // Price, delta and gamma from the grid; theta from the PDE itself,
    // theta = r V - r S delta - sigma^2 S^2 gamma / 2.
    Greeks solve(const Payoff& payoff, double spot) const;

    // Solution at tau = T on the whole grid, for inspection; spots are e^{x_i}.
    std::vector<double> solveGrid(const Payoff& payoff, double spot, std::vector<double>& spots) const;

private:
    double T_, r_, sigma_;
    GridParams grid_;

    // one step of the theta-scheme from tau to tau + dt, in place on V
    void step(std::vector<double>& V, const std::vector<double>& spots, const Payoff& payoff,
              double tau, double dt, double theta, double dx) const;
};

}  // namespace pricer
