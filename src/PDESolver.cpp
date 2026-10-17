#include "pricer/PDESolver.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

#include "pricer/Tridiagonal.hpp"

namespace pricer {

PDESolver::PDESolver(double maturity, double rate, double vol, GridParams grid)
    : T_(maturity), r_(rate), sigma_(vol), grid_(grid) {
    if (!(T_ > 0.0) || !(sigma_ > 0.0)) {
        throw std::invalid_argument("PDESolver: maturity and volatility must be positive");
    }
    if (grid_.spaceNodes < 5 || grid_.spaceNodes % 2 == 0) {
        throw std::invalid_argument("PDESolver: spaceNodes must be odd and at least 5");
    }
    if (grid_.timeSteps < 1) throw std::invalid_argument("PDESolver: timeSteps must be positive");
    if (grid_.theta < 0.0 || grid_.theta > 1.0) {
        throw std::invalid_argument("PDESolver: theta must lie in [0, 1]");
    }
    if (!(grid_.stdDevs > 0.0)) throw std::invalid_argument("PDESolver: stdDevs must be positive");
}

void PDESolver::step(std::vector<double>& V, const std::vector<double>& spots, const Payoff& payoff,
                     double tau, double dt, double theta, double dx) const {
    const std::size_t n = V.size();
    const double s2 = sigma_ * sigma_;
    const double nu = r_ - 0.5 * s2;

    // L V_i = l V_{i-1} + c V_i + u V_{i+1}, centred differences
    const double l = 0.5 * s2 / (dx * dx) - nu / (2.0 * dx);
    const double c = -s2 / (dx * dx) - r_;
    const double u = 0.5 * s2 / (dx * dx) + nu / (2.0 * dx);

    // (I - theta dt L) V^{new} = (I + (1 - theta) dt L) V^{old} on the interior nodes
    const std::size_t m = n - 2;
    const double e = (1.0 - theta) * dt;
    const double a = theta * dt;
    std::vector<double> lower(m, -a * l), diag(m, 1.0 - a * c), upper(m, -a * u), rhs(m);
    for (std::size_t i = 0; i < m; ++i) {
        rhs[i] = V[i + 1] + e * (l * V[i] + c * V[i + 1] + u * V[i + 2]);
    }

    const double left = payoff.lowerBoundary(spots.front(), tau + dt, r_);
    const double right = payoff.upperBoundary(spots.back(), tau + dt, r_);
    rhs.front() += a * l * left;
    rhs.back() += a * u * right;

    std::vector<double> interior;
    solveTridiagonal(lower, diag, upper, rhs, interior);
    std::copy(interior.begin(), interior.end(), V.begin() + 1);
    V.front() = left;
    V.back() = right;
}

std::vector<double> PDESolver::solveGrid(const Payoff& payoff, double spot, std::vector<double>& spots) const {
    if (!(spot > 0.0)) throw std::invalid_argument("PDESolver: spot must be positive");
    const std::size_t n = grid_.spaceNodes;
    const std::size_t mid = n / 2;

    // half-width: stdDevs standard deviations of ln S_T, and at least 1.5 times the distance to the strike
    const double halfWidth = std::max(grid_.stdDevs * sigma_ * std::sqrt(T_),
                                      1.5 * std::abs(std::log(payoff.strike() / spot)));
    const double dx = halfWidth / static_cast<double>(mid);

    spots.resize(n);
    std::vector<double> V(n);
    for (std::size_t i = 0; i < n; ++i) {
        spots[i] = spot * std::exp((static_cast<double>(i) - static_cast<double>(mid)) * dx);
        V[i] = payoff(spots[i]);
    }

    // Rannacher start: the payoff has a kink at the strike, which Crank-Nicolson propagates as
    // spurious oscillations; the first steps are replaced by two implicit half-steps each.
    const double dt = T_ / static_cast<double>(grid_.timeSteps);
    const std::size_t smoothed = grid_.theta < 1.0 ? std::min(grid_.rannacherSteps, grid_.timeSteps) : 0;
    double tau = 0.0;
    for (std::size_t k = 0; k < grid_.timeSteps; ++k) {
        if (k < smoothed) {
            step(V, spots, payoff, tau, 0.5 * dt, 1.0, dx);
            step(V, spots, payoff, tau + 0.5 * dt, 0.5 * dt, 1.0, dx);
        } else {
            step(V, spots, payoff, tau, dt, grid_.theta, dx);
        }
        tau += dt;
    }
    return V;
}

Greeks PDESolver::solve(const Payoff& payoff, double spot) const {
    std::vector<double> spots;
    const std::vector<double> V = solveGrid(payoff, spot, spots);
    const std::size_t i = V.size() / 2;
    const double dx = std::log(spots[i + 1] / spots[i]);

    const double Vx = (V[i + 1] - V[i - 1]) / (2.0 * dx);
    const double Vxx = (V[i + 1] - 2.0 * V[i] + V[i - 1]) / (dx * dx);
    const double delta = Vx / spot;
    const double gamma = (Vxx - Vx) / (spot * spot);
    const double theta = r_ * V[i] - r_ * spot * delta - 0.5 * sigma_ * sigma_ * spot * spot * gamma;
    return {V[i], delta, gamma, theta};
}

}  // namespace pricer
