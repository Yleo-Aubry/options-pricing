// Accuracy and convergence of the PDE solver against the closed form, printed as CSV.
//
// Usage:
//   benchmark cases        20 calls and 20 puts on the default grid
//   benchmark convergence  at-the-money call, grid refined by a factor 2 at each line

#include <chrono>
#include <cmath>
#include <cstdio>
#include <string>
#include <vector>

#include "pricer/BlackScholes.hpp"
#include "pricer/PDESolver.hpp"

using pricer::OptionType;

namespace {

struct Case {
    double S, K, T, r, sigma;
};

// the 20 parameter sets of the original test suite
const std::vector<Case> kCases = {
    {100, 100, 1.0, 0.05, 0.20}, {80, 100, 1.0, 0.05, 0.20},  {120, 100, 1.0, 0.05, 0.20},
    {100, 90, 1.0, 0.05, 0.20},  {100, 110, 1.0, 0.05, 0.20}, {100, 100, 0.5, 0.05, 0.20},
    {100, 100, 2.0, 0.05, 0.20}, {100, 100, 1.0, 0.01, 0.20}, {100, 100, 1.0, 0.10, 0.20},
    {100, 100, 1.0, 0.05, 0.10}, {100, 100, 1.0, 0.05, 0.40}, {95, 100, 1.0, 0.03, 0.25},
    {105, 100, 1.0, 0.03, 0.25}, {90, 100, 2.0, 0.02, 0.30},  {110, 100, 2.0, 0.02, 0.30},
    {70, 100, 1.5, 0.04, 0.35},  {130, 100, 1.5, 0.04, 0.35}, {150, 100, 0.75, 0.06, 0.18},
    {100, 120, 1.5, 0.04, 0.35}, {100, 80, 1.5, 0.04, 0.35},
};

void cases() {
    std::printf("type,S,K,T,r,sigma,price_bs,price_pde,err_price,err_delta,err_gamma,err_theta\n");
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        for (const Case& c : kCases) {
            const pricer::PDESolver solver(c.T, c.r, c.sigma);
            const pricer::Call call(c.K);
            const pricer::Put put(c.K);
            const pricer::Payoff& payoff = type == OptionType::Call ? static_cast<const pricer::Payoff&>(call) : put;
            const pricer::Greeks p = solver.solve(payoff, c.S);
            const pricer::Greeks b = pricer::blackScholes(type, c.S, c.K, c.T, c.r, c.sigma);
            std::printf("%s,%g,%g,%g,%g,%g,%.6f,%.6f,%.2e,%.2e,%.2e,%.2e\n",
                        type == OptionType::Call ? "call" : "put", c.S, c.K, c.T, c.r, c.sigma, b.price, p.price,
                        std::abs(p.price - b.price), std::abs(p.delta - b.delta), std::abs(p.gamma - b.gamma),
                        std::abs(p.theta - b.theta));
        }
    }
}

void convergence() {
    const double S = 100, K = 100, T = 1, r = 0.05, sigma = 0.2;
    const pricer::Call call(K);
    const pricer::Greeks b = pricer::blackScholes(OptionType::Call, S, K, T, r, sigma);
    std::printf("scheme,nodes,steps,err_price,ratio,err_delta,err_gamma,time_ms\n");
    for (double theta : {0.5, 1.0}) {
        double previous = 0.0;
        for (std::size_t m = 50; m <= 1600; m *= 2) {
            pricer::GridParams g;
            g.spaceNodes = 2 * m + 1;
            g.timeSteps = m;
            g.theta = theta;
            const auto start = std::chrono::steady_clock::now();
            const pricer::Greeks p = pricer::PDESolver(T, r, sigma, g).solve(call, S);
            const double ms = std::chrono::duration<double, std::milli>(std::chrono::steady_clock::now() - start).count();
            const double err = std::abs(p.price - b.price);
            std::printf("%s,%zu,%zu,%.3e,%s,%.3e,%.3e,%.1f\n", theta == 0.5 ? "crank-nicolson" : "implicit",
                        g.spaceNodes, g.timeSteps, err,
                        previous > 0 ? std::to_string(previous / err).substr(0, 4).c_str() : "", std::abs(p.delta - b.delta),
                        std::abs(p.gamma - b.gamma), ms);
            previous = err;
        }
    }
}

}  // namespace

int main(int argc, char** argv) {
    const std::string mode = argc > 1 ? argv[1] : "cases";
    if (mode == "cases") cases();
    else if (mode == "convergence") convergence();
    else {
        std::fprintf(stderr, "usage: benchmark [cases|convergence]\n");
        return 1;
    }
    return 0;
}
