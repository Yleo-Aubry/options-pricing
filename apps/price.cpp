// Prices a European call or put with the PDE solver and compares it with the closed form.
//
// Usage:
//   price [--type call|put] [--spot S] [--strike K] [--maturity T] [--rate r] [--vol sigma]
//         [--space N] [--time M] [--theta theta]

#include <cmath>
#include <cstdio>
#include <memory>
#include <stdexcept>
#include <string>

#include "pricer/BlackScholes.hpp"
#include "pricer/PDESolver.hpp"

namespace {

void usage() {
    std::puts("usage: price [--type call|put] [--spot S] [--strike K] [--maturity T] [--rate r]\n"
              "             [--vol sigma] [--space N] [--time M] [--theta theta]\n"
              "defaults: call, S = 100, K = 100, T = 1, r = 0.05, sigma = 0.2, N = 401, M = 400, theta = 0.5");
}

}  // namespace

int main(int argc, char** argv) {
    std::string type = "call";
    double S = 100.0, K = 100.0, T = 1.0, r = 0.05, sigma = 0.2;
    pricer::GridParams grid;

    try {
        for (int i = 1; i < argc; ++i) {
            const std::string arg = argv[i];
            if (arg == "--help" || arg == "-h") {
                usage();
                return 0;
            }
            if (i + 1 >= argc) throw std::invalid_argument("missing value after " + arg);
            const std::string value = argv[++i];
            if (arg == "--type") type = value;
            else if (arg == "--spot") S = std::stod(value);
            else if (arg == "--strike") K = std::stod(value);
            else if (arg == "--maturity") T = std::stod(value);
            else if (arg == "--rate") r = std::stod(value);
            else if (arg == "--vol") sigma = std::stod(value);
            else if (arg == "--space") grid.spaceNodes = std::stoul(value);
            else if (arg == "--time") grid.timeSteps = std::stoul(value);
            else if (arg == "--theta") grid.theta = std::stod(value);
            else throw std::invalid_argument("unknown option " + arg);
        }
        if (type != "call" && type != "put") throw std::invalid_argument("--type must be call or put");

        const bool isCall = type == "call";
        std::unique_ptr<pricer::Payoff> payoff;
        if (isCall) payoff = std::make_unique<pricer::Call>(K);
        else payoff = std::make_unique<pricer::Put>(K);

        const pricer::PDESolver solver(T, r, sigma, grid);
        const pricer::Greeks pde = solver.solve(*payoff, S);
        const pricer::Greeks bs = pricer::blackScholes(isCall ? pricer::OptionType::Call : pricer::OptionType::Put,
                                                       S, K, T, r, sigma);

        std::printf("European %s, S = %g, K = %g, T = %g, r = %g, sigma = %g\n", type.c_str(), S, K, T, r, sigma);
        std::printf("grid: %zu nodes in ln S, %zu time steps, theta = %g\n\n", grid.spaceNodes, grid.timeSteps,
                    grid.theta);
        std::printf("%-8s %14s %14s %12s\n", "", "PDE", "closed form", "abs. error");
        const char* names[] = {"price", "delta", "gamma", "theta"};
        const double p[] = {pde.price, pde.delta, pde.gamma, pde.theta};
        const double b[] = {bs.price, bs.delta, bs.gamma, bs.theta};
        for (int k = 0; k < 4; ++k) {
            std::printf("%-8s %14.6f %14.6f %12.2e\n", names[k], p[k], b[k], std::abs(p[k] - b[k]));
        }
    } catch (const std::exception& e) {
        std::fprintf(stderr, "error: %s\n", e.what());
        usage();
        return 1;
    }
    return 0;
}
