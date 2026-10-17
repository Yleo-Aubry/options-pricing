// Unit tests, run by CTest. Each test returns normally or records a failure with CHECK.

#include <cmath>
#include <cstdio>
#include <functional>
#include <random>
#include <stdexcept>
#include <string>
#include <vector>

#include "pricer/BlackScholes.hpp"
#include "pricer/PDESolver.hpp"
#include "pricer/Payoff.hpp"
#include "pricer/Tridiagonal.hpp"

using namespace pricer;

namespace {

int failures = 0;

#define CHECK(cond)                                                              \
    do {                                                                         \
        if (!(cond)) {                                                           \
            std::printf("  FAILED %s:%d: %s\n", __FILE__, __LINE__, #cond);      \
            ++failures;                                                          \
        }                                                                        \
    } while (0)

#define CHECK_NEAR(a, b, tol)                                                                    \
    do {                                                                                         \
        const double a_ = (a), b_ = (b);                                                         \
        if (!(std::abs(a_ - b_) <= (tol))) {                                                     \
            std::printf("  FAILED %s:%d: |%s - %s| = %.3e > %.1e\n", __FILE__, __LINE__, #a, #b, \
                        std::abs(a_ - b_), static_cast<double>(tol));                            \
            ++failures;                                                                          \
        }                                                                                        \
    } while (0)

template <class E>
bool throws(const std::function<void()>& f) {
    try {
        f();
    } catch (const E&) {
        return true;
    }
    return false;
}

struct Case {
    double S, K, T, r, sigma;
};

const std::vector<Case> kCases = {
    {100, 100, 1.0, 0.05, 0.20}, {80, 100, 1.0, 0.05, 0.20},  {120, 100, 1.0, 0.05, 0.20},
    {100, 90, 1.0, 0.05, 0.20},  {100, 110, 1.0, 0.05, 0.20}, {100, 100, 0.5, 0.05, 0.20},
    {100, 100, 2.0, 0.05, 0.20}, {100, 100, 1.0, 0.01, 0.20}, {100, 100, 1.0, 0.10, 0.20},
    {100, 100, 1.0, 0.05, 0.10}, {100, 100, 1.0, 0.05, 0.40}, {95, 100, 1.0, 0.03, 0.25},
    {105, 100, 1.0, 0.03, 0.25}, {90, 100, 2.0, 0.02, 0.30},  {110, 100, 2.0, 0.02, 0.30},
    {70, 100, 1.5, 0.04, 0.35},  {130, 100, 1.5, 0.04, 0.35}, {150, 100, 0.75, 0.06, 0.18},
    {100, 120, 1.5, 0.04, 0.35}, {100, 80, 1.5, 0.04, 0.35},
};

// ---------------------------------------------------------------------------
// tridiagonal solver
// ---------------------------------------------------------------------------
void testTridiagonalRandomSystems() {
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> unif(-1.0, 1.0);
    for (std::size_t n : {1, 2, 3, 10, 200}) {
        std::vector<double> a(n), b(n), c(n), x_true(n), d(n, 0.0), x;
        for (std::size_t i = 0; i < n; ++i) {
            a[i] = unif(gen);
            c[i] = unif(gen);
            b[i] = 2.5 + std::abs(unif(gen));  // diagonally dominant
            x_true[i] = unif(gen);
        }
        for (std::size_t i = 0; i < n; ++i) {
            d[i] = b[i] * x_true[i] + (i > 0 ? a[i] * x_true[i - 1] : 0.0) + (i + 1 < n ? c[i] * x_true[i + 1] : 0.0);
        }
        solveTridiagonal(a, b, c, d, x);
        CHECK(x.size() == n);
        for (std::size_t i = 0; i < n; ++i) CHECK_NEAR(x[i], x_true[i], 1e-12);
    }
}

void testTridiagonalErrors() {
    std::vector<double> x;
    CHECK(throws<std::invalid_argument>([&] { solveTridiagonal({0, 1}, {1, 1}, {1}, {1, 1}, x); }));
    CHECK(throws<std::invalid_argument>([&] { solveTridiagonal({}, {}, {}, {}, x); }));
    // [[1, 1], [1, 1]] is singular: the second pivot is 0
    CHECK(throws<std::runtime_error>([&] { solveTridiagonal({0, 1}, {1, 1}, {1, 0}, {1, 2}, x); }));
}

// ---------------------------------------------------------------------------
// payoffs and closed form
// ---------------------------------------------------------------------------
void testPayoffs() {
    const Call call(100);
    const Put put(100);
    CHECK(call(120) == 20 && call(80) == 0);
    CHECK(put(80) == 20 && put(120) == 0);
    // boundary values satisfy put-call parity: C - P = S - K exp(-r tau)
    const double tau = 0.7, r = 0.03;
    for (double S : {0.5, 1000.0}) {
        const double c = S < 100 ? call.lowerBoundary(S, tau, r) : call.upperBoundary(S, tau, r);
        const double p = S < 100 ? put.lowerBoundary(S, tau, r) : put.upperBoundary(S, tau, r);
        CHECK_NEAR(c - p, S - 100 * std::exp(-r * tau), 1e-12);
    }
    CHECK(throws<std::invalid_argument>([] { Call bad(0.0); }));
}

void testBlackScholes() {
    // reference values for S = K = 100, T = 1, r = 5%, sigma = 20%
    const Greeks c = blackScholes(OptionType::Call, 100, 100, 1, 0.05, 0.2);
    CHECK_NEAR(c.price, 10.450583572185565, 1e-12);
    CHECK_NEAR(c.delta, 0.636830651175619, 1e-12);
    const Greeks p = blackScholes(OptionType::Put, 100, 100, 1, 0.05, 0.2);
    CHECK_NEAR(p.price, 5.573526022256971, 1e-12);

    for (const Case& k : kCases) {
        const Greeks call = blackScholes(OptionType::Call, k.S, k.K, k.T, k.r, k.sigma);
        const Greeks put = blackScholes(OptionType::Put, k.S, k.K, k.T, k.r, k.sigma);
        CHECK_NEAR(call.price - put.price, k.S - k.K * std::exp(-k.r * k.T), 1e-10);
        CHECK_NEAR(call.delta - put.delta, 1.0, 1e-12);
        CHECK_NEAR(call.gamma, put.gamma, 1e-14);
        // the closed form solves the PDE: theta + r S delta + sigma^2 S^2 gamma / 2 = r V
        CHECK_NEAR(call.theta + k.r * k.S * call.delta + 0.5 * k.sigma * k.sigma * k.S * k.S * call.gamma,
                   k.r * call.price, 1e-10);
    }
}

// ---------------------------------------------------------------------------
// PDE solver
// ---------------------------------------------------------------------------
void testPDEAgainstClosedForm() {
    for (OptionType type : {OptionType::Call, OptionType::Put}) {
        for (const Case& k : kCases) {
            const PDESolver solver(k.T, k.r, k.sigma);
            const Call call(k.K);
            const Put put(k.K);
            const Payoff& payoff = type == OptionType::Call ? static_cast<const Payoff&>(call) : put;
            const Greeks pde = solver.solve(payoff, k.S);
            const Greeks bs = blackScholes(type, k.S, k.K, k.T, k.r, k.sigma);
            CHECK_NEAR(pde.price, bs.price, 3e-3);
            CHECK_NEAR(pde.delta, bs.delta, 1e-4);
            CHECK_NEAR(pde.gamma, bs.gamma, 1e-5);
            CHECK_NEAR(pde.theta, bs.theta, 2e-3);
        }
    }
}

void testPDEPutCallParity() {
    const double S = 93, K = 105, T = 1.3, r = 0.04, sigma = 0.27;
    const PDESolver solver(T, r, sigma);
    const Greeks c = solver.solve(Call(K), S);
    const Greeks p = solver.solve(Put(K), S);
    // C - P solves the scheme with the linear payoff S - K, which has no kink: parity holds up to
    // the discretisation error of a forward contract, far below the error on each option
    CHECK_NEAR(c.price - p.price, S - K * std::exp(-r * T), 1e-4);
    CHECK_NEAR(c.delta - p.delta, 1.0, 5e-5);
}

double priceError(double theta, std::size_t m) {
    GridParams g;
    g.spaceNodes = 2 * m + 1;
    g.timeSteps = m;
    g.theta = theta;
    const double pde = PDESolver(1, 0.05, 0.2, g).solve(Call(100), 100).price;
    return std::abs(pde - blackScholes(OptionType::Call, 100, 100, 1, 0.05, 0.2).price);
}

void testConvergenceOrders() {
    // Crank-Nicolson with Rannacher start: second order, the error is divided by 4
    for (std::size_t m : {100, 200, 400}) {
        const double ratio = priceError(0.5, m) / priceError(0.5, 2 * m);
        CHECK(ratio > 3.8 && ratio < 4.2);
    }
    // fully implicit: first order in time, the error is divided by about 2
    const double ratio = priceError(1.0, 800) / priceError(1.0, 1600);
    CHECK(ratio > 1.9 && ratio < 2.2);
}

void testRannacherRemovesOscillations() {
    // with large time steps relative to dx, plain Crank-Nicolson oscillates around the strike
    GridParams g;
    g.spaceNodes = 801;
    g.timeSteps = 20;
    const double T = 0.25, r = 0.05, sigma = 0.2;
    const double exact = blackScholes(OptionType::Call, 100, 100, T, r, sigma).gamma;
    g.rannacherSteps = 0;
    const double plain = PDESolver(T, r, sigma, g).solve(Call(100), 100).gamma;
    g.rannacherSteps = 2;
    const double smoothed = PDESolver(T, r, sigma, g).solve(Call(100), 100).gamma;
    CHECK(std::abs(smoothed - exact) < 0.1 * std::abs(plain - exact));
    CHECK_NEAR(smoothed, exact, 1e-4);
}

void testSolutionIsMonotoneAndBounded() {
    std::vector<double> spots;
    const std::vector<double> V = PDESolver(1, 0.05, 0.2).solveGrid(Call(100), 100, spots);
    for (std::size_t i = 0; i < V.size(); ++i) {
        CHECK(V[i] >= -1e-10);                                   // no negative price
        CHECK(V[i] <= spots[i] + 1e-10);                         // a call is worth less than the stock
        CHECK(V[i] >= spots[i] - 100 * std::exp(-0.05) - 1e-8);  // and more than S - K exp(-rT)
        if (i > 0) CHECK(V[i] >= V[i - 1] - 1e-10);              // increasing in the spot
    }
    CHECK_NEAR(spots[V.size() / 2], 100.0, 1e-12);
}

void testInvalidParameters() {
    GridParams even;
    even.spaceNodes = 400;
    CHECK(throws<std::invalid_argument>([&] { PDESolver(1, 0.05, 0.2, even); }));
    GridParams badTheta;
    badTheta.theta = 1.5;
    CHECK(throws<std::invalid_argument>([&] { PDESolver(1, 0.05, 0.2, badTheta); }));
    CHECK(throws<std::invalid_argument>([] { PDESolver(0, 0.05, 0.2); }));
    CHECK(throws<std::invalid_argument>([] { PDESolver(1, 0.05, -0.2); }));
    CHECK(throws<std::invalid_argument>([] { PDESolver(1, 0.05, 0.2).solve(Call(100), 0); }));
}

}  // namespace

int main() {
    const std::vector<std::pair<std::string, void (*)()>> tests = {
        {"tridiagonal random systems", testTridiagonalRandomSystems},
        {"tridiagonal errors", testTridiagonalErrors},
        {"payoffs and boundaries", testPayoffs},
        {"closed form", testBlackScholes},
        {"PDE against closed form, 40 cases", testPDEAgainstClosedForm},
        {"PDE put-call parity", testPDEPutCallParity},
        {"convergence orders", testConvergenceOrders},
        {"Rannacher smoothing", testRannacherRemovesOscillations},
        {"monotone and bounded solution", testSolutionIsMonotoneAndBounded},
        {"invalid parameters", testInvalidParameters},
    };
    for (const auto& [name, test] : tests) {
        const int before = failures;
        test();
        std::printf("%-40s %s\n", name.c_str(), failures == before ? "ok" : "FAILED");
    }
    std::printf("\n%d failure(s)\n", failures);
    return failures == 0 ? 0 : 1;
}
