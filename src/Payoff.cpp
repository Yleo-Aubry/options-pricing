#include "pricer/Payoff.hpp"

#include <algorithm>
#include <cmath>
#include <stdexcept>

namespace pricer {

namespace {
double checkedStrike(double K) {
    if (!(K > 0.0)) throw std::invalid_argument("strike must be positive");
    return K;
}
}  // namespace

Call::Call(double strike) : K_(checkedStrike(strike)) {}

double Call::operator()(double spot) const { return std::max(spot - K_, 0.0); }

double Call::lowerBoundary(double, double, double) const { return 0.0; }

double Call::upperBoundary(double spot, double tau, double r) const {
    return spot - K_ * std::exp(-r * tau);
}

Put::Put(double strike) : K_(checkedStrike(strike)) {}

double Put::operator()(double spot) const { return std::max(K_ - spot, 0.0); }

double Put::lowerBoundary(double spot, double tau, double r) const {
    return K_ * std::exp(-r * tau) - spot;
}

double Put::upperBoundary(double, double, double) const { return 0.0; }

}  // namespace pricer
