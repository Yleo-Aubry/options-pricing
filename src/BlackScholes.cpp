#include "pricer/BlackScholes.hpp"

#include <cmath>
#include <stdexcept>

namespace pricer {

namespace {
constexpr double kInvSqrt2Pi = 0.39894228040143267794;

double normCdf(double x) { return 0.5 * std::erfc(-x / std::sqrt(2.0)); }
double normPdf(double x) { return kInvSqrt2Pi * std::exp(-0.5 * x * x); }
}  // namespace

Greeks blackScholes(OptionType type, double S, double K, double T, double r, double sigma) {
    if (!(S > 0.0) || !(K > 0.0) || !(T > 0.0) || !(sigma > 0.0)) {
        throw std::invalid_argument("blackScholes: spot, strike, maturity and volatility must be positive");
    }
    const double sqrtT = std::sqrt(T);
    const double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * sqrtT);
    const double d2 = d1 - sigma * sqrtT;
    const double df = std::exp(-r * T);

    const double gamma = normPdf(d1) / (S * sigma * sqrtT);
    const double timeDecay = -S * normPdf(d1) * sigma / (2.0 * sqrtT);

    if (type == OptionType::Call) {
        return {S * normCdf(d1) - K * df * normCdf(d2), normCdf(d1), gamma,
                timeDecay - r * K * df * normCdf(d2)};
    }
    return {K * df * normCdf(-d2) - S * normCdf(-d1), normCdf(d1) - 1.0, gamma,
            timeDecay + r * K * df * normCdf(-d2)};
}

}  // namespace pricer
