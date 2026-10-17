#pragma once

namespace pricer {

enum class OptionType { Call, Put };

// Price and sensitivities of a European option.
struct Greeks {
    double price;
    double delta;  // dV/dS
    double gamma;  // d2V/dS2
    double theta;  // dV/dt, per year of calendar time (negative for a long call without dividends)
};

// Closed-form Black-Scholes values, used as the reference for the PDE solver.
Greeks blackScholes(OptionType type, double spot, double strike, double maturity, double rate,
                    double vol);

}  // namespace pricer
