#pragma once

namespace pricer {

// European payoff with the Dirichlet values the pricer imposes at the edges of the grid.
//
// tau is the time to maturity and r the risk-free rate. The boundary values are the
// asymptotic prices of the option when the spot tends to 0 or to infinity, evaluated at the
// finite edges of the grid.
class Payoff {
public:
    virtual ~Payoff() = default;

    virtual double operator()(double spot) const = 0;
    virtual double lowerBoundary(double spot, double tau, double r) const = 0;
    virtual double upperBoundary(double spot, double tau, double r) const = 0;

    // Strike, used to make sure the grid covers the region where the payoff is not linear.
    virtual double strike() const = 0;
};

class Call : public Payoff {
public:
    explicit Call(double strike);

    double operator()(double spot) const override;
    double lowerBoundary(double spot, double tau, double r) const override;  // 0
    double upperBoundary(double spot, double tau, double r) const override;  // S - K exp(-r tau)
    double strike() const override { return K_; }

private:
    double K_;
};

class Put : public Payoff {
public:
    explicit Put(double strike);

    double operator()(double spot) const override;
    double lowerBoundary(double spot, double tau, double r) const override;  // K exp(-r tau) - S
    double upperBoundary(double spot, double tau, double r) const override;  // 0
    double strike() const override { return K_; }

private:
    double K_;
};

}  // namespace pricer
