#include "pricer/Tridiagonal.hpp"

#include <cmath>
#include <stdexcept>
#include <string>

namespace pricer {

void solveTridiagonal(const std::vector<double>& lower, const std::vector<double>& diag,
                      const std::vector<double>& upper, const std::vector<double>& d,
                      std::vector<double>& x) {
    const std::size_t n = d.size();
    if (n == 0 || lower.size() != n || diag.size() != n || upper.size() != n) {
        throw std::invalid_argument("solveTridiagonal: the four vectors must have the same non-zero size");
    }

    // forward elimination: c[i] and x[i] hold the modified upper diagonal and right-hand side
    std::vector<double> c(n);
    x.assign(n, 0.0);
    double pivot = diag[0];
    for (std::size_t i = 0; i < n; ++i) {
        if (i > 0) pivot = diag[i] - lower[i] * c[i - 1];
        if (std::abs(pivot) < 1e-300) {
            throw std::runtime_error("solveTridiagonal: zero pivot at row " + std::to_string(i));
        }
        c[i] = (i + 1 < n) ? upper[i] / pivot : 0.0;
        x[i] = (d[i] - (i > 0 ? lower[i] * x[i - 1] : 0.0)) / pivot;
    }

    // back substitution
    for (std::size_t i = n - 1; i-- > 0;) {
        x[i] -= c[i] * x[i + 1];
    }
}

}  // namespace pricer
