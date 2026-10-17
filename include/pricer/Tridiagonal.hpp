#pragma once

#include <vector>

namespace pricer {

// Solves the tridiagonal system A x = d with the Thomas algorithm, in O(n) operations.
//
// lower, diag and upper are the sub-, main and super-diagonals of A, all of size n;
// lower[0] and upper[n - 1] are not used. x is resized to n.
//
// The algorithm performs Gaussian elimination without pivoting, which is stable when A is
// diagonally dominant; this is the case for the matrices of the pricer.
//
// Throws std::invalid_argument if the sizes differ and std::runtime_error on a zero pivot.
void solveTridiagonal(const std::vector<double>& lower, const std::vector<double>& diag,
                      const std::vector<double>& upper, const std::vector<double>& d,
                      std::vector<double>& x);

}  // namespace pricer
