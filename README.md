# European option pricing by finite differences

C++17 solver of the Black-Scholes partial differential equation for European calls and puts: a
theta-scheme on a log-spot grid, solved at each time step by the Thomas algorithm. It returns the price
and the delta, gamma and theta, and is validated against the closed-form formulas.

Course project in C++ (December 2025) by Yléo Aubry, Yepi Sylvian Jean-Mederick Bessi and Ayman
Graziana, revised since; the original report is in [`report/`](report).

On the default grid (401 nodes, 400 time steps, about 3 ms per price), the largest error over 40 test
options is $1.7 \cdot 10^{-3}$ on the price and $6 \cdot 10^{-5}$ on the delta, and the error is
divided by 4 each time the grid is refined by a factor 2.

## Problem

Under the Black-Scholes model, with constant rate $r$ and volatility $\sigma$, the price $V(S, t)$ of a
European option with payoff $h$ solves

```math
\frac{\partial V}{\partial t} + \frac{1}{2} \sigma^2 S^2 \frac{\partial^2 V}{\partial S^2} + r S \frac{\partial V}{\partial S} - r V = 0, \qquad V(S, T) = h(S),
```

with $h(S) = (S - K)^+$ for a call and $(K - S)^+$ for a put.

### Log-spot transform

The coefficients $\frac12 \sigma^2 S^2$ and $rS$ depend on $S$. With $x = \ln S$ and the time to maturity
$\tau = T - t$, the chain rule

```math
\frac{\partial V}{\partial S} = \frac{1}{S} \frac{\partial V}{\partial x}, \qquad \frac{\partial^2 V}{\partial S^2} = \frac{1}{S^2} \left( \frac{\partial^2 V}{\partial x^2} - \frac{\partial V}{\partial x} \right)
```

turns the equation into a convection-diffusion equation with constant coefficients, solved forward in
$\tau$:

```math
\frac{\partial V}{\partial \tau} = \frac{\sigma^2}{2} \frac{\partial^2 V}{\partial x^2} + \nu \frac{\partial V}{\partial x} - r V, \qquad \nu = r - \frac{\sigma^2}{2}, \qquad V(x, 0) = h(e^x).
```

The matrices of the scheme are then the same on every row of the grid, and a uniform step in $x$ is a
constant relative step in $S$, which matches the scale of the problem.

## Discretisation

### Grid

The grid has $N$ (odd) nodes $x_i = \ln S_0 + (i - i_0) \Delta x$, centred on the spot: $x_{i_0} = \ln S_0$.
Its half-width is

```math
L = \max\left(6 \, \sigma \sqrt{T}, \; 1.5 \, \lvert \ln(K / S_0) \rvert\right),
```

six standard deviations of $\ln S_T$, extended if needed so that the strike lies well inside the grid.
Placing the spot on a node means that the price and the Greeks are read on the grid, without
interpolation.

### Theta-scheme

With centred differences, the discrete operator is

```math
(L_h V)_i = \frac{\sigma^2}{2} \frac{V_{i+1} - 2 V_i + V_{i-1}}{\Delta x^2} + \nu \frac{V_{i+1} - V_{i-1}}{2 \Delta x} - r V_i,
```

and one step from $\tau_n$ to $\tau_{n+1} = \tau_n + \Delta\tau$ solves

```math
\left(I - \theta \, \Delta\tau \, L_h\right) V^{n+1} = \left(I + (1 - \theta) \, \Delta\tau \, L_h\right) V^n
```

on the interior nodes, a tridiagonal system. $\theta = 1/2$ is the Crank-Nicolson scheme, second order
in time and space and unconditionally stable; $\theta = 1$ is the fully implicit scheme, first order in
time; $\theta = 0$ is the explicit scheme, stable only for $\Delta\tau \lesssim \Delta x^2 / \sigma^2$.

Each system is solved by the Thomas algorithm, Gaussian elimination specialised to tridiagonal
matrices, in $O(N)$ operations instead of $O(N^3)$. The matrix $I - \theta \, \Delta\tau \, L_h$ is
diagonally dominant as soon as $\lvert \nu \rvert \Delta x < \sigma^2$, which makes the elimination
stable without pivoting.

### Boundary conditions

The values at the edges $S_{\min} = e^{x_0}$ and $S_{\max} = e^{x_{N-1}}$ are imposed (Dirichlet), using
the asymptotic prices of each option:

| | $S \to 0$ | $S \to \infty$ |
|---|---|---|
| call | $0$ | $S - K e^{-r\tau}$ |
| put | $K e^{-r\tau} - S$ | $0$ |

Each payoff class provides its own boundary values, so that the solver does not depend on the option
type.

### Rannacher start

The payoff has a kink at the strike. Crank-Nicolson damps the high-frequency components of the error
very little, and propagates this kink as oscillations that pollute the gamma near the money. The first
two time steps are therefore each replaced by two fully implicit half-steps (Rannacher smoothing),
which damp these components; the following steps use the chosen $\theta$. The global order of the
scheme is unchanged.

### Greeks

At the spot node $i_0$, with centred differences in $x$,

```math
\Delta = \frac{1}{S_0} \frac{V_{i_0+1} - V_{i_0-1}}{2 \Delta x}, \qquad \Gamma = \frac{1}{S_0^2} \left( \frac{V_{i_0+1} - 2 V_{i_0} + V_{i_0-1}}{\Delta x^2} - \frac{V_{i_0+1} - V_{i_0-1}}{2 \Delta x} \right),
```

both second order in $\Delta x$. The theta (sensitivity to calendar time, per year) is obtained from the
equation itself rather than from a time difference:

```math
\Theta = \frac{\partial V}{\partial t} = r V - r S_0 \Delta - \frac{1}{2} \sigma^2 S_0^2 \Gamma.
```

## Results

All errors are measured against the closed-form Black-Scholes price and Greeks.

### Convergence

At-the-money call ($S_0 = K = 100$, $T = 1$, $r = 5\%$, $\sigma = 20\%$), with $N = 2M + 1$ nodes and $M$
time steps:

| $M$ | Crank-Nicolson, error | ratio | implicit, error | ratio | delta error (CN) |
|---|---|---|---|---|---|
| 50 | $1.33 \cdot 10^{-2}$ | | $3.41 \cdot 10^{-2}$ | | $9.7 \cdot 10^{-5}$ |
| 100 | $3.33 \cdot 10^{-3}$ | 4.00 | $1.38 \cdot 10^{-2}$ | 2.47 | $2.4 \cdot 10^{-5}$ |
| 200 | $8.32 \cdot 10^{-4}$ | 4.00 | $6.06 \cdot 10^{-3}$ | 2.26 | $6.1 \cdot 10^{-6}$ |
| 400 | $2.08 \cdot 10^{-4}$ | 4.00 | $2.83 \cdot 10^{-3}$ | 2.14 | $1.5 \cdot 10^{-6}$ |
| 800 | $5.20 \cdot 10^{-5}$ | 4.00 | $1.36 \cdot 10^{-3}$ | 2.07 | $3.8 \cdot 10^{-7}$ |
| 1600 | $1.30 \cdot 10^{-5}$ | 4.00 | $6.68 \cdot 10^{-4}$ | 2.03 | $9.5 \cdot 10^{-8}$ |

Crank-Nicolson converges exactly at order 2: the error is divided by 4 when both steps are halved. The
implicit scheme tends to order 1, limited by its time error. The cost is linear in the number of nodes
per time step: the finest Crank-Nicolson grid (3,201 nodes, 1,600 steps) takes about 0.25 s.

### Test options

20 parameter sets ($S_0$ from 70 to 150, $K$ from 80 to 120, $T$ from 0.5 to 2 years, $r$ from 1% to 10%,
$\sigma$ from 10% to 40%), each priced as a call and as a put on the default grid. Largest absolute error:

| | price | delta | gamma | theta |
|---|---|---|---|---|
| calls | $1.7 \cdot 10^{-3}$ | $5.8 \cdot 10^{-5}$ | $5.3 \cdot 10^{-6}$ | $9.8 \cdot 10^{-4}$ |
| puts | $1.7 \cdot 10^{-3}$ | $3.0 \cdot 10^{-5}$ | $5.3 \cdot 10^{-6}$ | $9.6 \cdot 10^{-4}$ |

The largest errors occur for the highest volatility ($\sigma = 40\%$), where the grid is widest and
$\Delta x$ largest. The full table is printed by `benchmark cases`.

### Changes from the course version

The version described in the report used a fixed grid from $S_{\max}/3000$ to $S_{\max}$, the boundary
values of a call for every option, the Crank-Nicolson scheme only, linear interpolation at the spot and
a one-sided difference for the delta. On the same 20 calls, its largest price error was
$3.2 \cdot 10^{-2}$, and at the money its delta was 0.6008 against 0.6368. The present version adds:

- the general theta-scheme, with $\theta$ chosen by the user (the course version asked for it but did
  not use it);
- boundary values specific to each payoff, so that puts get their own asymptotic values;
- the grid centred on the spot and the centred Greeks, which remove the interpolation error;
- the Rannacher start and the theta;
- a CMake build with an executable for pricing, one for benchmarks and a test suite run by CTest.

## Repository structure

```
include/pricer/
    Payoff.hpp            payoffs and their boundary values (call, put)
    PDESolver.hpp         grid, theta-scheme, Rannacher start, Greeks
    Tridiagonal.hpp       Thomas algorithm
    BlackScholes.hpp      closed-form price and Greeks
src/                      implementations
apps/
    price.cpp             command-line pricer
    benchmark.cpp         accuracy on the test options and convergence study
tests/test_pricer.cpp     unit tests
report/                   original course report (in French) and its LaTeX source
```

## Usage

A C++17 compiler and CMake 3.20 or later.

```
cmake -S . -B build
cmake --build build
ctest --test-dir build
```

Pricing (all arguments are optional; the defaults are shown):

```
build/price --type call --spot 100 --strike 100 --maturity 1 --rate 0.05 --vol 0.2 --space 401 --time 400 --theta 0.5
```

```
European call, S = 100, K = 100, T = 1, r = 0.05, sigma = 0.2
grid: 401 nodes in ln S, 400 time steps, theta = 0.5

                    PDE    closed form   abs. error
price         10.449767      10.450584     8.17e-04
delta          0.636837       0.636831     6.32e-06
gamma          0.018764       0.018762     2.16e-06
theta         -6.414532      -6.414028     5.04e-04
```

`build/benchmark cases` and `build/benchmark convergence` print the two studies above as CSV.

The tests check the Thomas algorithm on random systems, the closed form against reference values and
put-call parity, the 40 test options within $3 \cdot 10^{-3}$ on the price and $10^{-4}$ on the delta,
the convergence orders of both schemes, the effect of the Rannacher start on the gamma, the
monotonicity and no-arbitrage bounds of the computed prices, and the rejection of invalid parameters.
