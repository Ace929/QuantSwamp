#include <vector>
#include <cmath>
#include <iostream>

// Crank-Nicolson solver for Black-Scholes PDE
double crankNicolsonPricing(double S0, double K, double r, double sigma, double T, int N, int M) {
    // Grid parameters
    double dt = T / M;            // Time step
    double dS = 2 * S0 / N;       // Price step (adjust grid size as needed)
    std::vector<double> S(N + 1); // Spot price grid
    std::vector<double> V(N + 1); // Option value grid

    // Initialize price grid
    for (int i = 0; i <= N; ++i) {
        S[i] = i * dS;
    }

    // Terminal condition (payoff at maturity)
    for (int i = 0; i <= N; ++i) {
        V[i] = std::max(S[i] - K, 0.0); // Call option payoff
    }

    // Coefficients for tridiagonal system (alpha, beta, gamma)
    std::vector<double> alpha(N + 1), beta(N + 1), gamma(N + 1);

    // Time-stepping loop
    for (int t = M - 1; t >= 0; --t) {
        std::vector<double> rhs(N + 1); // Right-hand side vector

        // Populate coefficients and RHS for Crank-Nicolson
        for (int i = 1; i < N; ++i) {
            double sig2 = sigma * sigma;
            double a = 0.25 * dt * (sig2 * i * i - r * i);
            double b = -dt * 0.5 * (sig2 * i * i + r);
            double c = 0.25 * dt * (sig2 * i * i + r * i);

            alpha[i] = -a;
            beta[i] = 1 - b;
            gamma[i] = -c;

            rhs[i] = a * V[i-1] + (1 + b) * V[i] + c * V[i+1];
        }

        // Boundary conditions (Dirichlet for call option)
        V[0] = 0.0;                // S=0: option value is 0
        V[N] = S[N] - K * exp(-r * (T - t * dt)); // S=infinity: discounted payoff

        // Solve tridiagonal system (Thomas algorithm)
        for (int i = 1; i < N; ++i) {
            beta[i] = beta[i] - alpha[i] * gamma[i-1] / beta[i-1];
            rhs[i] = rhs[i] - alpha[i] * rhs[i-1] / beta[i-1];
        }

        // Back-substitution
        V[N-1] = rhs[N-1] / beta[N-1];
        for (int i = N-2; i > 0; --i) {
            V[i] = (rhs[i] - gamma[i] * V[i+1]) / beta[i];
        }
    }

    // Linear interpolation to find V(S0)
    int idx = static_cast<int>(S0 / dS);
    double w = (S0 - S[idx]) / dS;
    return (1 - w) * V[idx] + w * V[idx + 1];
}

int main() {
    double S0 = 100.0;   // Spot price
    double K = 100.0;     // Strike
    double r = 0.05;      // Risk-free rate
    double sigma = 0.2;    // Volatility
    double T = 1.0;       // Time to maturity
    int N = 100;          // Price grid steps
    int M = 100;          // Time steps

    double price = crankNicolsonPricing(S0, K, r, sigma, T, N, M);
    std::cout << "European Call Price: " << price << std::endl;

    return 0;
}