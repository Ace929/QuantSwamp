#include <iostream>
#include <vector>
#include <algorithm> // For max()

double implicitFDMAmericanPut(double S0, double K, double r, double sigma, double T, int N, int M) {
    double dt = T / M;
    double dS = 2 * S0 / N;
    std::vector<double> S(N + 1), V_old(N + 1), V_new(N + 1);

    // Initialize grid and terminal payoff
    for (int i = 0; i <= N; ++i) {
        S[i] = i * dS;
        V_old[i] = std::max(K - S[i], 0.0);
    }

    // Coefficients for tridiagonal system
    std::vector<double> a(N + 1), b(N + 1), c(N + 1), rhs(N + 1);

    for (int t = M - 1; t >= 0; --t) {
        for (int i = 1; i < N; ++i) {
            double sig2 = sigma * sigma;
            double alpha = 0.5 * dt * (sig2 * i * i - r * i);
            double beta = 1 + dt * (sig2 * i * i + r);
            double gamma = 0.5 * dt * (sig2 * i * i + r * i);

            a[i] = -alpha;
            b[i] = beta;
            c[i] = -gamma;
            rhs[i] = V_old[i];
        }

        // Boundary conditions
        rhs[0] = K * exp(-r * (T - t * dt)); // S=0
        rhs[N] = 0.0;                        // S=maxS

        // Thomas algorithm
        for (int i = 1; i < N; ++i) {
            double m = a[i] / b[i - 1];
            b[i] -= m * c[i - 1];
            rhs[i] -= m * rhs[i - 1];
        }

        V_new[N] = rhs[N] / b[N];
        for (int i = N - 1; i >= 0; --i) {
            V_new[i] = (rhs[i] - c[i] * V_new[i + 1]) / b[i];
            // Early exercise check
            V_new[i] = std::max(V_new[i], K - S[i]);
        }

        V_old = V_new;
    }

    // Interpolate to S0
    int idx = static_cast<int>(S0 / dS);
    double w = (S0 - S[idx]) / dS;
    return (1 - w) * V_old[idx] + w * V_old[idx + 1];
}

int main() {
    double price = implicitFDMAmericanPut(100.0, 100.0, 0.05, 0.2, 1.0, 200, 100);
    std::cout << "American Put Price (Implicit FDM): " << price << std::endl;
    return 0;
}