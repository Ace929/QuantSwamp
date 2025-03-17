#include <iostream>
#include <vector>
#include <cmath>

double explicitFDMPut(double S0, double K, double r, double sigma, double T, int N, int M) {
    double dt = T / M;       // Time step
    double dS = 2 * S0 / N;  // Price step
    double maxS = N * dS;    // Maximum price

    // Stability check (CFL condition)
    double stability = dt * (sigma * sigma * N * N + r * N);
    if (stability > 1.0) {
        std::cerr << "Warning: Explicit scheme may be unstable. Reduce time steps." << std::endl;
    }

    std::vector<double> S(N + 1), V_old(N + 1), V_new(N + 1);

    // Initialize price grid and terminal payoff
    for (int i = 0; i <= N; ++i) {
        S[i] = i * dS;
        V_old[i] = std::max(K - S[i], 0.0); // Put payoff
    }

    // Explicit time-stepping
    for (int t = M - 1; t >= 0; --t) {
        for (int i = 1; i < N; ++i) {
            double delta = (V_old[i + 1] - V_old[i - 1]) / (2 * dS);
            double gamma = (V_old[i + 1] - 2 * V_old[i] + V_old[i - 1]) / (dS * dS);
            double theta = 0.5 * sigma * sigma * S[i] * S[i] * gamma 
                         + r * S[i] * delta - r * V_old[i];
            V_new[i] = V_old[i] - dt * theta;
        }

        // Boundary conditions
        V_new[0] = K * exp(-r * (T - t * dt)); // S=0: discounted strike
        V_new[N] = 0.0;                        // S=maxS: put value 0

        V_old = V_new;
    }

    // Interpolate to S0
    int idx = S0 / dS;
    double w = (S0 - S[idx]) / dS;
    return (1 - w) * V_old[idx] + w * V_old[idx + 1];
}

int main() {
    double price = explicitFDMPut(100.0, 100.0, 0.05, 0.2, 1.0, 100, 1000);
    std::cout << "European Put Price (Explicit FDM): " << price << std::endl;
    return 0;
}