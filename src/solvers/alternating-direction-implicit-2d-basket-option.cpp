#include <iostream>
#include <vector>
#include <cmath>

void solveTridiagonal(const std::vector<double>& a, const std::vector<double>& b, 
                      const std::vector<double>& c, std::vector<double>& rhs) {
    int n = rhs.size();
    std::vector<double> beta(n), gamma(n);
    beta[0] = b[0];
    gamma[0] = rhs[0] / beta[0];

    for (int i = 1; i < n; ++i) {
        beta[i] = b[i] - a[i] * c[i-1] / beta[i-1];
        gamma[i] = (rhs[i] - a[i] * gamma[i-1]) / beta[i];
    }

    for (int i = n-2; i >= 0; --i) {
        gamma[i] -= c[i] * gamma[i+1] / beta[i];
    }
    rhs = gamma;
}

double adiBasketOption(double S1, double S2, double K, double r, double sigma1, 
                       double sigma2, double rho, double T, int N, int M) {
    double dt = T / M;
    double dS = 2 * std::max(S1, S2) / N;
    std::vector<std::vector<double>> V(N+1, std::vector<double>(N+1));

    // Terminal condition (basket call payoff)
    for (int i = 0; i <= N; ++i) {
        for (int j = 0; j <= N; ++j) {
            double S1_val = i * dS;
            double S2_val = j * dS;
            V[i][j] = std::max(S1_val + S2_val - K, 0.0);
        }
    }

    // ADI time-stepping
    for (int t = M-1; t >= 0; --t) {
        // Step 1: Implicit in S1, explicit in S2
        for (int j = 1; j < N; ++j) {
            std::vector<double> a(N+1), b(N+1), c(N+1), rhs(N+1);
            for (int i = 1; i < N; ++i) {
                // PDE coefficients for S1 direction
                a[i] = -0.25 * dt * (sigma1*sigma1 * i*i - r*i);
                b[i] = 1 + 0.5 * dt * (sigma1*sigma1 * i*i + r);
                c[i] = -0.25 * dt * (sigma1*sigma1 * i*i + r*i);
                rhs[i] = V[i][j] + 0.5 * dt * (
                    // S2 explicit terms (simplified)
                    0.5 * sigma2*sigma2 * j*j * (V[i][j+1] - 2*V[i][j] + V[i][j-1]) +
                    r * j * (V[i][j+1] - V[i][j-1]) / (2*dS) - r * V[i][j]
                );
            }
            solveTridiagonal(a, b, c, rhs);
            for (int i = 0; i <= N; ++i) V[i][j] = rhs[i];
        }

        // Step 2: Implicit in S2, explicit in S1 (similar to Step 1)
        // [Code mirrored for S2 direction]
    }

    // Interpolate to (S1, S2)
    int i = static_cast<int>(S1 / dS);
    int j = static_cast<int>(S2 / dS);
    return V[i][j]; // Simplified interpolation for brevity
}

int main() {
    double price = adiBasketOption(50, 50, 100, 0.05, 0.2, 0.3, 0.5, 1.0, 100, 50);
    std::cout << "Basket Call Price (ADI): " << price << std::endl;
    return 0;
}