#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Generate GARCH-M(1,1) return data
vector<double> generate_garch_m_returns(int T, double mu, double lambda, double omega, double alpha, double beta) {
    vector<double> returns(T, 0.0);
    vector<double> h(T, 1.0); // Initialize variance to 1

    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal(0, 1);

    for (int t = 1; t < T; ++t) {
        double z_t = normal(gen); // Standard normal shock
        h[t] = omega + alpha * pow(returns[t - 1], 2) + beta * h[t - 1]; // GARCH variance update
        returns[t] = mu + lambda * h[t] + sqrt(h[t]) * z_t; // GARCH-M return
    }
    
    return returns;
}

// Run GARCH-M model and print results
void garch_m(int T, double mu, double lambda, double omega, double alpha, double beta) {
    vector<double> returns = generate_garch_m_returns(T, mu, lambda, omega, alpha, beta);

    cout << "Time | Return | Variance" << endl;
    for (int t = 0; t < T; ++t) {
        cout << t << " | " << returns[t] << " | " << (omega + alpha * pow(returns[max(0, t - 1)], 2) + beta * 1.0) << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Set GARCH-M parameters
    double mu = 0.02;    // Mean return
    double lambda = 0.1; // Risk premium (volatility effect on return)
    double omega = 0.01; // Long-term variance
    double alpha = 0.1;  // ARCH effect (sensitivity to shocks)
    double beta = 0.85;  // GARCH effect (volatility persistence)

    // Run GARCH-M model
    garch_m(T, mu, lambda, omega, alpha, beta);

    return 0;
}
