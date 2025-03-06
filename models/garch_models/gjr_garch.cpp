#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Generate synthetic return data
vector<double> generate_synthetic_returns(int T) {
    vector<double> returns;
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal(0, 1);

    for (int t = 0; t < T; ++t) {
        returns.push_back(normal(gen)); // Normally distributed returns
    }
    return returns;
}

// GJR-GARCH(1,1) Model Computation
void gjr_garch(const vector<double>& returns, double omega, double alpha, double beta, double gamma) {
    int T = returns.size();
    vector<double> h(T, 1.0); // Initialize variance

    for (int t = 1; t < T; ++t) {
        double r2 = returns[t - 1] * returns[t - 1]; // Squared return
        double I_t = (returns[t - 1] < 0) ? 1.0 : 0.0; // Indicator for negative returns

        // Compute variance update with asymmetry effect
        h[t] = omega + alpha * r2 + gamma * r2 * I_t + beta * h[t - 1];

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << h[t] << " | Return: " << returns[t - 1] << " | I_t: " << I_t << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set GJR-GARCH parameters
    double omega = 0.01; // Constant term
    double alpha = 0.1;  // ARCH effect (past shocks)
    double beta = 0.85;  // GARCH effect (volatility persistence)
    double gamma = 0.15; // Asymmetry effect (leverage)

    // Run GJR-GARCH model
    gjr_garch(returns, omega, alpha, beta, gamma);

    return 0;
}
