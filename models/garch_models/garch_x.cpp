#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Generate synthetic return and exogenous data
void generate_data(int T, vector<double>& returns, vector<double>& X) {
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal(0, 1);

    for (int t = 0; t < T; ++t) {
        returns.push_back(normal(gen)); // Normally distributed returns
        X.push_back(0.5 + 0.1 * normal(gen)); // Exogenous variable (e.g., market volatility index)
    }
}

// GARCH-X(1,1) Model Computation
void garch_x(const vector<double>& returns, const vector<double>& X, double omega, double alpha, double beta, double delta) {
    int T = returns.size();
    vector<double> h(T, 1.0); // Initialize variance

    for (int t = 1; t < T; ++t) {
        double r2 = returns[t - 1] * returns[t - 1]; // Squared return

        // Compute variance update incorporating exogenous variable
        h[t] = omega + alpha * r2 + beta * h[t - 1] + delta * X[t];

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << h[t] << " | Exogenous X_t: " << X[t] << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps
    vector<double> returns, X;

    // Generate synthetic returns and exogenous variable
    generate_data(T, returns, X);

    // Set GARCH-X parameters
    double omega = 0.01; // Constant term
    double alpha = 0.1;  // ARCH effect
    double beta = 0.85;  // GARCH effect
    double delta = 0.05; // Effect of exogenous variable

    // Run GARCH-X model
    garch_x(returns, X, omega, alpha, beta, delta);

    return 0;
}
