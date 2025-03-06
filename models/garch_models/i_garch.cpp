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

// I-GARCH(1,1) Model Computation
void i_garch(const vector<double>& returns, double omega, double alpha) {
    int T = returns.size();
    vector<double> h(T, 1.0); // Initialize variance

    for (int t = 1; t < T; ++t) {
        double r2 = returns[t - 1] * returns[t - 1]; // Squared return

        // Compute variance update ensuring alpha + beta = 1
        h[t] = omega + alpha * r2 + (1 - alpha) * h[t - 1];

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << h[t] << " | Return: " << returns[t - 1] << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set I-GARCH parameters
    double omega = 0.01; // Small constant term
    double alpha = 0.1;  // ARCH effect (ensuring alpha + beta = 1)

    // Run I-GARCH model
    i_garch(returns, omega, alpha);

    return 0;
}
