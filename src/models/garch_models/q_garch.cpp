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

// Q-GARCH(1,1) Model Computation
void q_garch(const vector<double>& returns, double omega, double alpha, double beta, double gamma) {
    int T = returns.size();
    vector<double> h(T, 1.0); // Initialize variance

    for (int t = 1; t < T; ++t) {
        double r2 = returns[t - 1] * returns[t - 1]; // Squared return
        double r4 = r2 * r2; // Fourth power return for Q-GARCH effect

        // Compute variance update with quadratic term
        h[t] = omega + alpha * r2 + beta * h[t - 1] + gamma * r4;

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << h[t] << " | Return: " << returns[t - 1] << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set Q-GARCH parameters
    double omega = 0.01; // Constant term
    double alpha = 0.1;  // ARCH effect
    double beta = 0.85;  // GARCH effect
    double gamma = 0.02; // Quadratic effect (for extreme returns)

    // Run Q-GARCH model
    q_garch(returns, omega, alpha, beta, gamma);

    return 0;
}
