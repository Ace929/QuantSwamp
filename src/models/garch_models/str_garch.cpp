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

// Smooth transition function
double G(double s, double gamma, double c) {
    return 1.0 / (1.0 + exp(-gamma * (s - c)));
}

// STR-GARCH(1,1) Model Computation
void str_garch(const vector<double>& returns, double omega, double alpha, double beta, double gamma, double c) {
    int T = returns.size();
    vector<double> h(T, 1.0); // Initialize variance

    for (int t = 1; t < T; ++t) {
        double r2 = returns[t - 1] * returns[t - 1]; // Squared return
        double G_t = G(returns[t - 1], gamma, c); // Compute smooth transition function

        // Compute variance update with regime-dependent impact
        h[t] = omega + alpha * G_t * r2 + beta * h[t - 1];

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << h[t] << " | Return: " << returns[t - 1] 
             << " | G_t: " << G_t << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set STR-GARCH parameters
    double omega = 0.01; // Constant term
    double alpha = 0.1;  // ARCH effect
    double beta = 0.85;  // GARCH effect
    double gamma = 5.0;  // Smoothness of transition
    double c = 0.0;      // Threshold level

    // Run STR-GARCH model
    str_garch(returns, omega, alpha, beta, gamma, c);

    return 0;
}
