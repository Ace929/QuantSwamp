#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Function to generate synthetic return data
vector<double> generate_synthetic_returns(int T) {
    vector<double> returns;
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> d(0, 1);

    for (int t = 0; t < T; ++t) {
        returns.push_back(d(gen)); // Normally distributed returns
    }
    return returns;
}

// E-GARCH(1,1) Model Computation
void e_garch(const vector<double>& returns, double omega, double alpha, double beta, double gamma) {
    int T = returns.size();
    vector<double> log_h(T, log(1.0));  // Log-variance initialized

    for (int t = 1; t < T; ++t) {
        double r = returns[t - 1];
        double std_r = r / exp(log_h[t - 1] / 2.0); // Standardized return

        // Compute log variance update
        log_h[t] = omega + beta * log_h[t - 1] + alpha * (fabs(std_r) - sqrt(2.0 / M_PI)) + gamma * std_r;

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << exp(log_h[t]) << " | Log Variance: " << log_h[t] << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set E-GARCH parameters
    double omega = -0.5; // Constant
    double alpha = 0.1;  // Volatility magnitude effect
    double beta = 0.9;   // Persistence
    double gamma = -0.2; // Asymmetry (negative means bad news increases volatility more)

    // Run E-GARCH model
    e_garch(returns, omega, alpha, beta, gamma);

    return 0;
}
