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

// T-GARCH(1,1) Model Computation
void t_garch(const vector<double>& returns, double omega, double alpha_pos, double alpha_neg, double beta) {
    int T = returns.size();
    vector<double> h(T, 1.0); // Initialize variance

    for (int t = 1; t < T; ++t) {
        double r = returns[t - 1];
        double r_pos = (r > 0) ? r * r : 0; // Positive return effect
        double r_neg = (r < 0) ? r * r : 0; // Negative return effect

        // Compute variance update with different impacts for positive/negative shocks
        h[t] = omega + alpha_pos * r_pos + alpha_neg * r_neg + beta * h[t - 1];

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << h[t] << " | Return: " << r 
             << " | Positive Effect: " << r_pos << " | Negative Effect: " << r_neg << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set T-GARCH parameters
    double omega = 0.01;    // Constant term
    double alpha_pos = 0.1; // ARCH effect for positive returns
    double alpha_neg = 0.2; // ARCH effect for negative returns (higher for leverage effect)
    double beta = 0.85;     // GARCH effect (volatility persistence)

    // Run T-GARCH model
    t_garch(returns, omega, alpha_pos, alpha_neg, beta);

    return 0;
}
