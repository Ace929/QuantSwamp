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
    normal_distribution<> d(0, 1);

    for (int t = 0; t < T; ++t) {
        returns.push_back(d(gen)); // Random normal returns
    }
    return returns;
}

// GARCH(1,1) Model Computation
void garch(const vector<double>& returns, double omega, double alpha, double beta) {
    int T = returns.size();
    vector<double> h(T, 1.0);  // Initial variance

    for (int t = 1; t < T; ++t) {
        double r2 = returns[t - 1] * returns[t - 1]; // Squared return

        // Compute variance update
        h[t] = omega + alpha * r2 + beta * h[t - 1];

        // Output variance at each step
        cout << "Time t=" << t << " | Variance: " << h[t] << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set GARCH parameters
    double omega = 0.01; // Constant term
    double alpha = 0.1;  // Impact of past squared returns (ARCH effect)
    double beta = 0.85;  // Persistence of variance (GARCH effect)

    // Run GARCH model
    garch(returns, omega, alpha, beta);

    return 0;
}
