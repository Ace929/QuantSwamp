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
        returns.push_back(d(gen)); // Normally distributed returns
    }
    return returns;
}

// C-GARCH(1,1) Model Computation
void c_garch(const vector<double>& returns, double alpha, double beta, double theta) {
    int T = returns.size();
    
    // Initialize variance components
    double sigma2 = 1.0;  // Initial variance
    double q = 1.0;       // Long-run variance

    // Iterate over time steps
    for (int t = 1; t < T; ++t) {
        double r2 = returns[t - 1] * returns[t - 1]; // Squared return
        
        // Update long-run component
        q = theta * q + (1 - theta) * sigma2;
        
        // Update conditional variance
        sigma2 = q + alpha * (r2 - q) + beta * (sigma2 - q);

        // Output variance at each step
        cout << "Time t=" << t << " | Conditional Variance: " << sigma2 << " | Long-run Variance: " << q << endl;
    }
}

int main() {
    int T = 100;  // Number of time steps

    // Generate synthetic returns
    vector<double> returns = generate_synthetic_returns(T);

    // Set C-GARCH parameters
    double alpha = 0.1;  // Short-run effect
    double beta = 0.85;  // GARCH persistence
    double theta = 0.98; // Long-run mean reversion

    // Run C-GARCH model
    c_garch(returns, alpha, beta, theta);

    return 0;
}
