#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Function to simulate Geometric Brownian Motion
vector<double> simulateGBM(int steps, double delta_t, double S0, double mu, double sigma, int seed = 42) {
    vector<double> trajectory;
    trajectory.push_back(S0); // Initial stock price

    // Random number generator
    random_device rd;
    mt19937 generator(seed); 
    normal_distribution<double> distribution(0.0, 1.0); 

    for (int i = 1; i <= steps; ++i) {
        double Z = distribution(generator); // Random normal variable
        double St = trajectory.back() * exp((mu - 0.5 * sigma * sigma) * delta_t + sigma * sqrt(delta_t) * Z);
        trajectory.push_back(St);
    }

    return trajectory;
}

int main() {
    int steps = 100;      // Number of time steps
    double delta_t = 1.0; // Time step size (e.g., days)
    double S0 = 100.0;    // Initial stock price
    double mu = 0.05;     // Expected return (5% annualized)
    double sigma = 0.2;   // Volatility (20% annualized)

    vector<double> trajectory = simulateGBM(steps, delta_t, S0, mu, sigma);

    // Print the simulated stock price path
    cout << "Geometric Brownian Motion (Stock Prices):\n";
    for (size_t i = 0; i < trajectory.size(); ++i) {
        cout << "Step " << i << ": $" << trajectory[i] << "\n";
    }

    return 0;
}
