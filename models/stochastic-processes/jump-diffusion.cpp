#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Function to simulate a Jump Diffusion Process (Merton Model)
vector<double> simulateJumpDiffusion(int steps, double delta_t, double S0, double mu, double sigma, double lambda, double jump_mu, double jump_sigma, int seed = 42) {
    vector<double> trajectory;
    trajectory.push_back(S0); // Initial stock price

    // Random number generators
    random_device rd;
    mt19937 generator(seed);
    normal_distribution<double> normal_dist(0.0, 1.0); // Standard normal for Brownian motion
    poisson_distribution<int> poisson_dist(lambda * delta_t); // Poisson for jump events
    normal_distribution<double> jump_dist(jump_mu, jump_sigma); // Log-normal jump size

    for (int i = 1; i <= steps; ++i) {
        double Z = normal_dist(generator); // Brownian motion noise
        int N = poisson_dist(generator); // Number of jumps
        double J = (N > 0) ? exp(jump_dist(generator)) - 1 : 0; // Jump impact

        double St = trajectory.back() * exp((mu - 0.5 * sigma * sigma) * delta_t + sigma * sqrt(delta_t) * Z + J * N);
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
    double lambda = 0.5;  // Jump frequency (0.5 jumps per unit time)
    double jump_mu = -0.1; // Average jump size (log-normal mean)
    double jump_sigma = 0.2; // Jump size volatility

    vector<double> trajectory = simulateJumpDiffusion(steps, delta_t, S0, mu, sigma, lambda, jump_mu, jump_sigma);

    // Print the simulated stock price path
    cout << "Jump Diffusion Model (Stock Prices):\n";
    for (size_t i = 0; i < trajectory.size(); ++i) {
        cout << "Step " << i << ": $" << trajectory[i] << "\n";
    }

    return 0;
}
