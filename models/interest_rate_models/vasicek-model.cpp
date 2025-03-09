#include <iostream>
#include <vector>
#include <cmath>
#include <random>

using namespace std;

// Vasicek Model Parameters
const double mu = 0.05;   // Long-term mean interest rate
const double theta = 0.15; // Mean reversion speed
const double sigma = 0.02; // Volatility
const double r0 = 0.03;   // Initial interest rate
const double T = 1.0;     // Time horizon (years)
const int steps = 100;    // Number of time steps

// Function to simulate Vasicek short rate path
vector<double> simulateVasicek(double dt, int seed = 42) {
    vector<double> r(steps, r0);

    // Random number generator
    random_device rd;
    mt19937 generator(seed);
    normal_distribution<double> normalDist(0.0, 1.0);

    // Euler-Maruyama method
    for (int i = 1; i < steps; ++i) {
        double dW = normalDist(generator) * sqrt(dt);
        r[i] = r[i - 1] + theta * (mu - r[i - 1]) * dt + sigma * dW;
    }

    return r;
}

int main() {
    double dt = T / steps; // Time step size
    vector<double> rates = simulateVasicek(dt);

    // Print simulated interest rate path
    cout << "Vasicek Interest Rate Path:\n";
    for (size_t i = 0; i < rates.size(); ++i) {
        cout << "Time " << i * dt << ": r = " << rates[i] << endl;
    }

    return 0;
}
