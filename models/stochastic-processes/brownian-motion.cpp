#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>

using namespace std;
using namespace Eigen;

// Function to simulate Brownian Motion
vector<VectorXd> simulateBrownianMotion(int steps, double delta_t, double sigma, VectorXd initial_position, int seed = 42) {
    int dimensions = initial_position.size();
    vector<VectorXd> trajectory;
    trajectory.push_back(initial_position);

    // Random number generator (normal distribution)
    random_device rd;
    mt19937 generator(seed); 
    normal_distribution<double> distribution(0.0, 1.0); 

    for (int i = 1; i <= steps; ++i) {
        VectorXd step(dimensions);
        for (int j = 0; j < dimensions; ++j) {
            step(j) = sigma * sqrt(delta_t) * distribution(generator);
        }
        initial_position += step;
        trajectory.push_back(initial_position);
    }

    return trajectory;
}

int main() {
    int steps = 100;             // Number of time steps
    double delta_t = 1.0;        // Time step size
    double sigma = 1.0;          // Volatility factor
    VectorXd initial_position(2); // 2D Brownian motion (x, y)
    initial_position << 0, 0;   // Start at (0,0)

    vector<VectorXd> trajectory = simulateBrownianMotion(steps, delta_t, sigma, initial_position);

    // Print the trajectory
    cout << "Brownian Motion Trajectory:\n";
    for (size_t i = 0; i < trajectory.size(); ++i) {
        cout << "Step " << i << ": (" << trajectory, " << trajectory  << 
    }

    return 0;
}
