#include <iostream>
#include <vector>
#include <cmath>
#include <random>

double monteCarloAsianCall(double S0, double K, double r, double sigma, double T, int M, int N) {
    double dt = T / M;
    double sumPayoffs = 0.0;

    // Random number setup
    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int sim = 0; sim < N; ++sim) {
        double S = S0;
        double pathSum = S;

        // Simulate one path
        for (int t = 0; t < M; ++t) {
            double Z = dist(gen);
            S *= exp((r - 0.5 * sigma * sigma) * dt + sigma * sqrt(dt) * Z);
            pathSum += S;
        }

        // Arithmetic average and payoff
        double avgPrice = pathSum / (M + 1);
        sumPayoffs += std::max(avgPrice - K, 0.0);
    }

    // Discounted average payoff
    return exp(-r * T) * sumPayoffs / N;
}

int main() {
    double price = monteCarloAsianCall(100.0, 100.0, 0.05, 0.2, 1.0, 252, 100000);
    std::cout << "Asian Call Price (Monte Carlo): " << price << std::endl;
    return 0;
}