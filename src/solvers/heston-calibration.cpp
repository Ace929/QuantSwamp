#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>
#include <cmath>

using namespace Eigen;

// Heston model price via COS method (simplified)
double hestonPriceCOS(double K, double T, double r, double kappa, double theta, 
                      double sigma, double rho, double v0, int N=200) {
    // Simplified COS method implementation (placeholder for brevity)
    // Assume complex integration replaced with closed-form for example
    return std::max(100.0 - K, 0.0); // Replace with actual COS code
}

double calibrationError(const std::vector<double>& params, const std::vector<double>& marketPrices,
                        const std::vector<double>& strikes, double r, double T) {
    double kappa = params[0], theta = params[1], sigma = params[2], rho = params[3], v0 = params[4];
    double error = 0.0;
    for (size_t i = 0; i < strikes.size(); ++i) {
        double modelPrice = hestonPriceCOS(strikes[i], T, r, kappa, theta, sigma, rho, v0);
        error += (modelPrice - marketPrices[i]) * (modelPrice - marketPrices[i]);
    }
    return error;
}

std::vector<double> differentialEvolution(int populationSize, int maxGen, 
                                          const std::vector<double>& marketPrices,
                                          const std::vector<double>& strikes,
                                          double r, double T) {
    std::mt19937 gen(42);
    std::uniform_real_distribution<double> dist(0.0, 1.0);
    const int dim = 5; // kappa, theta, sigma, rho, v0

    // Parameter bounds [min, max]
    std::vector<std::pair<double, double>> bounds = {
        {0.1, 5.0},   // kappa
        {0.01, 0.5},  // theta
        {0.1, 1.0},   // sigma (vol-of-vol)
        {-0.9, 0.0},  // rho
        {0.01, 0.5}   // v0
    };

    // Initialize population
    std::vector<std::vector<double>> population(populationSize);
    for (int i = 0; i < populationSize; ++i) {
        for (int j = 0; j < dim; ++j) {
            population[i].push_back(bounds[j].first + dist(gen) * 
                                   (bounds[j].second - bounds[j].first));
        }
    }

    // Evolution loop
    for (int gen = 0; gen < maxGen; ++gen) {
        for (int i = 0; i < populationSize; ++i) {
            // Mutation: Pick 3 random vectors (a, b, c)
            int a = i, b, c;
            do { b = dist(gen) * populationSize; } while (b == i);
            do { c = dist(gen) * populationSize; } while (c == i || c == b);

            // Generate trial vector
            std::vector<double> trial(dim);
            for (int j = 0; j < dim; ++j) {
                trial[j] = population[a][j] + 0.8 * (population[b][j] - population[c][j]);
                // Apply bounds
                trial[j] = std::max(bounds[j].first, std::min(bounds[j].second, trial[j]));
            }

            // Crossover
            for (int j = 0; j < dim; ++j) {
                if (dist(gen) > 0.5) trial[j] = population[i][j];
            }

            // Selection
            double errorOld = calibrationError(population[i], marketPrices, strikes, r, T);
            double errorNew = calibrationError(trial, marketPrices, strikes, r, T);
            if (errorNew < errorOld) {
                population[i] = trial;
            }
        }
    }

    // Find best solution
    auto best = population[0];
    double bestError = calibrationError(best, marketPrices, strikes, r, T);
    for (const auto& p : population) {
        double error = calibrationError(p, marketPrices, strikes, r, T);
        if (error < bestError) {
            best = p;
            bestError = error;
        }
    }
    return best;
}

int main() {
    // Example market data
    std::vector<double> strikes = {90, 100, 110};
    std::vector<double> marketPrices = {12.5, 8.2, 5.1}; // Placeholder
    double r = 0.05, T = 1.0;

    auto params = differentialEvolution(50, 100, marketPrices, strikes, r, T);
    std::cout << "Calibrated Heston params: kappa=" << params[0] 
              << ", theta=" << params[1] << ", sigma=" << params[2] 
              << ", rho=" << params[3] << ", v0=" << params[4] << std::endl;
    return 0;
}