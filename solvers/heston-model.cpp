#include <iostream>
#include <cmath>
#include <random>
#include <vector>

// Struct to hold Heston model parameters
struct HestonParams {
    double S0;       // Initial stock price
    double V0;       // Initial variance
    double K;        // Strike price
    double r;        // Risk-free rate
    double kappa;    // Mean reversion speed of variance
    double theta;    // Long-run variance
    double sigma;    // Volatility of volatility (vol-of-vol)
    double rho;      // Correlation between asset and volatility
    double T;        // Time to maturity
    int N;           // Number of time steps
    int M;           // Number of Monte Carlo simulations
};

// Standard normal cumulative distribution function
double normalCDF(double x) {
    return 0.5 * std::erfc(-x * std::sqrt(0.5));
}

// Generate correlated random variables using Cholesky decomposition
void generateCorrelatedGaussians(double rho, double &Z1, double &Z2) {
    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);

    Z1 = d(gen);
    Z2 = rho * Z1 + std::sqrt(1 - rho * rho) * d(gen);
}

// Monte Carlo simulation for the Heston model
double hestonMonteCarlo(const HestonParams &params) {
    double dt = params.T / params.N;
    std::vector<double> ST(params.M);

    std::random_device rd;
    std::mt19937 gen(rd());
    std::normal_distribution<> d(0.0, 1.0);

    for (int i = 0; i < params.M; ++i) {
        double S = params.S0;
        double V = params.V0;

        for (int j = 0; j < params.N; ++j) {
            double Z1, Z2;
            generateCorrelatedGaussians(params.rho, Z1, Z2);

            // Ensure variance remains positive using max(V, 0)
            double sqrtV = std::sqrt(std::max(V, 0.0));
            S *= std::exp((params.r - 0.5 * V) * dt + sqrtV * std::sqrt(dt) * Z1);
            V += params.kappa * (params.theta - V) * dt + params.sigma * sqrtV * std::sqrt(dt) * Z2;

            // Ensure non-negative variance using max(V, 0)
            V = std::max(V, 0.0);
        }
        ST[i] = std::max(S - params.K, 0.0); // Payoff for European call option
    }

    // Compute expected option price under risk-neutral measure
    double sumPayoff = 0.0;
    for (double payoff : ST) {
        sumPayoff += payoff;
    }

    return std::exp(-params.r * params.T) * (sumPayoff / params.M);
}

int main() {
    HestonParams params;
    
    // User input parameters
    std::cout << "Enter stock price (S0): ";
    std::cin >> params.S0;
    std::cout << "Enter initial variance (V0): ";
    std::cin >> params.V0;
    std::cout << "Enter strike price (K): ";
    std::cin >> params.K;
    std::cout << "Enter risk-free interest rate (r) in decimal (e.g., 0.05 for 5%): ";
    std::cin >> params.r;
    std::cout << "Enter mean reversion speed (kappa): ";
    std::cin >> params.kappa;
    std::cout << "Enter long-run variance (theta): ";
    std::cin >> params.theta;
    std::cout << "Enter volatility of volatility (sigma): ";
    std::cin >> params.sigma;
    std::cout << "Enter correlation (rho): ";
    std::cin >> params.rho;
    std::cout << "Enter time to maturity (T) in years: ";
    std::cin >> params.T;
    std::cout << "Enter number of time steps (N): ";
    std::cin >> params.N;
    std::cout << "Enter number of Monte Carlo simulations (M): ";
    std::cin >> params.M;

    // Compute option price using Heston model
    double optionPrice = hestonMonteCarlo(params);
    std::cout << "Heston Model Option Price: " << optionPrice << std::endl;

    return 0;
}
