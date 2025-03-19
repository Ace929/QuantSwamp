#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense>
using namespace Eigen;

class QuantumMonteCarlo {
public:
    QuantumMonteCarlo(int dim, int numSamples) : dim(dim), numSamples(numSamples), gen(42) {}

    double priceEuropeanCall(double S0, double K, double r, double sigma, double T) {
        VectorXd sobol = generateSobolSequence(numSamples, dim); // Use external Sobol library
        double payoffSum = 0.0;
        for (int i = 0; i < numSamples; ++i) {
            double Z = inverseGaussianCDF(sobol(i)); // Box-Muller transform for Sobol→Gaussian
            double ST = S0 * exp((r - 0.5*sigma*sigma)*T + sigma*sqrt(T)*Z);
            payoffSum += std::max(ST - K, 0.0);
        }
        return exp(-r*T) * payoffSum / numSamples;
    }

private:
    int dim, numSamples;
    std::mt19937 gen;

    double inverseGaussianCDF(double u) {
        return sqrt(2) * erfinv(2*u - 1); // Approximate transform
    }
};

int main() {
    QuantumMonteCarlo qmc(1, 10000);
    double price = qmc.priceEuropeanCall(100.0, 100.0, 0.05, 0.2, 1.0);
    std::cout << "QMC Call Price: " << price << std::endl;
    return 0;
}