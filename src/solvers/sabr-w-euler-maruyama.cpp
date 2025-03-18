#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>

using namespace Eigen;

double sabrMonteCarlo(double F0, double alpha, double beta, double rho, 
                      double sigma0, double T, int M, int N) {
    double dt = T / M;
    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);
    Matrix2d corrMatrix;
    corrMatrix << 1.0, rho, rho, 1.0;
    LLT<Matrix2d> llt(corrMatrix);
    Matrix2d L = llt.matrixL();

    double payoffSum = 0.0;
    for (int sim = 0; sim < N; ++sim) {
        double F = F0, sigma = sigma0;
        for (int t = 0; t < M; ++t) {
            Vector2d Z = L * Vector2d(dist(gen), dist(gen));
            F += sigma * std::pow(F, beta) * Z[0] * sqrt(dt);
            sigma += alpha * sigma * Z[1] * sqrt(dt);
        }
        payoffSum += std::max(F - F0, 0.0); // Call option
    }
    return exp(-alpha * T) * payoffSum / N; // Simplified discounting
}

int main() {
    double price = sabrMonteCarlo(100.0, 0.3, 0.5, -0.6, 0.2, 1.0, 252, 100000);
    std::cout << "SABR Call Price: " << price << std::endl;
    return 0;
}