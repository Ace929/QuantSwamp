#include <iostream>
#include <Eigen/Dense>
#include <random>
using namespace Eigen;

MatrixXd gaussianCopula(int nAssets, double rho, int nSims) {
    MatrixXd corr = MatrixXd::Constant(nAssets, nAssets, rho);
    for (int i = 0; i < nAssets; ++i) corr(i, i) = 1.0;

    LLT<MatrixXd> llt(corr);
    MatrixXd L = llt.matrixL(); // Cholesky decomposition

    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);

    MatrixXd copula(nSims, nAssets);
    for (int sim = 0; sim < nSims; ++sim) {
        VectorXd Z(nAssets);
        for (int i = 0; i < nAssets; ++i) Z(i) = dist(gen);
        VectorXd correlatedZ = L * Z;
        copula.row(sim) = correlatedZ;
    }
    return copula; // Returns matrix of correlated Gaussian variates
}

int main() {
    MatrixXd copula = gaussianCopula(5, 0.3, 1000);
    std::cout << "Copula samples (first 5):\n" << copula.topRows(5) << std::endl;
    return 0;
}