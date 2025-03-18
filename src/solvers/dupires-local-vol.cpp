#include <vector>
#include <Eigen/Dense>
using namespace Eigen;

VectorXd dupireLocalVol(const VectorXd& strikes, const VectorXd& maturities,
                        const MatrixXd& marketPrices, double r) {
    int N = strikes.size(), M = maturities.size();
    MatrixXd localVols(N, M);

    for (int j = 0; j < M; ++j) {
        double T = maturities(j);
        for (int i = 1; i < N-1; ++i) {
            // Compute ∂C/∂T, ∂C/∂K, ∂²C/∂K² using finite differences
            double dCdT = (marketPrices(i,j+1) - marketPrices(i,j)) / (maturities(j+1) - maturities(j));
            double dCdK = (marketPrices(i+1,j) - marketPrices(i-1,j)) / (strikes(i+1) - strikes(i-1));
            double d2CdK2 = (marketPrices(i+1,j) - 2*marketPrices(i,j) + marketPrices(i-1,j)) 
                            / pow(strikes(i+1) - strikes(i), 2);
            // Dupire's formula: σ² = [∂C/∂T + rK∂C/∂K] / (0.5 K² ∂²C/∂K²)
            localVols(i,j) = sqrt((dCdT + r * strikes(i) * dCdK) / (0.5 * strikes(i)*strikes(i) * d2CdK2));
        }
    }
    return localVols;
}

// PDE solver for local volatility (similar to Black-Scholes FDM)
// [Implement finite difference solver using localVols grid]