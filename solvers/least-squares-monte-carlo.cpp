#include <iostream>
#include <vector>
#include <cmath>
#include <random>
#include <Eigen/Dense> // For regression

using namespace Eigen;

double lsmAmericanPut(double S0, double K, double r, double sigma, double T, int M, int N) {
    double dt = T / M;
    std::vector<std::vector<double>> paths(N, std::vector<double>(M + 1));
    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);

    // Simulate paths
    for (int sim = 0; sim < N; ++sim) {
        paths[sim][0] = S0;
        for (int t = 1; t <= M; ++t) {
            double Z = dist(gen);
            paths[sim][t] = paths[sim][t - 1] * exp((r - 0.5 * sigma * sigma) * dt 
                           + sigma * sqrt(dt) * Z);
        }
    }

    // Backward induction
    std::vector<double> cashflow(N);
    for (int sim = 0; sim < N; ++sim) {
        cashflow[sim] = std::max(K - paths[sim][M], 0.0);
    }

    for (int t = M - 1; t > 0; --t) {
        MatrixXd X(N, 3); // Basis functions: 1, S, S^2
        VectorXd Y(N);
        int inMoneyCount = 0;

        for (int sim = 0; sim < N; ++sim) {
            double S_t = paths[sim][t];
            if (S_t < K) { // Only consider in-the-money paths
                X(inMoneyCount, 0) = 1.0;
                X(inMoneyCount, 1) = S_t;
                X(inMoneyCount, 2) = S_t * S_t;
                Y(inMoneyCount) = cashflow[sim] * exp(-r * dt);
                inMoneyCount++;
            }
        }

        if (inMoneyCount > 0) {
            X.conservativeResize(inMoneyCount, 3);
            Y.conservativeResize(inMoneyCount);
            VectorXd beta = X.bdcSvd(ComputeThinU | ComputeThinV).solve(Y);

            // Update cashflows
            int idx = 0;
            for (int sim = 0; sim < N; ++sim) {
                double S_t = paths[sim][t];
                if (S_t < K) {
                    double continuation = beta[0] + beta[1] * S_t + beta[2] * S_t * S_t;
                    double exercise = K - S_t;
                    if (exercise > continuation) {
                        cashflow[sim] = exercise;
                    }
                    idx++;
                }
            }
        }
    }

    // Discount final cashflows to present
    double price = 0.0;
    for (int sim = 0; sim < N; ++sim) {
        price += cashflow[sim] * exp(-r * dt);
    }
    return price / N;
}

int main() {
    double price = lsmAmericanPut(100.0, 100.0, 0.05, 0.2, 1.0, 252, 10000);
    std::cout << "American Put Price (LSM): " << price << std::endl;
    return 0;
}