#include <iostream>
#include <vector>
#include <Eigen/Dense>
using namespace Eigen;

void solveTridiagonal(const VectorXd& a, const VectorXd& b, const VectorXd& c, VectorXd& rhs) {
    int n = rhs.size();
    VectorXd beta(n), gamma(n);
    beta(0) = b(0);
    gamma(0) = rhs(0) / beta(0);
    for (int i = 1; i < n; ++i) {
        beta(i) = b(i) - a(i) * c(i-1) / beta(i-1);
        gamma(i) = (rhs(i) - a(i) * gamma(i-1)) / beta(i);
    }
    for (int i = n-2; i >= 0; --i) {
        gamma(i) -= c(i) * gamma(i+1) / beta(i);
    }
    rhs = gamma;
}

double hestonADI(double S0, double K, double r, double kappa, double theta, 
                double sigma, double rho, double v0, double T, int N, int M) {
    double dt = T / M;
    double dS = 2 * S0 / N, dv = 0.5 / N;
    MatrixXd V = MatrixXd::Zero(N+1, N+1); // Grid: S (rows), v (cols)

    // Terminal condition
    for (int i = 0; i <= N; ++i) {
        double S = i * dS;
        V.row(i) = VectorXd::Constant(N+1, std::max(S - K, 0.0));
    }

    // ADI steps
    for (int t = M-1; t >= 0; --t) {
        // Step 1: Implicit in S, explicit in v (handle mixed derivative explicitly)
        for (int j = 1; j < N; ++j) {
            double v = j * dv;
            VectorXd a(N+1), b(N+1), c(N+1), rhs(N+1);
            for (int i = 1; i < N; ++i) {
                double S = i * dS;
                // Coefficients for S-direction
                a(i) = -0.25 * dt * (v * S*S / (dS*dS) + r * S / dS);
                b(i) = 1 + dt * (v * S*S / (dS*dS) + r);
                c(i) = -0.25 * dt * (v * S*S / (dS*dS) - r * S / dS);
                // RHS with cross-term (explicit)
                rhs(i) = V(i,j) + 0.5 * dt * (
                    rho * sigma * v * S * (V(i+1,j+1) - V(i-1,j+1) - V(i+1,j-1) + V(i-1,j-1)) / (4*dS*dv)
                );
            }
            solveTridiagonal(a, b, c, rhs);
            V.col(j) = rhs;
        }

        // Step 2: Implicit in v, explicit in S (similar structure)
        // [Mirror Step 1 for v-direction]
    }

    // Interpolate to (S0, v0)
    return V((int)(S0/dS), (int)(v0/dv));
}

int main() {
    double price = hestonADI(100, 100, 0.05, 2.0, 0.04, 0.3, -0.7, 0.04, 1.0, 100, 50);
    std::cout << "Heston Call Price: " << price << std::endl;
    return 0;
}