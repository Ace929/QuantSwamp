#include <iostream>
#include <vector>
#include <algorithm>

void projectedSOR(std::vector<double>& V, const std::vector<double>& payoff, 
                  const std::vector<std::vector<double>>& A, double tol=1e-6, int maxIter=1000) {
    int n = V.size();
    double omega = 1.2; // Relaxation parameter

    for (int iter = 0; iter < maxIter; ++iter) {
        double error = 0.0;
        for (int i = 0; i < n; ++i) {
            double sigma = 0.0;
            for (int j = 0; j < n; ++j) {
                if (j != i) sigma += A[i][j] * V[j];
            }
            double V_new = (payoff[i] - sigma) / A[i][i];
            V_new = std::max(payoff[i], V[i] + omega * (V_new - V[i])); // Projection
            error += std::abs(V_new - V[i]);
            V[i] = V_new;
        }
        if (error < tol) break;
    }
}

int main() {
    // Example: Solve LCP for American option (simplified system)
    int n = 100;
    std::vector<double> V(n, 1.0), payoff(n);
    std::vector<std::vector<double>> A(n, std::vector<double>(n, 0.0));

    // Define A (e.g., tridiagonal matrix from PDE discretization)
    for (int i = 0; i < n; ++i) {
        A[i][i] = 1.5; // Diagonal dominance for convergence
        if (i > 0) A[i][i-1] = -0.5;
        if (i < n-1) A[i][i+1] = -0.5;
        payoff[i] = std::max(100.0 - i, 0.0); // Example payoff
    }

    projectedSOR(V, payoff, A);
    std::cout << "American Option Value: " << V[50] << std::endl;
    return 0;
}