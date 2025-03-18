#include <iostream>
#include <Eigen/Dense> // Requires Eigen library (header-only)

using namespace Eigen;

VectorXd markowitzOptimization(const MatrixXd& covariance, const VectorXd& returns, double targetReturn) {
    int n = returns.size();

    // Define QP matrices: 0.5 * w^T Σ w (quadratic objective)
    MatrixXd Q = covariance;
    VectorXd c = VectorXd::Zero(n);

    // Constraints: return = targetReturn, sum(w) = 1
    MatrixXd A(2, n);
    A << returns.transpose(), VectorXd::Ones(n).transpose();
    VectorXd b(2);
    b << targetReturn, 1.0;

    // Solve the system using Lagrange multipliers
    MatrixXd M(n + 2, n + 2);
    M << Q, A.transpose(), 
         A, MatrixXd::Zero(2, 2);

    VectorXd rhs(n + 2);
    rhs << -c, b;

    VectorXd solution = M.fullPivHouseholderQr().solve(rhs);
    return solution.head(n);
}

int main() {
    // Example: 3 assets
    MatrixXd covariance(3, 3);
    covariance << 0.1, 0.02, 0.01,
                  0.02, 0.2, 0.03,
                  0.01, 0.03, 0.15;
    VectorXd returns(3);
    returns << 0.08, 0.12, 0.06;

    VectorXd weights = markowitzOptimization(covariance, returns, 0.1);
    std::cout << "Optimal weights:\n" << weights << std::endl;
    return 0;
}