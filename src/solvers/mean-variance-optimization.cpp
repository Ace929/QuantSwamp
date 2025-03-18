#include <iostream>
#include <Eigen/Dense>
#include <vector>

using namespace Eigen;
using namespace std;

// Function to compute the efficient frontier
VectorXd meanVarianceOptimization(const MatrixXd &covMatrix, const VectorXd &expectedReturns, double targetReturn) {
    int n = covMatrix.rows();  // Number of assets

    // Construct the necessary matrices
    MatrixXd A(n + 2, n + 2);
    A.setZero();
    
    // Fill in covariance matrix in top-left
    A.block(0, 0, n, n) = covMatrix;

    // Constraint vectors
    VectorXd b(n + 2);
    b.setZero();

    // Enforce sum of weights = 1
    A.block(n, 0, 1, n) = RowVectorXd::Ones(n);
    A.block(0, n, n, 1) = VectorXd::Ones(n);
    b(n) = 1;

    // Enforce expected return = targetReturn
    A.block(n + 1, 0, 1, n) = expectedReturns.transpose();
    A.block(0, n + 1, n, 1) = expectedReturns;
    b(n + 1) = targetReturn;

    // Solve the system Ax = b
    VectorXd x = A.colPivHouseholderQr().solve(b);
    
    // Extract the optimal weights (first n values)
    return x.head(n);
}

int main() {
    int n;
    cout << "Enter the number of assets: ";
    cin >> n;

    MatrixXd covMatrix(n, n);
    VectorXd expectedReturns(n);
    double targetReturn;

    cout << "Enter the expected returns of each asset: ";
    for (int i = 0; i < n; ++i) {
        cin >> expectedReturns(i);
    }

    cout << "Enter the covariance matrix row-wise: ";
    for (int i = 0; i < n; ++i) {
        for (int j = 0; j < n; ++j) {
            cin >> covMatrix(i, j);
        }
    }

    cout << "Enter target return: ";
    cin >> targetReturn;

    VectorXd optimalWeights = meanVarianceOptimization(covMatrix, expectedReturns, targetReturn);

    cout << "Optimal portfolio weights: \n" << optimalWeights << endl;

    return 0;
}
