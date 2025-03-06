#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <random>

using namespace Eigen;
using namespace std;

// Generate synthetic multivariate return data
vector<VectorXd> generate_synthetic_returns(int T, int n) {
    vector<VectorXd> returns;
    random_device rd;
    mt19937 gen(rd());
    normal_distribution<> normal(0, 1);

    for (int t = 0; t < T; ++t) {
        VectorXd r(n);
        for (int i = 0; i < n; ++i) {
            r(i) = normal(gen); // Normally distributed returns
        }
        returns.push_back(r);
    }
    return returns;
}

// M-GARCH(1,1) Model Computation (VECH-GARCH)
void m_garch(const vector<VectorXd>& returns, MatrixXd C, MatrixXd A, MatrixXd B) {
    int T = returns.size();
    int n = returns[0].size();

    // Initialize covariance matrix (Identity for simplicity)
    MatrixXd H = MatrixXd::Identity(n, n);

    // Store conditional covariance matrices
    vector<MatrixXd> H_t(T, MatrixXd::Zero(n, n));

    // Iterate over time steps
    for (int t = 0; t < T; ++t) {
        if (t > 0) {
            MatrixXd outer_product = returns[t - 1] * returns[t - 1].transpose();
            H = C + A * outer_product + B * H;
        }
        H_t[t] = H;

        // Output covariance matrix at each step
        cout << "Time t=" << t << " | Covariance Matrix H_t:\n" << H << "\n\n";
    }
}

int main() {
    int T = 100; // Number of time steps
    int n = 2;   // Number of assets

    // Generate synthetic returns
    vector<VectorXd> returns = generate_synthetic_returns(T, n);

    // Initialize M-GARCH model parameters
    MatrixXd C = MatrixXd::Identity(n, n) * 0.01; // Constant covariance
    MatrixXd A = MatrixXd::Identity(n, n) * 0.2;  // ARCH effect
    MatrixXd B = MatrixXd::Identity(n, n) * 0.75; // GARCH effect

    // Run M-GARCH model
    m_garch(returns, C, A, B);

    return 0;
}
