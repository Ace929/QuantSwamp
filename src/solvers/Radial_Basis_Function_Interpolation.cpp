#include <iostream>
#include <Eigen/Dense>
#include <vector>
#include <cmath>
using namespace Eigen;

class RBFInterpolator {
public:
    RBFInterpolator(const std::vector<Vector2d>& points, const VectorXd& values, double epsilon=1.0)
        : points(points), values(values), epsilon(epsilon) {
        // Build interpolation matrix
        int n = points.size();
        MatrixXd A(n, n);
        for (int i = 0; i < n; ++i) {
            for (int j = 0; j < n; ++j) {
                double r = (points[i] - points[j]).norm();
                A(i,j) = exp(-pow(epsilon * r, 2)); // Gaussian RBF
            }
        }
        weights = A.colPivHouseholderQr().solve(values);
    }

    double interpolate(const Vector2d& x) const {
        double result = 0.0;
        for (int i = 0; i < points.size(); ++i) {
            double r = (x - points[i]).norm();
            result += weights[i] * exp(-pow(epsilon * r, 2));
        }
        return result;
    }

private:
    std::vector<Vector2d> points;
    VectorXd weights, values;
    double epsilon;
};

int main() {
    // Example: Volatility surface interpolation
    std::vector<Vector2d> points = {{1.0, 0.1}, {1.0, 0.5}, {2.0, 0.1}}; // (Strike, Maturity)
    VectorXd values(3);
    values << 0.2, 0.25, 0.18;

    RBFInterpolator rbf(points, values, 1.0);
    std::cout << "Volatility at (1.5, 0.3): " << rbf.interpolate({1.5, 0.3}) << std::endl;
    return 0;
}