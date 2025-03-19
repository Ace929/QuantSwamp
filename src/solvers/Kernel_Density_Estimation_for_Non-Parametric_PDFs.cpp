#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <Eigen/Dense>
using namespace Eigen;

class KDE {
public:
    KDE(const VectorXd& data) : data(data) {
        n = data.size();
        h = silvermanBandwidth();
    }

    double estimate(double x) const {
        double density = 0.0;
        for (int i = 0; i < n; ++i) {
            density += gaussianKernel((x - data[i]) / h);
        }
        return density / (n * h);
    }

private:
    VectorXd data;
    double h;
    int n;

    double silvermanBandwidth() const {
        double std = std::sqrt((data.array() - data.mean()).square().sum() / (n - 1));
        return 1.06 * std * pow(n, -0.2);
    }

    double gaussianKernel(double u) const {
        return exp(-0.5 * u * u) / sqrt(2 * M_PI);
    }
};

int main() {
    VectorXd data(1000);
    std::generate(data.data(), data.data() + data.size(), [](){ return 2.0 + random(); });
    KDE kde(data);
    std::cout << "Density at x=2.0: " << kde.estimate(2.0) << std::endl;
    return 0;
}