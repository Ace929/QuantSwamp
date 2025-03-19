#include <iostream>
#include <Eigen/Dense>
#include <cmath>
#include <random>
using namespace Eigen;

class HMC {
public:
    HMC(std::function<double(const VectorXd&)> U, // Potential energy (negative log-posterior)
        std::function<VectorXd(const VectorXd&)> gradU, // Gradient of potential
        double stepSize, int nLeapfrog)
        : U(U), gradU(gradU), stepSize(stepSize), nLeapfrog(nLeapfrog), gen(42), dist(0, 1) {}

    VectorXd sample(const VectorXd& initial) {
        VectorXd q = initial;
        VectorXd p = VectorXd::Random(initial.size()); // Momentum
        double currentU = U(q);
        VectorXd currentGrad = gradU(q);

        // Simulate Hamiltonian dynamics (leapfrog integrator)
        VectorXd qNew = q;
        VectorXd pNew = p - 0.5 * stepSize * currentGrad;
        for (int i = 0; i < nLeapfrog; ++i) {
            qNew += stepSize * pNew;
            if (i != nLeapfrog - 1) {
                VectorXd grad = gradU(qNew);
                pNew -= stepSize * grad;
            }
        }
        pNew -= 0.5 * stepSize * gradU(qNew);

        // Metropolis acceptance
        double proposedU = U(qNew);
        double deltaH = proposedU - currentU + 0.5 * (pNew.squaredNorm() - p.squaredNorm());
        if (log(dist(gen)) < -deltaH) {
            q = qNew;
        }
        return q;
    }

private:
    std::function<double(const VectorXd&)> U;
    std::function<VectorXd(const VectorXd&)> gradU;
    double stepSize;
    int nLeapfrog;
    std::mt19937 gen;
    std::uniform_real_distribution<double> dist;
};

// Example: Sampling from a 2D Gaussian
int main() {
    VectorXd mu(2); mu << 1.0, 2.0;
    MatrixXd Sigma(2,2); Sigma << 1.0, 0.5, 0.5, 1.0;
    MatrixXd SigmaInv = Sigma.inverse();

    auto U = [&](const VectorXd& x) { return 0.5 * (x - mu).dot(SigmaInv * (x - mu)); };
    auto gradU = [&](const VectorXd& x) { return SigmaInv * (x - mu); };

    HMC hmc(U, gradU, 0.1, 10);
    VectorXd sample = hmc.sample(VectorXd::Zero(2));
    std::cout << "HMC sample: " << sample.transpose() << std::endl;
    return 0;
}