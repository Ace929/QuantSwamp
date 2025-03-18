#include <iostream>
#include <Eigen/Dense>
#include <vector>
using namespace Eigen;

class NeuralNetwork {
public:
    NeuralNetwork(int inputSize, int hiddenSize, int outputSize) 
        : W1(MatrixXd::Random(inputSize, hiddenSize)), 
          W2(MatrixXd::Random(hiddenSize, outputSize)) {}

    MatrixXd forward(const MatrixXd& X) {
        Z = X * W1;
        A = Z.unaryExpr([](double x) { return 1.0 / (1.0 + exp(-x)); }); // Sigmoid
        return A * W2;
    }

    void sgdStep(const MatrixXd& X, const MatrixXd& y, double lr=0.01) {
        MatrixXd yPred = forward(X);
        MatrixXd loss = yPred - y;

        // Backpropagation
        MatrixXd dW2 = A.transpose() * loss;
        MatrixXd dW1 = X.transpose() * (loss * W2.transpose()).cwiseProduct(A.cwiseProduct(1 - A));

        W2 -= lr * dW2;
        W1 -= lr * dW1;
    }

private:
    MatrixXd W1, W2, Z, A;
};

int main() {
    // Example: Learn y = 2x + 1 with noise
    NeuralNetwork nn(1, 4, 1);
    MatrixXd X(100, 1), y(100, 1);
    for (int i = 0; i < 100; ++i) {
        X(i) = i / 10.0;
        y(i) = 2 * X(i) + 1 + 0.1 * rand()/RAND_MAX;
    }

    for (int epoch = 0; epoch < 1000; ++epoch) {
        nn.sgdStep(X, y, 0.01);
    }
    std::cout << "Prediction at x=5: " << nn.forward(MatrixXd::Constant(1,1,5.0)) << std::endl;
    return 0;
}