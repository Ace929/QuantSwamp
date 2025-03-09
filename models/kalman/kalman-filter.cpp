#include <iostream>
#include <Eigen/Dense>

using namespace Eigen;
using namespace std;

class KalmanFilter {
public:
    // State vector
    VectorXd x;
    // State covariance matrix
    MatrixXd P;
    // State transition matrix
    MatrixXd F;
    // Process covariance matrix
    MatrixXd Q;
    // Measurement matrix
    MatrixXd H;
    // Measurement covariance matrix
    MatrixXd R;

    KalmanFilter(int state_size, int measurement_size) {
        x = VectorXd::Zero(state_size);
        P = MatrixXd::Identity(state_size, state_size);
        F = MatrixXd::Identity(state_size, state_size);
        Q = MatrixXd::Identity(state_size, state_size);
        H = MatrixXd::Zero(measurement_size, state_size);
        R = MatrixXd::Identity(measurement_size, measurement_size);
    }

    void predict() {
        x = F * x; // State prediction
        P = F * P * F.transpose() + Q; // Covariance prediction
    }

    void update(const VectorXd &z) {
        VectorXd y = z - H * x; // Measurement residual
        MatrixXd S = H * P * H.transpose() + R; // Residual covariance
        MatrixXd K = P * H.transpose() * S.inverse(); // Kalman gain

        x = x + K * y; // Updated state estimate
        P = (MatrixXd::Identity(x.size(), x.size()) - K * H) * P; // Updated covariance
    }

    void printState() {
        cout << "State vector:\n" << x << endl;
        cout << "Covariance matrix:\n" << P << endl;
    }
};

int main() {
    // Define dimensions (2D position tracking example)
    int state_size = 4;       // [x, y, vx, vy]
    int measurement_size = 2; // [x, y] position measurements

    KalmanFilter kf(state_size, measurement_size);

    // Initialize state (e.g., initial position and velocity)
    kf.x << 0, 0, 1, 1; // Initial position (0,0) with velocity (1,1)

    // State transition matrix (assuming constant velocity model)
    kf.F << 1, 0, 1, 0,
            0, 1, 0, 1,
            0, 0, 1, 0,
            0, 0, 0, 1;

    // Measurement matrix (we only observe position)
    kf.H << 1, 0, 0, 0,
            0, 1, 0, 0;

    // Set measurement covariance (sensor noise)
    kf.R << 0.1, 0,
            0, 0.1;

    // Set process noise
    kf.Q << 0.01, 0, 0, 0,
            0, 0.01, 0, 0,
            0, 0, 0.01, 0,
            0, 0, 0, 0.01;

    // Simulate a sequence of measurements
    VectorXd z1(2);
    z1 << 1.1, 0.9;

    VectorXd z2(2);
    z2 << 2.0, 2.1;

    cout << "Initial State:\n";
    kf.printState();

    // Prediction step
    kf.predict();
    cout << "\nAfter Prediction:\n";
    kf.printState();

    // First update step
    kf.update(z1);
    cout << "\nAfter First Update:\n";
    kf.printState();

    // Second prediction step
    kf.predict();
    cout << "\nAfter Second Prediction:\n";
    kf.printState();

    // Second update step
    kf.update(z2);
    cout << "\nAfter Second Update:\n";
    kf.printState();

    return 0;
}
