#include <iostream>
#include <vector>
#include <cmath>
#include <random>
using namespace std;

double europeanOptionMC(double S0, double K, double r, double sigma, double T, int steps) {
    mt19937 gen(42);
    normal_distribution<double> dist(0, 1);
    double dt = T / steps;
    double ST = S0;
    for (int i = 0; i < steps; ++i) {
        ST *= exp((r - 0.5*sigma*sigma)*dt + sigma*sqrt(dt)*dist(gen));
    }
    return exp(-r*T) * max(ST - K, 0.0);
}

double mlmcEstimator(int L, int M0, double S0, double K, double r, double sigma, double T) {
    vector<int> N(L+1), steps(L+1);
    vector<double> sumY(L+1, 0.0), cost(L+1, 0.0);

    // Configure levels (geometric refinement)
    steps[0] = 1;
    for (int l = 1; l <= L; ++l) steps[l] = 2 * steps[l-1];
    N[0] = M0;
    for (int l = 1; l <= L; ++l) N[l] = N[l-1] / 4; // Heuristic sample count

    // Compute MLMC estimator
    double total = 0.0;
    for (int l = 0; l <= L; ++l) {
        double Pfine = 0.0, Pcoarse = 0.0;
        for (int i = 0; i < N[l]; ++i) {
            Pfine = europeanOptionMC(S0, K, r, sigma, T, steps[l]);
            if (l > 0) {
                Pcoarse = europeanOptionMC(S0, K, r, sigma, T, steps[l-1]);
                sumY[l] += Pfine - Pcoarse;
            } else {
                sumY[l] += Pfine;
            }
        }
        total += sumY[l] / N[l];
    }
    return total;
}

int main() {
    double price = mlmcEstimator(4, 1000, 100.0, 100.0, 0.05, 0.2, 1.0);
    cout << "MLMC European Call Price: " << price << endl;
    return 0;
}