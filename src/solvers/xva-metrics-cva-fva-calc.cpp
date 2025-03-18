#include <iostream>
#include <vector>
#include <random>
using namespace std;

struct Exposure {
    vector<double> epe;  // Expected Positive Exposure
    vector<double> eff;  // Expected Funding Exposure
};

Exposure computeExposure(int paths, int steps, double T, double r) {
    mt19937 gen(42);
    normal_distribution<double> dist(0, 1);
    Exposure exp;
    exp.epe.resize(steps);
    exp.eff.resize(steps);

    for (int sim = 0; sim < paths; ++sim) {
        double S = 100.0, dt = T / steps;
        vector<double> path(steps + 1);
        path[0] = S;
        for (int t = 1; t <= steps; ++t) {
            S *= exp((r - 0.5*0.2*0.2)*dt + 0.2*sqrt(dt)*dist(gen));
            path[t] = S;
        }
        // Compute exposures at each time step
        for (int t = 0; t < steps; ++t) {
            double payoff = max(path[t] - 100.0, 0.0); // Example: call option
            exp.epe[t] += max(payoff, 0.0) / paths;
            exp.eff[t] += payoff / paths; // Simplified funding exposure
        }
    }
    return exp;
}

double computeCVA(const Exposure& exp, double pd, double lgd) {
    double cva = 0.0;
    for (size_t t = 0; t < exp.epe.size(); ++t) {
        cva += lgd * pd * exp.epe[t]; // Assume constant PD for simplicity
    }
    return cva;
}

int main() {
    Exposure exp = computeExposure(10000, 12, 1.0, 0.05);
    double cva = computeCVA(exp, 0.02, 0.6); // PD=2%, LGD=60%
    cout << "CVA: " << cva << endl;
    return 0;
}