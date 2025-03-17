#include <iostream>
#include <vector>
#include <random>
#include <algorithm>

void computeRiskMetrics(const std::vector<double>& pnl, double confidence) {
    std::vector<double> sortedPnl = pnl;
    std::sort(sortedPnl.begin(), sortedPnl.end());
    int index = static_cast<int>((1 - confidence) * sortedPnl.size());
    double VaR = -sortedPnl[index]; // Loss is negative PnL
    double CVaR = 0.0;

    for (int i = 0; i <= index; ++i) {
        CVaR += -sortedPnl[i];
    }
    CVaR /= (index + 1);

    std::cout << "VaR (" << confidence * 100 << "%): " << VaR << std::endl;
    std::cout << "CVaR (" << confidence * 100 << "%): " << CVaR << std::endl;
}

int main() {
    double S0 = 100.0, mu = 0.1, sigma = 0.2, T = 1.0;
    int N = 100000;
    std::vector<double> pnl(N);
    std::mt19937 gen(42);
    std::normal_distribution<double> dist(0.0, 1.0);

    for (int sim = 0; sim < N; ++sim) {
        double Z = dist(gen);
        double ST = S0 * exp((mu - 0.5 * sigma*sigma) * T + sigma * sqrt(T) * Z);
        pnl[sim] = ST - S0; // PnL for long position
    }

    computeRiskMetrics(pnl, 0.95);
    return 0;
}