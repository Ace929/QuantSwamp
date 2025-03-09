#include <iostream>
#include <cmath>

// standard normal cumulative distribution function
double normalCDF(double x) {
    return 0.5 * std::erfc(-x * std::sqrt(0.5));
}

void blackScholes(double S, double K, double r, double sigma, double T) {
    double d1 = (std::log(S / K) + (r + 0.5 * sigma * sigma) * T) / (sigma * std::sqrt(T));
    double d2 = d1 - sigma * std::sqrt(T);

    double callPrice = S * normalCDF(d1) - K * std::exp(-r * T) * normalCDF(d2);
    double putPrice = K * std::exp(-r * T) * normalCDF(-d2) - S * normalCDF(-d1);

    std::cout << "Black-Scholes Option Pricing:\n";
    std::cout << "Call Option Price: " << callPrice << std::endl;
    std::cout << "Put Option Price: " << putPrice << std::endl;
}

int main() {
    double S, K, r, sigma, T;

    std::cout << "Enter stock price (S): ";
    std::cin >> S;
    std::cout << "Enter strike price (K): ";
    std::cin >> K;
    std::cout << "Enter risk-free interest rate (r) in decimal (e.g., 0.05 for 5%): ";
    std::cin >> r;
    std::cout << "Enter volatility (sigma) in decimal (e.g., 0.2 for 20%): ";
    std::cin >> sigma;
    std::cout << "Enter time to maturity (T) in years: ";
    std::cin >> T;

    blackScholes(S, K, r, sigma, T);

    return 0;
}
