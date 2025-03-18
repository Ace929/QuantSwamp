#include <iostream>
#include <vector>
#include <cmath>
#include <complex>
#include <fftw3.h> // Requires FFTW library

double fftOptionPrice(double S0, double K, double r, double T, double sigma, double alpha=0.25, int N=4096) {
    double eta = 0.25; // Grid spacing in log-strike space
    double lambda = 2 * M_PI / (N * eta);
    double delta = 1e-4; // Damping factor

    std::vector<std::complex<double>> phi(N), integrand(N);
    fftw_plan plan = fftw_plan_dft_1d(N, reinterpret_cast<fftw_complex*>(&integrand[0]),
                                      reinterpret_cast<fftw_complex*>(&phi[0]), FFTW_FORWARD, FFTW_ESTIMATE);

    // Characteristic function of log-price (Black-Scholes)
    for (int j = 0; j < N; ++j) {
        double v_j = (j - N/2) * eta;
        double k = log(K/S0);
        std::complex<double> i(0, 1);
        std::complex<double> exponent = i*(v_j - (alpha+1)*i)*k - 0.5*sigma*sigma*T*(v_j - alpha*i)*(v_j - alpha*i);
        integrand[j] = exp(exponent) * exp(-r*T) / (alpha*alpha + alpha - v_j*v_j + i*(2*alpha + 1)*v_j);
    }

    fftw_execute(plan);
    fftw_destroy_plan(plan);

    // Find the strike index and interpolate
    int index = static_cast<int>(log(K/S0)/lambda + N/2);
    double price = real(phi[index]) * exp(-alpha*log(K/S0)) / M_PI;
    return price;
}

int main() {
    double price = fftOptionPrice(100.0, 100.0, 0.05, 1.0, 0.2);
    std::cout << "FFT Call Price: " << price << std::endl;
    return 0;
}