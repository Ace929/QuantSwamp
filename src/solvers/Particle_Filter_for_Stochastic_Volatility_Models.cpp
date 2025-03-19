#include <iostream>
#include <vector>
#include <random>
#include <Eigen/Dense>
using namespace Eigen;

struct Particle {
    double vol;   // Hidden state (volatility)
    double weight;
};

class ParticleFilter {
public:
    ParticleFilter(int numParticles, double vol0, double processNoise, double obsNoise)
        : numParticles(numParticles), processNoise(processNoise), obsNoise(obsNoise), gen(42) {
        particles.resize(numParticles);
        for (auto& p : particles) {
            p.vol = vol0 + 0.1 * randn();
            p.weight = 1.0 / numParticles;
        }
    }

    void update(double observedReturn) {
        // Predict step (volatility follows mean-reverting process)
        for (auto& p : particles) {
            p.vol = 0.9 * p.vol + 0.1 * randn() * sqrt(processNoise);
        }

        // Update weights (likelihood of observed return)
        double totalWeight = 0.0;
        for (auto& p : particles) {
            p.weight *= exp(-0.5 * pow(observedReturn, 2) / (2 * p.vol * p.vol + obsNoise));
            totalWeight += p.weight;
        }

        // Normalize weights
        for (auto& p : particles) p.weight /= totalWeight;

        // Resample using systematic resampling
        resample();
    }

    double estimateVol() const {
        double vol = 0.0;
        for (const auto& p : particles) vol += p.vol * p.weight;
        return vol;
    }

private:
    std::vector<Particle> particles;
    int numParticles;
    double processNoise, obsNoise;
    std::mt19937 gen;

    double randn() {
        static std::normal_distribution<double> dist(0.0, 1.0);
        return dist(gen);
    }

    void resample() {
        std::vector<Particle> newParticles(numParticles);
        double step = 1.0 / numParticles;
        double u = rand() * step;
        double c = particles[0].weight;
        int i = 0;
        for (int j = 0; j < numParticles; ++j) {
            while (u > c) c += particles[++i].weight;
            newParticles[j] = particles[i];
            newParticles[j].weight = 1.0 / numParticles;
            u += step;
        }
        particles = newParticles;
    }
};

int main() {
    ParticleFilter pf(1000, 0.2, 0.01, 0.05);
    pf.update(0.03); // Simulated return observation
    std::cout << "Estimated Volatility: " << pf.estimateVol() << std::endl;
    return 0;
}