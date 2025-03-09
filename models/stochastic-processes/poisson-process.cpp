#include <iostream>
#include <vector>
#include <random>

using namespace std;

// Function to simulate a Poisson Process
vector<double> simulatePoissonProcess(double lambda, double T, int seed = 42) {
    vector<double> event_times;
    double time = 0.0;

    // Random number generator for exponential inter-arrival times
    random_device rd;
    mt19937 generator(seed);
    exponential_distribution<double> exp_dist(lambda);

    while (time < T) {
        double inter_arrival_time = exp_dist(generator);
        time += inter_arrival_time;
        if (time < T) {
            event_times.push_back(time);
        }
    }
    return event_times;
}

int main() {
    double lambda = 2.0; // Average events per time unit (e.g., 2 events per second)
    double T = 10.0;     // Total time to simulate

    vector<double> event_times = simulatePoissonProcess(lambda, T);

    // Print event times
    cout << "Poisson Process Event Times:\n";
    for (size_t i = 0; i < event_times.size(); ++i) {
        cout << "Event " << i + 1 << " at time " << event_times[i] << "\n";
    }

    return 0;
}
