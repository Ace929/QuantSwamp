#include <iostream>
#include <vector>
#include <cmath>
#include <limits>
#include <random>

using namespace std;

// Parameters
const int T = 5;  // Time periods
const int S_max = 10; // Maximum inventory
const int X_max = 5;  // Max order size
const double h = 1.0; // Holding cost per unit
const double c = 2.0; // Ordering cost per unit
const double p = 5.0; // Penalty for unmet demand

// Demand distribution (Normal)
random_device rd;
mt19937 generator(rd());
normal_distribution<double> demand_distribution(3.0, 1.0);

// Function to compute cost
double cost(int s, int x, int d) {
    int new_inventory = max(0, s + x - d);
    int unmet_demand = max(0, d - s - x);
    return h * new_inventory + c * x + p * unmet_demand;
}

// DP Table (state x time)
vector<vector<double>> V(S_max + 1, vector<double>(T + 1, numeric_limits<double>::infinity()));
vector<vector<int>> policy(S_max + 1, vector<int>(T, 0)); // Stores optimal order decisions

void solveDP() {
    // Base Case: Final period cost is zero
    for (int s = 0; s <= S_max; ++s) {
        V[s][T] = 0;
    }

    // Backward induction
    for (int t = T - 1; t >= 0; --t) {
        for (int s = 0; s <= S_max; ++s) {
            for (int x = 0; x <= X_max; ++x) {
                double expected_cost = 0;
                for (int i = 0; i < 10; ++i) { // Monte Carlo sampling for demand
                    int d = max(0, (int)round(demand_distribution(generator)));
                    d = min(d, S_max); // Cap demand
                    expected_cost += cost(s, x, d) + V[min(S_max, s + x - d)][t + 1];
                }
                expected_cost /= 10; // Average over scenarios

                // Choose minimum cost action
                if (expected_cost < V[s][t]) {
                    V[s][t] = expected_cost;
                    policy[s][t] = x;
                }
            }
        }
    }
}

void printResults() {
    cout << "Optimal Order Policy:\n";
    for (int t = 0; t < T; ++t) {
        cout << "Time " << t << ":\n";
        for (int s = 0; s <= S_max; ++s) {
            cout << "  Inventory " << s << " → Order " << policy[s][t] << "\n";
        }
    }
}

int main() {
    solveDP();
    printResults();
    return 0;
}
