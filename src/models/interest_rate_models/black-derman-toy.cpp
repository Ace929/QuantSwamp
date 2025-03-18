#include <iostream>
#include <vector>
#include <cmath>

using namespace std;

// BDT Parameters
const int T = 5; // Number of periods
const double r0 = 0.05; // Initial short rate
const double dt = 1.0; // Time step (years)
const vector<double> sigma = {0.2, 0.18, 0.15, 0.12, 0.1}; // Volatility for each period

// Function to generate BDT interest rate tree
vector<vector<double>> buildBDTTree(double r0, const vector<double>& sigma) {
    vector<vector<double>> tree(T);
    tree[0].push_back(r0);

    for (int t = 1; t < T; ++t) {
        int levels = t + 1;
        double prevRate = tree[t - 1][0];
        tree[t].resize(levels);

        for (int i = 0; i < levels; ++i) {
            tree[t][i] = prevRate * exp((2 * i - (t - 1)) * sigma[t - 1]);
        }
    }
    return tree;
}

// Function to print the BDT tree
void printBDTTree(const vector<vector<double>>& tree) {
    cout << "Black-Derman-Toy Interest Rate Tree:\n";
    for (size_t t = 0; t < tree.size(); ++t) {
        cout << "Time " << t << ": ";
        for (size_t i = 0; i < tree[t].size(); ++i) {
            cout << tree[t][i] << " ";
        }
        cout << endl;
    }
}

int main() {
    vector<vector<double>> bdtTree = buildBDTTree(r0, sigma);
    printBDTTree(bdtTree);
    return 0;
}
