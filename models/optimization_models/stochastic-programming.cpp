#include <iostream>
#include <vector>
#include <random>
#include <gurobi_c++.h> // Gurobi for stochastic programming

using namespace std;

// Function to simulate demand (normal distribution)
vector<pair<double, double>> generateDemandScenarios(int scenarios, double meanA, double stdA, double meanB, double stdB) {
    vector<pair<double, double>> demandScenarios;
    random_device rd;
    mt19937 generator(rd());
    normal_distribution<double> distA(meanA, stdA);
    normal_distribution<double> distB(meanB, stdB);

    for (int i = 0; i < scenarios; i++) {
        double d1 = max(0.0, distA(generator)); // Demand for A
        double d2 = max(0.0, distB(generator)); // Demand for B
        demandScenarios.push_back({d1, d2});
    }
    return demandScenarios;
}

int main() {
    try {
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        int numScenarios = 100; // Monte Carlo scenarios
        double meanDemandA = 5, stdDevA = 2;
        double meanDemandB = 4, stdDevB = 1.5;
        
        vector<pair<double, double>> demandScenarios = generateDemandScenarios(numScenarios, meanDemandA, stdDevA, meanDemandB, stdDevB);

        // Decision Variables
        GRBVar x1 = model.addVar(0, 10, 0, GRB_CONTINUOUS, "x1"); // Production of A
        GRBVar x2 = model.addVar(0, 8, 0, GRB_CONTINUOUS, "x2"); // Production of B

        vector<GRBVar> y1(numScenarios);
        vector<GRBVar> y2(numScenarios);
        
        for (int i = 0; i < numScenarios; i++) {
            y1[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y1_" + to_string(i));
            y2[i] = model.addVar(0, GRB_INFINITY, 0, GRB_CONTINUOUS, "y2_" + to_string(i));
        }

        // Objective: Maximize Expected Profit
        GRBQuadExpr objective;
        for (int i = 0; i < numScenarios; i++) {
            objective += (5 * y1[i] + 4 * y2[i]) / numScenarios;
        }
        objective -= 2 * x1 + 3 * x2
