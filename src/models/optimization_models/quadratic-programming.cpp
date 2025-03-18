#include <iostream>
#include <vector>
#include <glpk.h> // GLPK for LP
#include <gurobi_c++.h> // Gurobi for QP

using namespace std;

int main() {
    try {
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // Decision Variables: Investment in two assets (x1, x2)
        GRBVar x1 = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x1");
        GRBVar x2 = model.addVar(0.0, 1.0, 0.0, GRB_CONTINUOUS, "x2");

        // Quadratic Risk Matrix (Covariances)
        double Q[2][2] = {
            {0.02, 0.01},  // Covariance between x1, x1 and x1, x2
            {0.01, 0.03}   // Covariance between x2, x1 and x2, x2
        };

        // Set Quadratic Objective Function: Minimize (1/2) x^T Q x
        GRBQuadExpr objective;
        objective += 0.5 * Q[0][0] * x1 * x1;
        objective += 0.5 * Q[1][1] * x2 * x2;
        objective += Q[0][1] * x1 * x2; // Cross-term

        model.setObjective(objective, GRB_MINIMIZE);

        // Constraint: Total investment must sum to 1
        model.addConstr(x1 + x2 == 1, "budget");

        // Optimize model
        model.optimize();

        // Output results
        cout << "Optimal Portfolio Allocation:\n";
        cout << "x1 (Asset 1) = " << x1.get(GRB_DoubleAttr_X) << endl;
        cout << "x2 (Asset 2) = " << x2.get(GRB_DoubleAttr_X) << endl;
        cout << "Minimum Risk Value = " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    } catch (GRBException e) {
        cout << "Error: " << e.getMessage() << endl;
    } catch (...) {
        cout << "Unknown error occurred." << endl;
    }

    return 0;
}
