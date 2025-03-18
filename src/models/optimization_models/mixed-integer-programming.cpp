#include <iostream>
#include <gurobi_c++.h> // Gurobi MIP solver

using namespace std;

int main() {
    try {
        GRBEnv env = GRBEnv();
        GRBModel model = GRBModel(env);

        // Decision Variables
        GRBVar x1 = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER, "x1"); // Product A units
        GRBVar x2 = model.addVar(0, GRB_INFINITY, 0, GRB_INTEGER, "x2"); // Product B units
        GRBVar y1 = model.addVar(0, 1, 0, GRB_BINARY, "y1"); // 1 if A is produced
        GRBVar y2 = model.addVar(0, 1, 0, GRB_BINARY, "y2"); // 1 if B is produced

        // Objective Function: Maximize Profit
        GRBLinExpr profit = 5 * x1 + 4 * x2;
        model.setObjective(profit, GRB_MAXIMIZE);

        // Constraints
        model.addConstr(2 * x1 + 3 * x2 <= 8, "Labor");
        model.addConstr(4 * x1 + x2 <= 6, "Materials");

        // Linking Constraints: Production only if selected
        model.addConstr(x1 <= 10 * y1, "Enable_x1");
        model.addConstr(x2 <= 10 * y2, "Enable_x2");

        // Solve Model
        model.optimize();

        // Output Results
        cout << "Optimal Production Plan:\n";
        cout << "x1 (Product A) = " << x1.get(GRB_DoubleAttr_X) << endl;
        cout << "x2 (Product B) = " << x2.get(GRB_DoubleAttr_X) << endl;
        cout << "y1 (Produce A?) = " << y1.get(GRB_DoubleAttr_X) << endl;
        cout << "y2 (Produce B?) = " << y2.get(GRB_DoubleAttr_X) << endl;
        cout << "Maximum Profit = " << model.get(GRB_DoubleAttr_ObjVal) << endl;

    } catch (GRBException e) {
        cout << "Error: " << e.getMessage() << endl;
    } catch (...) {
        cout << "Unknown error occurred." << endl;
    }

    return 0;
}
