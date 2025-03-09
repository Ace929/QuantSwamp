#include <iostream>
#include <glpk.h> // GLPK library

using namespace std;

int main() {
    // Create a problem object
    glp_prob *lp;
    lp = glp_create_prob();
    glp_set_prob_name(lp, "Resource Allocation");
    glp_set_obj_dir(lp, GLP_MAX); // Maximization problem

    // Define variables (x1 and x2)
    glp_add_cols(lp, 2);
    glp_set_col_name(lp, 1, "x1");
    glp_set_col_kind(lp, 1, GLP_CV); // Continuous variable
    glp_set_col_bnds(lp, 1, GLP_LO, 0.0, 0.0); // x1 >= 0

    glp_set_col_name(lp, 2, "x2");
    glp_set_col_kind(lp, 2, GLP_CV); // Continuous variable
    glp_set_col_bnds(lp, 2, GLP_LO, 0.0, 0.0); // x2 >= 0

    // Define objective function: Max Z = 5x1 + 4x2
    glp_set_obj_coef(lp, 1, 5.0); // Coefficient for x1
    glp_set_obj_coef(lp, 2, 4.0); // Coefficient for x2

    // Define constraints (2 constraints)
    glp_add_rows(lp, 2);
    
    // Constraint 1: 2x1 + 3x2 <= 8 (Labor)
    glp_set_row_name(lp, 1, "Labor");
    glp_set_row_bnds(lp, 1, GLP_UP, 0.0, 8.0); // Upper bound <= 8
    int ind1[] = {0, 1, 2}; // Indexes (GLPK uses 1-based index)
    double val1[] = {0, 2.0, 3.0}; // Coefficients: 2x1 + 3x2
    glp_set_mat_row(lp, 1, 2, ind1, val1);

    // Constraint 2: 4x1 + x2 <= 6 (Materials)
    glp_set_row_name(lp, 2, "Materials");
    glp_set_row_bnds(lp, 2, GLP_UP, 0.0, 6.0); // Upper bound <= 6
    int ind2[] = {0, 1, 2}; // Indexes
    double val2[] = {0, 4.0, 1.0}; // Coefficients: 4x1 + x2
    glp_set_mat_row(lp, 2, 2, ind2, val2);

    // Solve LP problem
    glp_simplex(lp, NULL);

    // Get results
    double x1 = glp_get_col_prim(lp, 1);
    double x2 = glp_get_col_prim(lp, 2);
    double maxProfit = glp_get_obj_val(lp);

    // Display results
    cout << "Optimal solution found:" << endl;
    cout << "x1 (Product A) = " << x1 << endl;
    cout << "x2 (Product B) = " << x2 << endl;
    cout << "Maximum Profit = " << maxProfit << endl;

    // Free memory
    glp_delete_prob(lp);

    return 0;
}
