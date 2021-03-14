/**
    Noah Rubin
    2020-06-29

    The following model implements (IP1) from the paper "An LP based proof of the nonexistence of a pair of orthogonal latin squares of order 6"

    Note that the model is infeasible - implying that there is no pair of orthogonal latin squares of order 6

**/


#include <iostream>
#include "gurobi_c++.h"
#include <math.h>

#define n 6

using namespace std;

int main(int argc, char* argv[]) {
    try {   

        // Setup environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "LatinSquare.log");
        env.start();

        GRBModel model = GRBModel(env);

        // Declare variables
        GRBVar x[n][n][n][n];

        int i = 0, j = 0, k = 0, l = 0;

        // Instantiate variables
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                        x[i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                    }
                }
            }
        }

        // Implement constraint (1)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                GRBLinExpr expr = 0;
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        expr += x[i][j][k][l];
                    }
                }
                string s = "1_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }


        // Implement constraint (2)

        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++) {
                    for (l = 0; l < n; l++) {
                        expr += x[i][j][k][l];
                    }
                }
                string s = "2_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (3)

        for (i = 0; i < n; i++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++) {
                    for (k = 0; k < n; k++) {
                        expr += x[i][j][k][l];
                    }
                }
                string s = "3_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (4)

        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (l = 0; l < n; l++) {
                        expr += x[i][j][k][l];
                    }
                }
                string s = "4_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (5)

        for (j = 0; j < n; j++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (k = 0; k < n; k++) {
                        expr += x[i][j][k][l];
                    }
                }
                string s = "5_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (6)

        for (k = 0; k < n; k++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        expr += x[i][j][k][l];
                    }
                }
                string s = "6_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Generate objective function as sum x_ijkl forall i,j,k,l in {0, ..., n - 1}

        GRBLinExpr objective = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        objective += x[i][j][k][l];
                    }
                }
            }
        }

        model.setObjective(objective);

        //model.set(GRB_IntParam_PoolSolutions, 1024);

        model.set(GRB_DoubleParam_TimeLimit, 1000000);
        model.set(GRB_DoubleParam_NodeLimit, 1000000);
        model.set(GRB_DoubleParam_IterationLimit, 10000000);

        model.optimize();

        cout << endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        if (x[i][j][k][l].get(GRB_DoubleAttr_X) > 0.5)
                            cout << to_string(k) << " ";
                    }
                }
            }
            cout << endl;
        }

        cout << endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        if (x[i][j][k][l].get(GRB_DoubleAttr_X) > 0.5)
                            cout << to_string(l) << " ";
                    }
                }
            }
            cout << endl;
        }
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }

    return 0;
}



