/**
    Noah Rubin
    2020-06-29
    The following model implements (IP1) from the paper "An LP based proof of the nonexistence of a pair of orthogonal latin squares of order 6"
    Note that the model is infeasible - implying that there is no pair of orthogonal latin squares of order 6
**/
#include <iostream>
#include <math.h>
#include "gurobi_c++.h"
using namespace std;
int main(int argc, char* argv[]) {
    
    int n = 5;
    int maxslack = 25;
    if (argc >= 2) {
        n = atoi(argv[1]);
        cout << "n = " << n << "\n";
    }
    if (argc >= 3) {
        maxslack = atoi(argv[2]);
        cout << "maxslack = " << maxslack << "\n";
    }
    
    int nFifth = (int)pow(n, 5), nForth = (int)pow(n, 4), nCube = (int)pow(n, 3), nSqr = (int)pow(n, 2);

    
    try {
        // Setup environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "LatinSquare.log");
        env.start();
        GRBModel model = GRBModel(env);
        // Declare variables
        GRBVar *x = new GRBVar[nFifth];
        GRBVar *slack = new GRBVar[nSqr];

        int i = 0, j = 0, k = 0, l = 0, m = 0;
        // Instantiate variables
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        for (m = 0; m < n; m++) {
                            string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l) + "_" + to_string(m);
                            x[i * nForth + j * nCube + k * nSqr + l * n + m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                        }
                    }
                }
            }
        }

        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                string name = "s_" + to_string(i) + "_" + to_string(j);
                slack[i * n + j] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
            }
        }

        // Implement constraint (1)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                GRBLinExpr expr = 0;
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        for (m = 0; m < n; m++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "a_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }
        // Implement constraint (2)
        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++) {
                    for (l = 0; l < n; l++) {
                        for (m = 0; m < n; m++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "b_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }
        // Implement constraint (3)
        for (i = 0; i < n; i++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++) {
                    for (k = 0; k < n; k++) {
                        for (m = 0; m < n; m++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "c_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }
        // Implement constraint (4)
        for (i = 0; i < n; i++) {
            for (m = 0; m < n; m++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++) {
                    for (k = 0; k < n; k++) {
                        for (l = 0; l < n; l++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "d_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }
        // Implement constraint (5)
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (l = 0; l < n; l++) {
                        for (m = 0; m < n; m++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "e_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }
        // Implement constraint (6)
        for (j = 0; j < n; j++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (k = 0; k < n; k++) {
                        for (m = 0; m < n; m++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "f_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }
        // Implement constraint (7)
        for (j = 0; j < n; j++) {
            for (m = 0; m < n; m++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (k = 0; k < n; k++) {
                        for (l = 0; l < n; l++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "g_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }
        // Implement constraint (8)
        for (k = 0; k < n; k++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        for (m = 0; m < n; m++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "h_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr == 1.0, s);
            }
        }
        // Implement constraint (9)
        for (k = 0; k < n; k++) {
            for (m = 0; m < n; m++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        for (l = 0; l < n; l++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "i_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr == 1.0, s);
            }
        }


        // Implement constraint (10)
        for (l = 0; l < n; l++) {
            for (m = 0; m < n; m++) {
                GRBLinExpr expr = slack[l * n + m];
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        for (k = 0; k < n; k++) {
                            expr += x[i * nForth + j * nCube + k * nSqr + l * n + m];
                        }
                    }
                }
                string s = "j_" + to_string(i) + "_" + to_string(j);
                if(l == m)
                    model.addConstr(expr >= 1.0, s);
                else
                    model.addConstr(expr == 1.0, s);
            }
        }
        // Implement Symmetry Breaking

        /**
        // Fix first column of L1 and domains of first column of L2
        for (i = 1; i < n; i++) {
            for (j = i; j < n; j++) {
                for (k = 1; k < n; k++) {
                    for (m = 0; m < n; m++) {
                        //cant be i or > i + 1
                        if (j == i || j > i + 1 || k != i) {
                            GRBLinExpr expr = 0;
                            expr += x[i][0][k][j][m];
                            model.addConstr(expr == 0);
                        }
                    }
                }
            }
        }
        **/

        // Fix first row of both squares
        for (i = 0; i < n; i++) {
            GRBLinExpr expr = 0;
            expr += x[i * nCube + i * nSqr + i * n + i];
            model.addConstr(expr == 1.0);
        }

        GRBLinExpr sumSlack;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                sumSlack += slack[i * n + j];
            }
        }

        model.addConstr(sumSlack <= maxslack);

        // Maybe test that x_10411 = x_23111 - pair (1,1) appears twice accross L and M?
        //model.addConstr(x[nForth + 0 * nCube + 4 * nSqr + n + 1] == 1);
        //model.addConstr(x[2 * nForth + 1 * nCube + 1 * nSqr + n + 1] == 1);


        model.set(GRB_IntParam_Threads, 1);

        model.write("C:\\Users\\Noah\\3MOLS_n=" + to_string(n) + ".lp");
        model.optimize();
        cout << endl;
        

        int s = 0;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                if (slack[i * n + j].get(GRB_DoubleAttr_X) > 0.5) {
                    cout << "Slack " << i << ", " << j << " = 1\n";
                    s += 1;
                }
            }
        }
        cout << "\nsumslack = " << s << "\n" << endl;
        
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        for (m = 0; m < n; m++) {
                            if (x[i * nForth + j * nCube + k * nSqr + l * n + m].get(GRB_DoubleAttr_X) > 0.5)
                                cout << to_string(k) << " ";
                        }
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
                        for (m = 0; m < n; m++) {
                            if (x[i * nForth + j * nCube + k * nSqr + l * n + m].get(GRB_DoubleAttr_X) > 0.5)
                                cout << to_string(l) << " ";
                        }
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
                        for (m = 0; m < n; m++) {
                            if (x[i * nForth + j * nCube + k * nSqr + l * n + m].get(GRB_DoubleAttr_X) > 0.5)
                                cout << to_string(m) << " ";
                        }
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
