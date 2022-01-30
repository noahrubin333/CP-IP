/**
    Noah Rubin
    2020-06-29

    The following model implements (IP1) from the paper "An LP based proof of the nonexistence of a pair of orthogonal latin squares of order 6"
    except generalized to search for 3 mutually orthogonal Latin squares

    The order of the Latin squares to search for needs to be given as the macro definition N
**/


#include <iostream>
#include "gurobi_c++.h"
#include <string>

using namespace std;

int main(int argc, char* argv[]) {
    try {   

        // Setup environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "LatinSquare.log");
        env.start();

        GRBModel model = GRBModel(env);
        GRBModel TransversalFixer_Model = GRBModel(env);


        // Declare variables
        GRBVar x[N][N][N][N][N];
        GRBVar transversal_fixer_vars[N][N][N];

        //fix_transversals(transversal_fixer_vars, TransversalFixer_Model);

        int i = 0, j = 0, k = 0, l = 0, m = 0;

        // Instantiate variables
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        for (m = 0; m < N; m++) {
                            string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l) + "_" + to_string(m);
                            x[i][j][k][l][m] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                        }
                    }
                }
            }
        }

        // Implement constraint (1)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                GRBLinExpr expr = 0;
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        for (m = 0; m < N; m++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "a_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (2)
        for (i = 0; i < N; i++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    for (l = 0; l < N; l++) {
                        for (m = 0; m < N; m++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "b_" + to_string(i) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (3)
        for (i = 0; i < N; i++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    for (k = 0; k < N; k++) {
                        for (m = 0; m < N; m++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "c_" + to_string(i) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (4)
        for (i = 0; i < N; i++) {
            for (m = 0; m < N; m++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    for (k = 0; k < N; k++) {
                        for (l = 0; l < N; l++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "d_" + to_string(i) + "_" + to_string(m);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (5)
        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (l = 0; l < N; l++) {
                        for (m = 0; m < N; m++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "e_" + to_string(j) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (6)
        for (j = 0; j < N; j++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (k = 0; k < N; k++) {
                        for (m = 0; m < N; m++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "f_" + to_string(j) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (7)
        for (j = 0; j < N; j++) {
            for (m = 0; m < N; m++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (k = 0; k < N; k++) {
                        for (l = 0; l < N; l++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "g_" + to_string(j) + "_" + to_string(m);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (8)
        for (k = 0; k < N; k++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (j = 0; j < N; j++) {
                        for (m = 0; m < N; m++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "h_" + to_string(k) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (9)
        for (k = 0; k < N; k++) {
            for (m = 0; m < N; m++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (j = 0; j < N; j++) {
                        for (l = 0; l < N; l++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "i_" + to_string(k) + "_" + to_string(m);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (10)
        for (l = 0; l < N; l++) {
            for (m = 0; m < N; m++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (j = 0; j < N; j++) {
                        for (k = 0; k < N; k++) {
                            expr += x[i][j][k][l][m];
                        }
                    }
                }
                string s = "j_" + to_string(l) + "_" + to_string(m);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }


        // Implement Symmetry Breaking


        // Fix first column of L1
        for (i = 1; i < N; i++) {
            for (j = 0; j < N; j++) {
                if(i != j) {
                    for (k = 0; k < N; k++) {
                        for (m = 0; m < N; m++) {
                            GRBLinExpr expr = 0;
                            expr += x[i][0][j][k][m];
                            model.addConstr(expr == 0);
                        }
                    }
                }
            }
        }

#ifdef FIXCOL
        // Fix domains of first column of L2
        for (i = 1; i < N; i++) {
            for (j = i; j < N; j++) {
                for (k = 1; k < N; k++) {
                    for (m = 0; m < N; m++) {
                        if (i+1 != j) {
                            GRBLinExpr expr = 0;
                            expr += x[i][0][k][j][m];
                            model.addConstr(expr == 0);
                        }
                    }
                }
            }
        }
#else
        // Fix domains of first column of L2
        for (i = 1; i < N; i++) {
            for (j = i; j < N; j++) {
                for (k = 1; k < N; k++) {
                    for (m = 0; m < N; m++) {
                        //cant be i or > i + 1
                        if (j != i + 1) {
                            GRBLinExpr expr = 0;
                            expr += x[i][0][k][j][m];
                            model.addConstr(expr == 0);
                        }
                    }
                }
            }
        }
#endif

        // Fix first row of both squares
        for (i = 0; i < N; i++) {
            GRBLinExpr expr = 0;
            expr += x[0][i][i][i][i];
            model.addConstr(expr == 1.0);
        }

        // Fix values of entries at (1,0) to be (1,2)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (!(i == 1 && j == 2)) {
                    for (k = 0; k < N; k++) {
                        GRBLinExpr expr = 0;
                        expr += x[1][0][i][j][k];
                        model.addConstr(expr == 0.0);
                    }
                }
            }
        }

        // Fix values of entries at (N-2,0) to be (N-2,N-1)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if (!(i == N-2 && j == N-1)) {
                    for (k = 0; k < N; k++) {
                        GRBLinExpr expr = 0;
                        expr += x[N-2][0][i][j][k];
                        model.addConstr(expr == 0.0);
                    }
                }
            }
        }
        
        model.set(GRB_IntParam_SolutionLimit, 1);
        model.set(GRB_IntParam_MIPFocus, 1);
        model.set(GRB_IntParam_Threads, 1);
        if(argc >= 2) {
            int seed = stoi(argv[1]);
            model.set(GRB_IntParam_Seed, seed);
            cout << "Using random seed " << seed << endl;
        }
        model.set("GURO_PAR_BINHEUR", "0");
        
        model.write("MOLS_n=" + to_string(N) + ".lp");

        model.optimize();

        cout << endl;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        for (m = 0; m < N; m++) {
                            if (x[i][j][k][l][m].get(GRB_DoubleAttr_X) > 0.5)
                                cout << to_string(k) << " ";
                        }
                    }
                }
            }
            cout << endl;
        }

        cout << endl;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        for (m = 0; m < N; m++) {
                            if (x[i][j][k][l][m].get(GRB_DoubleAttr_X) > 0.5)
                                cout << to_string(l) << " ";
                        }
                    }
                }
            }
            cout << endl;
        }
        
        cout << endl;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        for (m = 0; m < N; m++) {
                            if (x[i][j][k][l][m].get(GRB_DoubleAttr_X) > 0.5)
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
