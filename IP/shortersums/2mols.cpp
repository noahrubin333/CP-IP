/**
    Noah Rubin
    2020-06-29

    The following model implements (IP1) from the paper "An LP based proof of the nonexistence of a pair of orthogonal latin squares of order 6"

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

        // Declare variables
        GRBVar x[N][N][N][N];
        GRBVar s1[N][N][N];
        GRBVar s2[N][N][N];
        GRBVar s3[N][N][N];

        int i = 0, j = 0, k = 0, l = 0;

        // Instantiate variables
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                        x[i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                    }
                    string name = "s1_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    s1[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                    name = "s2_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    s2[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                    name = "s3_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    s3[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                }
            }
        }

        // Implement constraint (1)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                GRBLinExpr expr = 0;
                for (k = 0; k < N; k++) {
                    GRBLinExpr expr2 = 0;
                    for (l = 0; l < N; l++) {
                        expr2 += x[i][j][k][l];
                    }
                    string s = "1a_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    model.addConstr(expr2, GRB_EQUAL, s1[i][j][k], s);
                    expr += s1[i][j][k];
                }
                string s = "1b_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }


        // Implement constraint (2)

        for (i = 0; i < N; i++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    expr += s1[i][j][k];
                }
                string s = "2_" + to_string(i) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (3)

        for (i = 0; i < N; i++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    GRBLinExpr expr2 = 0;
                    for (k = 0; k < N; k++) {
                        expr2 += x[i][j][k][l];
                    }
                    string s = "3a_" + to_string(i) + "_" + to_string(l) + "_" + to_string(j);
                    model.addConstr(expr2, GRB_EQUAL, s2[i][l][j], s);
                    expr += s2[i][l][j];
                }
                string s = "3b_" + to_string(i) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (4)

        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    expr += s1[i][j][k];
                }
                string s = "4_" + to_string(j) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (5)

        for (j = 0; j < N; j++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    expr += s2[i][j][l];
                }
                string s = "5_" + to_string(j) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (6)

        for (k = 0; k < N; k++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    GRBLinExpr expr2 = 0;
                    for (j = 0; j < N; j++) {
                        expr2 += x[i][j][k][l];
                    }
                    string s = "6a_" + to_string(k) + "_" + to_string(l) + "_" + to_string(i);
                    model.addConstr(expr2, GRB_EQUAL, s3[k][l][i], s);
                    expr += s3[k][l][i];
                }
                string s = "6b_" + to_string(k) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement Symmetry Breaking

        // Fix first row of both squares
        for (i = 0; i < N; i++) {
            GRBLinExpr expr = 0;
            expr += x[0][i][i][i];
            model.addConstr(expr == 1.0);
        }

#ifdef FIXCOL
        // Fix first column of both squares
        for (i = 1; i < N-1; i++) {
            GRBLinExpr expr = 0;
            expr += x[i][0][i][i+1];
            model.addConstr(expr == 1.0);
        }

        // Fix values of entries at (N-1,0) to be (N-1,1)
        {
            GRBLinExpr expr = 0;
            expr += x[N-1][0][N-1][1];
            model.addConstr(expr == 1.0);
        }
#else
        // Fix first column of L1
        for (i = 1; i < N; i++) {
            for (j = 0; j < N; j++) {
                if(i != j) {
                    for (k = 0; k < N; k++) {
                        GRBLinExpr expr = 0;
                        expr += x[i][0][j][k];
                        model.addConstr(expr == 0);
                    }
                }
            }
        }

        // Fix domains of first column of L2
        for (i = 1; i < N; i++) {
            for (j = i; j < N; j++) {
                for (k = 1; k < N; k++) {
                    //cant be i or > i + 1
                    if (j != i + 1) {
                        GRBLinExpr expr = 0;
                        expr += x[i][0][k][j];
                        model.addConstr(expr == 0);
                    }
                }
            }
        }

        // Fix values of entries at (1,0) to be (1,2)
        {
            GRBLinExpr expr = 0;
            expr += x[1][0][1][2];
            model.addConstr(expr == 1.0);
        }

        // Fix values of entries at (N-2,0) to be (N-2,N-1)
        {
            GRBLinExpr expr = 0;
            expr += x[N-2][0][N-2][N-1];
            model.addConstr(expr == 1.0);
        }
#endif

        // Generate objective function as sum x_ijkl forall i,j,k,l in {0, ..., n - 1}

        GRBLinExpr objective = 0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        objective += x[i][j][k][l];
                    }
                }
            }
        }

        model.set(GRB_IntParam_Threads, 1);
        if(argc >= 2) {
            int seed = stoi(argv[1]);
            model.set(GRB_IntParam_Seed, seed);
            cout << "Using random seed " << seed << endl;
        }

        model.optimize();

        cout << endl;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        if (x[i][j][k][l].get(GRB_DoubleAttr_X) > 0.5)
                            cout << to_string(k) << " ";
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



