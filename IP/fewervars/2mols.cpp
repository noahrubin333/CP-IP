/**
    The following implements model (4.9) from the paper "Modelling for feasibility - The case of mutually orthogonal Latin squares problem"
    except modified to use O(N^3) variables instead of O(N^4) variables

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
        GRBVar x[N][N][N];
        GRBVar y[N][N][N];
#ifndef NONEWVARS
        GRBVar z[N][N][N];
#ifdef EXTRA_VARS
        GRBVar w[N][N][N];
#endif
#endif

        int i = 0, j = 0, k = 0, l = 0;

        // Instantiate variables
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    x[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);

                    name = "y_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    y[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
#ifndef NONEWVARS
                    name = "z_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    z[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
#ifdef EXTRA_VARS
                    name = "w_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k);
                    w[i][j][k] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
#endif
#endif
                }
            }
        }

        // Implement constraint on values
        
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                GRBLinExpr expr = 0;
                for (k = 0; k < N; k++) {
                    expr += x[i][j][k];
                }
                string s = "1a_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);

                expr = 0;
                for (k = 0; k < N; k++) {
                    expr += y[i][j][k];
                }
                s = "1b_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }


        // Implement constraint on rows

        for (i = 0; i < N; i++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    expr += x[i][j][k];
                }
                string s = "2a_" + to_string(i) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);

                expr = 0;
                for (j = 0; j < N; j++) {
                    expr += y[i][j][k];
                }
                s = "2b_" + to_string(i) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint on columns

        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    expr += x[i][j][k];
                }
                string s = "3a_" + to_string(i) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);

                expr = 0;
                for (i = 0; i < N; i++) {
                    expr += y[i][j][k];
                }
                s = "3b_" + to_string(i) + "_" + to_string(k);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint on z and w variables

#ifndef NONEWVARS
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        {
                            GRBLinExpr expr = -1;
                            expr += x[i][j][k];
                            expr += y[i][j][l];
                            string s = "4a_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                            model.addConstr(expr, GRB_LESS_EQUAL, z[i][k][l], s);
                        }
#ifdef EXTRA_CONSTRAINTS
                        {
                            GRBLinExpr expr = -1;
                            expr += x[i][j][k];
                            expr += z[i][k][l];
                            string s = "4b_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                            model.addConstr(expr, GRB_LESS_EQUAL, y[i][j][l], s);
                        }
                        {
                            GRBLinExpr expr = -1;
                            expr += z[i][k][l];
                            expr += y[i][j][l];
                            string s = "4c_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                            model.addConstr(expr, GRB_LESS_EQUAL, x[i][j][k], s);
                        }
#endif
#ifdef EXTRA_VARS
                        {
                            GRBLinExpr expr = -1;
                            expr += x[i][j][k];
                            expr += y[i][j][l];
                            string s = "4d_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                            model.addConstr(expr, GRB_LESS_EQUAL, w[j][k][l], s);
                        }
#endif
                    }
                }
            }
        }

        // Implement orthogonality constraint

        for (k = 0; k < N; k++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    expr += z[i][k][l];
                }
                string s = "5a_" + to_string(k) + "_" + to_string(l);
                GRBConstr const5a = model.addConstr(expr, GRB_EQUAL, 1.0, s);

#ifdef EXTRA_ORTHOGONALITY
                expr = 0;
                for (i = 0; i < N; i++) {
                    expr += z[k][i][l];
                }
                s = "5b_" + to_string(k) + "_" + to_string(l);
                GRBConstr const5b = model.addConstr(expr, GRB_EQUAL, 1.0, s);

                expr = 0;
                for (i = 0; i < N; i++) {
                    expr += z[k][l][i];
                }
                s = "5c_" + to_string(k) + "_" + to_string(l);
                GRBConstr const5c = model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif

#ifdef LAZY_CONSTRAINTS
                const5a.set(GRB_IntAttr_Lazy, 1);
#ifdef EXTRA_ORTHOGONALITY
                const5b.set(GRB_IntAttr_Lazy, 1);
                const5c.set(GRB_IntAttr_Lazy, 1);
#endif
#endif

#ifdef EXTRA_VARS
                expr = 0;
                for (i = 0; i < N; i++) {
                    expr += w[i][k][l];
                }
                s = "5d_" + to_string(k) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif
            }
        }
#else
        int x1 = 0, x2 = 0, y1 = 0, y2 = 0;

        for (k = 0; k < N; k++)
            for (l = 0; l < N; l++)
                for (x1 = 0; x1 < N; x1++)
                    for (x2 = 0; x2 < N; x2++)
                        for (y1 = 0; y1 < N; y1++)
                            for (y2 = 0; y2 < N; y2++)
                            {
                                if(x1 != x2 || y1 != y2)
                                {
                                    GRBLinExpr expr = 0;
                                    expr += x[x1][y1][k];
                                    expr += y[x1][y1][l];
                                    expr += x[x2][y2][k];
                                    expr += y[x2][y2][l];
                                    GRBConstr const5 = model.addConstr(expr, GRB_LESS_EQUAL, 3.0);
                                    const5.set(GRB_IntAttr_Lazy, 1);
                                }
                            }
#endif

        // Implement Symmetry Breaking

        // Fix first row of both squares
        for (i = 0; i < N; i++) {
            GRBLinExpr expr = 0;
            expr += x[0][i][i];
            model.addConstr(expr == 1.0);

            expr = 0;
            expr += y[0][i][i];
            model.addConstr(expr == 1.0);
        }

#ifdef FIXCOL
        cout << "Fixing first column of second square to be (0, 2, 3, 4, ..., n-1, 1)" << endl;

        // Fix first column of both squares
        for (i = 1; i < N; i++) {
            GRBLinExpr expr = 0;
            expr += x[i][0][i];
            model.addConstr(expr == 1.0);

            if(i < N-1) {
                GRBLinExpr expr = 0;
                expr += y[i][0][i+1];
                model.addConstr(expr == 1.0);
            }
        }
        GRBLinExpr expr = 0;
        expr += y[N-1][0][1];
        model.addConstr(expr == 1.0);
#else
        // Fix first column of first square
        for (i = 1; i < N; i++) {
            GRBLinExpr expr = 0;
            expr += x[i][0][i];
            model.addConstr(expr == 1.0);

        }

        // Eliminate possibilities for first column of second square
        for (i = 1; i < N; i++) {
            for (j = i; j < N; j++) {
                if (j != i+1) {
                    GRBLinExpr expr = 0;
                    expr += y[i][0][j];
                    model.addConstr(expr == 0.0);
                }
            }
        }

        // Fix entries of second square
        {
            GRBLinExpr expr = 0;
            expr += y[1][0][2];
            model.addConstr(expr == 1.0);
        }

        {
            GRBLinExpr expr = 0;
            expr += y[N-2][0][N-1];
            model.addConstr(expr == 1.0);
        }
#endif

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
                    if (x[i][j][k].get(GRB_DoubleAttr_X) > 0.5)
                        cout << to_string(k) << " ";
                }
            }
            cout << endl;
        }

        cout << endl;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    if (y[i][j][k].get(GRB_DoubleAttr_X) > 0.5)
                        cout << to_string(k) << " ";
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



