/**

    Carleton University RTA Project - ain.cpp

    Noah Rubin
    
    Under supervision of Dr. Kevin Cheung and Dr. Brett Stevens.

    Let L1, L2 and L3 be Latin Squares of order n (as defined in OLS_HEADER.h).

    This program fills in n / 2 transversals of L3 (with Gurobi) and then attempts to determine whether L1 and L2 can
    be filled in to be orthogonal to each other and to the partially filled L3 (also w/ Gurobi). If the system is feasable
    then we may proceed to filling in the remainder of L3 with the CP solver (or-tools).

**/

// Relevant includes
#include "gurobi_c++.h"
#include <iostream>
#include <ctime>
#include <fstream>

// Use std
using namespace std;

// Define n and t

// Begin main
int main(int argc, char* argv[]) {

    int n = 0;
    if (argc == 2) {
        n = atoi(argv[1]);
    }
    else {
        n = 8;
    }

    std::cout << "n: " << n << "\n";

    // Try block is required for running Gurobi programs
    try {

        // Setup environment and begin logging
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "LatinSquare.log");
        env.start();
        GRBModel model = GRBModel(env);

        // Declare counters
        int i = 0, j = 0, k = 0, l = 0;

        GRBVar* x = new GRBVar[(int)pow(n, 4)];

        int nCube = (int)pow(n, 3);
        int nSqr = (int)pow(n, 2);

        // Instantiate variables for L1 and L2
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                        x[i * nCube + j * nSqr + k * n + l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
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
                        expr += x[i * nCube + j * nSqr + k * n + l];
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
                        expr += x[i * nCube + j * nSqr + k * n + l];
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
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "c_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (4)
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (l = 0; l < n; l++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "d_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (5)
        for (j = 0; j < n; j++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (k = 0; k < n; k++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "e_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (6)
        for (k = 0; k < n; k++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "f_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Fix first row of both squares
        
        for (i = 0; i < n; i++) {
            GRBLinExpr expr = 0;
            expr += x[i * nSqr + i * n + i];
            model.addConstr(expr == 1.0);
        }
        
        // Implement symmetry breaking - First column of each square
        for (i = 1; i < n; i++) {
            for (j = 1; j < n; j++) {
                for (k = 0; k < n; k++) {
                    //cant be i or > i + 1
                    if (j == i || j > i + 1 || k != i) {
                        GRBLinExpr expr = 0.0;
                        expr += x[i * nCube + k * n + j];
                        model.addConstr(expr == 0.0);
                    }
                }
            }
        }
        // Set model params
        model.set(GRB_IntParam_Threads, 1);

        // Write model to a .lp file so we can test w/ CPLEX
        model.write("C:\\Users\\Noah\\MOLS_n=" + to_string(n) + ".lp");

        std::srand(std::time(nullptr));
        int seed = abs(std::rand() * std::rand());
        model.set(GRB_IntParam_Seed, seed);

        std::cout << "seed: " << seed << endl;
        // Execute model
        model.optimize();
        
        // Print L1
        std::cout << endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
                            std::cout << to_string(k) << " ";
                    }
                }
            }
            std::cout << endl;
        }

        // Print L2
        std::cout << endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
                            std::cout << to_string(l) << " ";
                    }
                }
            }
            std::cout << endl;
        }
    }

    // Handle exceptions
    catch (GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << endl;
        std::cout << e.getMessage() << endl;
    }
    catch (...) {
        std::cout << "Exception during optimization" << endl;
    }
    return 0;
}
