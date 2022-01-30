/**

    Carleton University RTA Project - main.cpp

    Noah Rubin
    
    Under supervision of Dr. Kevin Cheung and Dr. Brett Stevens.

    Let L1, L2 and L3 be Latin Squares of order n (as defined in OLS_HEADER.h).

    This program fills in n / 2 transversals of L3 (with Gurobi) and then attempts to determine whether L1 and L2 can
    be filled in to be orthogonal to each other and to the partially filled L3 (also w/ Gurobi). If the system is feasable
    then we may proceed to filling in the remainder of L3 with the CP solver (or-tools).

**/

// Relevant includes
#include "gurobi_c++.h"
#define n N
#include <iostream>
#include <vector>
#include <tuple>
#include <ctime>

// Use std
using namespace std;

// Define t
#define t (n / 2)
#define x xvar

GRBVar x[n][n][n][n];

#ifdef CALLBACK_CUTS
#include "callback.cpp"
#endif

// Define function to fill in the first t entries in L3
void fix_transversals(GRBVar[n][n][n], GRBModel, vector<vector<tuple<int, int>>>&, vector<int>&, int, char**);
void fix_transversals(GRBVar x[n][n][n], GRBModel model, vector<vector<tuple<int,int>>>& Tvec, vector<int>& firstCol, int argc, char** argv) {
	try {

		// Declare Counters
		int i = 0, j = 0, v = 0;

		//Instantiate variables for L3
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				for (v = 0; v < n; v++) {
					string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(v);
					x[i][j][v] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
				}
			}
		}

		//Enforce cell uniqueness
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				GRBLinExpr expr = 0;
				for (v = 0; v < n; v++) {
					expr += x[i][j][v];
				}
				string s = "V_" + to_string(i) + "_" + to_string(j);
				model.addConstr(expr == 1.0);
			}
		}

		//Enforce row uniqueness upto t
		for (i = 0; i < n; i++) {
			for (v = 0; v < t; v++) {
				GRBLinExpr expr = 0;
				for (j = 0; j < n; j++) {
					expr += x[i][j][v];
				}
				string s = "R_" + to_string(i) + "_" + to_string(v);
				model.addConstr(expr == 1.0);
			}
		}

		//Enforce column uniqueness upto t
		for (j = 0; j < n; j++) {
			for (v = 0; v < t; v++) {
				GRBLinExpr expr = 0;
				for (i = 0; i < n; i++) {
					expr += x[i][j][v];
				}
				string s = "C_" + to_string(j) + "_" + to_string(v);
				model.addConstr(expr == 1.0);
			}
		}

		// Symmetry breaking
		for (i = 0; i < t; i++) {
			GRBLinExpr expr = 0;
			expr += x[0][i][i];
			model.addConstr(expr == 1.0);
		}


		// Need to take first col L2 into account

		for (i = 1; i < n; i++) {
			GRBLinExpr expr = 0;
			expr += x[i][0][i];
			model.addConstr(expr == 0.0);
		}

		for (i = 1; i < n; i++) {
			GRBLinExpr expr = 0;
			expr += x[i][0][firstCol[i]];
			model.addConstr(expr == 0.0);
		}

		//Solve

		srand(std::time(nullptr));
		model.set(GRB_IntParam_Seed, std::rand());
		model.set(GRB_IntParam_LogToConsole, 0);
		model.set(GRB_StringParam_LogFile, "test.log");
		if(argc >= 2) {
			int seed = stoi(argv[1]);
			model.set(GRB_IntParam_Seed, seed);
			cout << "Using random seed " << seed << endl;
		}
		model.optimize();

		// Print what is currently filled in L3

		
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				for (int k = 0; k < n; k++) {
					if (x[i][j][k].get(GRB_DoubleAttr_X) == 1.0) {
						if (k < t) {
							std::cout << to_string(k) + " ";
						}
						else {
							std::cout << ". ";
						}
					}
				}
			}
			std::cout << endl;
		}
		

		// Encode info into vector (called Tvec)
		for (int m = 0; m < t; m++) {
			// Fill in T_m0 for m0 = 0, ..., t - 1
			vector<std::tuple<int, int>> T_m;
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					// If x_ijk == 1 then add (i, j) to T_m0
					if (x[i][j][m].get(GRB_DoubleAttr_X) == 1) {
						T_m.push_back(make_tuple(i, j));
					}
				}
			}
			Tvec.push_back(T_m);
		}
	}

	// Handle exceptions
	catch (GRBException e) {
		cout << "Error code = " << e.getErrorCode() << endl;
		cout << e.getMessage() << endl;
	}
	catch (...) {
		cout << "Exception during optimization" << endl;
	}
}

// Begin main
int main(int argc, char* argv[]) {

    // Try block is required for running Gurobi programs
    try {

        // Setup environment and begin logging
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "LatinSquare.log");
        env.start();

        // Build models for fixing transversals and solving L1, L2, L3 system
        GRBModel model = GRBModel(env);
        GRBModel transversal_fixer_model = GRBModel(env);

        // Declare counters
        int i = 0, j = 0, k = 0, l = 0;

        // Declare variables for whole model and fixing L3 transversals
        GRBVar Ts[n][n][n];

        // Instantiate variables for L1 and L2
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
                string s = "2_" + to_string(i) + "_" + to_string(k);
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
                string s = "3_" + to_string(i) + "_" + to_string(l);
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
                string s = "4_" + to_string(j) + "_" + to_string(k);
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
                string s = "5_" + to_string(j) + "_" + to_string(l);
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
                string s = "6_" + to_string(k) + "_" + to_string(l);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Declare vector for getting L3 transversal info
        vector<vector<tuple<int, int>>> Tvec;

        vector<int> firstCol_L2 = {0};
        for(i = 2; i < n; i++)
            firstCol_L2.push_back(i);
        firstCol_L2.push_back(1);

        // Attempt to fix t transversals of L3
        fix_transversals(Ts, transversal_fixer_model, Tvec, firstCol_L2, argc, argv);

        // Implement constraint (13) for k in K
        for (k = 0; k < n; k++) {
            for (int m = 0; m < t; m++) {
                GRBLinExpr expr = 0;
                for (auto Tm0 : Tvec[m]) {
                    for (l = 0; l < n; l++) {
                        expr += x[std::get<0>(Tm0)][std::get<1>(Tm0)][k][l];
                    }
                }
                model.addConstr(expr == 1.0);
            }
        }

        // Implement constraint (13) for l in L
        for (l = 0; l < n; l++) {
            for (int m = 0; m < t; m++) {
                GRBLinExpr expr = 0;
                for (auto Tm0 : Tvec[m]) {
                    for (k = 0; k < n; k++) {
                        expr += x[std::get<0>(Tm0)][std::get<1>(Tm0)][k][l];
                    }
                }
                model.addConstr(expr == 1.0);
            }
        }

        // Last step is symmetry breaking

#ifdef FIXCOL
        cout << "Fixing first column of second square to be (0, 2, 3, 4, ..., n-1, 1)" << endl;

        // Fix first column of both squares
        for (i = 1; i < n-1; i++) {
            GRBLinExpr expr = 0;
            expr += x[i][0][i][i+1];
            model.addConstr(expr == 1.0);
        }

        // Fix values of entries at (N-1,0) to be (N-1,1)
        {
            GRBLinExpr expr = 0;
            expr += x[n-1][0][N-1][1];
            model.addConstr(expr == 1.0);
        }
#else
        // Fix first column of L1 and domains of first column of L2
        for (i = 1; i < n; i++) {
            for (j = 1; j < n; j++) {
                for (k = 0; k < n; k++) {
                    //cant be i or > i + 1
                    if (j == i || j > i + 1 || k != i) {
                        GRBLinExpr expr = 0.0;
                        expr += x[i][0][k][j];
                        model.addConstr(expr == 0.0);
                    }
                }
            }
        }
#endif

        // Fix first row of both squares
        for (i = 0; i < n; i++) {
            GRBLinExpr expr = 0;
            expr += x[0][i][i][i];
            model.addConstr(expr == 1.0);
        }

        for (i = 1; i < n; i++) {
            GRBLinExpr expr = 0;
            expr += x[i][0][i][firstCol_L2[i]];
            model.addConstr(expr == 1.0);
        }


        // Set model params
        model.set(GRB_IntParam_Threads, 1); // Enable multithreading accross 12 cores

        // Seed w/ random number
        //srand(std::time(nullptr));
        //model.set(GRB_IntParam_Seed, std::rand() * std::rand());

        // Write model to a .lp file so we can test w/ CPLEX
        //model.write("C:\\Users\\Noah\\MOLS_n=" + to_string(n) + ".lp");

#ifdef CALLBACK_CUTS
        cout << "Using custom clique cuts" << endl;
        model.setCallback(new MyCallback);
#endif

        // Execute model and catch exceptions
        try {
            model.optimize();
        }
        catch (GRBException e) {
            cout << "Error code = " << e.getErrorCode() << endl;
            cout << e.getMessage() << endl;
            srand(std::time(nullptr));
            model.set(GRB_IntParam_Seed, std::rand() * std::rand());
        }
        catch (...) {
            cout << "Exception during optimization" << endl;
        }
        
        if (model.get(GRB_IntAttr_Status) == 3) {
            exit(0);
        }

        // Print L1
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

        // Print L2
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
