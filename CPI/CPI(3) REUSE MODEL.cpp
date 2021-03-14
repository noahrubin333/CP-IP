/**
	Noah Rubin
	2020-09-27
	This program uses constraint programming to generate partially filled (upto t transversals) 1MOLS(n).
	It then attempts to solve the 2.5MOLS(n) system as follows:

	Generate k Gurobi Models (corresponding to k possible fillings for column 1 of L2)

	WHILE TRUE DO:
		 i)	Generate next filling of L3
		ii) Create constraints corresponding to filling of L3
	   iii) Encode constraints into all k models
		iv) Optimize each model in batch of k threads
		 v) FOR EACH MODEL DO: 
				IF MODEL->STATUS == OPTIMAL DO:
					EXIT(0)

**/

// Memcheck includes
#define _CRTDBG_MAP_ALLOC
#include <stdlib.h>
#include <crtdbg.h>

// OR-TOOLS and Gurobi includes
#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "gurobi_c++.h"

// Windows API
#include <Windows.h>

// Required elements from std
#include <vector>
#include <tuple>
#include <string>

// Namespaces
using namespace std;
using namespace operations_research;
using namespace sat;


// Define PRINT and ENEUMERATE_ALL flags
#define PRINT 0
#define ENEUMERATE_ALL 1

// Enter main
int main(int argc, char* argv[]) {

	// Set memcheck flags to auto
	_CrtSetDbgFlag(_CRTDBG_ALLOC_MEM_DF | _CRTDBG_LEAK_CHECK_DF);

	int n = 0, t = 0;
	if (argc == 3) {
		n = atoi(argv[1]);
		t = atoi(argv[2]);
	}
	else {
		n = 5;
		t = n / 2;
	}

	n = 7;
	t = 3;

	chrono::steady_clock::time_point tic, toc;

	// =================================================================
	// BUILD CP MODEL
	// =================================================================

	// CP Model Builder
	CpModelBuilder cp_model;

	// CP Domain (general case)
	Domain domain(0, n - 1);

	// Declare variables T_ij for CP Model L3
	IntVar * T = new IntVar [n * t];

	//Counters
	int i = 0, j = 0, k = 0, l = 0;

	Domain* restricted_domains_first_col_L3 = new Domain[n];

	for (i = 0; i < n; i++) {
		vector<int64> Domain_Temp;
		for (j = 0; j < n; j++) {
			if (j != i) {
				Domain_Temp.push_back(j);
			}
		}
		restricted_domains_first_col_L3[i] = Domain::FromValues(Domain_Temp);
	}

	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < t; j++) {
			// First row
			if (i == 0) {
				T[i * t + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("x_" + to_string(i) + to_string(j));
			}
			// General case
			else {
				T[i * t + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
			}
		}
	}

	// Column uniqueness constraints
	for (j = 0; j < t; j++) {
		vector<IntVar> col_i;
		for (i = 0; i < n; i++) {
			col_i.push_back(T[i * t + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}

	// Row uniqueness constraints
	for (i = 0; i < n; i++) {
		vector<IntVar> entry_ij;
		for (j = 0; j < t; j++) {
			entry_ij.push_back(T[i * t + j]);
		}
		cp_model.AddAllDifferent(entry_ij);
	}

	// =================================================================
	// BUILD GUROBI MODEL
	//
	//	-> We can build a single model for all generic constraints, then
	//	-> copy this model for each possible column fixing in L2
	// =================================================================

	// Setup environment
	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "LatinSquare.log");
	env.start();

	// Declare base model
	GRBModel model = GRBModel(env);

	// Declare variables for MIP model
	// Need to make dynamic as follows:

	// x = new GRBVar *[n ^ 4]
	// x[i][j][k][l] -> x[i * n^3 + j * n^2 + k * n + l]
	GRBVar *x = new GRBVar[(int)pow(n, 4)];

	int nCube = (int)pow(n, 3);
	int nSqr  = (int)pow(n, 2);

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
					expr += x[i * nCube + j * nSqr + k * n + l];
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
					expr += x[i * nCube + j * nSqr + k * n + l];
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
					expr += x[i * nCube + j * nSqr + k * n + l];
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
					expr += x[i * nCube + j * nSqr + k * n + l];
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
					expr += x[i * nCube + j * nSqr + k * n + l];
				}
			}
			string s = "6_" + to_string(i) + "_" + to_string(j);
			model.addConstr(expr, GRB_EQUAL, 1.0, s);
		}
	}
	/**
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
	**/
	// Implement symmetry breaking - First row of each square
	
	for (i = 0; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += x[i * nSqr + i * n + i];
		model.addConstr(expr == 1.0);
	}
	/*
	for (i = 0; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += x[i * nCube + i * nSqr + i * n];
		model.addConstr(expr == 1.0);
	}
	*/
	// Force constraints to update before we clone the model
	model.update();
	model.set(GRB_IntParam_LogToConsole, 0);
	model.set(GRB_IntParam_Threads, 1);
	

	// =========================================================
	// Begin loop
	// =========================================================

	// CP Model - to be built
	Model CPmodel;

	// Keep track of trial count
	int num_solutions = 0;


	model.optimize();
	

	// Implement behavior for each fixing of L3 transversals
	CPmodel.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {

		GRBModel *BASE_MODEL = new GRBModel(model);

		// Vector which represents indexes of filled in values in L3
		vector<vector<tuple<int, int>>> Tvec;

		// Grab L3 info from the CP Solver and encode into vector
		for (j = 0; j < t; j++) {
			vector<tuple<int, int>> T_j;
			for (i = 0; i < n; i++) {
				T_j.push_back(make_tuple(i, SolutionIntegerValue(r, T[i * t + j])));
			}
			Tvec.push_back(T_j);
		}

		// Begin MIP solve
		try {
			// Implement constraint (13) for k in K
			for (k = 0; k < n; k++) {
				for (int m = 0; m < t; m++) {
					GRBLinExpr expr = 0;
					for (auto Tm0 : Tvec[m]) {
						for (l = 0; l < n; l++) {
							expr += x[std::get<0>(Tm0) * nCube + std::get<1>(Tm0) * nSqr + k * n + l];
						}
					}

					BASE_MODEL->addConstr(expr == 1.0, "c_13_k");
				}
			}

			// Implement constraint (13) for l in L
			for (l = 0; l < n; l++) {
				for (int m = 0; m < t; m++) {
					GRBLinExpr expr = 0;
					for (auto Tm0 : Tvec[m]) {
						for (k = 0; k < n; k++) {
							expr += x[std::get<0>(Tm0) * nCube + std::get<1>(Tm0) * nSqr + k * n + l];
						}
					}

					BASE_MODEL->addConstr(expr == 1.0, "c_13_k");

				}
			}

			BASE_MODEL->optimize();
			
			if (!(num_solutions % 1000)) {
				cout << "Trial #" + to_string(num_solutions) + "\n";
			}
			// If clone returned feasible then print solution and exit
			if (BASE_MODEL->get(GRB_IntAttr_Status) == 2) {
				toc = chrono::high_resolution_clock::now();
				
				cout << endl;
				// Print L1
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						for (k = 0; k < n; k++) {
							for (l = 0; l < n; l++) {
								if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
										cout << to_string(k) << " ";
							}
						}
					}
					cout << endl;
				}
				// Print L2
				cout << "\n===============" << endl;
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						for (k = 0; k < n; k++) {
							for (l = 0; l < n; l++) {
								if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
										cout << to_string(l) << " ";
							}
						}
					}
					cout << endl;
				}
				cout << "\n===============" << endl;
				for (i = 0; i < n; i++) {
					for (j = 0; j < t; j++) {
						cout << "x_" + to_string(i) + to_string(j) << " = " << SolutionIntegerValue(r, T[i * t + j]) << " ";
					}
					cout << "\n";
				}
				cout << "\n===============" << endl;

				cout << "Executed " + to_string(num_solutions) + " trials in " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms\n";
				exit(0);
				
			}
		}
		catch (GRBException e) {
			cout << "Error code = " << e.getErrorCode() << endl;
			cout << e.getMessage() << endl;
		}
		catch (...) {
			cout << "Exception during optimization" << endl;
		}

		delete BASE_MODEL;
		num_solutions++;
		}));

	// Params for CP Model
	SatParameters parameters;
	if (ENEUMERATE_ALL) {
		parameters.set_enumerate_all_solutions(true);
	}
	parameters.set_enumerate_all_solutions(true);
	CPmodel.Add(NewSatParameters(parameters));

	// Start time
	tic = chrono::high_resolution_clock::now();

	// Build and execute CP Model
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &CPmodel);

	// End time, should never be reached
	toc = chrono::high_resolution_clock::now();

	// Report time
	cout << "Time taken to exhaust solution space: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms";

	// Fin
	return EXIT_SUCCESS;
}