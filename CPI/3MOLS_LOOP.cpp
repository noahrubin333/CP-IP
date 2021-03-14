/**
    Noah Rubin
    2020-06-24
    This program uses OR-TOOLS and constraint programming to generate pairs of OLS(n).
**/


//FIX FIRST ROW AND COLUMN
#include "gurobi_c++.h"
#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"

#include <vector>
#include <tuple>
#include <string>
#include <iostream>

using namespace std;
using namespace operations_research;
using namespace sat;

chrono::steady_clock::time_point tic, toc;

int main(int argc, char* argv[]) {
	int n;
	if (argc == 2) {
		n = atoi(argv[1]);
	}
	else {
		cout << "Enter n\n";
		exit(0);
	}


	// =================================================================
	// BUILD CP MODEL
	// =================================================================
	// CP Model Builder
	CpModelBuilder cp_model;
	// CP Domain (general case)
	Domain domain(0, n - 1);
	Domain *firstCol = new Domain[n];
	for (int i = 0; i < n; i++) {
		vector<int64> domain_temp;
		for (int j = 0; j < n; j++) {
			if (j != i) {
				domain_temp.push_back(j);
			}
		}
		firstCol[i] = Domain::FromValues(domain_temp);
	}


	// Declare variables T_ij for CP Model L3
	IntVar* T = new IntVar[n * n];
	//Counters
	int i = 0, j = 0, k = 0, l = 0;
	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// First row
			if (i == 0) {
				T[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("x_" + to_string(i) + to_string(j));
			}
			// First col
			else if (j == 0 && i > 0) {
				T[i * n + j] = cp_model.NewIntVar(firstCol[i]).WithName("x_" + to_string(i) + to_string(j));
			}
			// General case
			else {
				T[i * n + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
			}
		}
	}

	// Column uniqueness constraints
	for (j = 0; j < n; j++) {
		vector<IntVar> col_i;
		for (i = 0; i < n; i++) {
			col_i.push_back(T[i * n + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}

	// Row uniqueness constraints
	for (i = 0; i < n; i++) {
		vector<IntVar> entry_ij;
		for (j = 0; j < n; j++) {
			entry_ij.push_back(T[i * n + j]);
		}
		cp_model.AddAllDifferent(entry_ij);
	}

	// =================================================================
	// BUILD MIP MODEL
	// =================================================================
	// Setup environment
	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "TimeTesting_3MOLS_n=" + to_string(n) + ".log");
	env.start();
	// Declare base model
	GRBModel ip_model = GRBModel(env);
	// Declare variables for MIP model
	// Need to make dynamic as follows:
	// x = new GRBVar *[n ^ 4]
	// x[i][j][k][l] -> x[i * n^3 + j * n^2 + k * n + l]
	GRBVar* x = new GRBVar[(int)pow(n, 4)];
	int nCube = (int)pow(n, 3);
	int nSqr = (int)pow(n, 2);
	// Instantiate variables for L1 and L2
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				for (l = 0; l < n; l++) {
					string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
					x[i * nCube + j * nSqr + k * n + l] = ip_model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
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
			ip_model.addConstr(expr, GRB_EQUAL, 1.0, s);
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
			ip_model.addConstr(expr, GRB_EQUAL, 1.0, s);
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
			ip_model.addConstr(expr, GRB_EQUAL, 1.0, s);
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
			ip_model.addConstr(expr, GRB_EQUAL, 1.0, s);
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
			ip_model.addConstr(expr, GRB_EQUAL, 1.0, s);
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
			ip_model.addConstr(expr, GRB_EQUAL, 1.0, s);
		}
	}

	// Implement symmetry breaking - First row of each square
	for (i = 0; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += x[i * nSqr + i * n + i];
		ip_model.addConstr(expr == 1.0);
	}

	// Force constraints to update before we clone the model
	ip_model.update();
	ip_model.set(GRB_IntParam_LogToConsole, 0);
	ip_model.set(GRB_IntParam_Threads, 1);
	ip_model.set(GRB_DoubleParam_TimeLimit, 50000);

	Model model;

	// Sets a time limit of 10 seconds.
	SatParameters parameters;
	parameters.set_max_time_in_seconds(50000.0);
	parameters.set_enumerate_all_solutions(true);
	parameters.set_randomize_search(true);
	model.Add(NewSatParameters(parameters));
	
	int num_solutions = 0;

	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {

		num_solutions++;

		// Vector which represents indexes of filled in values in L3
		vector<vector<tuple<int, int>>> Tvec;
		// Grab L3 info from the CP Solver and encode into vector
		for (j = 0; j < n; j++) {
			vector<tuple<int, int>> T_j;
			for (i = 0; i < n; i++) {
				T_j.push_back(make_tuple(i, SolutionIntegerValue(r, T[i * n + j])));
			}
			Tvec.push_back(T_j);
		}

		// Implement constraint (13) for k in K
		for (k = 0; k < n; k++) {
			for (int m = 0; m < n; m++) {
				GRBLinExpr expr = 0;
				for (auto Tm0 : Tvec[m]) {
					for (l = 0; l < n; l++) {
						expr += x[std::get<0>(Tm0) * nCube + std::get<1>(Tm0) * nSqr + k * n + l];
					}
				}
				ip_model.addConstr(expr == 1.0, "c_13_k");
			}
		}

		// Implement constraint (13) for l in L
		for (l = 0; l < n; l++) {
			for (int m = 0; m < n; m++) {
				GRBLinExpr expr = 0;
				for (auto Tm0 : Tvec[m]) {
					for (k = 0; k < n; k++) {
						expr += x[std::get<0>(Tm0) * nCube + std::get<1>(Tm0) * nSqr + k * n + l];
					}
				}
				// Add constraint for L2 to each clone
				ip_model.addConstr(expr == 1.0, "c_13_l");

			}
		}

		ip_model.optimize();

		// Print
		if (ip_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL ) {
			cout << "Found solution after " << num_solutions << " squares tested\n\n";

			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					cout << SolutionIntegerValue(r, T[i * n + j]) << " ";
				}
				cout << "\n";
			}

			cout << "\n\n";

			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					for (k = 0; k < n; k++) {
						for (l = 0; l < n; l++) {
							if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5) {
								cout << k << " ";
							}
						}
					}
				}
				cout << "\n";
			}

			cout << "\n\n";

			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					for (k = 0; k < n; k++) {
						for (l = 0; l < n; l++) {
							if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5) {
								cout << l << " ";
							}
						}
					}
				}
				cout << "\n";
			}

			cout << "\nTotal time:" << chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() << "ms\n";

			exit(0);
		}

		if (ip_model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
			cout << "IP model timed out after evaluating " << num_solutions << "squares";
		}


		// Remove variable fixing from this iteration
		for (i = 0; i < nSqr; i++) {
			GRBConstr ck = ip_model.getConstrByName("c_13_k");
			GRBConstr cl = ip_model.getConstrByName("c_13_l");

			ip_model.remove(ck);
			ip_model.remove(cl);
			ip_model.update();

		}


		// Tally total iterations so far
		if (!(num_solutions % 1000)) {
			cout << "Evaluated " << num_solutions << " possible squares so far\n";
		}

		if (chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() >= 50000000) {
			cout << "Timed out and evaluated " << num_solutions << "squares";
		}
	}));

	tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
	cout << "Timed out and evaluated " << num_solutions << "squares";
    return 0;
}