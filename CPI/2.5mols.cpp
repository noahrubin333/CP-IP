/**
	Noah Rubin
	2020-09-15
	This program uses constraint programming to generate partially filled (upto t transversals) 1MOLS(n).
	It then attempts to solve the 2.5MOLS(n) system sequentially until a feasible solution is found.

**/

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include "gurobi_c++.h"

#include <vector>
#include <tuple>
#include <string>

using namespace std;
using namespace operations_research;
using namespace sat;

#define n ORDER
#define t ( n )

#define PRINT 0
#define ENEUMERATE_ALL 1
#define VERBOSE 0
#define NO_CUTS 0
#define RANDOMIZE_SEARCH 0
#define TIMEOUT 60000		// Timeout in seconds

chrono::time_point<chrono::high_resolution_clock> tic, toc;

int main(int argc, char* argv[]) {
	// Model
	CpModelBuilder cp_model;

	Domain domain(0, n - 1);
	Domain reduced_domain(1, n - 1);

	// Declare variables T_ij for CP Model L3
	IntVar T[n][t];

	//Counters
	int i = 0, j = 0, k = 0, l = 0;

	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < t; j++) {
			// x_0k for k = 0,...,t
			if (i == 0) {
				T[i][j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("x_" + to_string(i) + to_string(j));
			}
#ifdef FIXCOL
			else if (i == j || (i+1 == j) || (i==n-1 && j==1)) {
#else
			else if (i == j || (i==1 && j==2) || (i==n-2 && j==n-1)) {
#endif
				T[i][j] = cp_model.NewIntVar(reduced_domain).WithName("x_" + to_string(i) + to_string(j));
			}
			// General case
			else {
				T[i][j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
			}
		}
	}

	// Column uniqueness constraints
	for (j = 0; j < t; j++) {
		vector<IntVar> col_i;
		for (i = 0; i < n; i++) {
			col_i.push_back(T[i][j]);
		}
		cp_model.AddAllDifferent(col_i);
	}

	// Row uniqueness
	for (i = 0; i < n; i++) {
		vector<IntVar> entry_ij;
		for (j = 0; j < t; j++) {
			entry_ij.push_back(T[i][j]);
		}
		cp_model.AddAllDifferent(entry_ij);
	}

	// =================================================================
	// IMPLEMENT GUROBI MODEL HERE
	// =================================================================

	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "LatinSquare.log");
	env.start();

	// Build models for fixing transversals and solving L1, L2, L3 system
	GRBModel model = GRBModel(env);

	// Declare variables for MIP model and fixing L3 transversals
	GRBVar x[n][n][n][n];

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

	// Implement symmetry breaking
	
	// Fix first row of both squares
	for (i = 0; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += x[0][i][i][i];
		model.addConstr(expr == 1.0);
	}

#ifdef FIXCOL
	// Fix first column of L1 and L2
	for (i = 1; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += x[i][0][i][i==n-1 ? 1 : i+1];
		model.addConstr(expr == 1.0);
	}
#else
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

	model.set(GRB_IntParam_LogToConsole, VERBOSE > 1);
	model.set(GRB_IntParam_Threads, 1);
#if NO_CUTS
	model.set(GRB_IntParam_Cuts, 0);
#endif
	model.update();
	model.set(GRB_DoubleParam_NodeLimit, 1);
	model.optimize();
	model.set(GRB_DoubleParam_NodeLimit, INFINITY);

	// =========================================================
	// Tell model how to deal with solutions
	// =========================================================

	Model CPmodel;
	int num_solutions = 0;
	CPmodel.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {

		
		int XX[n][n];
		for(int i=0; i<n; i++)
			for(int j=0; j<n; j++)
				XX[i][j] = -1;

		for (i = 0; i < n; i++) {
			for (j = 0; j < t; j++) {
				XX[i][SolutionIntegerValue(r, T[i][j])] = j;
			}
		}

#if VERBOSE > 1
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
			    if(XX[i][j]==-1)
					cout << '.';
				else
					cout << XX[i][j];
			}
			cout << endl;
		}
		cout << endl;
#endif

		vector<vector<tuple<int, int>>> Tvec;
		for (j = 0; j < t; j++) {
			vector<tuple<int, int>> T_j;
			//cout << "T_" + to_string(j) + " = ";
			for (i = 0; i < n; i++) {
				T_j.push_back(make_tuple(i, SolutionIntegerValue(r, T[i][j])));
			  	//cout << "(" + to_string(i) + ", " + to_string(SolutionIntegerValue(r, T[i][j])) + ") ";
			}
			//cout << endl;
			Tvec.push_back(T_j);
		}
		try {
			
			GRBModel newmodel = GRBModel(model);

			// Implement constraint (13) for k in K
			for (k = 0; k < n; k++) {
				for (int m = 0; m < t; m++) {
					GRBLinExpr expr = 0;
					for (auto Tm0 : Tvec[m]) {
						for (l = 0; l < n; l++) {
							expr += x[std::get<0>(Tm0)][std::get<1>(Tm0)][k][l];
						}
					}
					string s = "C13_k_" + to_string(k) + "_" + to_string(m);
					newmodel.addConstr(expr == 1.0, "c_13a_" + to_string(k) + "_" + to_string(m));
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
					string s = "C13_l_" + to_string(l) + "_" + to_string(m);
					newmodel.addConstr(expr == 1.0, "c_13b_" + to_string(l) + "_" + to_string(m));
				}
			}

			newmodel.optimize();
			num_solutions++;
			toc = chrono::high_resolution_clock::now();

			if (newmodel.get(GRB_IntAttr_Status) == 2) {
				cout << endl;
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						for (k = 0; k < n; k++) {
							for (l = 0; l < n; l++) {
								if (newmodel.getVarByName("x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l)).get(GRB_DoubleAttr_X) > 0.5)
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
								if (newmodel.getVarByName("x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l)).get(GRB_DoubleAttr_X) > 0.5)
									cout << to_string(l) << " ";
							}
						}
					}
					cout << endl;
				}
				
				// Print L3
				cout << endl;
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						cout << to_string(XX[i][j]) << " ";
					}
					cout << endl;
				}

				cout << "Executed " + to_string(num_solutions) + " trials.\n";
				exit(0);
			}
			else {
				if (VERBOSE) {
					cout << "Trial #" << to_string(num_solutions) << " Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl << flush;
				}
			}
			if (chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() > TIMEOUT*1000) {
				cout << "Stopping after trial # " + to_string(num_solutions) + "\n";
				cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl;
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

		if (PRINT) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < t; j++) {
					cout << "x_" + to_string(i) + to_string(j) << " = " << SolutionIntegerValue(r, T[i][j]) << " ";
				}
				cout << "\n";
			}
		}
		}));

	// Params

	SatParameters parameters;
	if (ENEUMERATE_ALL) {
		parameters.set_enumerate_all_solutions(true);
	}
	parameters.set_enumerate_all_solutions(true);
#if RANDOMIZE_SEARCH == 1
	parameters.set_randomize_search(true);
#endif
	CPmodel.Add(NewSatParameters(parameters));

	// Execute
	tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &CPmodel);
	toc = chrono::high_resolution_clock::now();

	// Report time
	cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms" << endl;

	// Fin
	return EXIT_SUCCESS;
}
