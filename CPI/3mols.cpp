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
#include <iomanip>

using namespace std;
using namespace operations_research;
using namespace sat;

chrono::time_point<chrono::high_resolution_clock> tic, toc, tic2, toc2;

#if defined(CYCLETYPE)
	#define NO_DUPLICATE_CYCLE_TYPES 1
#else
	#define NO_DUPLICATE_CYCLE_TYPES 0
#endif

#define TIMEOUT 60000
#define n ORDER
#define STOP_ON_SOL 1
#define VERBOSE 0
#if defined(NOFIXCOL) || (defined(FIXCOLTHIRD) && !defined(DOMAINRED))
	#define DOMAIN_REDUCTION 0
#else
	#define DOMAIN_REDUCTION 1
#endif
#define RANDOMIZE_SEARCH 1
#define LOGGING 0

#if NO_DUPLICATE_CYCLE_TYPES==1 && !defined(NOFIXCOL) && !defined(FIXCOL)
#include "cycle_types.h"
#endif

int main(int argc, char* argv[]) {
	cout << std::setprecision(1) << fixed;
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
#if defined(FIXCOLTHIRD)
		domain_temp.push_back(i);
#else
		for (int j = 0; j < n; j++) {
			if (j != i
#if DOMAIN_REDUCTION == 1 && !defined(FIXCOL)
						&& !(i == 1 && j == 2) && !(i==n-2 && j==n-1)
#endif
			) {
				domain_temp.push_back(j);
			}
		}
#endif
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
#if DOMAIN_REDUCTION == 1 && !defined(FIXCOL)
				T[i * n + j] = cp_model.NewIntVar(firstCol[i]).WithName("x_" + to_string(i) + to_string(j));
#elif defined(FIXCOL)
				Domain d = Domain::FromValues({ i }).Complement();
				Domain e = Domain::FromValues({ i == n-1 ? 1 : i+1 }).Complement();
				T[i * n + j] = cp_model.NewIntVar(domain.IntersectionWith(d).IntersectionWith(e)).WithName("x_" + to_string(i) + to_string(j));
#else
				T[i * n + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
#endif
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
#if LOGGING == 0
	env.set("LogFile", "TimeTesting_3MOLS_n=" + to_string(n) + ".log");
#endif
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
			string s = "c1_" + to_string(i) + "_" + to_string(j);
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
			string s = "c2_" + to_string(i) + "_" + to_string(k);
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
			string s = "c3_" + to_string(i) + "_" + to_string(l);
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
			string s = "c4_" + to_string(j) + "_" + to_string(k);
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
			string s = "c5_" + to_string(j) + "_" + to_string(l);
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
			string s = "c6_" + to_string(k) + "_" + to_string(l);
			ip_model.addConstr(expr, GRB_EQUAL, 1.0, s);
		}
	}

	// Implement symmetry breaking - First row of each square
	for (i = 0; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += x[i * nSqr + i * n + i];
		ip_model.addConstr(expr == 1.0);
	}

#if DOMAIN_REDUCTION == 1 && !defined(FIXCOL)
	// Implement symmetry breaking - First column
	for (i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			for (k = 0; k < n; k++) {
				//cant be i or > i + 1
#if defined(FIXCOLTHIRD)
				if (j == i || j > i + 1) {
#else
				if (j == i || j > i + 1 || k != i) {
#endif
					GRBLinExpr expr = 0.0;
					expr += x[i*nCube + 0*nSqr + k*n+ j];
					ip_model.addConstr(expr == 0.0);
				}
			}
		}
	}
#endif

#if defined(FIXCOL)
	// Fix first column of L1 and L2
	for (i = 1; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += x[i*n*n*n+0*n*n+i*n+(i==n-1 ? 1 : i+1)];
		ip_model.addConstr(expr == 1.0);
	}
#endif

#if NO_DUPLICATE_CYCLE_TYPES==1 && !defined(NOFIXCOL) && !defined(FIXCOL)
	if (n==6) {
		GRBLinExpr expr1 = 0.0;
		for(int i = 0; i < n; i++) {
			ip_model.addConstr(x[5*n*n*n + 0*n*n + i*n + 4] == 0);
			expr1 += x[3*n*n*n + 0*n*n + i*n + 4];
		}
		ip_model.addConstr(expr1 == 1.0);
	} else if (n==7) {
		for(int i = 0; i < n; i++) {
			ip_model.addConstr(x[4*n*n*n + 0*n*n + i*n + 1] == 0);
		}
	} else if (n==8) {
		GRBLinExpr expr1 = 0.0;
		for(int i = 0; i < n; i++) {
			ip_model.addConstr(x[4*n*n*n + 0*n*n + i*n + 1] == 0);
			expr1 += x[5*n*n*n + 0*n*n + i*n + 6];
		}
		ip_model.addConstr(expr1 == 1.0);
	} else if (n==9) {
		GRBLinExpr expr1 = 0.0;
		GRBLinExpr expr2 = 0.0;
		for(int i = 0; i < n; i++) {
			ip_model.addConstr(x[6*n*n*n + 0*n*n + i*n + 3] == 0);
			ip_model.addConstr(x[6*n*n*n + 0*n*n + i*n + 1] == 0);
			ip_model.addConstr(x[6*n*n*n + 0*n*n + i*n + 4] == 0);
			expr1 += x[8*n*n*n + 0*n*n + i*n + 6];
			expr1 += x[2*n*n*n + 0*n*n + i*n + 3];
			expr2 += x[8*n*n*n + 0*n*n + i*n + 7];
			expr2 += x[2*n*n*n + 0*n*n + i*n + 3];
		}
		ip_model.addConstr(expr1 <= 1.0);
		ip_model.addConstr(expr2 <= 1.0);
	} else if (n==10) {
		GRBLinExpr expr1 = 0.0;
		GRBLinExpr expr2 = 0.0;
		for(int i = 0; i < n; i++) {
			ip_model.addConstr(x[5*n*n*n + 0*n*n + i*n + 1] == 0);
			ip_model.addConstr(x[5*n*n*n + 0*n*n + i*n + 4] == 0);
			ip_model.addConstr(x[6*n*n*n + 0*n*n + i*n + 1] == 0);
			ip_model.addConstr(x[6*n*n*n + 0*n*n + i*n + 3] == 0);
			expr1 += x[9*n*n*n + 0*n*n + i*n + 7];
			expr1 += x[4*n*n*n + 0*n*n + i*n + 1];
			expr2 += x[7*n*n*n + 0*n*n + i*n + 8];
		}
		ip_model.addConstr(expr1 <= 1.0);
		ip_model.addConstr(expr2 == 1.0);
	}
	else{
		for(int i=0; i<blocked_cycle_types.size(); i++) {
			GRBLinExpr expr = 0.0;
			for(int j=2; j<n; j++) {
				if(j != n-2) {
					int k = blocked_cycle_types[i][j];
					for (int l=0; l<n; l++)
						expr += x[j*nCube + 0*nSqr + l*n + k];
				}
			}
			ip_model.addConstr(expr <= n-5);
		}
	}
#endif

	// Force constraints to update before we clone the model
	ip_model.update();
	ip_model.set(GRB_IntParam_LogToConsole, 0);
	ip_model.set(GRB_IntParam_Threads, 1);
	ip_model.set(GRB_DoubleParam_TimeLimit, TIMEOUT);
	//ip_model.set(GRB_IntParam_Presolve, 0);
	//ip_model.set(GRB_IntParam_Cuts, 0);
	//ip_model.set(GRB_IntParam_CliqueCuts, 0);
	//ip_model.set(GRB_IntParam_ZeroHalfCuts, 1);
	//ip_model.set(GRB_DoubleParam_Heuristics, 0);

	Model model;

	// Sets a time limit of 10 seconds.
	SatParameters parameters;
	parameters.set_max_time_in_seconds(TIMEOUT);
	parameters.set_enumerate_all_solutions(true);
#if RANDOMIZE_SEARCH == 1
	parameters.set_randomize_search(true);
#endif
	model.Add(NewSatParameters(parameters));
	
	int num_solutions = 0;

	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {

		num_solutions++;

#if VERBOSE > 1
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				printf("%ld ", SolutionIntegerValue(r, T[i*n+j]));
			}
			cout << endl;
		}
#endif

		// Vector which represents indexes of filled in values in L3
		vector<vector<tuple<int, int>>> Tvec;
		// Grab L3 info from the CP Solver and encode into vector
		for (int z = 0; z < n; z++) {
			vector<tuple<int, int>> T_z;
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					if (SolutionIntegerValue(r, T[i * n + j]) == z)
					{
						T_z.push_back(make_tuple(i, j));
					}
				}
			}
			Tvec.push_back(T_z);
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
				ip_model.addConstr(expr == 1.0, "c_13_k_" + to_string(k) + "_" + to_string(m));
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
				ip_model.addConstr(expr == 1.0, "c_13_l_" + to_string(l) + "_" + to_string(m));

			}
		}

		tic2 = chrono::high_resolution_clock::now();
		ip_model.optimize();
		toc2 = chrono::high_resolution_clock::now();

		// Print
		if (ip_model.get(GRB_IntAttr_Status) == GRB_OPTIMAL ) {
			cout << "Found solution after " << num_solutions << " squares tested\n\n";

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

			cout << "\n";

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

			cout << "\n";

			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					cout << SolutionIntegerValue(r, T[i * n + j]) << " ";
				}
				cout << "\n";
			}

			cout << "\n";

			cout << "Total time:" << chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() << "ms" << endl;

			if (STOP_ON_SOL) {
				exit(0);
			}
		} else if(VERBOSE > 1) {
			cout << "No solutions" << endl;
		}

		if (ip_model.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT) {
			cout << "IP model timed out after evaluating " << num_solutions << " squares" << endl;
		}


		// Remove variable fixing from this iteration
		for (int k = 0; k < n; k++) {
			for (int m = 0; m < n; m++) {
				GRBConstr ck = ip_model.getConstrByName("c_13_k_" + to_string(k) + "_" + to_string(m));
				GRBConstr cl = ip_model.getConstrByName("c_13_l_" + to_string(k) + "_" + to_string(m));
				ip_model.remove(ck);
				ip_model.remove(cl);
			}
		}
		ip_model.update();

		// Tally total iterations so far
		if (VERBOSE || !(num_solutions % 1000))
		{
			toc = chrono::high_resolution_clock::now();
			cout << "Trial # " << num_solutions << " Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc2 - tic2).count() << " ms" << " Total time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms Avg. time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()/((double)1000*num_solutions) << " s" << endl << flush;
		}

		if (chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() >= TIMEOUT*1000) {
			cout << "Timed out after " << chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() << " ms and evaluated " << num_solutions << " squares" << endl;
			exit(0);
		}
	}));

	tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
	cout << "Timed out and evaluated " << num_solutions << " squares" << endl;
#if VERBOSE > 0
	cout << CpSolverResponseStats(response);
#endif
    return 0;
}
