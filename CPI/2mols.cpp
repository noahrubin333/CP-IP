/**
	Noah Rubin
	2020-09-02
	(modified by Curtis Bright)
	This program uses OR-Tools to fix the entries 0,1...,t-1 in square X and the corresponding t transversals of square Y.
	It then passes the model to Gurobi which attempts to fill in the remaining entries.
**/

#include "gurobi_c++.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"

using namespace std;
using namespace operations_research;
using namespace sat;

#define LOGGING 0
#define VERBOSE 1
#define SET_NODE_LIMIT 0	// Set to 1 to limit the # of nodes solved
#define NODE_LIMIT 10000
#define NO_CUTS 0

GRBVar* Gx;
int order;

#ifdef CALLBACK_CUTS
#include "callback.cpp"
#endif

int num_solutions_1, num_solutions_2;

void PRINT(GRBVar* x, int n, int t);

int main(int argc, char* argv[]) {

	int n = 9, t = 1; // Default

	if (argc != 3) {
		cout << "INVALID NUMBER OF ARGUMENTS, NEED n AND t\n";
		exit(0);
	}
	else {
		n = atoi(argv[1]);
		t = atoi(argv[2]);
	}

	order = n;
	int nCube = (int)pow(n, 3);
	int nSqr = (int)pow(n, 2);

	// CP Model Builder
	CpModelBuilder cp_model;

	// CP Domain (general case)
	Domain domain(0, n - 1);

	// Declare variables for the two models
	IntVar* x = new IntVar[n * n];
	IntVar* y = new IntVar[n * n];
	IntVar* z = new IntVar[n * t];

	//Counters
	int i = 0, j = 0, k = 0, l = 0;

	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// First row
			if (i == 0) {
				x[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("x_" + to_string(i) + to_string(j));
				y[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("y_" + to_string(i) + to_string(j));
				if (j < t)
					z[i * t + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("z_" + to_string(i) + to_string(j));
			}
			// General case
			else {
				x[i * n + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
				if(j==0)
					y[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ i })).WithName("y_" + to_string(i) + to_string(j));
				else
					y[i * n + j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
				if (j < t)
					z[i * t + j] = cp_model.NewIntVar(domain).WithName("z_" + to_string(i) + to_string(j));
			}
		}
	}

	//=====================================================================
	// Now build model for fixing transversals in X
	//====================================================================

	// Column uniqueness constraints
	for (j = 0; j < n; j++) {
		vector<IntVar> col_i;
		for (i = 0; i < n; i++) {
			col_i.push_back(x[i * n + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}

	for (j = 0; j < t; j++) {
		vector<IntVar> col_i;
		for (i = 0; i < n; i++) {
			col_i.push_back(z[i * t + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}

	for (j = 0; j < n; j++) {
		vector<IntVar> Ys;
		for (i = 0; i < n; i++) {
			Ys.push_back(y[i * n + j]);
		}
		cp_model.AddAllDifferent(Ys);
	}

	// Row uniqueness constraints
	for (i = 0; i < n; i++) {
		vector<IntVar> entry_ij;
		for (j = 0; j < n; j++) {
			entry_ij.push_back(x[i * n + j]);
		}
		cp_model.AddAllDifferent(entry_ij);
	}

	for (i = 0; i < n; i++) {
		vector<IntVar> entry_ij;
		for (j = 0; j < t; j++) {
			entry_ij.push_back(z[i * t + j]);
		}
		cp_model.AddAllDifferent(entry_ij);
	}

	// Orthogonality constraints
	for (i = 0; i < n; i++) {
		vector<IntVar> Ys;
		for (j = 0; j < n; j++) {
			Ys.push_back(y[i * n + j]);
		}
		cp_model.AddAllDifferent(Ys);

		for (j = 0; j < t; j++) {
			cp_model.AddVariableElement(x[i * n + j], Ys, z[i * t + j]);
		}
	}


	//==================================================================================
	// GUROBI MODEL
	//=================================================================================
	GRBEnv env = GRBEnv(true);
#if LOGGING
	env.set("LogFile", "LatinSquare.log");
#endif
	env.start();
	GRBModel model = GRBModel(env);

	Gx = new GRBVar[(int)pow(n, 4)];


	// Instantiate variables for L1 and L2
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				for (l = 0; l < n; l++) {
					string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
					Gx[i * nCube + j * nSqr + k * n + l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
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
					expr += Gx[i * nCube + j * nSqr + k * n + l];
				}
			}
			string s = "a_" + to_string(i) + "_" + to_string(j);
			model.addConstr(expr == 1.0, s);
		}
	}

	// Implement constraint (2)
	for (i = 0; i < n; i++) {
		for (k = 0; k < n; k++) {
			GRBLinExpr expr = 0;
			for (j = 0; j < n; j++) {
				for (l = 0; l < n; l++) {
					expr += Gx[i * nCube + j * nSqr + k * n + l];
				}
			}
			string s = "b_" + to_string(i) + "_" + to_string(k);
			model.addConstr(expr == 1.0, s);
		}
	}

	// Implement constraint (3)
	for (i = 0; i < n; i++) {
		for (l = 0; l < n; l++) {
			GRBLinExpr expr = 0;
			for (j = 0; j < n; j++) {
				for (k = 0; k < n; k++) {
					expr += Gx[i * nCube + j * nSqr + k * n + l];
				}
			}
			string s = "c_" + to_string(i) + "_" + to_string(l);
			model.addConstr(expr == 1.0, s);
		}
	}

	// Implement constraint (4)
	for (j = 0; j < n; j++) {
		for (k = 0; k < n; k++) {
			GRBLinExpr expr = 0;
			for (i = 0; i < n; i++) {
				for (l = 0; l < n; l++) {
					expr += Gx[i * nCube + j * nSqr + k * n + l];
				}
			}
			string s = "d_" + to_string(j) + "_" + to_string(k);
			model.addConstr(expr == 1.0, s);
		}
	}

	// Implement constraint (5)
	for (j = 0; j < n; j++) {
		for (l = 0; l < n; l++) {
			GRBLinExpr expr = 0;
			for (i = 0; i < n; i++) {
				for (k = 0; k < n; k++) {
					expr += Gx[i * nCube + j * nSqr + k * n + l];
				}
			}
			string s = "e_" + to_string(j) + "_" + to_string(l);
			model.addConstr(expr == 1.0, s);
		}
	}

	// Implement constraint (6)
	for (k = 0; k < n; k++) {
		for (l = 0; l < n; l++) {
			GRBLinExpr expr = 0;
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					expr += Gx[i * nCube + j * nSqr + k * n + l];
				}
			}
			string s = "f_" + to_string(k) + "_" + to_string(l);
			model.addConstr(expr == 1.0, s);
		}
	}

	// Fix first row of both squares

	for (i = 0; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += Gx[i * nSqr + i * n + i];
		model.addConstr(expr == 1.0);
	}

	// Fix first column of first square

	/*for (i = 1; i < n; i++) {
		for (j = 0; j < n; j++) {
			if(i != j)
			{
				for (k = 0; k < n; k++) {
					GRBLinExpr expr = 0;
					expr += Gx[i * nCube + 0 * nSqr + j * n + k];
					model.addConstr(expr == 0.0);
				}
			}
		}
	}*/

#if !defined(NOFIXCOL)	
	// Implement symmetry breaking - First column of each square
	for (i = 1; i < n; i++) {
		for (j = 1; j < n; j++) {
			for (k = 0; k < n; k++) {
				//cant be i or > i + 1
				if (j == i || j > i + 1 || k != i) {
					GRBLinExpr expr = 0.0;
					expr += Gx[i * nCube + k * n + j];
					model.addConstr(expr == 0.0);
				}
			}
		}
	}
#endif

	//================================================================================
	// BEGIN CP SOLVE FOR TRANSVERSALS IN X
	//================================================================================

	// First CP model keeps track of transversals in X
	Model model1;

	// Keep track of how many times we fix t transversals in x
	int num_calls = 0, infeasible_calls = 0;
	model1.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {

#if VERBOSE
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << SolutionIntegerValue(r, y[i * n + j]);
			}
			cout << endl;
		}
		cout << endl;
		
		for (i = 0; i < n; i++) {
			for (j = 0; j < t; j++) {
				cout << SolutionIntegerValue(r, x[i * n + j]);
			}
			cout << endl;
		}
		cout << endl;
#endif

		// Fix variables in Gurobi model to 1 based on variables fixed in X and Y 
		for (i = 0; i < n; i++) {
			for (j = 0; j < t; j++) {
				// Looks really ugly, but the logic is as follows:
				// i = row, x[i][j] = col, j = value in Y, y[i][x[i][j]] = y[i][col] = value in X

				// Set gurobi variables Gx[i][x[i][j]][j][y[i][x[i][j]] = 1
				GRBLinExpr expr = Gx[i * nCube + SolutionIntegerValue(r, x[i * n + j]) * nSqr + j + n * SolutionIntegerValue(r, y[i * n + SolutionIntegerValue(r, x[i * n + j])])];
				model.addConstr(expr == 1.0, "expr_" + to_string(i) + "_" + to_string(j));
			}
		}

#ifdef CALLBACK_CUTS
		// Don't allow presolve to translate constraints
		model.set(GRB_IntParam_PreCrush, 1);

		cout << "Using custom clique cuts" << endl;
		model.setCallback(new MyCallback);
#endif

#if SET_NODE_LIMIT
		model.set(GRB_DoubleParam_NodeLimit, NODE_LIMIT);
#endif
#if NO_CUTS
		model.set(GRB_IntParam_Cuts, 0);
#endif
		//model.set(GRB_IntParam_LogToConsole, 0);
		std::srand(std::time(nullptr));
		int seed = abs(std::rand() * std::rand());
		model.set(GRB_IntParam_Threads, 1);
		model.set(GRB_IntParam_Seed, seed);

		//==============================================================================
		// EXECUTE IP MODEL
		//==============================================================================
		try {
			model.optimize();
		}
		catch (GRBException e) {
			std::cout << "Error code = " << e.getErrorCode() << endl;
			std::cout << e.getMessage() << endl;
		}
		catch (...) {
			std::cout << "Unexpected exception during optimization" << endl;
		}

		num_calls++;

		// Model was proven to be infeasible
		if (model.get(GRB_IntAttr_Status) == 3) {
			infeasible_calls++;
		}
		
#if SET_NODE_LIMIT && NODELIMIT == 0
		cout << "Calls infeasible at root: " << infeasible_calls << "/" << num_calls << " = " << to_string(infeasible_calls/(float)num_calls) << endl;
#endif
		
		// If we have a feasible solution then print and exit
		if (model.get(GRB_IntAttr_Status) == 2) {
			printf("Optimal solution found\n");
			PRINT(Gx, n, t);
			exit(0);
		}

		// Remove variable fixing from this iteration
		for (i = 0; i < n; i++) {
			for (j = 0; j < t; j++) {
				GRBConstr c = model.getConstrByName("expr_" + to_string(i) + "_" + to_string(j));
				model.remove(c);
			}
		}

		// Update model
		model.update();

		}));

	// Execute model
	SatParameters param;
	param.set_enumerate_all_solutions(true);
	std::srand(std::time(nullptr));
	int seed = abs(std::rand() * std::rand());
	param.set_random_seed(seed);
	param.set_randomize_search(true);
	model1.Add(NewSatParameters(param));
	auto tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model1);
	auto toc = chrono::high_resolution_clock::now();
	// Report time
	LOG(INFO) << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms";
	// Fin
	return EXIT_SUCCESS;
}

// Helper function to print X and Y from the Gurobi Model
void PRINT(GRBVar* Gx, int n, int t) {
	int nCube = (int)pow(n, 3);
	int nSqr = (int)pow(n, 2);

	int i = 0, j = 0, k = 0, l = 0;
	// Print L1
	std::cout << endl;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			for (k = 0; k < n; k++) {
				for (l = 0; l < n; l++) {
					if (Gx[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
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
					if (Gx[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
						std::cout << to_string(l) << " ";
				}
			}
		}
		std::cout << endl;
	}
}
