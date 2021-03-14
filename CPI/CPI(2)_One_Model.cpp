/**
	Noah Rubin
	2020-09-02
	This program uses OR-TOOLS to generate t transversals of square x, it then fills in the same entries in y and
	passes the model to gurobi, which attempts to solve
**/
#pragma once
#include "gurobi_c++.h"

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"

using namespace std;
using namespace operations_research;
using namespace sat;

void PRINT(GRBVar* x, int n, int t);

// Start time, end time
std::chrono::steady_clock::time_point tic, toc;

int main(int argc, char* argv[]) {

	int n = -1, t = -1; // Default

	if (argc != 3) {
		cout << "INVALID NUMBER OF ARGUMENTS, NEED n AND t\n";
		exit(0);
	}
	else {
		n = atoi(argv[1]);
		t = atoi(argv[2]);
	}

	// Keep track of n^3 and n^2 for indexing
	int nCube = (int)pow(n, 3);
	int nSqr = (int)pow(n, 2);

	// CP Model Builder
	CpModelBuilder cp_model;

	// CP Domain (general case)
	Domain domain(0, n - 1);

	// Declare variables for the two models
	IntVar* x = new IntVar[n * t];
	IntVar* z = new IntVar[n * t];
	IntVar* y = new IntVar[n * n];

	//Counters
	int i = 0, j = 0, k = 0, l = 0;

	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// First row
			if (i == 0) {
				if (j < t) {
					x[i * t + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("x_" + to_string(i) + to_string(j));
					z[i * t + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("z_" + to_string(i) + to_string(j));
				}
				y[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("y_" + to_string(i) + to_string(j));
			}
			// General case
			else {
				if (j < t) {
					x[i * t + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
					z[i * t + j] = cp_model.NewIntVar(domain).WithName("z_" + to_string(i) + to_string(j));
				}
				y[i * n + j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
			}
		}
	}

	//=====================================================================
	// Now build model for fixing transversals in X
	//====================================================================

	// Column uniqueness constraints
	for (j = 0; j < t; j++) {
		vector<IntVar> col_i;
		for (i = 0; i < n; i++) {
			col_i.push_back(x[i * t + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}

	// Row uniqueness constraints
	for (i = 0; i < n; i++) {
		vector<IntVar> entry_ij;
		for (j = 0; j < t; j++) {
			entry_ij.push_back(x[i * t + j]);
		}
		cp_model.AddAllDifferent(entry_ij);
	}

	// Enforce z[i][k] = y[i, x[i][k]]
	for (i = 0; i < n; i++) {
		vector<IntVar> Ys;
		for (j = 0; j < n; j++) {
			Ys.push_back(y[i * n + j]);
		}

		for (j = 0; j < t; j++) {
			cp_model.AddVariableElement(x[i * t + j], Ys, z[i * t + j]);
		}
	}

	// Column uniqueness for Y
	for (j = 0; j < n; j++) {
		vector<IntVar> Ys;
		for (i = 0; i < n; i++) {
			Ys.push_back(y[i * n + j]);
		}
		cp_model.AddAllDifferent(Ys);
	}

	// All differents for Z
	for (j = 0; j < t; j++) {
		vector<IntVar> Zs;
		for (i = 0; i < n; i++) {
			Zs.push_back(z[i * t + j]);
		}
		cp_model.AddAllDifferent(Zs);
	}
	for (i = 0; i < n; i++) {
		vector<IntVar> Zs;
		for (j = 0; j < t; j++) {
			Zs.push_back(z[i * t + j]);
		}
		cp_model.AddAllDifferent(Zs);
	}


	//==================================================================================
	// GUROBI MODEL
	//=================================================================================
	GRBEnv env = GRBEnv(true);
	env.set("LogFile", "LatinSquare.log");
	env.start();
	GRBModel model = GRBModel(env);

	GRBVar* Gx = new GRBVar[(int)pow(n, 4)];


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
			string s = "b_" + to_string(i) + "_" + to_string(j);
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
			string s = "c_" + to_string(i) + "_" + to_string(j);
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
			string s = "d_" + to_string(i) + "_" + to_string(j);
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
			string s = "e_" + to_string(i) + "_" + to_string(j);
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
			string s = "f_" + to_string(i) + "_" + to_string(j);
			model.addConstr(expr == 1.0, s);
		}
	}

	// Fix first row of both squares

	for (i = 0; i < n; i++) {
		GRBLinExpr expr = 0;
		expr += Gx[i * nSqr + i * n + i];
		model.addConstr(expr == 1.0);
	}

	// Fix model params
	//model.set(GRB_DoubleParam_NodeLimit, 1);
	model.set(GRB_IntParam_LogToConsole, 0);

	/**
	// Implement symmetry breaking - First column of each square
	// Only seems to slow down model
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
	**/
	//================================================================================
	// BEGIN CP SOLVE FOR TRANSVERSALS IN X
	//================================================================================

	// First CP model keeps track of transversals in X
	Model model1;

	// Keep track of how many times we fix t transversals in x
	int num_solutions = 0;
	model1.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
		// Fix variables in Gurobi model to 1 based on variables fixed in X and Y 
		for (i = 0; i < n; i++) {
			for (j = 0; j < t; j++) {
				// Looks really ugly, but the logic is as follows:
				// i = row, x[i][j] = col, j = value in X, y[i][x[i][j]] = y[i][col] = value in Y

				// Set gurobi variables Gx[i][x[i][j]][j][y[i][x[i][j]] = 1
				GRBLinExpr expr = Gx[i * nCube + SolutionIntegerValue(r, x[i * t + j]) * nSqr + j * n + SolutionIntegerValue(r, y[i * n + SolutionIntegerValue(r, x[i * t + j])])];
				string name = "expr_" + to_string(i) + to_string(j);
				model.addConstr(expr == 1.0, name);
			}
		}

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

		// If we have a feasible solution then print and exit
		if (model.get(GRB_IntAttr_Status) == 2) {
			PRINT(Gx, n, t);
			toc = chrono::high_resolution_clock::now();
			cout << "\nTime elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms\n";
			cout << "Evaluated " << num_solutions << " trials before finding a solution\n";
			exit(EXIT_SUCCESS);
		}

		// Remove variable fixing from this iteration
		for (i = 0; i < n; i++) {
			for (j = 0; j < t; j++) {
				string name = "expr_" + to_string(i) + to_string(j);
				GRBConstr c = model.getConstrByName(name);
				model.remove(c);
			}
		}

		// Force model to update
		model.update();

		// Increment solution count
		num_solutions++;
		}));

	// Execute CP model
	SatParameters param;
	param.set_enumerate_all_solutions(true);
	model1.Add(NewSatParameters(param));

	// Begin timing
	tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model1);
	
	// Never reached
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