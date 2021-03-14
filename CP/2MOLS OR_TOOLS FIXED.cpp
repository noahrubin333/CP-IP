/**
	Noah Rubin
	2020-09-02
	This program uses OR-TOOLS and constraint programming to generate pairs of OLS(n).
	- Added symmetry breaking to improve runtime
**/
#pragma once
#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
using namespace std;
using namespace operations_research;
using namespace sat;
#define n 8
int main(int argc, char* argv[]) {
	// Model
	CpModelBuilder cp_model;
	// Fix general domains for x_ij, y_ij, z_ij
	Domain domain(0, n - 1);
	Domain domain2(0, n * n - 1);
	// Declare constrained domains for symmetry breaking
	Domain firstRow[n];
	Domain firstCol_L1[n];
	Domain firstCol_L2[n];
	// Declare variables x_ij forall i,j in {0, ..., n - 1}
	IntVar x[n][n], y[n][n], z[n][n];
	//Counters
	int i = 0, j = 0;
	// Fix first row domains
	for (i = 0; i < n; i++) {
		firstRow[i] = Domain::FromValues({ i });
		firstCol_L1[i] = Domain::FromValues({ i });
	}

	// Fix entry (0,1) to 2 in L2 
	firstCol_L2[1] = Domain::FromValues({ 2 });

	/**
	// Fix remaining constrained domains for column 1 in L2
	for (i = 2; i < n; i++) {
		vector<int64> domainTemp;
		if (i == (n - 2)) {
			domainTemp.push_back((n - 1));
			firstCol_L2[i] = Domain::FromValues(domainTemp);
			continue;
		}
		if (i == (n - 1)) {
			domainTemp.push_back(1);
			for (j = 3; j <= (n - 2); j++) {
				domainTemp.push_back(j);
			}
			firstCol_L2[i] = Domain::FromValues(domainTemp);
			break;
		}
		for (j = 1; j <= i + 1; j++) {
			if (j != i && j != 2) {
				domainTemp.push_back(j);
			}
		}
		firstCol_L2[i] = Domain::FromValues(domainTemp);
	}
	**/

	int64 firstCol[] = { 3,4,5,6,7,1 };

	for (i = 2; i < n; i++) {
		firstCol_L2[i] = Domain::FromValues({ firstCol[i - 2] });
	}

	//Declare variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == 0) { // First row
				x[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("x_" + to_string(i) + to_string(j));
				y[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("y_" + to_string(i) + to_string(j));
			}
			else if (j == 0 && i > 0) { // First column
				x[i][j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
				y[i][j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
			}
			else { // Generic
				x[i][j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
				y[i][j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
			}
			// z_ij remains the same always
			z[i][j] = cp_model.NewIntVar(domain2).WithName("z_" + to_string(i) + to_string(j));
		}
	}
	// Column uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(x[i][j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	// Row uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(x[j][i]);
		}
		cp_model.AddAllDifferent(row_i);
	}
	// Column uniqueness constraints L2
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(y[i][j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	// Row uniqueness constraints L2
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(y[j][i]);
		}
		cp_model.AddAllDifferent(row_i);
	}
	// RHS of Linear Expression is 0
	LinearExpr expr_0 = LinearExpr(0);

	// Build Z_ij using the logic Z_ij = X_ij + n * Y_ij  <==>  X_ij + n *  Y_ij - Z_ij = 0
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			LinearExpr expr = LinearExpr(0);
			expr.AddTerm(x[i][j], 1);  // X_ij
			expr.AddTerm(y[i][j], n);  //      + n * Y_ij
			expr.AddTerm(z[i][j], -1); //                 - Z_ij
			cp_model.AddEquality(expr, expr_0); //               = 0
		}
	}
	// Enforce all_different predicate on Z_ij forall i,j in {0,...,n}
	vector<IntVar> Z_VARS;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			Z_VARS.push_back(z[i][j]);
		}
	}
	cp_model.AddAllDifferent(Z_VARS);

	// Tell model how to count solutions
	Model model;
	int num_solutions = 0;
	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
		cout << "\n";
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << SolutionIntegerValue(r, x[i][j]) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n";
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << SolutionIntegerValue(r, y[i][j]) << " ";
			}
			cout << "\n";
		}
		cout << "\n";
		num_solutions++;
		}));

	// Execute model
	auto tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
	auto toc = chrono::high_resolution_clock::now();
	// Report time
	LOG(INFO) << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms";
	// Fin
	return EXIT_SUCCESS;
}