/**
	Curtis Bright | Noah Rubin
	2021-03-12

	This program implements the pure CP model for solving 2MOLS(n).

	This program allows for the following options:

	i)		n	 					-> run the model in default state, full symm breaking
	ii)		n	S					-> run the model with specified level of symm breaking
										-> S = 'N' for none, S = 'R' for first rows Lexocographic, 'D' for dominance detection
	iii)	n	[l1,l2,...,ln]	-> run the model with full symm breaking and column 1 of Y fixed to [l1, l2, ..., ln]
										-> pass list without spaces, otherwise compiler sees multiple arguments
	**/

#pragma once

using namespace std;
using namespace operations_research;
using namespace sat;

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include <string>

// Here configure macros
//	Z_VARS_LAST -> apply ordering to 
#define Z_VARS_LAST 1

int main(int argc, char* argv[]) {


	Domain* firstCol_passed = NULL; // Expect possibility of column fixing
	int passed_list = 0; // Flag set if column fixing is passed at runtime


	int n = 0; // Store n

	int i = 0, j = 0, k = 0, l = 0; // General purpose counters

	char choice = 'D'; // Default to full symm breaking (w/ dominance detection)

	if (argc < 2) { // 
		cout << "Enter dimension of squares";
		exit(0);
	}
	if (argc == 2) {
		n = atoi(argv[1]);
	}

	// Allow passing of first column of Y as a list [Y_00, Y_10, ..., Y_(n-1)0]
	else if (strlen(argv[2]) > 1) {
		n = atoi(argv[1]); // Grab n
		firstCol_passed = new Domain[n];
		passed_list = 1; // Tell program a list was passed

		// Parse list from arguments
		int count = 0;
		for (i = 0; i < strlen(argv[2]) - 1; i++) {
			if (48 <= (int)argv[2][i] && (int)argv[2][i] <= 57) {
				firstCol_passed[count] = Domain::FromValues({ (int)argv[2][i] - 48 });
				count++;
			}
		}

		// Make sure a valid list was passed
		if (n != count) {
			cout << "n = " << n << " but you entered a list of length " << count << "\n";
			exit(0);
		}
	}

	// Handle case of symmetry breaking strategy being passed
	else if (argc == 3) {
		n = atoi(argv[1]);
		choice = argv[2][0];
	}

	// Handle any weird params being passed
	else {
		cout << "Unexpected parameters passsed, see code for possible values\n";
		exit(0);
	}

	// Model
	CpModelBuilder cp_model;
	// Fix general domains for x_ij, y_ij, z_ij
	Domain domain(0, n - 1);
#if !defined(BOOL) && !defined(INDEX)
	Domain domain2(0, n * n - 1);
#endif
	// Declare constrained domains for symmetry breaking
	Domain* firstRow = new Domain[n];
	Domain* firstCol_L1 = new Domain[n];
	Domain* firstCol_L2 = new Domain[n];
	// Declare variables x_ij forall i,j in {0, ..., n - 1}
	IntVar* x = new IntVar[n * n];
	IntVar* y = new IntVar[n * n];
	IntVar* z = new IntVar[n * n];

	// Fix first row domains
	for (i = 0; i < n; i++) {
		firstRow[i] = Domain::FromValues({ i });
		firstCol_L1[i] = Domain::FromValues({ i });
	}
#ifdef FIXCOL
	for (i = 1; i < n - 1; i++) {
		firstCol_L2[i] = Domain::FromValues({ i + 1 });
	}
	firstCol_L2[n - 1] = Domain::FromValues({ 1 });
#else
	// Implement dominance detection by propogating the following fixings:

	// Fix entry (0,1) to 2 in L2
	firstCol_L2[1] = Domain::FromValues({ 2 });
	// Fix entry (0,N-2) to N-1 in L2
	firstCol_L2[n - 2] = Domain::FromValues({ n - 1 });


	// Now fix remaining constrained domains for column 1 in L2
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
#endif
	//Declare variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == 0) { // First row
				if (choice == 'R' || choice == 'D' || passed_list) {
					x[i * n + j] = cp_model.NewIntVar(firstRow[j]).WithName("x_" + to_string(i) + to_string(j));
					y[i * n + j] = cp_model.NewIntVar(firstRow[j]).WithName("y_" + to_string(i) + to_string(j));
				}
				else {
					x[i * n + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
					y[i * n + j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
				}
			}
			else if (j == 0 && i > 0) { // First column
				if (passed_list) { // Set domains according to passed col 1 of Y
					x[i * n + j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
					y[i * n + j] = cp_model.NewIntVar(firstCol_passed[i]).WithName("y_" + to_string(i) + to_string(j));
				}
				else if (choice == 'D') { // Constrain domains of col 1 of Y
					x[i * n + j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
					y[i * n + j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
				}
				else { // Do not apply symm breaking to col 1 of Y
					x[i * n + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
					y[i * n + j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
				}
			}
			else { // Generic case
#ifdef FIXDIAG
				if (i == j)
					x[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ i == n - 1 ? 1 : i + 1 })).WithName("y_" + to_string(i) + to_string(j));
				else
#elif defined(FIXANTIDIAG)
				if (i == j)
					x[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ 0 })).WithName("y_" + to_string(i) + to_string(j));
				else if (i == n - j - 1)
					x[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ n - 1 })).WithName("y_" + to_string(i) + to_string(j));
				else
#endif
					// Generic domain for all entries not in first row / col of X
					x[i * n + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));

#ifdef FIXDIAG
				if (j == i - 1)
					y[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ 0 })).WithName("x_" + to_string(i) + to_string(j));
				else if (j == i + 1)
					y[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ 1 })).WithName("x_" + to_string(i) + to_string(j));
				else
#endif
					// Generic domain for all entries not in first row / col of Y
					y[i * n + j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
			}
#if Z_VARS_LAST == 1
			// We can add all Z_ij to the model last which forces OR-tools to branch on them after choosing values for X, Y
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
#endif
			// z_ij remains the same always
#if defined(BOOL) || defined(INDEX)
			// Instantiate Z_ij for using indexing constraints to encode orthogonality
			z[i * n + j] = cp_model.NewIntVar(domain).WithName("z_" + to_string(i) + to_string(j));
#else
			// Instantiate Z_ij for using linear constraints to encode orthogonality - much slower
			z[i * n + j] = cp_model.NewIntVar(domain2).WithName("z_" + to_string(i) + to_string(j));
#endif
		}
	}

	// Latin Property: Column uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(x[i * n + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}

	// Latin Property: Row uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(x[j * n + i]);
		}
		cp_model.AddAllDifferent(row_i);
	}


	// Latin Property: Column uniqueness constraints L2
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(y[i * n + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	// Latin Property: Row uniqueness constraints L2
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(y[j * n + i]);
		}
		cp_model.AddAllDifferent(row_i);
	}

	// What does this do?
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			i += 1;
		}
	}

#if defined(BOOL) || defined(INDEX)
	// Column uniqueness constraints on Z
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(z[i * n + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	// Row uniqueness constraints on Z
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(z[j * n + i]);
		}
		cp_model.AddAllDifferent(row_i);
	}
	// Constraints defining Z variables
#ifdef INDEX
	for (int i = 0; i < n; i++)
	{
		vector<IntVar> Zvec_i;
		for (int j = 0; j < n; j++)
			Zvec_i.push_back(z[i * n + j]);
		for (int j = 0; j < n; j++)
			cp_model.AddVariableElement(x[i * n + j], Zvec_i, y[i * n + j]);
	}
#else
	BoolVar a[n][n][n];
	BoolVar b[n][n][n];
	BoolVar c[n][n][n];
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
			{
				a[i * n + j][k] = cp_model.NewBoolVar();
				b[i * n + j][k] = cp_model.NewBoolVar();
				c[i * n + j][k] = cp_model.NewBoolVar();
				cp_model.AddNotEqual(x[i * n + j], k).OnlyEnforceIf(Not(a[i * n + j][k]));
				cp_model.AddNotEqual(y[i * n + j], k).OnlyEnforceIf(Not(b[i * n + j][k]));
				cp_model.AddEquality(z[i * n + j], k).OnlyEnforceIf(c[i * n + j][k]);
				/* Constraints that don't seem to help
								cp_model.AddEquality(x[i * n + j], k).OnlyEnforceIf(a[i * n + j][k]);
								cp_model.AddEquality(y[i * n + j], k).OnlyEnforceIf(b[i * n + j][k]);
								cp_model.AddNotEqual(z[i * n + j], k).OnlyEnforceIf(Not(c[i * n + j][k]));
				*/
			}
	for (int i = 0; i < n; i++)
		for (int j = 0; j < n; j++)
			for (int k = 0; k < n; k++)
				for (int l = 0; l < n; l++)
				{
					cp_model.AddBoolOr({ Not(a[i * n + j][k]),Not(b[i * n + j][l]),c[i][k][l] });
#ifdef EXTRA
					cp_model.AddBoolOr({ Not(a[i * n + j][k]),b[i * n + j][l],Not(c[i][k][l]) });
					cp_model.AddBoolOr({ a[i * n + j][k],Not(b[i * n + j][l]),Not(c[i][k][l]) });
#endif
				}
#endif

#else
	// RHS of Linear Expression is 0
	LinearExpr expr_0 = LinearExpr(0);
	// Build Z_ij using the logic Z_ij = X_ij + n * Y_ij  <==>  X_ij + n *  Y_ij - Z_ij = 0
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			LinearExpr expr = LinearExpr(0);
			expr.AddTerm(x[i * n + j], 1);  // X_ij
			expr.AddTerm(y[i * n + j], n);  //      + n * Y_ij
			expr.AddTerm(z[i * n + j], -1); //                 - Z_ij
			cp_model.AddEquality(expr, expr_0); //               = 0
		}
	}
	// Enforce all_different predicate on Z_ij forall i,j in {0,...,n}
	vector<IntVar> Z_VARS;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			Z_VARS.push_back(z[i * n + j]);
		}
	}
	cp_model.AddAllDifferent(Z_VARS);
#endif
	// Tell model how to count solutions
	Model model;

	// Add param to dictate the max time allowed
	SatParameters param;
	param.set_max_time_in_seconds(60000);
	// param.set_enumerate_all_solutions(true);
	model.Add(NewSatParameters(param));;

	int num_solutions = 0; // Can only ever be 1 unless solver is configured to find all 2MOLS(n)
	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
		// This callback function allows the solver to report solutions via stdout
		cout << "\n";
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << SolutionIntegerValue(r, x[i * n + j]) << " ";
			}
			cout << "\n";
		}
		cout << "\n\n";
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << SolutionIntegerValue(r, y[i * n + j]) << " ";
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
	if (num_solutions > 0)
		cout << "Solution found" << endl;

	// Report time, branches and conflict count

	cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << endl;
	cout << "Number of branches explored: " << response.num_branches() << "\n";
	cout << "Number of conflicts: " << response.num_conflicts() << "\n";

	// Fin
	return EXIT_SUCCESS;
}