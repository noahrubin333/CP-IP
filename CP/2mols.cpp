/**
	Noah Rubin
	2020-09-02
	This program uses OR-TOOLS and constraint programming to generate pairs of OLS(n).
	- Added symmetry breaking to improve runtime
**/
#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
using namespace std;
using namespace operations_research;
using namespace sat;
#define n ORDER
#define Z_VARS_LAST 1
#define APPA_BRANCHING 0
int main(int argc, char* argv[]) {
	// Model
	CpModelBuilder cp_model;
	// Fix general domains for x_ij, y_ij, z_ij
	Domain domain(0, n - 1);
#if !defined(BOOL) && !defined(INDEX)
	Domain domain2(0, n * n - 1);
#endif
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

#ifdef FIXCOL
	for (i = 1; i < n-1; i++) {
		firstCol_L2[i] = Domain::FromValues({ i + 1 });
	}
	firstCol_L2[n-1] = Domain::FromValues({ 1 });
#elif defined(NOSYM)
	// Do not fix first row or column
	for (i = 0; i < n; i++) {
		firstRow[i] = domain;
		firstCol_L1[i] = domain;
		firstCol_L2[i] = domain;
	}
#else
	// Fix entry (0,1) to 2 in L2
	firstCol_L2[1] = Domain::FromValues({ 2 });
	// Fix entry (0,N-2) to N-1 in L2
	firstCol_L2[n-2] = Domain::FromValues({ n - 1 });

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
#endif

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
#ifdef FIXDIAG
				if(i==j)
					x[i][j] = cp_model.NewIntVar(Domain::FromValues({ i==n-1 ? 1 : i+1 })).WithName("y_" + to_string(i) + to_string(j));
				else
#elif defined(FIXANTIDIAG)
				if(i==j)
					x[i][j] = cp_model.NewIntVar(Domain::FromValues({ 0 })).WithName("y_" + to_string(i) + to_string(j));
				else if(i==n-j-1)
					x[i][j] = cp_model.NewIntVar(Domain::FromValues({ n-1 })).WithName("y_" + to_string(i) + to_string(j));
				else
#endif
					x[i][j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
					
#ifdef FIXDIAG
				if(j==i-1)
					y[i][j] = cp_model.NewIntVar(Domain::FromValues({ 0 })).WithName("x_" + to_string(i) + to_string(j));
				else if(j==i+1)
					y[i][j] = cp_model.NewIntVar(Domain::FromValues({ 1 })).WithName("x_" + to_string(i) + to_string(j));
				else
#endif
					y[i][j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
			}
#if Z_VARS_LAST == 1
		}
	}
	for(i = 0; i < n; i++) {
		for(j = 0; j < n; j++) {
#endif
			// z_ij remains the same always
#if defined(BOOL) || defined(INDEX)
			z[i][j] = cp_model.NewIntVar(domain).WithName("z_" + to_string(i) + to_string(j));
#else
			z[i][j] = cp_model.NewIntVar(domain2).WithName("z_" + to_string(i) + to_string(j));
#endif
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
#if defined(BOOL) || defined(INDEX)
	// Column uniqueness constraints on Z
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(z[i][j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	// Row uniqueness constraints on Z
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(z[j][i]);
		}
		cp_model.AddAllDifferent(row_i);
	}

	// Constraints defining Z variables
#ifdef INDEX
	for(int i=0; i<n; i++)
	{	vector<IntVar> Zvec_i;
		for(int j=0; j<n; j++)
			Zvec_i.push_back(z[i][j]);
		for(int j=0; j<n; j++)
			cp_model.AddVariableElement(x[i][j], Zvec_i, y[i][j]);
	}
#else
	BoolVar a[n][n][n];
	BoolVar b[n][n][n];
	BoolVar c[n][n][n];
	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			for(int k=0; k<n; k++)
			{
				a[i][j][k] = cp_model.NewBoolVar();
				b[i][j][k] = cp_model.NewBoolVar();
				c[i][j][k] = cp_model.NewBoolVar();
				cp_model.AddNotEqual(x[i][j], k).OnlyEnforceIf(Not(a[i][j][k]));
				cp_model.AddNotEqual(y[i][j], k).OnlyEnforceIf(Not(b[i][j][k]));
				cp_model.AddEquality(z[i][j], k).OnlyEnforceIf(c[i][j][k]);
/* Constraints that don't seem to help
				cp_model.AddEquality(x[i][j], k).OnlyEnforceIf(a[i][j][k]);
				cp_model.AddEquality(y[i][j], k).OnlyEnforceIf(b[i][j][k]);
				cp_model.AddNotEqual(z[i][j], k).OnlyEnforceIf(Not(c[i][j][k]));
*/
			}

	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			for(int k=0; k<n; k++)
				for(int l=0; l<n; l++)
				{
					cp_model.AddBoolOr({Not(a[i][j][k]),Not(b[i][j][l]),c[i][k][l]});
#ifdef EXTRA
					cp_model.AddBoolOr({Not(a[i][j][k]),b[i][j][l],Not(c[i][k][l])});
					cp_model.AddBoolOr({a[i][j][k],Not(b[i][j][l]),Not(c[i][k][l])});
#endif
				}
#endif
				
#else // Use linear equations to encode orthogonality

#ifdef LINEAR_ELEM // Encode linear equality by element constraints
	int64_t Z_modulo_n[n*n], Z_division_n[n*n];
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			Z_modulo_n[i*n+j] = j;
			Z_division_n[i*n+j] = i;
		}
	}
	// Encode z[i][j] = x[i][j] + n * y[i][j]
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			cp_model.AddElement(z[i][j], Z_modulo_n, x[i][j]);
			cp_model.AddElement(z[i][j], Z_division_n, y[i][j]);
		}
	}
#else
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
#endif

	// Enforce all_different predicate on Z_ij forall i,j in {0,...,n}
	vector<IntVar> Z_VARS;
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			Z_VARS.push_back(z[i][j]);
		}
	}
	cp_model.AddAllDifferent(Z_VARS);
#endif

	// Tell model how to count solutions
	Model model;
	// Model parameters
	SatParameters param;
	// Use a single CPU core
	param.set_num_search_workers(1);

	if(argc >= 2) {
		int seed = stoi(argv[1]);
		param.set_random_seed(seed);
		// Randomize preferred variable order
		param.set_preferred_variable_order(SatParameters::IN_RANDOM_ORDER);
		cout << "Using random seed " << seed << endl;
	}

#if APPA_BRANCHING == 1
	vector <IntVar> VARS;
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			VARS.push_back(x[i][j]);
			VARS.push_back(y[i][j]);
			VARS.push_back(z[i][j]);
		}
	}
	cp_model.AddDecisionStrategy(VARS, DecisionStrategyProto::CHOOSE_FIRST, DecisionStrategyProto::SELECT_MIN_VALUE);
	param.set_search_branching(SatParameters::FIXED_SEARCH);
#endif

	model.Add(NewSatParameters(param));

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
	if(num_solutions>0)
		cout << "Solution found" << endl;
   	cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl;
	cout << CpSolverResponseStats(response);
	// Fin
	return EXIT_SUCCESS;
}
