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
#define Z_VARS_LAST 1
#define INDEX 1
#define NOSYM 1
int main(int argc, char* argv[]) {

	char choice = 'D'; // Default to full symm breaking (w/ dominance detection)
	int passed_list = 0; // Flag set if column fixing is passed at runtime

	vector<Domain> firstCol_passed; // Expect possibility of column fixing
	if (argc < 2) { // 
		cout << "Enter dimension of squares";
		exit(0);
	}

	int n = stoi(argv[1]);

	// Allow passing of first column of Y as a list [Y_00, Y_10, ..., Y_(n-1)0]
	if (argc == 3) {
		if (strlen(argv[2]) > 1) {
			passed_list = 1;
			string list = string(argv[2]);
			istringstream argvStream(list.substr(1, list.size() - 2));
			string token;
			int count = 0;
			while (getline(argvStream, token, ',')) {
				firstCol_passed.push_back(Domain::FromValues({ stoi(token) }));
				count++;
			}

			// Make sure a valid list was passed
			if (n != count) {
				cout << "n = " << n << " but you entered a list of length " << count << "\n";
				exit(0);
			}
		}
		else {
			n = stoi(argv[1]);
			choice = argv[2][0];
		}
	}

	// Model
	CpModelBuilder cp_model;
	// Fix general domains for x_ij, y_ij, z_ij
	Domain domain(0, n - 1);
#if !defined(BOOL) && !defined(INDEX)
	Domain domain2(0, n * n - 1);
#endif
	// Declare constrained domains for symmetry breaking
	Domain *firstRow = new Domain[n];
	Domain *firstCol_L1 = new Domain[n];
	Domain *firstCol_L2 = new Domain[n];
	// Declare variables x_ij forall i,j in {0, ..., n - 1}
	IntVar *x = new IntVar[n * n];
	IntVar* y = new IntVar[n * n];
	IntVar* z = new IntVar[n * n];

	//Counters
	int i = 0, j = 0;
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
	firstCol_L2[n - 2] = Domain::FromValues({ n - 1 });

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
				x[i * n + j] = cp_model.NewIntVar(firstRow[j]).WithName("x_" + to_string(i) + to_string(j));
				y[i * n + j] = cp_model.NewIntVar(firstRow[j]).WithName("y_" + to_string(i) + to_string(j));
			}
			else if (j == 0 && i > 0) { // First column
				if (passed_list) { // Set domains according to passed col 1 of Y
					x[i * n + j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
					y[i * n + j] = cp_model.NewIntVar(firstCol_passed[i]).WithName("y_" + to_string(i) + to_string(j));
				}
				else {
					x[i * n + j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
					y[i * n + j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
				}
			}
			else { // Generic
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
					x[i * n + j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));

#ifdef FIXDIAG
				if (j == i - 1)
					y[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ 0 })).WithName("x_" + to_string(i) + to_string(j));
				else if (j == i + 1)
					y[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ 1 })).WithName("x_" + to_string(i) + to_string(j));
				else
#endif
					y[i * n + j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
			}
#if Z_VARS_LAST == 1
		}
	}
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
#endif
			// z_ij remains the same always
#if defined(BOOL) || defined(INDEX)
			z[i * n + j] = cp_model.NewIntVar(domain).WithName("z_" + to_string(i) + to_string(j));
#else
			z[i * n + j] = cp_model.NewIntVar(domain2).WithName("z_" + to_string(i) + to_string(j));
#endif
		}
	}
	// Column uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(x[i * n + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	// Row uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(x[j * n + i]);
		}
		cp_model.AddAllDifferent(row_i);
	}
	// Column uniqueness constraints L2
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(y[i * n + j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	// Row uniqueness constraints L2
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(y[j * n + i]);
		}
		cp_model.AddAllDifferent(row_i);
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

#else // Use linear equations to encode orthogonality

#ifdef LINEAR_ELEM // Encode linear equality by element constraints
	vector<int64> Z_modulo_n(n* n);
	vector<int64> Z_division_n(n* n);
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			Z_modulo_n[i * n + j] = j;
			Z_division_n[i * n + j] = i;
		}
	}
	// Encode z[i * n + j] = x[i * n + j] + n * y[i * n + j]
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			cp_model.AddElement(z[i * n + j], Z_modulo_n, x[i * n + j]);
			cp_model.AddElement(z[i * n + j], Z_division_n, y[i * n + j]);
		}
	}
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
#endif

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

	if (argc >= 2) {
		SatParameters param;
		int seed = stoi(argv[1]);
		param.set_random_seed(seed);
		param.set_max_time_in_seconds(60000);
		param.set_randomize_search(true);
		model.Add(NewSatParameters(param));
		cout << "Using random seed " << seed << endl;
	}

#if APPA_BRANCHING == 1
	vector <IntVar> VARS;
	for (int i = 0; i < n; i++) {
		for (int j = 0; j < n; j++) {
			VARS.push_back(x[i * n + j]);
			VARS.push_back(y[i * n + j]);
			VARS.push_back(z[i * n + j]);
		}
	}
	cp_model.AddDecisionStrategy(VARS, DecisionStrategyProto::CHOOSE_FIRST, DecisionStrategyProto::SELECT_MIN_VALUE);
#endif

	int num_solutions = 0;
	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
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
	// Report time
	if (num_solutions > 0)
		cout << "Solution found" << endl;
	cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl;
	cout << CpSolverResponseStats(response);
	// Fin
	return EXIT_SUCCESS;
}