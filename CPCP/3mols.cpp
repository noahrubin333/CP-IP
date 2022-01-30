/**
	Noah Rubin
	2020-09-15
	Modified by Curtis Bright
	Solves 3MOLS(n) by using a separate CP model to initialize the third square
**/

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"

#include <vector>
#include <tuple>
#include <string>
#include <iomanip>

using namespace std;
using namespace operations_research;
using namespace sat;

#if defined(CYCLETYPE)
	#define NO_DUPLICATE_CYCLE_TYPES 1
#else
	#define NO_DUPLICATE_CYCLE_TYPES 0
#endif

#define VERBOSE 0

#define n ORDER

#define STOP_ON_SOL 1

#define ENUMERATE_ALL 1

#define RANDOMIZE_SEARCH 1

#define EXTRA_SQUARE 0

#define TIMEOUT 60000			// Timeout in seconds

#if NO_DUPLICATE_CYCLE_TYPES==1 && !defined(NOFIXCOL) && !defined(FIXCOL)
#include "cycle_types.h"
#endif

chrono::time_point<chrono::high_resolution_clock> tic, toc, tic2, toc2;

int main(int argc, char* argv[]) {
	std::cout << std::setprecision(1);
	// Fix general domains for x_ij, y_ij, z_ij
	Domain domain(0, n - 1);

	//Counters
	int i = 0, j = 0;

	// L3 Model
	CpModelBuilder cp_model_1;

	// Declare variables L_ij for CP Model L3
	IntVar L[n][n];
	
	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			// z_0k for k = 0,...,n
			if (i == 0) {
				L[i][j] = cp_model_1.NewIntVar(Domain::FromValues({ j })).WithName("l_" + to_string(i) + to_string(j));
			}
			else if (j == 0) {
				Domain d = Domain::FromValues({ i }).Complement();
#if defined(FIXCOL)
				Domain e = Domain::FromValues({ i == n-1 ? 1 : i+1 }).Complement();
				L[i][j] = cp_model_1.NewIntVar(domain.IntersectionWith(d).IntersectionWith(e)).WithName("l_" + to_string(i) + to_string(j));
#elif defined(NOFIXCOL)
				L[i][j] = cp_model_1.NewIntVar(domain).WithName("l_" + to_string(i) + to_string(j));
#elif defined(FIXCOLTHIRD)
				L[i][j] = cp_model_1.NewIntVar(Domain(i)).WithName("l_" + to_string(i) + to_string(j));
#else
				Domain e = domain;
				if(i == 1)
					e = Domain::FromValues({ 2 }).Complement();
				if(i == n-2)
					e = Domain::FromValues({ n-1 }).Complement();
				L[i][j] = cp_model_1.NewIntVar(domain.IntersectionWith(d).IntersectionWith(e)).WithName("l_" + to_string(i) + to_string(j));
#endif
			}
			// General case
			else {
				L[i][j] = cp_model_1.NewIntVar(domain).WithName("l_" + to_string(i) + to_string(j));
			}
		}
	}

	// Column uniqueness constraints
	for (j = 0; j < n; j++) {
		vector<IntVar> col_j;
		for (i = 0; i < n; i++) {
			col_j.push_back(L[i][j]);
		}
		cp_model_1.AddAllDifferent(col_j);
	}

	// Row uniqueness
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(L[i][j]);
		}
		cp_model_1.AddAllDifferent(row_i);
	}

	// =========================================================
	// Tell model how to deal with solutions
	// =========================================================

	Model CPmodel;
	int num_solutions = 0;
	CPmodel.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
		num_solutions++;

		toc = chrono::high_resolution_clock::now();

		if (chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() > TIMEOUT*1000) {
			cout << "Completed " << to_string(num_solutions-1) << " trials" << endl;
			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl;
			exit(0);
		}

		int L3[n][n];

		for (i = 0; i < n; i++) {
		    for (j = 0; j < n; j++) {
		        L3[i][j] = SolutionIntegerValue(r, L[i][j]);
		    }
        }

#if VERBOSE > 1
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << L3[i][j] << " ";
			}
			cout << endl;
		}
#endif

		// Model
		CpModelBuilder cp_model;

		// Declare constrained domains for symmetry breaking
		Domain firstRow[n];
		Domain firstCol_L1[n];
		Domain firstCol_L2[n];
		// Declare variables x_ij forall i,j in {0, ..., n - 1}
		IntVar x[n][n], y[n][n], z[n][n], a[n][n], b[n][n], c[n][n], d[n][n];
		// Fix first row domains
		for (i = 0; i < n; i++) {
			firstRow[i] = Domain::FromValues({ i });
#if defined(NOFIXCOL) || defined(FIXCOLTHIRD)
			firstCol_L1[i] = domain;
#else
			firstCol_L1[i] = Domain::FromValues({ i });
#endif
		}

#ifdef FIXCOL
		for (i = 1; i < n-1; i++) {
			firstCol_L2[i] = Domain::FromValues({ i + 1 });
		}
		firstCol_L2[n-1] = Domain::FromValues({ 1 });
#elif defined(NOFIXCOL) || (defined(FIXCOLTHIRD) && !defined(DOMAINRED))
		for (i = 1; i < n; i++) {
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

		// Declare variables and provide their associated domains
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == 0) { // First row
					x[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("x_" + to_string(i) + to_string(j));
					y[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("y_" + to_string(i) + to_string(j));
				}
#ifndef NOFIXCOL
				else if (j == 0 && i > 0) { // First column
					x[i][j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
					y[i][j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
				}
#endif
				else { // Generic
					x[i][j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
					y[i][j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
				}
				z[i][j] = cp_model.NewIntVar(Domain::FromValues({ L3[i][j] })).WithName("z_" + to_string(i) + to_string(j));
			}
		}
		
		// Declare auxiliary variables to enforce orthogonality
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				a[i][j] = cp_model.NewIntVar(domain).WithName("a_" + to_string(i) + to_string(j));
				b[i][j] = cp_model.NewIntVar(domain).WithName("b_" + to_string(i) + to_string(j));
				c[i][j] = cp_model.NewIntVar(domain).WithName("c_" + to_string(i) + to_string(j));
#if EXTRA_SQUARE == 1
				d[i][j] = cp_model.NewIntVar(domain).WithName("d_" + to_string(i) + to_string(j));
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
		// Column uniqueness constraints L3
		for (i = 0; i < n; i++) {
			vector<IntVar> col_i;
			for (j = 0; j < n; j++) {
				col_i.push_back(z[i][j]);
			}
			cp_model.AddAllDifferent(col_i);
		}
		// Row uniqueness constraints L3
		for (i = 0; i < n; i++) {
			vector<IntVar> row_i;
			for (j = 0; j < n; j++) {
				row_i.push_back(z[j][i]);
			}
			cp_model.AddAllDifferent(row_i);
		}
		// Column uniqueness constraints on A
		for (i = 0; i < n; i++) {
			vector<IntVar> col_i;
			for (j = 0; j < n; j++) {
				col_i.push_back(a[i][j]);
			}
			cp_model.AddAllDifferent(col_i);
		}
		// Row uniqueness constraints on A
		for (i = 0; i < n; i++) {
			vector<IntVar> row_i;
			for (j = 0; j < n; j++) {
				row_i.push_back(a[j][i]);
			}
			cp_model.AddAllDifferent(row_i);
		}
		// Column uniqueness constraints on B
		for (i = 0; i < n; i++) {
			vector<IntVar> col_i;
			for (j = 0; j < n; j++) {
				col_i.push_back(b[i][j]);
			}
			cp_model.AddAllDifferent(col_i);
		}
		// Row uniqueness constraints on B
		for (i = 0; i < n; i++) {
			vector<IntVar> row_i;
			for (j = 0; j < n; j++) {
				row_i.push_back(b[j][i]);
			}
			cp_model.AddAllDifferent(row_i);
		}
		// Column uniqueness constraints on C
		for (i = 0; i < n; i++) {
			vector<IntVar> col_i;
			for (j = 0; j < n; j++) {
				col_i.push_back(c[i][j]);
			}
			cp_model.AddAllDifferent(col_i);
		}
		// Row uniqueness constraints on C
		for (i = 0; i < n; i++) {
			vector<IntVar> row_i;
			for (j = 0; j < n; j++) {
				row_i.push_back(c[j][i]);
			}
			cp_model.AddAllDifferent(row_i);
		}
#if EXTRA_SQUARE == 1
		// Column uniqueness constraints on D
		for (i = 0; i < n; i++) {
			vector<IntVar> col_i;
			for (j = 0; j < n; j++) {
				col_i.push_back(d[i][j]);
			}
			cp_model.AddAllDifferent(col_i);
		}
		// Row uniqueness constraints on D
		for (i = 0; i < n; i++) {
			vector<IntVar> row_i;
			for (j = 0; j < n; j++) {
				row_i.push_back(d[j][i]);
			}
			cp_model.AddAllDifferent(row_i);
		}
#endif

		// Constraints defining A variables
		for(int i=0; i<n; i++)
		{	vector<IntVar> Avec_i;
			for(int j=0; j<n; j++)
				Avec_i.push_back(a[i][j]);
			for(int j=0; j<n; j++)
				cp_model.AddVariableElement(x[i][j], Avec_i, y[i][j]);
		}
		// Constraints defining B variables
		for(int i=0; i<n; i++)
		{	vector<IntVar> Bvec_i;
			for(int j=0; j<n; j++)
				Bvec_i.push_back(b[i][j]);
			for(int j=0; j<n; j++)
				cp_model.AddVariableElement(x[i][j], Bvec_i, z[i][j]);
		}
		// Constraints defining C variables
		for(int i=0; i<n; i++)
		{	vector<IntVar> Cvec_i;
			for(int j=0; j<n; j++)
				Cvec_i.push_back(c[i][j]);
			for(int j=0; j<n; j++)
				cp_model.AddVariableElement(y[i][j], Cvec_i, z[i][j]);
		}
#if EXTRA_SQUARE == 1
		// Constraints defining D variables
		for(int i=0; i<n; i++)
		{	vector<IntVar> Dvec_i;
			for(int j=0; j<n; j++)
				Dvec_i.push_back(d[j][i]);
			for(int j=0; j<n; j++)
				cp_model.AddVariableElement(a[j][i], Dvec_i, b[j][i]);
		}
#endif

#if NO_DUPLICATE_CYCLE_TYPES==1 && !defined(NOFIXCOL) && !defined(FIXCOL)
	BoolVar by[n][n];
	for(int i=2; i<n; i++) {
		if((n==9 && (i==8||i==2)) || (n==10 && (i==9||i==4)) || (n==11 && (i==10||i==4)) || n >= 12) {
			for(int j=1; j<n; j++) {
				by[i][j] = cp_model.NewBoolVar();
				cp_model.AddNotEqual(y[i][0], j).OnlyEnforceIf(Not(by[i][j]));
			}
		}
	}
	if (n==6) {
		cp_model.AddEquality(y[3][0], 4);
	} else if (n==7) {
		cp_model.AddNotEqual(y[4][0], 1);
	} else if (n==8) {
		cp_model.AddEquality(y[5][0], 6);
		cp_model.AddNotEqual(y[4][0], 1);
	} else if (n==9) {
		cp_model.AddNotEqual(y[6][0], 3);
		cp_model.AddNotEqual(y[6][0], 1);
		cp_model.AddNotEqual(y[6][0], 4);
		cp_model.AddBoolOr({Not(by[8][6]), Not(by[2][3])});
		cp_model.AddBoolOr({Not(by[8][7]), Not(by[2][3])});
	} else if (n==10) {
		cp_model.AddEquality(y[7][0], 8);
		cp_model.AddNotEqual(y[5][0], 1);
		cp_model.AddNotEqual(y[5][0], 4);
		cp_model.AddNotEqual(y[6][0], 1);
		cp_model.AddNotEqual(y[6][0], 3);
		cp_model.AddBoolOr({Not(by[9][7]), Not(by[4][1])});
	} else if (n==11) {
		cp_model.AddNotEqual(y[8][0], 1);
		cp_model.AddNotEqual(y[8][0], 3);
		cp_model.AddNotEqual(y[8][0], 4);
		cp_model.AddNotEqual(y[8][0], 5);
		cp_model.AddNotEqual(y[8][0], 6);
		cp_model.AddNotEqual(y[7][0], 1);
		cp_model.AddNotEqual(y[7][0], 3);
		cp_model.AddNotEqual(y[7][0], 4);
		cp_model.AddNotEqual(y[7][0], 6);
		cp_model.AddNotEqual(y[6][0], 1);
		cp_model.AddNotEqual(y[5][0], 4);
		cp_model.AddBoolOr({Not(by[10][7]), Not(by[4][1])});
		cp_model.AddBoolOr({Not(by[10][8]), Not(by[4][1])});
		cp_model.AddBoolOr({Not(by[10][9]), Not(by[4][1])});
		cp_model.AddBoolOr({Not(by[10][9]), Not(by[4][5])});
	} else {
		for(int i=0; i<blocked_cycle_types.size(); i++) {
			vector<BoolVar> clause;
			for(int j=2; j<n; j++) {
				if(j != n-2) {
					int k = blocked_cycle_types[i][j];
					clause.push_back(Not(by[j][k]));
				}
			}
			cp_model.AddBoolOr(clause);
		}
	}
#endif

		// Tell model how to count solutions
		Model model_new;

		/*if(argc >= 2) {
			SatParameters param;
			int seed = stoi(argv[1]);
			param.set_random_seed(seed);
			param.set_randomize_search(true);
			model.Add(NewSatParameters(param));;
			cout << "Using random seed " << seed << endl;
		}*/

		model_new.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
			cout << "Solution found in trial # " << num_solutions << endl;
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
			cout << "\n\n";
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					cout << SolutionIntegerValue(r, z[i][j]) << " ";
				}
				cout << "\n";
			}
			cout << "\n";

			toc = chrono::high_resolution_clock::now();

			if (STOP_ON_SOL) {
				cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms" << endl;
				exit(0);
			}

			if (ENUMERATE_ALL) {
				// Report time
				cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms" << endl;
			}
		}));

		// Execute model
		tic2 = chrono::high_resolution_clock::now();
		const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model_new);
		toc2 = chrono::high_resolution_clock::now();
#if VERBOSE > 1
		cout << "No solutions" << endl;
#endif
		toc = chrono::high_resolution_clock::now();

#if VERBOSE == 0
		if(num_solutions % 1000 == 0)
#endif
		cout << fixed << "Trial # " << num_solutions << " Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc2 - tic2).count() << " ms" << " Total time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms Avg. time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count()/((double)1000*num_solutions) << " s" << endl;

		}));

	// Params

	SatParameters parameters;
	if (ENUMERATE_ALL) {
		parameters.set_enumerate_all_solutions(true);
	}
#if RANDOMIZE_SEARCH == 1
	parameters.set_randomize_search(true);
#endif
	CPmodel.Add(NewSatParameters(parameters));

	// Execute
	tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model_1.Build(), &CPmodel);
	toc = chrono::high_resolution_clock::now();

	// Report time
	cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms" << endl;

	// Fin
	return EXIT_SUCCESS;
}
