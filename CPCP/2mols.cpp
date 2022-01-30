/**
	Noah Rubin
	2020-09-15
	Modified by Curtis Bright
	Solves 2MOLS(n) by using a separate CP model to initialize the first square
**/

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"

#include <vector>
#include <tuple>
#include <string>

using namespace std;
using namespace operations_research;
using namespace sat;

#define VERBOSE 1

#define n ORDER

#define STOP_ON_SOL 1

#define ENUMERATE_ALL 1

#define TIMEOUT 60000			// Timeout in seconds

#define MAX_TRIALS 100000

#define ROW_WISE 1			// Order variables row-wise instead of column-wise

#define MIN_VALUE_DOMAIN_REDUCTION 1	// Select minimum value in the branching variable's domain

chrono::time_point<chrono::high_resolution_clock> tic, toc;

int main(int argc, char* argv[]) {
	// Model
	CpModelBuilder cp_model_1;
	
	// Fix general domains for x_ij
	Domain domain(0, n - 1);

	// Declare constrained domains for symmetry breaking
	Domain firstRow[n];
	Domain firstCol[n];

	// Declare variables x_ij forall i,j in {0, ..., n - 1}
	IntVar L[n][n];
	//Counters
	int i = 0, j = 0;
	// Fix first row domains
	for (i = 0; i < n; i++) {
		firstRow[i] = Domain::FromValues({ i });
		firstCol[i] = Domain::FromValues({ i });
	}

#if defined(NOSYM)
	// Do not fix first row or column
	for (i = 0; i < n; i++) {
		firstRow[i] = domain;
		firstCol[i] = domain;
	}
#endif

	//Declare variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == 0) { // First row
				L[i][j] = cp_model_1.NewIntVar(firstRow[j]).WithName("L_" + to_string(i) + to_string(j));
			}
			else if (j == 0 && i > 0) { // First column
				L[i][j] = cp_model_1.NewIntVar(firstCol[i]).WithName("L_" + to_string(i) + to_string(j));
			}
			else { // Generic
				L[i][j] = cp_model_1.NewIntVar(domain).WithName("L_" + to_string(i) + to_string(j));
			}
		}
	}
	// Column uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> col_i;
		for (j = 0; j < n; j++) {
			col_i.push_back(L[i][j]);
		}
		cp_model_1.AddAllDifferent(col_i);
	}
	// Row uniqueness constraints L1
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(L[j][i]);
		}
		cp_model_1.AddAllDifferent(row_i);
	}

	// =========================================================
	// Tell model how to deal with solutions
	// =========================================================

	Model CPmodel;
	int num_solutions = 0;
	CPmodel.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {

		toc = chrono::high_resolution_clock::now();

		if (chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() > TIMEOUT*1000 || num_solutions >= MAX_TRIALS) {
#if VERBOSE == 1
			cout << "Stopping on trial # " << to_string(num_solutions) << endl;
#endif
			cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl;
			exit(0);
		}

#if VERBOSE == 0
		if(num_solutions % 1000 == 0)
#endif
		cout << "Trial # " << num_solutions << " Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl;
		int L1[n][n];

		for (i = 0; i < n; i++) {
		    for (j = 0; j < n; j++) {
		        L1[i][j] = SolutionIntegerValue(r, L[i][j]);
		    }
        }

#if VERBOSE >= 1
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << L1[i][j] << " ";
			}
			cout << endl;
		}
#endif

		// Model
		CpModelBuilder cp_model;

		// Declare constrained domains for symmetry breaking
		Domain firstRow[n];
		Domain firstCol_L2[n];
		// Declare variables x_ij forall i,j in {0, ..., n - 1}
		IntVar x[n][n], y[n][n], a[n][n];
		// Fix first row domains
		for (i = 0; i < n; i++) {
			firstRow[i] = Domain::FromValues({ i });
		}

#ifdef FIXCOL
		for (i = 1; i < n-1; i++) {
			firstCol_L2[i] = Domain::FromValues({ i + 1 });
		}
		firstCol_L2[n-1] = Domain::FromValues({ 1 });
#elif defined(NOFIXCOL)
		for (i = 1; i < n; i++) {
			firstCol_L2[i] = domain;
		}
#else
		vector<Domain> firstCol_passed;
		// Allow passing of first column of Y as a list [Y_00,Y_10,...,Y_(n-1)0]
		if(argc >= 2) {
			string list = string(argv[1]);
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
		}
#endif

		// Declare variables and provide their associated domains
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				if (i == 0) { // First row
					y[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("y_" + to_string(i) + to_string(j));
				}
#if !defined(FIXCOL) && !defined(NOFIXCOL)
				else if (j == 0 && i > 0 && argc >= 2) {
					y[i][j] = cp_model.NewIntVar(firstCol_passed[i]).WithName("y_" + to_string(i) + to_string(j));
				}
				else if (j == 0 && i > 0) {
					y[i][j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
				}
#elif defined(FIXCOL)
				else if (j == 0 && i > 0) {
					y[i][j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
				}
#elif defined(NOFIXCOL)
				else if (j == 0 && i > 0) {
					y[i][j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
				}
#endif
				else { // Generic
					y[i][j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
				}
				x[i][j] = cp_model.NewIntVar(Domain::FromValues({ L1[i][j] })).WithName("z_" + to_string(i) + to_string(j));
			}
		}
		
		// Declare auxiliary variables to enforce orthogonality
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				a[i][j] = cp_model.NewIntVar(domain).WithName("a_" + to_string(i) + to_string(j));
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

		// Constraints defining A variables
		for(int i=0; i<n; i++)
		{	vector<IntVar> Avec_i;
			for(int j=0; j<n; j++)
				Avec_i.push_back(a[i][j]);
			for(int j=0; j<n; j++)
				cp_model.AddVariableElement(x[i][j], Avec_i, y[i][j]);
		}

		vector<IntVar> myvec;
#if ROW_WISE == 1
		for(int i=0; i < n; i++)
			for(int j=0; j < n; j++)
#else
		for(int j=0; j < n; j++)
			for(int i=0; i < n; i++)
#endif
				myvec.push_back(x[i][j]);
#if ROW_WISE == 1
		for(int i=0; i < n; i++)
			for(int j=0; j < n; j++)
#else
		for(int j=0; j < n; j++)
			for(int i=0; i < n; i++)
#endif
				myvec.push_back(y[i][j]);
#if ROW_WISE == 1
		for(int i=0; i < n; i++)
			for(int j=0; j < n; j++)
#else
		for(int j=0; j < n; j++)
			for(int i=0; i < n; i++)
#endif
				myvec.push_back(a[i][j]);

		absl::Span<IntVar> variables(myvec);		

#if MIN_VALUE_DOMAIN_REDUCTION == 1		
		cp_model.AddDecisionStrategy(variables, DecisionStrategyProto::CHOOSE_FIRST, DecisionStrategyProto::SELECT_MIN_VALUE);
#else
		cp_model.AddDecisionStrategy(variables, DecisionStrategyProto::CHOOSE_FIRST, DecisionStrategyProto::SELECT_MAX_VALUE);
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
#ifdef LATEX
			cout << "\\begin{tabular}{|cccccccc|}\n";
#endif
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					cout << SolutionIntegerValue(r, x[i][j]) << " ";
#ifdef LATEX
					cout << "& ";
#endif
				}
#ifdef LATEX
				cout << "\\\\";
#endif
				cout << "\n";
			}
#ifdef LATEX
			cout << "\\end{tabular}\n";
#endif
			cout << "\n";
			cout << "\n";
#ifdef LATEX
			cout << "\\begin{tabular}{|cccccccc|}\n";
#endif
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					cout << SolutionIntegerValue(r, y[i][j]) << " ";
#ifdef LATEX
					cout << "& ";
#endif
				}
#ifdef LATEX
				cout << "\\\\";
#endif
				cout << "\n";
			}
#ifdef LATEX
			cout << "\\end{tabular}\n";
#endif
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
		const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model_new);

#if VERBOSE > 1
		cout << "No solutions" << endl;
#endif

		num_solutions++;

		}));

	// Params

	SatParameters parameters;
	if (ENUMERATE_ALL) {
		parameters.set_enumerate_all_solutions(true);
	}

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
