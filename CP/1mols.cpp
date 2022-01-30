/**
	Noah Rubin
	2020-09-02
	This program uses OR-TOOLS and constraint programming to generate pairs of OLS(n).
	- Added symmetry breaking to improve runtime
	Modified by Curtis to search for a single Latin square
**/
#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
using namespace std;
using namespace operations_research;
using namespace sat;
#define n ORDER
#define APPA_BRANCHING 1
int main(int argc, char* argv[]) {
	
	// Model
	CpModelBuilder cp_model;
	
	// Fix general domains for x_ij
	Domain domain(0, n - 1);

	// Declare constrained domains for symmetry breaking
	Domain firstRow[n];
	Domain firstCol[n];

	// Declare variables x_ij forall i,j in {0, ..., n - 1}
	IntVar x[n][n];
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
				x[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("x_" + to_string(i) + to_string(j));
			}
			else if (j == 0 && i > 0) { // First column
				x[i][j] = cp_model.NewIntVar(firstCol[i]).WithName("x_" + to_string(i) + to_string(j));
			}
			else { // Generic
				x[i][j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
			}
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

	// Tell model how to count solutions
	Model model;

	if(argc >= 2) {
		SatParameters param;
		int seed = stoi(argv[1]);
		param.set_random_seed(seed);
		param.set_randomize_search(true);
		model.Add(NewSatParameters(param));
		cout << "Using random seed " << seed << endl;
	}

#if APPA_BRANCHING == 1
	vector <IntVar> VARS;
	for(int i=0; i<n; i++) {
		for(int j=0; j<n; j++) {
			VARS.push_back(x[i][j]);
		}
	}
	cp_model.AddDecisionStrategy(VARS, DecisionStrategyProto::CHOOSE_FIRST, DecisionStrategyProto::SELECT_MIN_VALUE);
#endif

	int num_solutions = 0;
	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
		cout << "\n";
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << SolutionIntegerValue(r, x[i][j]) << " ";
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
