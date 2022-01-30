/**
	Noah Rubin
	2020-09-15
	This program uses constraint programming to generate partially filled (upto t transversals) MOLS(1)
	- Added symmetry breaking to improve runtime

	Model:

	(i) Let x_ij denote the column index of value j in row i

	(ii) Variables x_ij for i = 0,...,n-1 and j = 0,...,t-1

	(iii) All_Different{x_ij for i = 0,...,n-1} for j = 0,...,t-1 -> Uniqueness of entries in each columns
	(iv)  All_Different{x_ij for j = 0,...,t-1} for i = 0,...,n-1 -> Entries in same row cannot share same column
**/

#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
using namespace std;
using namespace operations_research;
using namespace sat;

#define n ORDER
#define t (n / 2)

#define PRINT 1
#define ENEUMERATE_ALL 0

int main(int argc, char* argv[]) {
	// Model
	CpModelBuilder cp_model;

	// 

	/**
		Fix domains for general entries, first row and x_12 to obey structure

		0 1 2 ... t ...
		2 ...
		:
	   n-1
		:
	**/

	Domain domain(0, n - 1);
	Domain domain_x_12 = Domain::FromValues({ 0 });
	Domain domain_x_nMinus2_nMinus1 = Domain::FromValues({ 0 });
	Domain firstRow[t];

	// Declare variables x_ij
	IntVar x[n][t];
	
	//Counters
	int i = 0, j = 0;

	for (i = 0; i < t; i++) {
		firstRow[i] = Domain::FromValues({ i });
	}

	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < t; j++) {
			// x_0k for k = 0,...,t
			if (i == 0) {
				x[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("x_" + to_string(i) + to_string(j));
			}
			// x_12
			else if (i == 1 && j == 2) {
				x[i][j] = cp_model.NewIntVar(domain_x_12).WithName("x_" + to_string(i) + to_string(j));
			}
			// Won't happen unless t >= n - 2
			else if (i == n - 1 && j == n - 2) { 
				x[i][j] = cp_model.NewIntVar(domain_x_nMinus2_nMinus1).WithName("x_" + to_string(i) + to_string(j));
			}
			// General case
			else {
				x[i][j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
			}
		}
	}

	// Column uniqueness constraints
	for (j = 0; j < t; j++) {
		vector<IntVar> col_i;
		for (i = 0; i < n; i++) {
			col_i.push_back(x[i][j]);
		}
		cp_model.AddAllDifferent(col_i);
	}
	
	// Enforce uniqueness for individual entries
	// ie. for n = 7, x_11 != x_12 != x_13 since they cannot share the same column AND row
	for (i = 0; i < n; i++) {
		vector<IntVar> entry_ij;
		for (j = 0; j < t; j++) {
			entry_ij.push_back(x[i][j]);
		}
		cp_model.AddAllDifferent(entry_ij);
	}

	// Matrix to hold partial Latin square
	int X[n][n];

	for(int i=0; i<n; i++)
		for(int j=0; j<n; j++)
			X[i][j] = -1;

	// Tell model how to observe and print solutions
	Model model;
	int num_solutions = 0;
	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
		if (PRINT) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < t; j++) {
					cout << "x_" + to_string(i) + to_string(j) << " = " << SolutionIntegerValue(r, x[i][j]) << " ";
					X[i][SolutionIntegerValue(r, x[i][j])] = j;
				}
				cout << "\n";
			}
		}
		num_solutions++;
	}));

	// Params

	SatParameters parameters;
	if (ENEUMERATE_ALL) {
		parameters.set_enumerate_all_solutions(true);
	}
	model.Add(NewSatParameters(parameters));

	// Execute
	auto tic = chrono::high_resolution_clock::now();
	const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model);
	auto toc = chrono::high_resolution_clock::now();

	// Print partially completed Latin square
	for(int i=0; i<n; i++)
	{	
		for(int j=0; j<n; j++)
			cout << (X[i][j] == -1 ? "." : to_string(X[i][j])) << " ";
		cout << "\n";
	}

	// Report time
	cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms\n";
		
	// Fin
	return EXIT_SUCCESS;
}
