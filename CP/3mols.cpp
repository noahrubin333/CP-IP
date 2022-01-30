/**
	Noah Rubin
	2020-09-02
	Modification by Curtis Bright to generate triples of MOLS(n)
**/
#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
using namespace std;
using namespace operations_research;
using namespace sat;
#define n ORDER
int main(int argc, char* argv[]) {
	// Model
	CpModelBuilder cp_model;
	// Fix general domains for x_ij, y_ij, z_ij
	Domain domain(0, n - 1);
	// Declare constrained domains for symmetry breaking
	Domain firstRow[n];
	Domain firstCol_L1[n];
	Domain firstCol_L2[n];
	// Declare variables x_ij forall i,j in {0, ..., n - 1}
	IntVar x[n][n], y[n][n], z[n][n], a[n][n], b[n][n], c[n][n];
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
				z[i][j] = cp_model.NewIntVar(firstRow[j]).WithName("z_" + to_string(i) + to_string(j));
			}
			else if (j == 0 && i > 0) { // First column
				x[i][j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
				y[i][j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
				z[i][j] = cp_model.NewIntVar(domain).WithName("z_" + to_string(i) + to_string(j));
			}
			else { // Generic
				x[i][j] = cp_model.NewIntVar(domain).WithName("x_" + to_string(i) + to_string(j));
				y[i][j] = cp_model.NewIntVar(domain).WithName("y_" + to_string(i) + to_string(j));
				z[i][j] = cp_model.NewIntVar(domain).WithName("z_" + to_string(i) + to_string(j));
			}
			a[i][j] = cp_model.NewIntVar(domain).WithName("a_" + to_string(i) + to_string(j));
			b[i][j] = cp_model.NewIntVar(domain).WithName("b_" + to_string(i) + to_string(j));
			c[i][j] = cp_model.NewIntVar(domain).WithName("c_" + to_string(i) + to_string(j));
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
				
	// Tell model how to count solutions
	Model model;

	if(argc >= 2) {
		SatParameters param;
		int seed = stoi(argv[1]);
		param.set_random_seed(seed);
		param.set_randomize_search(true);
		model.Add(NewSatParameters(param));;
		cout << "Using random seed " << seed << endl;
	}

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
		cout << "\n\n";
		for (i = 0; i < n; i++) {
			for (j = 0; j < n; j++) {
				cout << SolutionIntegerValue(r, z[i][j]) << " ";
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
