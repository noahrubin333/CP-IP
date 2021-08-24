#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include <iostream>
#include <fstream>
using namespace std;
using namespace operations_research;
using namespace sat;
#include <vector>
#include "Cycle_Types.h"
#define PRINT 1
#define ENEUMERATE_ALL 1
int main(int argc, char* argv[]) {

	int passed_list = 0; // Flag set if column fixing is passed at runtime

	vector<int64> firstCol_passed; // Expect possibility of column fixing

	int i = 0, j = 0, k = 0, l = 0; // General purpose counters

	int n = stoi(argv[1]);

	// Allow passing of first column of Y as a list [Y_00, Y_10, ..., Y_(n-1)0]
	if (argc > 2 && strlen(argv[2]) > 1) {
		passed_list = 1;
		string list = string(argv[2]);
		istringstream argvStream(list.substr(1, list.size() - 2));
		string token;
		int count = 0;
		while (getline(argvStream, token, ',')) {
			firstCol_passed.push_back(stoi(token));
			count++;
		}
		if (n != count) {
			std::cout << "n = " << n << " but you entered a list of length " << count << "\n";
			exit(0);
		}
	}

	// Model
	CpModelBuilder cp_model;

	// Fix general domains for x_ij, y_ij, z_ij
	Domain domain(0, n - 1);
	// Declare variables L_ij for CP Model L3
	IntVar *L = new IntVar[n * n];

	// Instantiate variables and provide their associated domains
	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			if (i == 0) {
				L[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ j })).WithName("l_" + to_string(i) + to_string(j));
			}
			else if ((j == 0 && i > 0) && passed_list) {
				// First col
				Domain d = Domain::FromValues({0, firstCol_passed[i], i});
				L[i * n + j] = cp_model.NewIntVar(domain.IntersectionWith(d.Complement())).WithName("l_" + to_string(i) + to_string(j));
			}
			// General case
			else {
				L[i * n + j] = cp_model.NewIntVar(domain).WithName("l_" + to_string(i) + to_string(j));
			}

		}
	}
	// Column uniqueness constraints
	for (j = 0; j < n; j++) {
		vector<IntVar> col_j;
		for (i = 0; i < n; i++) {
			col_j.push_back(L[i * n + j]);
		}
		cp_model.AddAllDifferent(col_j);
	}
	// Row uniqueness
	for (i = 0; i < n; i++) {
		vector<IntVar> row_i;
		for (j = 0; j < n; j++) {
			row_i.push_back(L[i * n + j]);
		}
		cp_model.AddAllDifferent(row_i);
	}

	Model model;
	int num_solutions = 0;

	string colFixing = "";
	if (passed_list) {	
		for (i = 0; i < firstCol_passed.size() - 1; i++) {
			colFixing += to_string(firstCol_passed[i]) + ",";
		}
		colFixing += to_string(firstCol_passed.back());
	}
	ofstream out("C:\\Users\\Noah\\Desktop\\1MOLS_n=" + to_string(n) + (passed_list ? "_fixing=[" + colFixing + "]" : "") + ".txt");
	
	model.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
		if (PRINT) {
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					out << SolutionIntegerValue(r, L[i * n + j]) << " ";
				}
				out << "\n";
			}
			out << "===============\n";
		}
		num_solutions++;
		if (num_solutions == 1000000) { exit(0); }
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
	// Report time
	std::cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms\n";
	std::cout << num_solutions << " solutions";
	// Fin
	out.close();
	return EXIT_SUCCESS;
}