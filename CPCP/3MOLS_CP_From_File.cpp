/**
	Curtis Bright | Noah Rubin
	2021-08-14
	Solves 3MOLS(n) using Curtis' CPCP model, and reads squares as input from a file.

	File naming convention is "1MOLS_n=<n>.txt".


**/
#include "ortools/sat/cp_model.h"
#include "ortools/sat/model.h"
#include "ortools/sat/sat_parameters.pb.h"
#include <vector>
#include <tuple>
#include <string>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <fstream>
#include <assert.h>
using namespace std;
using namespace operations_research;
using namespace sat;

#define VERBOSE 0
#define STOP_ON_SOL 1
#define ENUMERATE_ALL 1
#define RANDOMIZE_SEARCH 1
#define EXTRA_SQUARE 0
#define TIMEOUT 60000

#if NO_DUPLICATE_CYCLE_TYPES == 1 && !defined(NOFIXCOL) && !defined(FIXCOL)
#include "cycle_types.h"
#endif
chrono::time_point<chrono::high_resolution_clock> tic, toc, tic2, toc2;
int main(int argc, char* argv[]) {
	int i, j, num_solutions, passed_list;
	i = j = num_solutions = passed_list = 0;

	int n = stoi(argv[1]);
	vector<int64_t> firstCol_passed; // Expect possibility of column fixing
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

	string colFixing = "";
	if (passed_list) {
		for (i = 0; i < firstCol_passed.size() - 1; i++) {
			colFixing += to_string(firstCol_passed[i]) + ",";
		}
		colFixing += to_string(firstCol_passed.back());
	}

	ifstream input_squares;
	input_squares.open("C:\\Users\\Noah\\Desktop\\1MOLS_n=" + to_string(n) + (passed_list ? "_fixing=[" + colFixing + "]" : "") + ".txt");
	if (!input_squares.is_open()) {
		std::cout << "Failed to open square file.\n";
		exit(-1);
	}
	string line;
	string individual_entry;
	istringstream lineStream;
	int* square = new int[n * n]();
	Domain default_domain(0, n - 1);

	tic = chrono::high_resolution_clock::now();

	vector<int> lineValues;
	int lineNum = 1;
	while (getline(input_squares, line)) {
		if (!(lineNum % (n + 1))) {
			num_solutions++;
			if (chrono::duration_cast<std::chrono::milliseconds>(chrono::high_resolution_clock::now() - tic).count() > TIMEOUT * 1000) {
				std::cout << "Completed " << to_string(num_solutions - 1) << " trials" << endl;
				std::cout << "Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms" << endl;
				exit(0);
			}

			int *L3 = new int[n * n]();
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					L3[i * n + j] = square[i * n + j];
				}
			}			
			
#if VERBOSE > 1
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					std::cout << L3[i * n + j] << " ";
				}
				std::cout << endl;
			}
#endif
			// Model
			CpModelBuilder cp_model;
			// Declare constrained domains for symmetry breaking
			Domain *firstRow = new Domain[n];
			Domain *firstCol_L1 = new Domain[n];
			Domain *firstCol_L2 = new Domain[n];
			// Declare variables x_ij forall i,j in {0, ..., n - 1}
			IntVar* x = new IntVar[n * n], * y = new IntVar[n * n], * z = new IntVar[n * n],
				* a = new IntVar[n * n], * b = new IntVar[n * n], * c = new IntVar[n * n], * d = new IntVar[n * n];

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
			for (i = 1; i < n - 1; i++) {
				firstCol_L2[i] = Domain::FromValues({ i + 1 });
			}
			firstCol_L2[n - 1] = Domain::FromValues({ 1 });
#elif defined(NOFIXCOL) || (defined(FIXCOLTHIRD) && !defined(DOMAINRED))
			for (i = 1; i < n; i++) {
				firstCol_L2[i] = domain;
			}
#else
			// Fix entry (0,1) to 2 in L2
			firstCol_L2[1] = Domain::FromValues({ 2 });
			// Fix entry (0,N-2) to N-1 in L2
			firstCol_L2[n - 2] = Domain::FromValues({ n - 1 });
			// Fix remaining constrained domains for column 1 in L2
			for (i = 2; i < n; i++) {
				vector<int64_t> domainTemp;
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
						x[i * n + j] = cp_model.NewIntVar(firstRow[j]).WithName("x_" + to_string(i) + to_string(j));
						y[i * n + j] = cp_model.NewIntVar(firstRow[j]).WithName("y_" + to_string(i) + to_string(j));
					}
#ifndef NOFIXCOL
					else if (j == 0 && i > 0) { // First column
						x[i * n + j] = cp_model.NewIntVar(firstCol_L1[i]).WithName("x_" + to_string(i) + to_string(j));
						y[i * n + j] = cp_model.NewIntVar(firstCol_L2[i]).WithName("y_" + to_string(i) + to_string(j));
					}
#endif
					else { // Generic
						x[i * n + j] = cp_model.NewIntVar(default_domain).WithName("x_" + to_string(i) + to_string(j));
						y[i * n + j] = cp_model.NewIntVar(default_domain).WithName("y_" + to_string(i) + to_string(j));
					}
					z[i * n + j] = cp_model.NewIntVar(Domain::FromValues({ L3[i * n + j] })).WithName("z_" + to_string(i) + to_string(j));
				}
			}

			// Declare auxiliary variables to enforce orthogonality
			for (i = 0; i < n; i++) {
				for (j = 0; j < n; j++) {
					a[i * n + j] = cp_model.NewIntVar(default_domain).WithName("a_" + to_string(i) + to_string(j));
					b[i * n + j] = cp_model.NewIntVar(default_domain).WithName("b_" + to_string(i) + to_string(j));
					c[i * n + j] = cp_model.NewIntVar(default_domain).WithName("c_" + to_string(i) + to_string(j));
#if EXTRA_SQUARE == 1
					d[i * n + j] = cp_model.NewIntVar(default_domain).WithName("d_" + to_string(i) + to_string(j));
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
			// Column uniqueness constraints L3
			for (i = 0; i < n; i++) {
				vector<IntVar> col_i;
				for (j = 0; j < n; j++) {
					col_i.push_back(z[i * n + j]);
				}
				cp_model.AddAllDifferent(col_i);
			}
			// Row uniqueness constraints L3
			for (i = 0; i < n; i++) {
				vector<IntVar> row_i;
				for (j = 0; j < n; j++) {
					row_i.push_back(z[j * n + i]);
				}
				cp_model.AddAllDifferent(row_i);
			}
			// Column uniqueness constraints on A
			for (i = 0; i < n; i++) {
				vector<IntVar> col_i;
				for (j = 0; j < n; j++) {
					col_i.push_back(a[i * n + j]);
				}
				cp_model.AddAllDifferent(col_i);
			}
			// Row uniqueness constraints on A
			for (i = 0; i < n; i++) {
				vector<IntVar> row_i;
				for (j = 0; j < n; j++) {
					row_i.push_back(a[j * n + i]);
				}
				cp_model.AddAllDifferent(row_i);
			}
			// Column uniqueness constraints on B
			for (i = 0; i < n; i++) {
				vector<IntVar> col_i;
				for (j = 0; j < n; j++) {
					col_i.push_back(b[i * n + j]);
				}
				cp_model.AddAllDifferent(col_i);
			}
			// Row uniqueness constraints on B
			for (i = 0; i < n; i++) {
				vector<IntVar> row_i;
				for (j = 0; j < n; j++) {
					row_i.push_back(b[j * n + i]);
				}
				cp_model.AddAllDifferent(row_i);
			}
			// Column uniqueness constraints on C
			for (i = 0; i < n; i++) {
				vector<IntVar> col_i;
				for (j = 0; j < n; j++) {
					col_i.push_back(c[i * n + j]);
				}
				cp_model.AddAllDifferent(col_i);
			}
			// Row uniqueness constraints on C
			for (i = 0; i < n; i++) {
				vector<IntVar> row_i;
				for (j = 0; j < n; j++) {
					row_i.push_back(c[j * n + i]);
				}
				cp_model.AddAllDifferent(row_i);
			}
#if EXTRA_SQUARE == 1
			// Column uniqueness constraints on D
			for (i = 0; i < n; i++) {
				vector<IntVar> col_i;
				for (j = 0; j < n; j++) {
					col_i.push_back(d[i * n + j]);
				}
				cp_model.AddAllDifferent(col_i);
			}
			// Row uniqueness constraints on D
			for (i = 0; i < n; i++) {
				vector<IntVar> row_i;
				for (j = 0; j < n; j++) {
					row_i.push_back(d[j * n + i]);
				}
				cp_model.AddAllDifferent(row_i);
			}
#endif
			// Constraints defining A variables
			for (int i = 0; i < n; i++)
			{
				vector<IntVar> Avec_i;
				for (int j = 0; j < n; j++)
					Avec_i.push_back(a[i * n + j]);
				for (int j = 0; j < n; j++)
					cp_model.AddVariableElement(x[i * n + j], Avec_i, y[i * n + j]);
			}
			// Constraints defining B variables
			for (int i = 0; i < n; i++)
			{
				vector<IntVar> Bvec_i;
				for (int j = 0; j < n; j++)
					Bvec_i.push_back(b[i * n + j]);
				for (int j = 0; j < n; j++)
					cp_model.AddVariableElement(x[i * n + j], Bvec_i, z[i * n + j]);
			}
			// Constraints defining C variables
			for (int i = 0; i < n; i++)
			{
				vector<IntVar> Cvec_i;
				for (int j = 0; j < n; j++)
					Cvec_i.push_back(c[i * n + j]);
				for (int j = 0; j < n; j++)
					cp_model.AddVariableElement(y[i * n + j], Cvec_i, z[i * n + j]);
			}
#if EXTRA_SQUARE == 1
			// Constraints defining D variables
			for (int i = 0; i < n; i++)
			{
				vector<IntVar> Dvec_i;
				for (int j = 0; j < n; j++)
					Dvec_i.push_back(d[j * n + i]);
				for (int j = 0; j < n; j++)
					cp_model.AddVariableElement(a[j * n + i], Dvec_i, b[j * n + i]);
			}
#endif
#if NO_DUPLICATE_CYCLE_TYPES == 1 && !defined(NOFIXCOL) && !defined(FIXCOL)
			BoolVar by[n][n];
			for (int i = 2; i < n; i++) {
				if (i != n - 2) {
					for (int j = 1; j < n; j++) {
						by[i * n + j] = cp_model.NewBoolVar();
						cp_model.AddNotEqual(y[i][0], j).OnlyEnforceIf(Not(by[i * n + j]));
					}
				}
			}
			for (int i = 0; i < blocked_cycle_types.size(); i++) {
				vector<BoolVar> clause;
				for (int j = 2; j < n; j++) {
					if (j != n - 2) {
						int k = blocked_cycle_types[i * n + j];
						clause.push_back(Not(by[j][k]));
					}
				}
				cp_model.AddBoolOr(clause);
			}
#endif
			Model model_new;
			model_new.Add(NewFeasibleSolutionObserver([&](const CpSolverResponse& r) {
				std::cout << "Solution found in trial # " << num_solutions << endl;
				std::cout << "\n";
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						std::cout << SolutionIntegerValue(r, x[i * n + j]) << " ";
					}
					std::cout << "\n";
				}
				std::cout << "\n\n";
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						std::cout << SolutionIntegerValue(r, y[i * n + j]) << " ";
					}
					std::cout << "\n";
				}
				std::cout << "\n\n";
				for (i = 0; i < n; i++) {
					for (j = 0; j < n; j++) {
						std::cout << SolutionIntegerValue(r, z[i * n + j]) << " ";
					}
					std::cout << "\n";
				}
				std::cout << "\n";
				toc = chrono::high_resolution_clock::now();
				if (STOP_ON_SOL) {
					std::cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms" << endl;
					exit(0);
				}
				if (ENUMERATE_ALL) {
					// Report time
					std::cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms" << endl;
				}
				}));
			// Execute model
			tic2 = chrono::high_resolution_clock::now();
			const CpSolverResponse response = SolveCpModel(cp_model.Build(), &model_new);
			toc2 = chrono::high_resolution_clock::now();
#if VERBOSE > 1
			std::cout << "No solutions" << endl;
#endif
			toc = chrono::high_resolution_clock::now();
#if VERBOSE == 0
			if (num_solutions % 100 == 0)
#endif
				std::cout << fixed << "Trial # " << num_solutions << " Time: " << chrono::duration_cast<std::chrono::milliseconds>(toc2 - tic2).count() << " ms" << " Total time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << " ms Avg. time: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() / ((double)1000 * num_solutions) << " s" << endl;
}
		else {
			lineStream.str(line);
			while (getline(lineStream, individual_entry, ' ')) {
				lineValues.push_back(stoi(individual_entry));
			}
			lineStream.clear();
			for (int i = 0; i < n; i++) {
				square[((lineNum - 1) % (n + 1)) * n + i] = lineValues[i];
			}
			lineValues.clear();
		}

		lineNum++;
	}

	// Report time
	std::cout << "Time elapsed: " << chrono::duration_cast<std::chrono::milliseconds>(toc - tic).count() << "ms" << endl;
	// Fin
	return EXIT_SUCCESS;
}