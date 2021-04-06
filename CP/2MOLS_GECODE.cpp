/*
	Noah Rubin
	2020-03-30

	This program implements the 2MOLS(n) CP model using GECODE.

	Notes:
		- Only provide path to ~/GECODE/lib and ~/GECODE/include (header automatically links)
		- Options are vague and not always honored by the search engine
		- Branching is specified manually, and we have a lot of control over the internals
		- The Space is copied many times during the search
		- Using LDS (Limited Discrepency Search) speeds up the time massively for n = 7
*/



#include "gecode/driver.hh"
#include "gecode/int.hh"
#include "gecode/search.hh"
#include <iostream>
#include <math.h>
#include <vector>
#include <chrono>
using namespace Gecode;
using namespace std;


int copyCount = 0;

// Override the default standard Tracer for integer variables
// We want to tell the user which variable is being pruned (stored as a char)
class StdIntTracer_VariableIDs : public StdIntTracer {
protected:
	char Type;
public:
	StdIntTracer_VariableIDs() {
		Type = 'N';
	}
	StdIntTracer_VariableIDs(char S) : Type(S) {}
	void prune(const Space& home, const IntTraceRecorder& trace, const ViewTraceInfo& vti, int i, IntTraceDelta& delta) {
		os << "trace<Int>::prune(id:" << trace.id();
		if (trace.group().in())
			os << ",g:" << trace.group().id();

		// Output variable current value (after pruning)
		os << "): " << Type << "[" << i << "] = " << trace[i] << " - {";

		// Output minimum prune value
		os << delta.min();

		// Output maximum prune value
		if (delta.width() > 1)
			os << ".." << delta.max();

		// Iterate over pruned values (delta) and output them
		++delta;
		while (delta()) {
			os << ',' << delta.min();
			if (delta.width() > 1)
				os << ".." << delta.max();
			++delta;
		}
		os << "} by " << vti << std::endl;
	}

	// Declare static tracers for X, Y, Z variables
	static StdIntTracer_VariableIDs Sx;
	static StdIntTracer_VariableIDs Sy;
	static StdIntTracer_VariableIDs Sz;
};

// Provide variable names (to be output in the prune() call, helps us figure out what is happening)
StdIntTracer_VariableIDs StdIntTracer_VariableIDs::Sx('X');
StdIntTracer_VariableIDs StdIntTracer_VariableIDs::Sy('Y');
StdIntTracer_VariableIDs StdIntTracer_VariableIDs::Sz('Z');


// We define the model as a sub-space
class TWO_MOLS : public Space {

// Store dimension and variables
protected:
	int n;
	IntVarArray x, y, z;

// Provide constructor
public:
	TWO_MOLS(int DIMENSION) : x(*this, (int)pow(DIMENSION, 2), 0, DIMENSION - 1), y(*this, (int)pow(DIMENSION, 2), 0, DIMENSION - 1), z(*this, (int)pow(DIMENSION, 2), 0, DIMENSION - 1) {
		n = DIMENSION;
		cout << "Calling construction\n";

		// Declare domains for lexocographic ordering
		vector<int> domain_lex(n);
		for (int i = 0; i < n; i++) {
			domain_lex[i] = i;
		}
		
		// Declare domains for Dominance Detection
		vector<vector<int>> domain_constrained(n);

		// Fix entry (0,1) to 2 in L2 
		domain_constrained[1].push_back(2);
		
		// Fix remaining constrained domains for column 1 in L2
		for (int i = 2; i < n; i++) {
			vector<int> domainTemp;
			if (i == (n - 2)) {
				domainTemp.push_back((n - 1));
				domain_constrained[i] = domainTemp;
				continue;
			}
			if (i == (n - 1)) {
				domainTemp.push_back(1);
				for (int j = 3; j <= (n - 2); j++) {
					domainTemp.push_back(j);
				}
				domain_constrained[i] = domainTemp;
				break;
			}
			for (int j = 1; j <= i + 1; j++) {
				if (j != i && j != 2) {
					domainTemp.push_back(j);
				}
			}
			domain_constrained[i] = domainTemp;
		}
		

		// Apply domain constraints to the first rows and columns
		for (int i = 0; i < n; i++) {
			dom(*this, x[i], domain_lex[i]);
			dom(*this, y[i], domain_lex[i]);
			dom(*this, x[i * n], domain_lex[i]);
			if (i > 0) {
				dom(*this, y[i * n], IntSet(IntArgs(domain_constrained[i])));
			}
		}

		// Latin property in rows of X
		for (int i = 0; i < n; i++) {
			vector<IntVar> rows_x;
			for (int j = 0; j < n; j++) {
				rows_x.push_back(x[i * n + j]);
			}
			distinct(*this, rows_x);
		}

		// Latin property in cols of X
		for (int i = 0; i < n; i++) {
			vector<IntVar> cols_x;
			for (int j = 0; j < n; j++) {
				cols_x.push_back(x[j * n + i]);
			}
			distinct(*this, cols_x);
		}

		// Latin property in rows of Y
		for (int i = 0; i < n; i++) {
			vector<IntVar> rows_y;
			for (int j = 0; j < n; j++) {
				rows_y.push_back(y[i * n + j]);
			}
			distinct(*this, rows_y);
		}

		// Latin property in cols of Y
		for (int i = 0; i < n; i++) {
			vector<IntVar> cols_y;
			for (int j = 0; j < n; j++) {
				cols_y.push_back(y[j * n + i]);
			}
			distinct(*this, cols_y);
		}

		// Row uniqueness of Z
		for (int i = 0; i < n; i++) {
			vector<IntVar> rows_z;
			for (int j = 0; j < n; j++) {
				rows_z.push_back(z[i * n + j]);
			}
			distinct(*this, rows_z);
		}

		// Column uniqueness of Z
		for (int i = 0; i < n; i++) {
			vector<IntVar> cols_z;
			for (int j = 0; j < n; j++) {
				cols_z.push_back(z[j * n + i]);
			}
			distinct(*this, cols_z);
		}

		// Enforce element constraints on Z, X, Y
		for (int i = 0; i < n; i++)
		{
			vector<IntVar> Zvec_i;
			for (int j = 0; j < n; j++) {
				Zvec_i.push_back(z[i * n + j]);

			}
			for (int j = 0; j < n; j++) {
				element(*this, IntVarArgs(Zvec_i), x[i * n + j], y[i * n + j]);
			}
		}
		

		// Branch first on Xs, then Ys, then Zs
		
		// Branch strategy: select variable w/ smallest domain size --> select its minimum value
		//
		//								|										|
		//								v										v
		branch(*this, x,		INT_VAR_SIZE_MIN(),						INT_VAL_MIN());
		branch(*this, y,		INT_VAR_SIZE_MIN(),						INT_VAL_MIN());
		branch(*this, z,		INT_VAR_SIZE_MIN(),						INT_VAL_MIN());

		// Apply tracers to search to print out variable domain reductions
		trace(*this, x, TE_COMMIT | TE_INIT | TE_DONE | TE_FAIL | TE_POST | TE_PROPAGATE | TE_PRUNE, StdIntTracer_VariableIDs::Sx);
		trace(*this, y, TE_COMMIT | TE_INIT | TE_DONE | TE_FAIL | TE_POST | TE_PROPAGATE | TE_PRUNE, StdIntTracer_VariableIDs::Sy);
		trace(*this, z, TE_COMMIT | TE_INIT | TE_DONE | TE_FAIL | TE_POST | TE_PROPAGATE | TE_PRUNE, StdIntTracer_VariableIDs::Sz);

	}


	// Provide copy constructor for committing
	TWO_MOLS(TWO_MOLS& T) : Space(T){
		x.update(*this, T.x);
		y.update(*this, T.y);
		this->n = T.n;
	}


	// Provide instruction for copy within search
	virtual Space* copy(void) {
		copyCount++;
		if (this->stable()) {
			const Choice* ch = this->choice();
			this->print(*ch, 0, std::cout);
			cout << "\n";
		}
		return new TWO_MOLS(*this);
	}

	// Provide function to print solutions
	void Print() {
		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				std::cout << x[i * n + j] << " ";
			}
			std::cout << "\n";
		}
		
		std::cout << "\n";

		for (int i = 0; i < n; i++) {
			for (int j = 0; j < n; j++) {
				std::cout << y[i * n + j] << " ";
			}
			std::cout << "\n";
		}
		
	}

	static void printVars(const Space& home, const Brancher& b, unsigned int a, IntVar, int i, const int& n, ostream & o) {
		const TWO_MOLS& T = static_cast<const TWO_MOLS&>(home);
		cout << static_cast<char>(n);
	}

};


// Driver
int main(int argc, char *argv[]) {
	// Parse dimension
	int n = 0;
	if (argc == 2) {
		n = atoi(argv[1]);
	}
	else {
		cout << "Invalid argument for n";
		exit(0);
	}

	// Declare new search space
	TWO_MOLS *T = new TWO_MOLS(n);

	// Declare options - the Tracer will tell us when the search branches and on what value
	Search::Options *O = new Search::Options();
	O->threads = 1;
	O->tracer = &StdSearchTracer::def;

	auto tic = chrono::high_resolution_clock::now();
	
	// Execute Depth-First Search on the space
	DFS<TWO_MOLS> e(T, *O);

	// Free up old space
	delete T;

	// Grab first solution from DFS or else we have no solutions
	if (TWO_MOLS * S = e.next()) {
		S->Print();
	}
	else {
		cout << "No solutions found\n";
	}

	auto toc = chrono::high_resolution_clock::now();

	cout << "Total time (s) = " << chrono::duration_cast<chrono::nanoseconds>(toc - tic).count() / 1000000000.0 << "\n";
	cout << "Called copy " << copyCount << " times\n";

	return 0;
}