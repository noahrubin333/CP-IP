// Noah Rubin, updated by Curtis Bright
// 2020-12-15, updated 2021-07-20

using namespace std;

#include <iostream>
#include <vector>
#include <algorithm>
#include <sstream>
#include <string>

// Forward declaration
std::vector<int> categorize(vector<int>, int);

// Driver
int main(int argc, char** argv) {

	if(argc == 1) {
		cout << "Need to provide permutation on the command-line; e.g., [0,2,3,4,5,6,7,1]" << endl;
		exit(0);
	}

	vector<int> perm;
	string list = string(argv[1]);
	istringstream argvStream(list.substr(1, list.size() - 2));
	string token;
	while (getline(argvStream, token, ',')) {
		perm.push_back(stoi(token));
	}

	vector<int> cycleTypes = categorize(perm, perm.size());

	for(int i = 0; i < cycleTypes.size(); i++) {
		cout << cycleTypes[i];
		if(i < cycleTypes.size()-1)
			cout << ",";
	}
	cout << endl;
}


/**
	This function will categorize a permutation s into a list of integers k1, k2, ..., kp
	where ki is the length of the ith cycle of s for i = 1,...,p

	@param  s - an array of integers which represents a permutation rule on the numbers 0, 1, ..., n - 1
		len - the length of the permutation n

	@return vector<int> - an ordered vector containing k1, k2, ..., kp

**/
std::vector<int> categorize(vector<int> s, int len) {
	
	// Let s be a permutation on (0, 1, ..., len - 1) defined as 
	//	s = ( 0    1	 ...  len - 1  )
	//	    (s[0] s[1]		s[len - 1] )

	// Holds {k_i} for all i cycles in sigma (s)
	vector<int> cycleTypes;

	// Binary array, visited[i] = 1 if we have visited i, else 0
	int *visited = new int[len];

	// We define a cycle length as k such that s^(k + 1)(x) = x
	// ie. x -> s(x) -> s(s(x)) -> ... -> s^k(x) -> x
	int i = 0, j = 0; // Counters
	int x = 0, s_of_x = 0; // Holds x and s(x), which are cycled as defined above
	int count = 0; // Holds k as defined above
	int start = 0; // Holds x to check when we cycle back to first element

	// Innitialize visited array to all unvisited
	for (i = 0; i < len; i++) {
		visited[i] = 0;
	}

	i = 0;

	// As long as we haven't checked all elements in the permutation
	while (j != len) {
		start = x = i; // Start current cycle at x
		s_of_x = s[i]; // Calculate s(x)
		count = 1; // Start current cycle length at 1
		visited[x] = 1; // Visit x

		// Start Depth First Search
		while (s_of_x != start) { // While we have not completed the cycle
			visited[s_of_x] = 1;	// Visit s(x)
			x = s_of_x;				// let x <- s(x)
			s_of_x = s[x];			// let s(x) <- s^(count + 1)(x), ie. visit next element
			count += 1;				// Increment count
		}

		// Finished DFS, push the current cycle length to vector
		cycleTypes.push_back(count);

		// Reset count
		count = 1;

		// Scan for next unchecked element
		for (j = 0; j < len; j++) {
			// If we have not yet visited an element and are still in the permutation
			if (visited[j] == 0) {
				i = j; // Assign next start to the element we have yet to check and break
				break;
			}
		}
	}

	// Sort and return entire vector from least to greatest length
	std::sort(cycleTypes.begin(), cycleTypes.end());
	return cycleTypes;
}
