# Searching for Mutually Orthogonal Latin Squares using Constraint and Integer Programming
Implementations by Noah Rubin; updates and Linux driver scripts by Curtis Bright; based on work with Kevin Cheung and Brett Stevens.

This repository contains the implementations and logs for programs which generate sets of Mutually Orthogonal Latin Squares.

Each directory contains a script "runtime-test.sh" that accepts the parameters of the search and then performs the following actions in sequence:

* Compiles the implementation if necessary
* Runs the compiled executable
* Saves the output to a log file
* Verifies that any solution found is in fact a k-MOLS(n)

Models can be modified using command-line flags; run "runtime-test.sh" without any parameters to see the options that are available.

The script "summary-table.sh" summarizes the running times for a given k in a table.

CP - Pure constraint programming model using OR-Tools
CPCP - Constraint programming model using OR-Tools and an additional call to OR-Tools to fill in one square
CPI - Hybrid constraint/integer programming model using OR-Tools and Gurobi
IP - Pure integer programming model using Gurobi

For running in Windows using Microsoft Visual C++ ensure the following are set in project properties:
* Additional Include Directories (C/C++): 	<OR_TOOLS_PATH>\include;<GUROBI_PATH>\win64\include
* Additional Library Directories (Linker): 	<OR_TOOLS_PATH>\lib;<GUROBI_PATH>\win64\lib
* Additional Dependencies (Linker): 		gurobi_c++md2019.lib;gurobi.lib;ortools.lib
