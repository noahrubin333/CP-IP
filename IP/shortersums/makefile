2mols_%: 2mols.cpp
	g++ -m64 -g -o 2mols_$* 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$*
2mols_fixcol_%: 2mols.cpp
	g++ -m64 -g -o 2mols_fixcol_$* 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXCOL
