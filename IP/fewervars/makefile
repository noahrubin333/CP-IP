2mols_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$*
2mols_fixcol_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXCOL
2mols_nonewvars_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DNONEWVARS
2mols_fixcol_nonewvars_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXCOL -DNONEWVARS
2mols_extravars_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DEXTRA_VARS
2mols_fixcol_extravars_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -I/opt/gurobi903/linux64/include/ -L/opt/gurobi903/linux64/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXCOL -DEXTRA_VARS
