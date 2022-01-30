GUROBIDIR = /opt/gurobi903/linux64
2mols_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$*
2mols_cuts_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DCALLBACK_CUTS
2mols_randobj_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DRANDOM_OBJECTIVE
2mols_nosym_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DNOSYM
2mols_fixcol_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXCOL
2mols_fixcol_cuts_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DCALLBACK_CUTS -DFIXCOL
2mols_fixcol_randobj_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DRANDOM_OBJECTIVE -DFIXCOL
2mols_nofixcol_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DNOFIXCOL
2mols_nofixcol_cuts_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DCALLBACK_CUTS -DNOFIXCOL
2mols_nofixcol_randobj_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DRANDOM_OBJECTIVE -DNOFIXCOL
2mols_backcirc_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DBACKCIRC
2mols_backcirc_cuts_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DCALLBACK_CUTS -DBACKCIRC
2mols_backcirc_randobj_%: 2mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DRANDOM_OBJECTIVE -DBACKCIRC
2.5mols_%: 2.5mols.cpp
	g++ -m64 -O4 -g -o $@ 2.5mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$*
2.5mols_fixcol_%: 2.5mols.cpp
	g++ -m64 -O4 -g -o $@ 2.5mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXCOL
2.5mols_cuts_%: 2.5mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2.5mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DCALLBACK_CUTS
2.5mols_fixcol_cuts_%: 2.5mols.cpp callback.cpp
	g++ -m64 -O4 -g -o $@ 2.5mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DCALLBACK_CUTS -DFIXCOL
3mols_%: 3mols.cpp
	g++ -m64 -O4 -g -o $@ 3mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$*
3mols_fixcol_%: 3mols.cpp
	g++ -m64 -O4 -g -o $@ 3mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXCOL
2mols_fixdiag1_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXDIAG1
2mols_fixdiag2_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXDIAG2
2mols_fixantidiag_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DFIXANTIDIAG
2mols_hints_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DHINTS
2mols_fixcol_hints_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DHINTS -DFIXCOL
2mols_nofixcol_hints_%: 2mols.cpp
	g++ -m64 -O4 -g -o $@ 2mols.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm -DN=$* -DHINTS -DNOFIXCOL
2MOLS_SLACK: 2MOLS_SLACK.cpp
	g++ -m64 -O4 -g -o $@ 2MOLS_SLACK.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm
3MOLS_SLACK: 3MOLS_SLACK.cpp
	g++ -m64 -O4 -g -o $@ 3MOLS_SLACK.cpp -I$(GUROBIDIR)/include/ -L$(GUROBIDIR)/lib -lgurobi_c++ -lgurobi90 -lm
