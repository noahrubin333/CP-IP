#ORDIR = /home/cbright/or-tools_Debian-10-64bit_v7.8.7959/
ORDIR = /home/cbright/or-tools_Ubuntu-20.04-64bit_v7.8.7959
2mols_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$*
2mols_fixcol_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$* -DFIXCOL
2mols_nofixcol_%: 2mols.cpp
	g++ -m64 -g -o $@ 2mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$* -DNOFIXCOL
3mols_%: 3mols.cpp cycle_types.h
	g++ -m64 -g -o $@ 3mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$*
3mols_fixcol_%: 3mols.cpp cycle_types.h
	g++ -m64 -g -o $@ 3mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$* -DFIXCOL
3mols_fixcolthird_%: 3mols.cpp cycle_types.h
	g++ -m64 -g -o $@ 3mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$* -DFIXCOLTHIRD
3mols_domainred_%: 3mols.cpp cycle_types.h
	g++ -m64 -g -o $@ 3mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$* -DFIXCOLTHIRD -DDOMAINRED
3mols_cycletype_%: 3mols.cpp cycle_types.h
	g++ -m64 -g -o $@ 3mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$* -DFIXCOLTHIRD -DDOMAINRED -DCYCLETYPE
3mols_nofixcol_%: 3mols.cpp cycle_types.h
	g++ -m64 -g -o $@ 3mols.cpp -lm -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -fPIC -std=c++11 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP -O4 -DORDER=$* -DNOFIXCOL
