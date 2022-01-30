ORDIR = /home/cbright/or-tools_Ubuntu-20.04-64bit_v7.8.7959
0.5mols_%: 0.5mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 0.5mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$*
1mols_%: 1mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 1mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$*
2mols_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$*
3mols_%: 3mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 3mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$*
3mols_fixcol_%: 3mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 3mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DFIXCOL
2mols_fixcol_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DFIXCOL
2mols_nosym_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DNOSYM
2mols_nosym_index_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DNOSYM -DINDEX
2mols_fixdiag_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DFIXCOL -DFIXDIAG
2mols_fixantidiag_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DFIXCOL -DFIXANTIDIAG
2mols_bool_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL
2mols_fixcol_bool_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL -DFIXCOL
2mols_fixdiag_bool_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL -DFIXCOL -DFIXDIAG
2mols_fixantidiag_bool_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL -DFIXCOL -DFIXANTIDIAG
2mols_linearelem_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DLINEAR_ELEM
2mols_fixcol_linearelem_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DLINEAR_ELEM -DFIXCOL
2mols_fixdiag_linearelem_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DLINEAR_ELEM -DFIXCOL -DFIXDIAG
2mols_fixantidiag_linearelem_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DLINEAR_ELEM -DFIXCOL -DFIXANTIDIAG
2mols_extra_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL -DEXTRA
2mols_fixcol_extra_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL -DEXTRA -DFIXCOL
2mols_fixdiag_extra_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL -DEXTRA -DFIXCOL -DFIXDIAG
2mols_fixantidiag_extra_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DBOOL -DEXTRA -DFIXCOL -DFIXANTIDIAG
2mols_index_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DINDEX
2mols_fixcol_index_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DINDEX -DFIXCOL
2mols_fixdiag_index_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DINDEX -DFIXCOL -DFIXDIAG
2mols_fixantidiag_index_%: 2mols.cpp
	g++ -fPIC -std=c++11 -O4 -DNDEBUG -I$(ORDIR)/include -DARCH_K8 -Wno-deprecated -DUSE_CBC -DUSE_CLP -DUSE_BOP -DUSE_GLOP 2mols.cpp -L$(ORDIR)/lib -lprotobuf -lglog -lgflags -lCbcSolver -lCbc -lOsiCbc -lCgl -lClpSolver -lClp -lOsiClp -lOsi -lCoinUtils -lortools -lz -lrt -lpthread -o $@ -DORDER=$* -DINDEX -DFIXCOL -DFIXANTIDIAG
