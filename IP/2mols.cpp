/**
    Noah Rubin
    2020-06-29

    The following model implements (IP1) from the paper "An LP based proof of the nonexistence of a pair of orthogonal latin squares of order 6"

    The order of the Latin squares to search for needs to be given as the macro definition N
**/


#include <iostream>
#include "gurobi_c++.h"
#include <string>

using namespace std;

// Declare variables
GRBVar xvar[N][N][N][N];

#ifdef CALLBACK_CUTS
#include "callback.cpp"
#endif

#define WRITE_MODEL 1 // Write model to lp file
#define USE_SOS 0 // Include SOS-1 constraints
#define USE_LINEAR 1 // Include linear constraints
#define NO_GUROBI_CUTS 1 // Disable Gurobi's cuts
#define SHADOW_PRICE 1 // Shadow price variable branching
#define MIP_FOCUS 1 // Focus on finding feasible solutions
#define RINS 0 // Use Relaxation Induced Neighborhood Search (RINS)
#define LP_METHOD 0 // Set LP solving to use
#define BRANCH_DIR 0 // Set branch direction
#define NODE_LIMIT 0 // Limit on nodes to solve
#ifdef NOSYM
#define SYMMETRY_DETECTION 0 // Gurobi symmetry detection disabled
#else
#define SYMMETRY_DETECTION 1 // Gurobi symmetry detection enabled
#endif

#ifdef RANDOM_OBJECTIVE
#include <algorithm>
#include <vector>
#endif

int main(int argc, char* argv[]) {
    try {   

        // Setup environment
        GRBEnv env = GRBEnv(true);
        env.set("LogFile", "LatinSquare.log");
        env.start();

        GRBModel model = GRBModel(env);

        int i = 0, j = 0, k = 0, l = 0;

        // Instantiate variables
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                        xvar[i][j][k][l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                    }
                }
            }
        }
        
        // Implement constraint (1)

        for (k = 0; k < N; k++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (j = 0; j < N; j++) {
                        expr += xvar[i][j][k][l];
                    }
                }
                string s = "C1_" + to_string(k) + "_" + to_string(l);
#if USE_LINEAR
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif
#if USE_SOS
                GRBVar vars[N*N];
                double weights[N*N];
                for (i = 0; i < N; i++) {
                    for (j = 0; j < N; j++) {
                        weights[i+N*j] = 1;
                        vars[i+N*j] = xvar[i][j][k][l];
                    }
                }
                model.addSOS(vars, weights, N*N, GRB_SOS_TYPE1);
#endif
            }
        }

        // Implement constraint (2)

        for (j = 0; j < N; j++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (k = 0; k < N; k++) {
                        expr += xvar[i][j][k][l];
                    }
                }
                string s = "C2_" + to_string(j) + "_" + to_string(l);
#if USE_LINEAR
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif
#if USE_SOS
                GRBVar vars[N*N];
                double weights[N*N];
                for (i = 0; i < N; i++) {
                    for (k = 0; k < N; k++) {
                        weights[i+N*k] = 1;
                        vars[i+N*k] = xvar[i][j][k][l];
                    }
                }
                model.addSOS(vars, weights, N*N, GRB_SOS_TYPE1);
#endif
            }
        }

        // Implement constraint (3)

        for (j = 0; j < N; j++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < N; i++) {
                    for (l = 0; l < N; l++) {
                        expr += xvar[i][j][k][l];
                    }
                }
                string s = "C3_" + to_string(j) + "_" + to_string(k);
#if USE_LINEAR
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif
#if USE_SOS
                GRBVar vars[N*N];
                double weights[N*N];
                for (i = 0; i < N; i++) {
                    for (l = 0; l < N; l++) {
                        weights[i+N*l] = 1;
                        vars[i+N*l] = xvar[i][j][k][l];
                    }
                }
                model.addSOS(vars, weights, N*N, GRB_SOS_TYPE1);
#endif
            }
        }

        // Implement constraint (4)

        for (i = 0; i < N; i++) {
            for (l = 0; l < N; l++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    for (k = 0; k < N; k++) {
                        expr += xvar[i][j][k][l];
                    }
                }
                string s = "C4_" + to_string(i) + "_" + to_string(l);
#if USE_LINEAR
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif
#if USE_SOS
                GRBVar vars[N*N];
                double weights[N*N];
                for (j = 0; j < N; j++) {
                    for (k = 0; k < N; k++) {
                        weights[j+N*k] = 1;
                        vars[j+N*k] = xvar[i][j][k][l];
                    }
                }
                model.addSOS(vars, weights, N*N, GRB_SOS_TYPE1);
#endif
            }
        }

        // Implement constraint (5)

        for (i = 0; i < N; i++) {
            for (k = 0; k < N; k++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < N; j++) {
                    for (l = 0; l < N; l++) {
                        expr += xvar[i][j][k][l];
                    }
                }
                string s = "C5_" + to_string(i) + "_" + to_string(k);
#if USE_LINEAR
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif
#if USE_SOS
                GRBVar vars[N*N];
                double weights[N*N];
                for (j = 0; j < N; j++) {
                    for (l = 0; l < N; l++) {
                        weights[j+N*l] = 1;
                        vars[j+N*l] = xvar[i][j][k][l];
                    }
                }
                model.addSOS(vars, weights, N*N, GRB_SOS_TYPE1);
#endif
            }
        }

        // Implement constraint (6)
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                GRBLinExpr expr = 0;
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        expr += xvar[i][j][k][l];
                    }
                }
                string s = "C6_" + to_string(i) + "_" + to_string(j);
#if USE_LINEAR
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
#endif
#if USE_SOS
                GRBVar vars[N*N];
                double weights[N*N];
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        weights[k+N*l] = 1;
                        vars[k+N*l] = xvar[i][j][k][l];
                    }
                }
                model.addSOS(vars, weights, N*N, GRB_SOS_TYPE1);
#endif
            }
        }

        // Implement Symmetry Breaking

#if !defined(NOSYM)
        // Fix first row of both squares
        for (i = 0; i < N; i++)
            xvar[0][i][i][i].set(GRB_DoubleAttr_LB, 1);

#ifdef FIXCOL
        cout << "Fixing first column of second square to be (0, 2, 3, 4, ..., n-1, 1)" << endl;

        // Fix first column of both squares
        for (i = 1; i < N-1; i++)
            xvar[i][0][i][i+1].set(GRB_DoubleAttr_LB, 1);

        // Fix values of entries at (N-1,0) to be (N-1,1)
        xvar[N-1][0][N-1][1].set(GRB_DoubleAttr_LB, 1);
#elif !defined(NOFIXCOL)
        // Fix first column of L1
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                if(i != j) {
                    for (k = 0; k < N; k++) {
                        xvar[i][0][j][k].set(GRB_DoubleAttr_UB, 0);
                    }
                }
            }
        }

        // Pairs (i,i) do not appear after first row
        for (i = 1; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    xvar[i][j][k][k].set(GRB_DoubleAttr_UB, 0);
                }
            }
        }

#if !defined(BACKCIRC)
        // Fix domains of first column of L2
        for (i = 0; i < N; i++) {
            for (j = i+2; j < N; j++) {
                xvar[i][0][i][j].set(GRB_DoubleAttr_UB, 0);
            }
        }
        
        // Fix values of entries at (1,0) to be (1,2)
        xvar[1][0][1][2].set(GRB_DoubleAttr_LB, 1);

        // Fix values of entries at (N-2,0) to be (N-2,N-1)
        xvar[N-2][0][N-2][N-1].set(GRB_DoubleAttr_LB, 1);
#endif
#endif
#endif

#if defined(FIXDIAG1)
        cout << "Fixing diagonal of first square to be zeros" << endl;
        // Fix diagonal of the first square to be all zeros
        for(i = 1; i < N; i++)
            for(j = 1; j < N; j++)
                for(k = 0; k < N; k++)
                    xvar[i][i][j][k].set(GRB_DoubleAttr_UB, 0);
#elif defined(FIXDIAG2)
        cout << "Fixing diagonal of second square to be zeros" << endl;
        // Fix diagonal of the second square to be all zeros
        for(i = 1; i < N; i++)
            for(j = 0; j < N; j++)
                for(k = 1; k < N; k++)
                    xvar[i][i][j][k].set(GRB_DoubleAttr_UB, 0);
#elif defined(BACKCIRC)
        cout << "Fixing first square to be back-circulant" << endl;
        // Fix first square to be back-circulant
        for(i = 1; i < N; i++)
            for(j = 0; j < N-1; j++)
                for(k = 0; k < N; k++)
                    for(l = 0; l < N; l++)
                        if((i+j)%N != k)
                            xvar[i][j][k][l].set(GRB_DoubleAttr_UB, 0);
#elif defined(FIXANTIDIAG)
        cout << "Fixing anti-diagonal of first square to be (n-1)s" << endl;
        // Fix anti-diagonal of the first square to be all (n-1)s
        for(i = 1; i < N; i++)
            for(j = 0; j < N-1; j++)
                for(k = 0; k < N; k++)
                    xvar[i][N-(i+1)][j][k].set(GRB_DoubleAttr_UB, 0);
#endif

#ifdef HINTS
        cout << "Using hints to suggest that the first square is back-circulant" << endl;
        for(i = 1; i < N; i++)
            for(j = 0; j < N-1; j++)
                for(k = 0; k < N; k++)
                    for(l = 0; l < N; l++)
                        if((i+j)%N != k)
                            xvar[i][j][k][l].set(GRB_DoubleAttr_VarHintVal, 0);
#endif

        // Generate objective function as sum x_ijkl forall i,j,k,l in {0, ..., n - 1}

#ifdef RANDOM_OBJECTIVE
        vector<int> weights;
        for(i = 1; i <= N*N*N*N; i++)
            weights.push_back(i);
        if(argc >= 2) {
            int seed = stoi(argv[1]);
            srand(seed);
        }
        random_shuffle(weights.begin(), weights.end());
#endif

        GRBLinExpr objective = 0;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
#ifdef RANDOM_OBJECTIVE
                        objective += weights[i+j*N+k*N*N+l*N*N*N]*xvar[i][j][k][l];
#else
                        objective += xvar[i][j][k][l];
#endif
                    }
                }
            }
        }

        //model.addConstr(objective == N*N);
        model.setObjective(objective, GRB_MAXIMIZE);

#if NODE_LIMIT
        // Solve 10000 nodes
        model.set(GRB_DoubleParam_NodeLimit, 10000);
#endif

#if SYMMETRY_DETECTION == 0
        // Disable symmetry detection
        model.set(GRB_IntParam_Symmetry, 0);
#endif

#if BRANCH_DIR
        // Force branch direction
        model.set(GRB_IntParam_BranchDir, 1);
#endif

#if LP_METHOD
        // Set LP solving method to match Appa et al.
        model.set(GRB_IntParam_Method, 2);
        model.set(GRB_IntParam_NodeMethod, 1);
#endif

#if SHADOW_PRICE
        model.set(GRB_IntParam_VarBranch, 1);
#endif
#if MIP_FOCUS
        model.set(GRB_IntParam_MIPFocus, 1);
#endif
#if RINS
        model.set(GRB_IntParam_RINS, 1);
        model.set(GRB_IntParam_SubMIPNodes, 1000);
#endif
        model.set(GRB_IntParam_SolutionLimit, 1);
        model.set(GRB_IntParam_Threads, 1);
        if(argc >= 2) {
            int seed = stoi(argv[1]);
            model.set(GRB_IntParam_Seed, seed);
            cout << "Using random seed " << seed << endl;
        }

#if NO_GUROBI_CUTS
        // Turn off all the cuts used by Gurobi
        model.set(GRB_IntParam_Cuts, 0);
#endif

#ifdef CALLBACK_CUTS
        // Don't allow presolve to translate constraints
        model.set(GRB_IntParam_PreCrush, 1);
        
        cout << "Using custom clique cuts" << endl;
        model.setCallback(new MyCallback);
#endif

#if WRITE_MODEL
#ifdef FIXCOL
        model.write("model/2mols.fixcol." + to_string(N) + ".lp");
#else
        model.write("model/2mols." + to_string(N) + ".lp");
#endif
#endif

        model.optimize();

#ifdef CALLBACK_CUTS
        cout << "Generated " << totalcuts << " custom clique cuts.\n";
#endif

        cout << endl;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        if (xvar[i][j][k][l].get(GRB_DoubleAttr_X) > 0.5)
                            cout << to_string(k) << " ";
                    }
                }
            }
            cout << endl;
        }

        cout << endl;
        for (i = 0; i < N; i++) {
            for (j = 0; j < N; j++) {
                for (k = 0; k < N; k++) {
                    for (l = 0; l < N; l++) {
                        if (xvar[i][j][k][l].get(GRB_DoubleAttr_X) > 0.5)
                            cout << to_string(l) << " ";
                    }
                }
            }
            cout << endl;
        }
    }
    catch (GRBException e) {
        cout << "Error code = " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
    }
    catch (...) {
        cout << "Exception during optimization" << endl;
    }

    return 0;
}



