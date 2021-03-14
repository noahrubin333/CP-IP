/**

    Carleton University RTA Project - ain.cpp

    Noah Rubin
    
    Under supervision of Dr. Kevin Cheung and Dr. Brett Stevens.

    Let L1, L2 and L3 be Latin Squares of order n (as defined in OLS_HEADER.h).

    This program fills in n / 2 transversals of L3 (with Gurobi) and then attempts to determine whether L1 and L2 can
    be filled in to be orthogonal to each other and to the partially filled L3 (also w/ Gurobi). If the system is feasable
    then we may proceed to filling in the remainder of L3 with the CP solver (or-tools).

**/

// Relevant includes
#include "OLSHeader.h"

// Begin main
int main(int argc, char* argv[]) {

    std::chrono::high_resolution_clock::time_point tic;

    int n = 0;
    int colPossibility = -1;

    // Declare counters
    int i = 0, j = 0, k = 0, l = 0;

    int *firstCol_L2 = NULL;
    
    if (argc == 2) {
        n = atoi(argv[1]);
    }
    else if (strlen(argv[2]) > 1) {
        n = atoi(argv[1]);
        firstCol_L2 = new int[n];
        int count = 0;
        for (i = 0; i < strlen(argv[2]) - 1; i++) {
            if (48 <= (int)argv[2][i] && (int)argv[2][i] <= 57) {
                firstCol_L2[count] = (int)argv[2][i] - 48;
                count++;
            }
        }

        if (n != count) {
            cout << "n = " << n << " but you entered a list of length " << count << "\n";
            exit(0);
        }

        //exit(0);
    }
    else if (argc == 3) {
        colPossibility = atoi(argv[2]);
        n = atoi(argv[1]);
    }
    else {
        cout << "Enter dimension of squares\n";
        exit(0);
    }

    std::cout << "n: " << n << "\n";

    // Try block is required for running Gurobi programs
    try {
        // Setup environment and begin logging
        GRBEnv env = GRBEnv(true);
        if (colPossibility != -1) {
            env.set("LogFile", "C:\\Users\\Noah\\Desktop\\Logs\\log_n = " + to_string(n) + "_" + to_string(colPossibility) + ".log");
        }
        else {
            env.set("LogFile", "C:\\Users\\Noah\\Desktop\\Logs\\log_n = " + to_string(n) + ".log");
        }
        env.start();
        GRBModel model = GRBModel(env);

        GRBVar* x = new GRBVar[(int)pow(n, 4)];

        int nCube = (int)pow(n, 3);
        int nSqr = (int)pow(n, 2);

        // Instantiate variables for L1 and L2
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        string name = "x_" + to_string(i) + "_" + to_string(j) + "_" + to_string(k) + "_" + to_string(l);
                        x[i * nCube + j * nSqr + k * n + l] = model.addVar(0.0, 1.0, 0.0, GRB_BINARY, name);
                    }
                }
            }
        }

        // Implement constraint (1)
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                GRBLinExpr expr = 0;
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "a_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (2)
        for (i = 0; i < n; i++) {
            for (k = 0; k < n; k++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++) {
                    for (l = 0; l < n; l++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "b_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (3)
        for (i = 0; i < n; i++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (j = 0; j < n; j++) {
                    for (k = 0; k < n; k++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "c_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (4)
        for (j = 0; j < n; j++) {
            for (k = 0; k < n; k++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (l = 0; l < n; l++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "d_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (5)
        for (j = 0; j < n; j++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (k = 0; k < n; k++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "e_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Implement constraint (6)
        for (k = 0; k < n; k++) {
            for (l = 0; l < n; l++) {
                GRBLinExpr expr = 0;
                for (i = 0; i < n; i++) {
                    for (j = 0; j < n; j++) {
                        expr += x[i * nCube + j * nSqr + k * n + l];
                    }
                }
                string s = "f_" + to_string(i) + "_" + to_string(j);
                model.addConstr(expr, GRB_EQUAL, 1.0, s);
            }
        }

        // Fix first row of both squares

        for (i = 0; i < n; i++) {
            GRBLinExpr expr = 0;
            expr += x[i * nSqr + i * n + i];
            model.addConstr(expr == 1.0);
        }

        /*
        model.addConstr(x[nCube + n + 2] == 1.0);
        model.addConstr(x[(n - 2) * nCube + (n - 2) * n + (n - 1)] == 1.0);

        // Implement symmetry breaking - First column of each square
        for (i = 2; i < n; i++) { // Row index
            if (i == n - 2) { continue; }
            for (j = 1; j < n; j++) { // Value in L2
                for (k = 0; k < n; k++) { // Value in L1
                    //cant be i or > i + 1
                    if (j == i || j > i + 1 || k != i) {
                        GRBLinExpr expr = 0.0;
                        expr += x[i * nCube + k * n + j];
                        model.addConstr(expr == 0.0);
                    }
                }
            }
        }
        */

        int cols_7[][7] = {
            {0, 2, 1, 4, 3, 6, 5},
            {0, 2, 3, 4, 1, 6, 5},
            {0, 2, 3, 1, 5, 6, 4},
            {0, 2, 1, 4, 5, 6, 3},
            {0, 2, 3, 4, 5, 6, 1},
             };


        
        int cols_8[][8] = { {0, 2, 3, 1, 5, 4, 7, 6},
            {0, 2, 1, 4, 5, 3, 7, 6 },
            {0, 2, 3, 4, 5, 1, 7, 6 },
            {0, 2, 1, 4, 3, 6, 7, 5 },
            {0, 2, 3, 4, 1, 6, 7, 5 },
            {0, 2, 3, 1, 5, 6, 7, 4 },
            {0, 2, 1, 4, 5, 6, 7, 3 },
            {0, 2, 3, 4, 5, 6, 7, 1 } };

         
        int cols_9[][9] = { { 0, 2, 1, 4, 3, 6, 5, 8, 7 },
            {0, 2, 3, 4, 1, 6, 5, 8, 7},
        {0, 2, 3, 1, 5, 6, 4, 8, 7},
        {0, 2, 1, 4, 5, 6, 3, 8, 7},
        {0, 2, 3, 4, 5, 6, 1, 8, 7},
        {0, 2, 3, 1, 5, 4, 7, 8, 6},
        {0, 2, 1, 4, 5, 3, 7, 8, 6},
        {0, 2, 3, 4, 5, 1, 7, 8, 6},
        {0, 2, 1, 4, 3, 6, 7, 8, 5},
        {0, 2, 3, 4, 1, 6, 7, 8, 5},
        {0, 2, 3, 1, 5, 6, 7, 8, 4},
        {0, 2, 1, 4, 5, 6, 7, 8, 3},
        {0, 2, 3, 4, 5, 6, 7, 8, 1} };
        
        
    int cols_10[][10] = { { 0, 2, 3, 1, 5, 4, 7, 6, 9, 8 },
        {0, 2, 1, 4, 5, 3, 7, 6, 9, 8 },
        {0, 2, 3, 4, 5, 1, 7, 6, 9, 8 },
        {0, 2, 1, 4, 3, 6, 7, 5, 9, 8 },
        {0, 2, 3, 4, 1, 6, 7, 5, 9, 8 },
        {0, 2, 3, 1, 5, 6, 7, 4, 9, 8 },
        {0, 2, 1, 4, 5, 6, 7, 3, 9, 8 },
        {0, 2, 3, 4, 5, 6, 7, 1, 9, 8 },
        {0, 2, 1, 4, 3, 6, 5, 8, 9, 7 },
        {0, 2, 3, 4, 1, 6, 5, 8, 9, 7 },
        {0, 2, 3, 1, 5, 6, 4, 8, 9, 7 },
        {0, 2, 1, 4, 5, 6, 3, 8, 9, 7 },
        {0, 2, 3, 4, 5, 6, 1, 8, 9, 7 },
        {0, 2, 3, 1, 5, 4, 7, 8, 9, 6 },
        {0, 2, 1, 4, 5, 3, 7, 8, 9, 6 },
        {0, 2, 3, 4, 5, 1, 7, 8, 9, 6 },
        {0, 2, 1, 4, 3, 6, 7, 8, 9, 5 },
        {0, 2, 3, 4, 1, 6, 7, 8, 9, 5 },
        {0, 2, 3, 1, 5, 6, 7, 8, 9, 4 },
        {0, 2, 1, 4, 5, 6, 7, 8, 9, 3 },
        {0, 2, 3, 4, 5, 6, 7, 8, 9, 1 } };
        
        if (argc > 2) {
        if (strlen(argv[2]) > 2) {
            for (j = 1; j < n; j++) {
                GRBLinExpr expr = 0;
                expr += x[j * nCube + j * n + firstCol_L2[j]];
                model.addConstr(expr == 1.0);
            }
        }
        }
        

        if (colPossibility != -1) {
            if (n == 7) {
                for (j = 1; j < n; j++) {
                    GRBLinExpr expr = 0;
                    expr += x[j * nCube + j * n + cols_7[colPossibility][j]];
                    model.addConstr(expr == 1.0);
                }
            }
            else if (n == 8) {
                for (j = 1; j < n; j++) {
                    GRBLinExpr expr = 0;
                    expr += x[j * nCube + j * n + cols_8[colPossibility][j]];
                    model.addConstr(expr == 1.0);
                }
            }
            else if (n == 9) {
                for (j = 1; j < n; j++) {
                    GRBLinExpr expr = 0;
                    expr += x[j * nCube + j * n + cols_9[colPossibility][j]];
                    model.addConstr(expr == 1.0);
                }
            }
            else if (n == 10) {
                for (j = 1; j < n; j++) {
                    GRBLinExpr expr = 0;
                    expr += x[j * nCube + j * n + cols_10[colPossibility][j]];
                    model.addConstr(expr == 1.0);
                }
            }
        }

        // Set model params
        model.set(GRB_IntParam_Threads, 1);
        model.set(GRB_IntParam_LogToConsole, 0);
        //model.set(GRB_DoubleParam_NodeLimit, 0);


        // Write model to a .lp file so we can test w/ CPLEX
        model.write("C:\\Users\\Noah\\MOLS_n=" + to_string(n) + ".lp");

        srand(time(nullptr));
        int seed = abs(rand() * std::rand());
        model.set(GRB_IntParam_Seed, seed);

        std::cout << "seed: " << seed << endl;

        tic = std::chrono::high_resolution_clock::now();

        // Execute model
        model.optimize();
        
        
        /*
        ofstream TIME_FILE;
        TIME_FILE.open("C:\\Users\\Noah\\Desktop\\timings_n=" + to_string(n) + "_" + to_string(colPossibility) + ".txt", std::ofstream::out | std::ofstream::app);
        TIME_FILE << to_string(totalTime) << "\n";
        TIME_FILE.close();
        */
        

        // Print L1
        std::cout << endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
                            std::cout << to_string(k) << " ";
                    }
                }
            }
            std::cout << endl;
        }

        // Print L2
        std::cout << endl;
        for (i = 0; i < n; i++) {
            for (j = 0; j < n; j++) {
                for (k = 0; k < n; k++) {
                    for (l = 0; l < n; l++) {
                        if (x[i * nCube + j * nSqr + k * n + l].get(GRB_DoubleAttr_X) > 0.5)
                            std::cout << to_string(l) << " ";
                    }
                }
            }
            std::cout << endl;
        }

        cout << "Time (s): " << model.get(GRB_DoubleAttr_Runtime) << endl;
        cout << "Nodes Processed: " << model.get(GRB_DoubleAttr_NodeCount) << endl;
    }
    // Handle exceptions
    catch (GRBException e) {
        std::cout << "Error code = " << e.getErrorCode() << endl;
        std::cout << e.getMessage() << endl;
        auto toc = std::chrono::high_resolution_clock::now();
        auto totalTime = chrono::duration_cast<chrono::milliseconds>(toc - tic).count();
        cout << "Time (s): " << to_string(totalTime) << endl;
        cout << "Nodes Processed: " << 0.0 << endl;
    }
    catch (...) {
        std::cout << "Exception during optimization" << endl;
    }
    return 0;
}

int check_fixing(int* fixing, int x) {
    for (int i = 0; i < x; i++) {
        if (fixing[i] == i) {
            return 0;
        }
    }
    return 1;
}