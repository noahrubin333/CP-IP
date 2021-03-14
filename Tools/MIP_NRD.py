import random
import os
import datetime
import sys

PATH_TO_IP_EXECUTABLE = "C:\\Users\\Noah\\source\\repos\\Gurobi Project1\\x64\\Debug"
os.chdir(PATH_TO_IP_EXECUTABLE)


log = open("C:\\Users\\Noah\\Desktop\\Random_Assignment_Log_MIP_NRD.log", "a")

log.write("LOGGING BEGIN " + str(datetime.datetime.now()) + "\n=================================\n")

log.close()


for n in range(5, 11):
    for choice in ('N', 'R', 'D'):
        print('"Gurobi Project1.exe" ' + str(n) + str(choice))
        result = os.popen('"Gurobi Project1.exe" ' + str(n) + " " + str(choice)).readlines()[4:]
        print(result)
        # Check for validity of solution - if invalid do not include
        if 'Error' in result[1] and n != 6:
            continue

        # Parse running time and nodes processed
        time = float(result[len(result) - 2][10:].replace("\n", ""))
        nodes = float(result[len(result) - 1][17:].replace("\n", ""))

        # Write actual squares, time data and node data to file
        log = open("C:\\Users\\Noah\\Desktop\\Random_Assignment_Log_MIP_NRD.log", "a")
        log.write("n = " + str(n) + "\nSymmetry Breaking: " + str(choice) + "\n\n")
        log.writelines(result)
        log.write("---------------------------------\n")
        log.close()

log = open("C:\\Users\\Noah\\Desktop\\Random_Assignment_Log_MIP_NRD.log", "a")
log.write("LOGGING END " + str(datetime.datetime.now()) + "\n=================================\n")
log.close()