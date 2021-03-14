import random
import os
import datetime
import sys

PATH_TO_IP_EXECUTABLE = "C:\\Users\\Noah\\source\\repos\\OR_TOOLS_TEST_1\\x64\\Release"
os.chdir(PATH_TO_IP_EXECUTABLE)


log = open("C:\\Users\\Noah\\Desktop\\NRD_Log_CP.log", "a")

log.write("LOGGING BEGIN " + str(datetime.datetime.now()) + "\n=================================\n")

log.close()


for n in range(5, 11):
    for choice in ('N', 'R', 'D'):
        print('"OR_TOOLS_TEST_1.exe" ' + str(n) + " " + str(choice))
        result = os.popen('"OR_TOOLS_TEST_1.exe" ' + str(n) + " " + str(choice)).readlines()
        print(result)
        # Check for validity of solution - if invalid do not include
        if 'Error' in result[1] and n != 6:
            continue

        # Parse running time and nodes processed
        time = float(result[len(result) - 3][14:].replace("\n", ""))
        nodes = float(result[len(result) - 1][20:].replace("\n", ""))

        # Write actual squares, time data and node data to file
        log = open("C:\\Users\\Noah\\Desktop\\NRD_Log_CP.log", "a")
        log.write("n = " + str(n) + "\nSymmetry Breaking: " + str(choice) + "\n")
        log.writelines(result)
        log.write("---------------------------------\n")
        log.close()

log = open("C:\\Users\\Noah\\Desktop\\NRD_Log_CP.log", "a")
log.write("LOGGING END " + str(datetime.datetime.now()) + "\n=================================\n")
log.close()