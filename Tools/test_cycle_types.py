import random
import os
import sys

CycleTypes_Checked_Already = []

def check_l(l):
    #if ([0] + l) in CycleTypes_Checked_Already:
    #   return False

    for i in range(len(l)):
        if l[i] == i + 1:
            return False
    return True

def categorize(s, len):
    cycleType = []
    visited = [0] * len
    i, j, count, start, x, s_of_x = 0, 0, 0, 0, 0, 0

    while j != len:
        count, start, x, s_of_x = 1, i, i, s[i]
        visited[x] = 1
        while s_of_x != start:
            visited[s_of_x] = 1
            x = s_of_x
            s_of_x = s[x]
            count += 1
        cycleType.append(count)
        #count = 1
        for j in range(len + 1):
            if j != len:
                if visited[j] == 0:
                    i = j
                    break
    cycleType.sort()
    return cycleType

n = 9

li = list(range(1, n))
PATH_TO_SOURCE = "C:\\Users\\Noah\\source\\repos"
os.chdir(PATH_TO_SOURCE)

cycleTypeTimings = {}
OEIS_A002865 = [1, 0, 1, 1, 2, 2, 4, 4, 7, 8, 12, 14, 21, 24, 34, 41, 55]
number_of_canonical_fixings = OEIS_A002865[n - 1]

CP_LOG_FILE_NAME = "C:\\Users\\Noah\\Desktop\\Cycle_Type_Log_CP.log"
IP_LOG_FILE_NAME = "C:\\Users\\Noah\\Desktop\\Cycle_Type_Log_MIP.log"
CPL_LOG_FILE_NAME = "C:\\Users\\Noah\\Desktop\\Cycle_Type_Log_CPL.log"
#cycles_checked = {}
"""
for i in range(10000):
    # Generate valid column fixing
    li = list(range(1, n))
    while not check_l(li) or li == list(range(1, n)):
        random.shuffle(li)

    li = [0] + li
    CycleTypes_Checked_Already.append(li)

    # Determine the cycle type of the column fixing
    cycle_type = categorize(li, len(li))

    if tuple(cycle_type) not in cycles_checked:
        cycles_checked.update({tuple(cycle_type) : [li]})
    elif len(cycles_checked[tuple(cycle_type)]) < 3:
        cycles_checked[tuple(cycle_type)].append(li)

f = open("C:\\Users\\Noah\\Desktop\\Check_Log.log", "a")
for i in cycles_checked.keys():
    f.write(str(i) + " : " + str(cycles_checked[i]) + "\n")
f.write("===================\n")
f.close()
"""
cycles_checked = {
(1, 2, 2, 4) : [[0, 5, 8, 7, 6, 1, 3, 4, 2], [0, 2, 8, 6, 7, 1, 3, 4, 5], [0, 6, 5, 4, 3, 7, 1, 8, 2]],
(1, 8) : [[0, 5, 1, 8, 6, 3, 2, 4, 7], [0, 2, 4, 8, 6, 1, 3, 5, 7], [0, 6, 8, 7, 3, 4, 2, 1, 5]],
(1, 2, 6) : [[0, 3, 4, 7, 6, 8, 1, 2, 5], [0, 5, 4, 7, 1, 8, 2, 3, 6], [0, 7, 5, 1, 6, 2, 8, 4, 3]],
(1, 3, 5) : [[0, 2, 3, 6, 5, 8, 7, 1, 4], [0, 4, 3, 6, 5, 7, 2, 8, 1], [0, 4, 7, 1, 3, 6, 8, 5, 2]],
(1, 2, 3, 3) : [[0, 6, 4, 5, 2, 8, 7, 1, 3], [0, 2, 1, 4, 7, 8, 5, 3, 6], [0, 2, 4, 7, 1, 3, 8, 5, 6]],
(1, 2, 2, 2, 2) : [[0, 3, 5, 1, 6, 2, 4, 8, 7], [0, 8, 7, 6, 5, 4, 3, 2, 1], [0, 6, 3, 2, 5, 4, 1, 8, 7]]}


for i in cycles_checked.keys():
    for j in cycles_checked[i]:
        result = os.popen('"OR_TOOLS_TEST_2\\x64\\Release\\OR_TOOLS_TEST_2.exe" ' + str(n) + " " + str(j).replace(" ", "")).readlines()

        time = float(result[len(result) - 3][14:].replace("\n", ""))
        nodes = float(result[len(result) - 1][20:].replace("\n", ""))

        if tuple(i) in cycleTypeTimings:
            cycleTypeTimings[tuple(i)].append([j, time, nodes])
        else:
            cycleTypeTimings.update({tuple(i): [[j, time, nodes]]})

        log = open(CPL_LOG_FILE_NAME, "a")
        log.write("n = " + str(n) + "\nFixing: " + str(tuple(j)) + "\nCycle Type: " + str(i) + "\n")
        log.writelines(result)
        log.write("---------------------------------\n")
        log.close()
"""
for i in cycles_checked.keys():
    for j in cycles_checked[i]:
        result = os.popen('"Gurobi Project1\\x64\\Debug\\Gurobi Project1.exe" ' + str(n) + " " + str(j).replace(" ", "")).readlines()

        time = float(result[len(result) - 2][10:].replace("\n", ""))
        nodes = float(result[len(result) - 1][17:].replace("\n", ""))

        if tuple(i) in cycleTypeTimings:
            cycleTypeTimings[tuple(i)].append([j, time, nodes])
        else:
            cycleTypeTimings.update({tuple(i): [[j, time, nodes]]})

        log = open(IP_LOG_FILE_NAME, "a")
        log.write("n = " + str(n) + "\nFixing: " + str(tuple(j)) + "\nCycle Type: " + str(i) + "\n")
        log.writelines(result)
        log.write("---------------------------------\n")
        log.close()
"""