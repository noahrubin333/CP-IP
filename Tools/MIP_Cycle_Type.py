import random
import os
import datetime
import sys

checked_already = []

def check_l(l):
    if ([0] + li) in checked_already:
        return False

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


n = int(sys.argv[1])
print("n = " + str(n))
li = list(range(1, n))
PATH_TO_IP_EXECUTABLE = "C:\\Users\\Noah\\source\\repos\\Gurobi Project1\\x64\\Debug"
os.chdir(PATH_TO_IP_EXECUTABLE)




cycleTypeTimings = {}
OEIS_A002865 = [1, 0, 1, 1, 2, 2, 4, 4, 7, 8, 12, 14, 21, 24, 34, 41, 55]
number_of_canonical_fixings = OEIS_A002865[n - 1]

log = open("C:\\Users\\Noah\\Desktop\\Random_Assignment_Log_MIP.log", "a")

log.write("LOGGING BEGIN " + str(datetime.datetime.now()) + "\n=================================\n")

log.close()
valid = float(0)
invalid = float(0)

for i in range(100):
    # Generate valid column fixing
    li = list(range(1, n))
    while not check_l(li) or li == list(range(1, n)):
        random.shuffle(li)

    li = [0] + li
    checked_already.append(li)

    # Determine the cycle type of the column fixing
    cycle_type = categorize(li, len(li))

    # Execute 2-MOLS(n) program using command prompt and record output
    result = os.popen('"Gurobi Project1.exe" ' + str(n) + " " + str(li).replace(" ", "")).readlines()[4:]

    # Check for validity of solution - if invalid do not include
    if 'Error' in result[1] and n != 6:
        invalid += 1
        print(str(tuple(li)) + "is invalid")
        continue
    else:
        valid += 1

    # Parse running time and nodes processed
    time = float(result[len(result) - 2][10:].replace("\n", ""))
    nodes = float(result[len(result) - 1][17:].replace("\n", ""))

    # If we already have seen fixings of this cycle type then append the column fixing and runtime
    if tuple(cycle_type) in cycleTypeTimings:
        cycleTypeTimings[tuple(cycle_type)].append([li, time, nodes])
    # Else create new entry in list of cycle types
    else:
        cycleTypeTimings.update({tuple(cycle_type): [[li, time, nodes]]})

    # Write actual squares, time data and node data to file
    log = open("C:\\Users\\Noah\\Desktop\\Random_Assignment_Log_MIP.log", "a")
    log.write("n = " + str(n) + "\nFixing: " + str(tuple(li)) + "\nCycle Type: " + str(cycle_type))
    log.writelines(result)
    log.write("---------------------------------\n")
    log.close()

log = open("C:\\Users\\Noah\\Desktop\\Random_Assignment_Log_MIP.log", "a")
log.write("LOGGING END " + str(datetime.datetime.now()) + "\n=================================\n")
log.close()

print(cycleTypeTimings.keys())
print(cycleTypeTimings)
if invalid != 0:
    print("Ratio of valid to invalid is ~" + str(valid / invalid))
