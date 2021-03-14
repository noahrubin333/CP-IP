import sys
import numpy as np

PATH_TO_CP_LOG = "C:\\Users\\Noah\\Desktop\\Log_CP_NRD.log"
PATH_TO_IP_LOG = "C:\\Users\\Noah\\Desktop\\Log_MIP_NRD.log"

if "-linux" in sys.argv:
    PATH_TO_CP_LOG = "../Logs/NRD_Log_CP_Z_LAST.log"
    PATH_TO_IP_LOG = "../Logs/NRD_Log_MIP_Z_LAST.log"

print("""\\begin{tabular}{c@{\\quad}c@{\\qquad}c@{\\quad}r@{\\quad}r@{\\quad}r@{\\quad}r@{\\quad}r}
& Symmetry & \multicolumn{6}{c}{Order} \\\\
Model & Breaking & """, end="")

for n in range(5, 11):
    print("\\df{%d}" % n, end="")
    if (n != 10):
        print(" & ", end="")
print(" \\\\ \\hline")

#=========================================================================

LOG_FILE_MIP = open(PATH_TO_IP_LOG, "r")

DATA_RAW = LOG_FILE_MIP.readlines()

LOG_FILE_MIP.close()

# Parse out logging timestamps and separators
for i in range(len(DATA_RAW)):
    if "LOGGING" in DATA_RAW[i] or "==" in DATA_RAW[i]:
        DATA_RAW[i] = ""

DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")[:-1]

TIMINGS = {'N' : {}, 'R' : {}, 'D' : {}}

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)

    n = int(DATA_SPLIT[0][3:])
    SYMM_BREAKING = DATA_SPLIT[1][19]
    TIME = float(DATA_SPLIT[LEN - 3][9:].strip("\n"))
    COUNT_NODES = float(DATA_SPLIT[LEN - 2][17:].strip("\n"))

    if n in TIMINGS[SYMM_BREAKING]:
        TIMINGS[SYMM_BREAKING][n][0].append(TIME)
        TIMINGS[SYMM_BREAKING][n][1].append(COUNT_NODES)
    else:
        TIMINGS[SYMM_BREAKING].update({n: [[TIME], [COUNT_NODES]]})

AVG_TIME = STD_DEV = 0.0

for S in ('N', 'R', 'D'):
    print("IP", end=" & ")
    if S == 'N':
        print("None", end=" & ")
    elif S == 'R':
        print("Row Only", end=" & ")
    elif S == 'D':
        print("Row \& Col", end=" & ")
    for n in range(5, 11):
        if not S in TIMINGS or not n in TIMINGS[S]:
            print("T", end=" ")
        else:
            AVG_TIME = np.round(np.average(TIMINGS[S][n][0]), 1)
            STD_DEV_TIME = np.round(np.std(TIMINGS[S][n][0]), 1)

            AVG_NODES = np.round(np.average(TIMINGS[S][n][1]), 1)
            STD_DEV_NODES = np.round(np.std(TIMINGS[S][n][1]), 1)

            if AVG_TIME == 60000:
                print("T", end=" ")
            else:
                print(str(AVG_TIME), end=" ")
        if (n != 10):
            print(" & ", end="")

    print(" \\\\")

LOG_FILE_CP = open(PATH_TO_CP_LOG, "r")

DATA_RAW = LOG_FILE_CP.readlines()

LOG_FILE_CP.close()

# Parse out logging timestamps and separators
for i in range(len(DATA_RAW)):
    if "LOGGING" in DATA_RAW[i] or "==" in DATA_RAW[i]:
        DATA_RAW[i] = ""

DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")[:-1]

#=========================================================================

TIMINGS = {'N' : {}, 'R' : {}, 'D' : {}}

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)

    n = int(DATA_SPLIT[0][3:])
    SYMM_BREAKING = DATA_SPLIT[1][19]
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n")) / 1000.0 #ms to s
    COUNT_BRANCHES = float(DATA_SPLIT[LEN - 3][29:].strip("\n"))
    if n in TIMINGS[SYMM_BREAKING]:
        TIMINGS[SYMM_BREAKING][n][0].append(TIME)
        TIMINGS[SYMM_BREAKING][n][1].append(COUNT_BRANCHES)
    else:
        TIMINGS[SYMM_BREAKING].update({n: [[TIME], [COUNT_BRANCHES]]})

AVG_TIME = STD_DEV = 0.0

for S in ('N', 'R', 'D'):
    print("CP-index", end=" & ")
    if S == 'N':
        print("None", end=" & ")
    elif S == 'R':
        print("Row Only", end=" & ")
    elif S == 'D':
        print("Row \& Col", end=" & ")
    for n in range(5, 11):

        AVG_TIME = np.round(np.average(TIMINGS[S][n][0]), 1)
        STD_DEV_TIME = np.round(np.std(TIMINGS[S][n][0]), 1)

        AVG_BRANCHES = np.round(np.average(TIMINGS[S][n][1]), 1)
        STD_DEV_BRANCHES = np.round(np.std(TIMINGS[S][n][1]), 1)

        print(str(AVG_TIME), end=" ")
        if (n != 10):
            print(" & ", end="")
    print(" \\\\" if S != 'D' else "")
print("""\end{tabular}""")

