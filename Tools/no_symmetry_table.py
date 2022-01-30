import sys
import numpy as np

if "-windows" in sys.argv:
    PATH_TO_CPL_LOG = "C:\\Users\\Noah\\Desktop\\CPL_No_Sym.txt"
    PATH_TO_CPM_LOG = "C:\\Users\\Noah\\Desktop\\CPL_INDEX_NOSYM.txt"
    PATH_TO_CPI_LOG = "C:\\Users\\Noah\\Desktop\\CP_INDEX_NOSYM.txt"
    PATH_TO_IP_LOG = "C:\\Users\\Noah\\Desktop\\NRD_Log_IP_No_Symm.log"
else:
    PATH_TO_CPL_LOG = "../Logs/NRD_Log_CP_lin.log"
    PATH_TO_CPM_LOG = "../Logs/CPL_INDEX_NOSYM.txt"
    PATH_TO_CPI_LOG = "../Logs/NRD_Log_CP.log"
    PATH_TO_IP_LOG = "../Logs/NRD_Log_IP_No_Symm.log"

print("""\\begin{tabular}{@{\,}c@{\\quad}r@{\\quad}r@{\\quad}r@{\\quad}r@{\\quad}r@{\\quad}r@{\,}}
Model & """, end="")

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

for S in ('N'):
    print("IP", end=" & ")
    for n in range(5, 11):
        if not n in TIMINGS[S]:
            print("", end="")
            continue

        if TIMINGS[S][n][0][0] > 59999.9:
            print("Timeout", end="")
        else:
            print("{:,.1f}".format(TIMINGS[S][n][0][0]), end="")
        if (n != 10):
            print(" & ", end="")

    print(" \\\\")

LOG_FILE_CP = open(PATH_TO_CPL_LOG, "r")

DATA_RAW = LOG_FILE_CP.readlines()

LOG_FILE_CP.close()

# Parse out logging timestamps and separators
for i in range(len(DATA_RAW)):
    if "LOGGING" in DATA_RAW[i] or "==" in DATA_RAW[i]:
        DATA_RAW[i] = ""

DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")[:-1]

#=========================================================================

timings = [-1 for i in range(13)]

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)

    n = int(DATA_SPLIT[0][3:])
    SYMM_BREAKING = DATA_SPLIT[1][19]
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n")) / 1000.0

    if timings[n] == -1 and SYMM_BREAKING == 'N':
        timings[n] = TIME

print("CP-\\rlap{linear}\\hspace{3em}", end=" & ")
for n in range(5, 11):
    if timings[n] == -1:
        print("", end="")
    elif timings[n] > 59999.9:
        print("Timeout", end="")
    else:
        print("{:,.1f}".format(timings[n]), end="")
    if (n != 10):
        print(" & ", end="")
print(" \\\\")

LOG_FILE_CP = open(PATH_TO_CPM_LOG, "r")

DATA_RAW = LOG_FILE_CP.readlines()

LOG_FILE_CP.close()

# Parse out logging timestamps and separators
for i in range(len(DATA_RAW)):
    if "LOGGING" in DATA_RAW[i] or "==" in DATA_RAW[i]:
        DATA_RAW[i] = ""

DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")[:-1]

#=========================================================================

TIMINGS = {'N' : {}}

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)

    n = int(DATA_SPLIT[0][18:])
    SYMM_BREAKING = 'N'
    TIME = float(DATA_SPLIT[LEN - 4][10:].strip("\n"))
    COUNT_BRANCHES = float(DATA_SPLIT[LEN - 8][10:].strip("\n"))
    if n in TIMINGS[SYMM_BREAKING]:
        TIMINGS[SYMM_BREAKING][n][0].append(TIME)
        TIMINGS[SYMM_BREAKING][n][1].append(COUNT_BRANCHES)
    else:
        TIMINGS[SYMM_BREAKING].update({n: [[TIME], [COUNT_BRANCHES]]})

AVG_TIME = STD_DEV = 0.0

print("CP-\\rlap{moddiv}\\hspace{3em}", end=" & ")
for n in range(5, 11):

    AVG_TIME = np.round(np.average(TIMINGS['N'][n][0]), 1)

    if AVG_TIME > 59999.9:
        print("Timeout", end="")
    else:
        print("{:,.1f}".format(AVG_TIME), end="")
    if (n != 10):
        print(" & ", end="")
print(" \\\\")

LOG_FILE_CP = open(PATH_TO_CPI_LOG, "r")

DATA_RAW = LOG_FILE_CP.readlines()

LOG_FILE_CP.close()

# Parse out logging timestamps and separators
for i in range(len(DATA_RAW)):
    if "LOGGING" in DATA_RAW[i] or "==" in DATA_RAW[i]:
        DATA_RAW[i] = ""

DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")[:-1]

#=========================================================================

timings = [-1 for i in range(13)]

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)

    n = int(DATA_SPLIT[0][3:])
    SYMM_BREAKING = DATA_SPLIT[1][19]
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n")) / 1000.0

    if timings[n] == -1 and SYMM_BREAKING == 'N':
        timings[n] = TIME

print("CP-\\rlap{index}\\hspace{3em}", end=" & ")
for n in range(5, 11):
    if timings[n] == -1:
        print("& ", end="")
    elif timings[n] > 59999.9 or timings[n] == 50000:
        print("Timeout", end="")
    else:
        print("{:,.1f}".format(timings[n]), end="")
    if (n != 10):
        print(" & ", end="")
print("")

print("""\end{tabular}""")

