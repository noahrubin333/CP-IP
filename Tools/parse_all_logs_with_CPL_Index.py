import sys
import numpy as np
CycleTypeTimings_CP = [{} for i in range(7)]
CycleTypeTimings_CPL = [{} for i in range(7)]
CycleTypeTimings_CPM = [{} for i in range(7)]
CycleTypeTimings_IP = [{} for i in range(7)]
NOSYM_Timings_CP = [[] for i in range(9)]
DOMRED_Timings_CP = [[] for i in range(9)]
NOSYM_Timings_IP = [[] for i in range(9)]
DOMRED_Timings_IP = [[] for i in range(9)]
NOSYM_Timings_CPL = [[] for i in range(9)]
DOMRED_Timings_CPL = [[] for i in range(9)]
NOSYM_Timings_CPM = [[] for i in range(9)]
DOMRED_Timings_CPM = [[] for i in range(9)]

PATH_TO_CP_LOG = "../Logs/Cycle_Type_Log_CP_INDEX.log"
PATH_TO_IP_LOG = "../Logs/Cycle_Type_Log_MIP.log"
PATH_TO_CPL_LOG = "../Logs/Cycle_Type_Log_CP_LINEAR.log"
PATH_TO_CPM_LOG = "../Logs/Cycle_Type_Log_CP_MODDIV.log"
PATH_TO_NOSYM_LOG_CP = "../Logs/CP_INDEX_NOSYM.txt"
PATH_TO_DOMRED_LOG_CP = "../Logs/CP_INDEX_DOMRED.txt"
PATH_TO_NOSYM_LOG_CPM = "../Logs/CP_MODDIV_NOSYM.txt"
PATH_TO_DOMRED_LOG_CPM = "../Logs/CP_MODDIV_DOMRED.txt"
PATH_TO_NOSYM_LOG_IP = "../Logs/MIP_NOSYM.txt"
PATH_TO_DOMRED_LOG_IP = "../Logs/MIP_DOMRED.txt"
PATH_TO_NOSYM_LOG_CPL = "../Logs/MIP_NOSYM.txt"
PATH_TO_DOMRED_LOG_CPL = "../Logs/CP_LINEAR_DOMRED.txt"
PATH_TO_NOSYM_LOG_CPL = "../Logs/CP_LINEAR_NOSYM.txt"

#=======================================================
# Parse Cycle Type logs for CP_Index
#=======================================================
LOG_FILE_CP = open(PATH_TO_CP_LOG, "r")
DATA_RAW = LOG_FILE_CP.readlines()
LOG_FILE_CP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    CYCLE_TYPE = tuple(map(int, DATA_SPLIT[2][13:].strip(")\n").split(",")))
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n"))  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 3][29:].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    if tuple(CYCLE_TYPE) in CycleTypeTimings_CP[n - 5]:
        CycleTypeTimings_CP[n - 5][CYCLE_TYPE][0].append(TIME)
        CycleTypeTimings_CP[n - 5][CYCLE_TYPE][1].append(COUNT_BRANCHES)
    else:
        CycleTypeTimings_CP[n - 5].update({tuple(CYCLE_TYPE): [[TIME], [COUNT_BRANCHES]]})


#=======================================================
# Parse Cycle Type logs for MIP
#=======================================================
LOG_FILE_IP = open(PATH_TO_IP_LOG, "r")
DATA_RAW = LOG_FILE_IP.readlines()
LOG_FILE_IP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    CYCLE_TYPE = tuple(map(int, DATA_SPLIT[2][13:].strip(")\n").split(",")))
    TIME = float(DATA_SPLIT[LEN - 3][10:].strip("\n"))  # ms to s

    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 2][17:].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    if tuple(CYCLE_TYPE) in CycleTypeTimings_IP[n - 5]:
        CycleTypeTimings_IP[n - 5][CYCLE_TYPE][0].append(TIME)
        CycleTypeTimings_IP[n - 5][CYCLE_TYPE][1].append(COUNT_BRANCHES)
    else:
        CycleTypeTimings_IP[n - 5].update({tuple(CYCLE_TYPE): [[TIME], [COUNT_BRANCHES]]})




#=======================================================
# Parse Cycle Type logs for CP_Linear
#=======================================================
LOG_FILE_CPL = open(PATH_TO_CPL_LOG, "r")
DATA_RAW = LOG_FILE_CPL.readlines()
LOG_FILE_CPL.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    CYCLE_TYPE = tuple(map(int, DATA_SPLIT[2][13:].strip(")\n").split(",")))
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n")) / 1000.0 # ms to s

    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 3][29:].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    if tuple(CYCLE_TYPE) in CycleTypeTimings_CPL[n - 5]:
        CycleTypeTimings_CPL[n - 5][CYCLE_TYPE][0].append(TIME)
        CycleTypeTimings_CPL[n - 5][CYCLE_TYPE][1].append(COUNT_BRANCHES)
    else:
        CycleTypeTimings_CPL[n - 5].update({tuple(CYCLE_TYPE): [[TIME], [COUNT_BRANCHES]]})



#=======================================================
# Parse Cycle Type logs for CP_Moddiv
#=======================================================
LOG_FILE_CPM = open(PATH_TO_CPM_LOG, "r")
DATA_RAW = LOG_FILE_CPM.readlines()
LOG_FILE_CPM.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    if LEN < 10: continue
    n = int(DATA_SPLIT[0][3:])
    CYCLE_TYPE = tuple(map(int, DATA_SPLIT[2][13:].strip(")\n").split(",")))
    TIME = float(DATA_SPLIT[LEN - 15][len("Time elapsed: "):].replace("ms", "").strip("\n")) / 1000.0  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 8][len('branches:'):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    if tuple(CYCLE_TYPE) in CycleTypeTimings_CPM[n - 5]:
        CycleTypeTimings_CPM[n - 5][CYCLE_TYPE][0].append(TIME)
        CycleTypeTimings_CPM[n - 5][CYCLE_TYPE][1].append(COUNT_BRANCHES)
    else:
        CycleTypeTimings_CPM[n - 5].update({tuple(CYCLE_TYPE): [[TIME], [COUNT_BRANCHES]]})



#=======================================================
# Parse NoSym logs for MIP
#=======================================================
NOSYM_LOG_FILE_IP = open(PATH_TO_NOSYM_LOG_IP, "r")
DATA_RAW = NOSYM_LOG_FILE_IP.readlines()
NOSYM_LOG_FILE_IP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 3][10:].strip("\n"))  # ms to s

    COUNT_NODES_RAW = float(DATA_SPLIT[LEN - 2][len('Nodes Processed: '):].strip("\n"))
    COUNT_NODES = np.round(np.log2(COUNT_NODES_RAW), 2) if COUNT_NODES_RAW != 0 else np.round(0.00, 2)

    NOSYM_Timings_IP[n - 5] = [TIME, COUNT_NODES]
#=======================================================


#=======================================================
# Parse DomRed logs for MIP
#=======================================================
DOMRED_LOG_FILE_IP = open(PATH_TO_DOMRED_LOG_IP, "r")
DATA_RAW = DOMRED_LOG_FILE_IP.readlines()
DOMRED_LOG_FILE_IP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 3][10:].strip("\n"))  # ms to s

    COUNT_NODES_RAW = float(DATA_SPLIT[LEN - 2][len('Nodes Processed: '):].strip("\n"))
    COUNT_NODES = np.round(np.log2(COUNT_NODES_RAW), 2) if COUNT_NODES_RAW != 0 else np.round(0.00, 2)

    DOMRED_Timings_IP[n - 5] = [TIME, COUNT_NODES]
#=======================================================



#=======================================================
# Parse NoSym logs for CP_Index
#=======================================================
NOSYM_LOG_FILE_CP = open(PATH_TO_NOSYM_LOG_CP, "r")
DATA_RAW = NOSYM_LOG_FILE_CP.readlines()
NOSYM_LOG_FILE_CP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 15][len("Time elapsed: "):].replace("ms", "").strip("\n")) / 1000.0  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 8][len('branches: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    NOSYM_Timings_CP[n - 5] = [TIME, COUNT_BRANCHES]
#=======================================================



#=======================================================
# Parse DomRed logs for CP_Index
#=======================================================
DOMRED_LOG_FILE_CP = open(PATH_TO_DOMRED_LOG_CP, "r")
DATA_RAW = DOMRED_LOG_FILE_CP.readlines()
DOMRED_LOG_FILE_CP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 15][len("Time elapsed: "):].replace("ms", "").strip("\n")) / 1000.0  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 8][len('branches: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    DOMRED_Timings_CP[n - 5] = [TIME, COUNT_BRANCHES]
#=======================================================



#=======================================================
# Parse NoSym logs for CP_Linear
#=======================================================
NOSYM_LOG_FILE_CPL = open(PATH_TO_NOSYM_LOG_CPL, "r")
DATA_RAW = NOSYM_LOG_FILE_CPL.readlines()
NOSYM_LOG_FILE_CPL.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 15][len("Time elapsed: "):].replace("ms", "").strip("\n")) / 1000.0  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 8][len('branches: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    NOSYM_Timings_CPL[n - 5] = [TIME, COUNT_BRANCHES]
#=======================================================



#=======================================================
# Parse DomRed logs for CP_Linear
#=======================================================
DOMRED_LOG_FILE_CPL = open(PATH_TO_DOMRED_LOG_CPL, "r")
DATA_RAW = DOMRED_LOG_FILE_CPL.readlines()
DOMRED_LOG_FILE_CPL.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 15][len("Time elapsed: "):].replace("ms", "").strip("\n")) / 1000.0  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 8][len('branches: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    DOMRED_Timings_CPL[n - 5] = [TIME, COUNT_BRANCHES]
#=======================================================



#=======================================================
# Parse NoSym logs for CP_Moddiv
#=======================================================
NOSYM_LOG_FILE_CPM = open(PATH_TO_NOSYM_LOG_CPM, "r")
DATA_RAW = NOSYM_LOG_FILE_CPM.readlines()
NOSYM_LOG_FILE_CPM.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 15][len("Time elapsed: "):].replace("ms", "").strip("\n")) / 1000.0  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 8][len('branches: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    NOSYM_Timings_CPM[n - 5] = [TIME, COUNT_BRANCHES]
#=======================================================



#=======================================================
# Parse DomRed logs for CP_Moddiv
#=======================================================
DOMRED_LOG_FILE_CPM = open(PATH_TO_DOMRED_LOG_CPM, "r")
DATA_RAW = DOMRED_LOG_FILE_CPM.readlines()
DOMRED_LOG_FILE_CPM.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 15][len("Time elapsed: "):].replace("ms", "").strip("\n")) / 1000.0  # ms to s
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 8][len('branches: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    DOMRED_Timings_CPM[n - 5] = [TIME, COUNT_BRANCHES]


endl = "\\\\"
if "-nodes" in sys.argv:
    print("\\begin{tabular}{c@{\;\,}c@{\quad}r@{\quad}r@{\quad}r@{\quad}r@{\quad}r@{\quad}r@{\quad}r}\n$n$ & \\df{Sym.~Breaking} & \\df{IP} & \\df{CP-l} & \\df{CP-m} & \\df{CP-i} & \\df{IP} & \\df{CP-l} & \\df{CP-i} \\\\ \\hline")
else:
    print("\\begin{tabular}{c@{\;\,}c@{\quad}r@{\quad}r@{\quad}r@{\quad}r}\n$n$ & \\df{Sym.~Breaking} & \\df{IP} & \\df{CP-linear} & \\df{CP-moddiv} & \\df{CP-index} \\\\ \\hline")

for i in range(7):
    nodecounts = " & " + format(NRDTimings_IP[i][1], ".2f") + " & " + format(NRDTimings_CP[i][1], ".2f") + " & " + format(NRDTimings_CPL[i][1], ".2f") if "-nodes" in sys.argv else ""

    if "-min" in sys.argv:
        print((str(i + 5) + " & Domain red. & {:,}".format((NRDTimings_IP[i][0]) if NRDTimings_IP[i][0] < 60000 else -1)
            + " & {:,}".format((NRDTimings_CPL[i][0]) if NRDTimings_CPL[i][0] < 60000 else -1)
            + " & {:,}".format((NRDTimings_CPM[i][0]) if NRDTimings_CPM[i][0] < 60000 else -1)
            + " & {:,}".format((NRDTimings_CP[i][0]) if NRDTimings_CP[i][0] < 60000 else -1)
            + nodecounts
            + " \\\\").replace("-1", "Timeout"))
        min_time_IP = np.min([np.min(CycleTypeTimings_IP[i][k][0]) for k in CycleTypeTimings_IP[i].keys()]+[60000])
        min_time_CP = np.min([np.min(CycleTypeTimings_CP[i][k][0]) for k in CycleTypeTimings_CP[i].keys()]+[60000])
        min_time_CPL = np.min([np.min(CycleTypeTimings_CPL[i][k][0]) for k in CycleTypeTimings_CPL[i].keys()]+[60000])
        min_time_CPM = np.min([np.min(CycleTypeTimings_CPM[i][k][0]) for k in CycleTypeTimings_CPM[i].keys()]+[60000])
        if min_time_IP >= 60000: min_time_IP = -1
        if min_time_CP >= 60000: min_time_CP = -1
        if min_time_CPL >= 60000: min_time_CPL = -1
        if min_time_CPM >= 60000: min_time_CPM = -1
        print("{} & Cycle type min. & {:,} & {:,} & {:,} & {:,} \\\\".format(i+5, min_time_IP, min_time_CPL, min_time_CPM, min_time_CP).replace("-1", "Timeout"))
    else:
        print((str(i + 5) + " & Domain red. & {:,}".format(int(NRDTimings_IP[i][0]) if NRDTimings_IP[i][0] < 60000 else -1)
            + " & {:,}".format(int(NRDTimings_CPL[i][0]) if NRDTimings_CPL[i][0] < 60000 else -1)
            + " & {:,}".format(int(NRDTimings_CPM[i][0]) if NRDTimings_CPM[i][0] < 60000 else -1)
            + " & {:,}".format(int(NRDTimings_CP[i][0]) if NRDTimings_CP[i][0] < 60000 else -1)
            + nodecounts
            + " \\\\").replace("-1", "Timeout"))

    # for CURR_CYCLE_TYPE in CYCLE_TYPE_TIMINGS[i].keys():
    for CURR_CYCLE_TYPE in CycleTypeTimings_CP[i].keys():

        if CURR_CYCLE_TYPE in CycleTypeTimings_IP[i].keys():
            time_IP = int(np.median(CycleTypeTimings_IP[i][CURR_CYCLE_TYPE][0]))
            nodes_IP = np.median(CycleTypeTimings_IP[i][CURR_CYCLE_TYPE][1])
        else:
            time_IP = -2
            nodes_IP = -2

        if CURR_CYCLE_TYPE in CycleTypeTimings_CP[i].keys():
            time_CP = int(np.median(CycleTypeTimings_CP[i][CURR_CYCLE_TYPE][0]))
            nodes_CP = np.median(CycleTypeTimings_CP[i][CURR_CYCLE_TYPE][1])
        else:
            time_CP = -2
            nodes_CP = -2

        if CURR_CYCLE_TYPE in CycleTypeTimings_CPL[i].keys():
            time_CPL = int(np.median(CycleTypeTimings_CPL[i][CURR_CYCLE_TYPE][0]))
            nodes_CPL = np.median(CycleTypeTimings_CPL[i][CURR_CYCLE_TYPE][1])
        else:
            time_CPL = -2
            nodes_CPL = -2

        if CURR_CYCLE_TYPE in CycleTypeTimings_CPM[i].keys():
            time_CPM = int(np.median(CycleTypeTimings_CPM[i][CURR_CYCLE_TYPE][0]))
            nodes_CPM = np.median(CycleTypeTimings_CPM[i][CURR_CYCLE_TYPE][1])
        else:
            time_CPM = -2
            nodes_CPM = -2

        nodecounts = " & " + format(nodes_IP, ".2f") + " & " + format(nodes_CP, ".2f") + " & " if "-nodes" in sys.argv else ""
        print((str(i + 5) + " & " + str(CURR_CYCLE_TYPE) +
            " & {:,}".format(time_IP if time_IP < 60000 else -1) +
            " & {:,}".format(time_CPL if time_CPL < 60000 else -1) +
            " & {:,}".format(time_CPM if time_CPM < 60000 else -1) +
            " & {:,}".format(time_CP if time_CP < 60000 else -1) +
            nodecounts + " " + endl)
              .replace("-1", "Timeout")
              .replace("-2", "")
              .replace(" .00", ""))
print("\\end{tabular}")
