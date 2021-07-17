import numpy as np
CycleTypeTimings_CP = [{} for i in range(7)]
CycleTypeTimings_CPL = [{} for i in range(7)]
CycleTypeTimings_IP = [{} for i in range(7)]
NRDTimings_CP = [[] for i in range(8)]
NRDTimings_IP = [[] for i in range(8)]
NRDTimings_CPL = [[] for i in range(8)]

PATH_TO_CP_LOG = "C:\\Users\\Noah\\Desktop\\CP-IP\\Logs\\Cycle_Type_Log_CP.log"
PATH_TO_IP_LOG = "C:\\Users\\Noah\\Desktop\\CP-IP\\Logs\\Cycle_Type_Log_MIP.log"
PATH_TO_CPL_LOG = "C:\\Users\\Noah\\Desktop\\Cycle_Type_Log_CPL.log"
PATH_TO_NRD_LOG_CP = "C:\\Users\\Noah\\Desktop\\CP-IP\\Logs\\NRD_Log_CP.log"
PATH_TO_NRD_LOG_IP = "C:\\Users\\Noah\\Desktop\\CP-IP\\Logs\\NRD_Log_MIP.log"
PATH_TO_NRD_LOG_CPL = "C:\\Users\\Noah\\Desktop\\CP-IP\\Logs\\NRD_Log_CP_lin.log"

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

LOG_FILE_CPL = open(PATH_TO_CPL_LOG, "r")
DATA_RAW = LOG_FILE_CPL.readlines()
LOG_FILE_CPL.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    CYCLE_TYPE = tuple(map(int, DATA_SPLIT[2][13:].strip(")\n").split(",")))
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n"))  # ms to s

    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 3][29:].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    if tuple(CYCLE_TYPE) in CycleTypeTimings_CPL[n - 5]:
        CycleTypeTimings_CPL[n - 5][CYCLE_TYPE][0].append(TIME)
        CycleTypeTimings_CPL[n - 5][CYCLE_TYPE][1].append(COUNT_BRANCHES)
    else:
        CycleTypeTimings_CPL[n - 5].update({tuple(CYCLE_TYPE): [[TIME], [COUNT_BRANCHES]]})


NRD_LOG_FILE_IP = open(PATH_TO_NRD_LOG_IP, "r")
DATA_RAW = NRD_LOG_FILE_IP.readlines()
NRD_LOG_FILE_IP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 3][10:].strip("\n"))  # ms to s
    SYMM_METHOD = DATA_SPLIT[1][19:].strip("\n")
    if SYMM_METHOD != "D": continue
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 2][len('Nodes Processed: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    NRDTimings_IP[n - 5] = [TIME, COUNT_BRANCHES]

NRD_LOG_FILE_CP = open(PATH_TO_NRD_LOG_CP, "r")
DATA_RAW = NRD_LOG_FILE_CP.readlines()
NRD_LOG_FILE_CP.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n")) / 1000.0  # ms to s
    SYMM_METHOD = DATA_SPLIT[1][19:].strip("\n")
    if SYMM_METHOD != "D": continue
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 3][len('Number of branches explored: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)


    NRDTimings_CP[n - 5] = [TIME, COUNT_BRANCHES]

NRD_LOG_FILE_CPL = open(PATH_TO_NRD_LOG_CPL, "r")
DATA_RAW = NRD_LOG_FILE_CPL.readlines()
NRD_LOG_FILE_CPL.close()
DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")

for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)
    n = int(DATA_SPLIT[0][3:])
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n")) / 1000.0  # ms to s
    SYMM_METHOD = DATA_SPLIT[1][19:].strip("\n")
    if SYMM_METHOD != "D": continue
    COUNT_BRANCHES_RAW = float(DATA_SPLIT[LEN - 3][len('Number of branches explored: '):].strip("\n"))
    COUNT_BRANCHES = np.round(np.log2(COUNT_BRANCHES_RAW), 2) if COUNT_BRANCHES_RAW != 0 else np.round(0.00, 2)

    NRDTimings_CPL[n - 5] = [TIME, COUNT_BRANCHES]



endl = "\\\\"
print(
    "\\begin{tabular}{c@{\;\,}c@{\quad}r@{\quad}r@{\quad}r@{\quad}r@{\quad}r@{\quad}r}\n$n$ & \\df{Symm. Breaking} & \\df{IP Time} & \\df{CP Time} & \\df{CPL Time} & \\df{IP Nodes} & \\df{CP Nodes} & \\df{CPL Branches} \\\\ \\hline")

for i in range(7):
    print((str(i + 5) + " & Dominance & " + str(int(NRDTimings_IP[i][0]) if NRDTimings_IP[i][0] < 60000 else "-1")
           + " & " + str(int(NRDTimings_CP[i][0]) if NRDTimings_CP[i][0] < 60000 else "-1")
           + " & " + str(int(NRDTimings_CPL[i][0]) if NRDTimings_CPL[i][0] < 60000 else "-1")
           + " & " + format(NRDTimings_IP[i][1], ".2f")
           + " & " + format(NRDTimings_CP[i][1], ".2f")
           + " & " + format(NRDTimings_CPL[i][1], ".2f") + "\\\\").replace("-1", "T"))

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
            time_CPL = int(np.median(CycleTypeTimings_CPL[i][CURR_CYCLE_TYPE][0]) / 1000.0)
            nodes_CPL = np.median(CycleTypeTimings_CPL[i][CURR_CYCLE_TYPE][1])
        else:
            time_CPL = -2
            nodes_CPL = -2

        print((str(i + 5) + " & " + str(CURR_CYCLE_TYPE) + " & {:,}".format(time_IP if time_IP < 60000 else -1) + " & {:,}".format(
            time_CP if time_CP < 60000 else -1) + " & " + "{:,}".format(time_CPL if time_CPL < 60000 else -1) + " & " + format(nodes_IP, ".2f") + " & " + format(nodes_CP, ".2f")
               + " & " + format(nodes_CPL, ".2f") + endl)
              .replace("-1", "T")
              .replace("-2", "")
              .replace(" .00", ""))
print("\\end{tabular}")