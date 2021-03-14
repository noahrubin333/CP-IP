import numpy as np

PATH_TO_CP_LOG = "../Logs/Cycle_Type_Log_CP_Z_LAST.log"

LOG_FILE = open(PATH_TO_CP_LOG, "r")

DATA_RAW = LOG_FILE.readlines()

LOG_FILE.close()

# Parse out logging timestamps and separators
for i in range(len(DATA_RAW)):
    if "LOGGING" in DATA_RAW[i] or "==" in DATA_RAW[i]:
        DATA_RAW[i] = ""

DATA_CLEANSED = "".join(DATA_RAW).split("---------------------------------\n")[:-1]



CYCLE_TYPE_TIMINGS = [{} for i in range(6)]


for LOG_ENTRY in DATA_CLEANSED:
    DATA_SPLIT = LOG_ENTRY.split("\n")
    LEN = len(DATA_SPLIT)

    n = int(DATA_SPLIT[0][3:])
    CYCLE_TYPE = tuple(map(int, DATA_SPLIT[2][13:].strip("]\n").split(",")))
    TIME = float(DATA_SPLIT[LEN - 4][14:].strip("\n")) / 1000.0 #ms to s
    COUNT_BRANCHES = float(DATA_SPLIT[LEN - 3][29:].strip("\n"))

    if tuple(CYCLE_TYPE) in CYCLE_TYPE_TIMINGS[n - 5]:
        CYCLE_TYPE_TIMINGS[n - 5][CYCLE_TYPE][0].append(TIME)
        CYCLE_TYPE_TIMINGS[n - 5][CYCLE_TYPE][1].append(COUNT_BRANCHES)
    else:
        CYCLE_TYPE_TIMINGS[n - 5].update({tuple(CYCLE_TYPE): [[TIME], [COUNT_BRANCHES]]})

PATH_TO_IP_LOG = "../Logs/Cycle_Type_Log_MIP_Z_LAST.log"

f = open(PATH_TO_IP_LOG, "r")

s = f.read().split("---------------------------------\n")[:-1]

f.close()

CycleTypeTimings = [{} for i in range(6)]


for data in s:
    if data.startswith("n = "):
        data_split = data.split("\n")
        n = int(data_split[0][3:])
        CycleType = list(map(int, data_split[2][13:].strip("]\n").split(",")))
        time = float(data_split[len(data_split) - 3][10:].strip("\n"))
        nodes = int(float(data_split[len(data_split) - 2][17:].strip("\n")))

        if n != 6:
            if data_split[len(data_split) - 5].split(" ")[0] == 'Error':
                continue


        if tuple(CycleType) in CycleTypeTimings[n - 5]:
            CycleTypeTimings[n - 5][tuple(CycleType)][0].append(time)
            CycleTypeTimings[n - 5][tuple(CycleType)][1].append(nodes)
        else:
            CycleTypeTimings[n - 5].update({tuple(CycleType) : [[time], [nodes]]})



print("\\begin{tabular}{c@{\;\,}l@{\quad}r@{\quad}r@{\quad}r@{\quad}r}\n$n$ & \\df{Cycle Type} & \\df{IP Time} & \\df{CP Time} & \\df{IP Nodes} & \\df{CP Branches} \\\\ \\hline")

for i in range(6):

    #for CURR_CYCLE_TYPE in CYCLE_TYPE_TIMINGS[i].keys():
    for CURR_CYCLE_TYPE in CycleTypeTimings[i].keys():

        if "1, 1" in str(CURR_CYCLE_TYPE):
            continue

        if CURR_CYCLE_TYPE in CYCLE_TYPE_TIMINGS[i].keys():
            TIME_LIST = CYCLE_TYPE_TIMINGS[i][CURR_CYCLE_TYPE][0]
            COUNT_BRANCH_LIST = CYCLE_TYPE_TIMINGS[i][CURR_CYCLE_TYPE][1]

            AVG_TIME = int(np.median(TIME_LIST))
            AVG_NODE = int(np.median(COUNT_BRANCH_LIST))
        else:
            AVG_TIME = -1
            AVG_NODE = -1

        time = int(np.median(CycleTypeTimings[i][CURR_CYCLE_TYPE][0]))
        nodes = int(np.median(CycleTypeTimings[i][CURR_CYCLE_TYPE][1]))
        
        if i == 1:
            time = int(time/1000)
       
        endl = " \\\\" if CURR_CYCLE_TYPE != list(CycleTypeTimings[5])[-1] else ""

        print(str(i+5) + " & " + str(CURR_CYCLE_TYPE) + " & {:,}".format(time) + " & {:,}".format(AVG_TIME) + " & {:,}".format(nodes) + " & {:,}".format(AVG_NODE) + endl)

print("\\end{tabular}")
