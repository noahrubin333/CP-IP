import sys
import numpy as np

PATH_TO_CP_LOG = "C:\\Users\\Noah\\Desktop\\Random_Assignment_Log_CP.log"

if "-linux" in sys.argv:
    PATH_TO_CP_LOG = "../Logs/Random_Assignment_Log_CP.log"

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


AVG_TIME = STD_DEV = 0.0

for i in range(6):
    print("\n\n\\begin{table}\n\\begin{center}\n\\begin{tabular}{c@{\quad}r@{\quad}r@{\quad}r@{\quad}r}\n\\df{Cycle Type} & \\df{Avg. Time} & \\df{Std Dev} & \\df{Avg. Branches} & \df{Runs} \\\\ \\hline")

    for CURR_CYCLE_TYPE in CYCLE_TYPE_TIMINGS[i].keys():

        TIME_LIST = CYCLE_TYPE_TIMINGS[i][CURR_CYCLE_TYPE][0]
        COUNT_BRANCH_LIST = CYCLE_TYPE_TIMINGS[i][CURR_CYCLE_TYPE][1]

        AVG_TIME = np.round(np.average(TIME_LIST), 8 - i)
        STD_DEV = np.round(np.std(TIME_LIST), 8 - i)

        print(str(CURR_CYCLE_TYPE) + " & " + str(AVG_TIME) + " & " + str(STD_DEV) + " & " + str(np.round(np.average(COUNT_BRANCH_LIST), 2)) + " & " + str(len(TIME_LIST)) + " \\\\")

    print("\\end{tabular}\n\\end{center}\n\\caption{CP timings for all cycle types for $n=" + str(i + 5) + "$.\label{cp_times_" + str(i + 5) + "}}\n\\end{table}")
