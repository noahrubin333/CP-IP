import sys
import numpy as np

PATH_TO_MIP_LOG = "C:\\Users\\Noah\\Desktop\\Random_Assignment_Log.log"

if "-linux" in sys.argv:
    PATH_TO_MIP_LOG = "../Logs/Random_Assignment_Log_MIP.log"

f = open(PATH_TO_MIP_LOG, "r")

s = f.read().split("---------------------------------\n")[:-1]

f.close()

CycleTypeTimings = [{} for i in range(6)]


for data in s:
    if data.startswith("n = "):
        data_split = data.split("\n")
        n = int(data_split[0][3:])
        CycleType = list(map(int, data_split[2][13:].strip("]\n").split(",")))
        time = float(data_split[len(data_split) - 3][10:].strip("\n"))

        if n != 6:
            if data_split[len(data_split) - 5].split(" ")[0] == 'Error':
                continue


        if tuple(CycleType) in CycleTypeTimings[n - 5]:
            CycleTypeTimings[n - 5][tuple(CycleType)].append(time)
        else:
            CycleTypeTimings[n - 5].update({tuple(CycleType) : [time]})


for i in range(6):
    print("\n\n\\begin{table}\n\\begin{center}\n\\begin{tabular}{c@{\quad}r@{\quad}r@{\quad}r}\n\\df{Cycle Type} & \\df{Avg Time} & \\df{Std Dev.} & \\df{Runs} \\\\ \\hline")
    for key in CycleTypeTimings[i].keys():
        if (i + 5) == 5:
            time = np.round(np.average(CycleTypeTimings[i][key]), 10)
            variance = np.round(np.std(CycleTypeTimings[i][key]), 10)
        elif (i + 5) == 6:
            time = np.round(np.average(CycleTypeTimings[i][key]), 7)
            variance = np.round(np.std(CycleTypeTimings[i][key]), 7)
        elif (i + 5) == 7:
            time = np.round(np.average(CycleTypeTimings[i][key]), 4)
            variance = np.round(np.std(CycleTypeTimings[i][key]), 4)
        else:
            time = np.round(np.average(CycleTypeTimings[i][key]), 2)
            variance = np.round(np.std(CycleTypeTimings[i][key]), 2)
        print(str(key) + " & " + str(time) + " & " + str(variance) + " & " + str(len(CycleTypeTimings[i][key])) + " \\\\")
    print("\\end{tabular}\n\\end{center}\n\\caption{IP timings for all cycle types for $n=" + str(i + 5) + "$.\label{ip_times_" + str(i + 5) + "}}\n\\end{table}")
