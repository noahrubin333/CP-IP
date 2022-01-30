# Script to generate all possibilities for first column of second square (up to symmetry breaking)
import sys

if len(sys.argv) < 2:
	print "Need order of square"
	quit()

n = int(sys.argv[1])

P = [[] for i in range(max(4,n))]
P[0] = [[0]]
P[1] = [[0,1]]
P[2] = [[0,2,1]]
P[3] = [[0,2,3,1]]

for i in range(4,n):
	P[i] = map(lambda x: x + [i, i-1], P[i-2]) + map(lambda x: x[:-1] + [i] + x[-1:], P[i-1])

for x in P[n-1]:
	print x
