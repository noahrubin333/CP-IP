# Check that outputted squares are actually k-MOLS(n)

import sys

fixcol = False
nosym = False

if "-f" in sys.argv:
	fixcol = True
	sys.argv.remove("-f")

if "-n" in sys.argv:
	nosym = True
	sys.argv.remove("-n")

if len(sys.argv) < 3:
	print "Need order of MOLS and # of squares"
	quit()

n = int(sys.argv[1])
k = int(sys.argv[2])

A = [[[] for i in range(n)] for p in range(k)]

linecount = 0
for line in sys.stdin:
	if line == "\n":
		continue

	A[linecount/n][linecount%n] = map(int, line.strip().split(" "))

	linecount += 1

	if linecount == n*k:
		break

# Check squares are Latin
for p in range(k):
	for i in range(n):
		assert(sorted([A[p][i][j] for j in range(n)])==range(n))

	for j in range(n):
		assert(sorted([A[p][i][j] for i in range(n)])==range(n))

# Check orthogonality
for p in range(k):
	for q in range(p+1, k):
		assert(sorted([[A[p][i][j], A[q][i][j]] for i in range(n) for j in range(n)])==[[x,y] for x in range(n) for y in range(n)])

# Check entries of first row and column
if nosym == False:
	for i in range(n):
		for p in range(k):
			assert(A[p][0][i] == i)
		assert(A[0][i][0] == i)
		if k >= 2 and (i in range(2,n-2) or i==n-1):
			assert(A[1][i][0] in range(3,i) or A[1][i][0] in {1,i+1})
	if k >= 2:
		assert(A[1][1][0] == 2)
		assert(A[1][n-2][0] == n-1)
	if fixcol and k >= 2:
		for i in range(1,n-1):
			assert(A[1][i][0] == i+1)
		assert(A[1][n-1][0] == 1)
