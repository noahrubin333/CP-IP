# Script to generate all possibilities for first column of second square (up to symmetry breaking)
import sys

if len(sys.argv) < 2:
	print "Need order of square"
	quit()

n = int(sys.argv[1])

def partitions(n):
	if n==0:
		return []
	r = []
	r.append([n])
	for i in range(1,n):
		for x in partitions(n-i):
			r.append(sorted(x+[i],reverse=True))
	
	s = set([])
	for x in r:
		s.add(tuple(x))
	
	t = []
	for x in s:
		t.append(list(x))

	return sorted(t)

def partitions_with_exactly_one_1(n):
	r = []
	for x in partitions(n):
		if x[-1] == 1 and (len(x) == 1 or x[-2] != 1):
			r += [x]
	return r

def min_undef(x):
	for i in range(n):
		if x[i] == -1:
			return i
	assert(False)

def representative(x):
	w = [-1 for i in range(n)]
	for i in range(len(x),0,-1):
		l = x[i-1]
		li = []
		while len(li) < l:
			m = min_undef(w)
			w[m] = 0
			li.append(m)

		w[li[-1]] = li[0]
		for i in range(len(li)-1):
			w[li[i]] = li[i+1]
		
	return w

for x in map(representative, partitions_with_exactly_one_1(n)):
	print x

