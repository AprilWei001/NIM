import sys, random

if ( len(sys.argv) == 2 ): random.seed( sys.argv[1] )
try:
	lines = sys.stdin.readlines()
	random.shuffle(lines)
	for line in lines: print (line, end=''),
except:
	lines
