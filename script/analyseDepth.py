from collections import defaultdict
import sys

covered = defaultdict(int)

with open(sys.argv[1], 'r') as depthfile:
	for line in depthfile:
		line = line.strip().split('\t')
		covered[line[0]] += 1

lens = defaultdict(int)
with open(sys.argv[2], 'r') as samfile:
	for line in samfile:
		if line[0:3] == "@SQ":
			line = line.split()
			curID = line[1][3:]
			lens[curID] = int(line[2][3:])

outfile = open(sys.argv[3], 'w')
outfile.write("Name\tLength\tPercent Covered\n")
for x in covered:
	outfile.write(x + '\t' + str(lens[x]) + '\t' + str(float(covered[x])/lens[x]*100) + '\n')	
