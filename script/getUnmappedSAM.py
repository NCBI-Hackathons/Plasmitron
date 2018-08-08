from collections import defaultdict
import sys

mapped = defaultdict(int)
with open(sys.argv[1], 'r') as samfile:
	for line in samfile:
		if (line[0] != '@') & (line.split('\t')[2] != '*'):
			mapped[line.split('\t')[0]] = 1

print(len(mapped))
seqs = defaultdict(str)
outfile = open(sys.argv[3], 'w')

with open(sys.argv[2], 'r') as fastafile:
	for line in fastafile:
		if line[0] == '>':
			curRead = line.split(' ')[0][1:]
		elif (mapped[curRead] != 1):
			outfile.write(">" + curRead + '\n' + line + '\n')
outfile.close()	
