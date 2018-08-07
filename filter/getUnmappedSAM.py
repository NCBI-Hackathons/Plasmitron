from collections import defaultdict

mapped = defaultdict(int)
with open("chromAlign.sam", 'r') as samfile:
	for line in samfile:
		if line[0] != '@':
			mapped[line.split('\t')[0]] = 1

seqs = defaultdict(str)
outfile = open("testout.fasta", 'w')

with open("SRR7445584.fasta", 'r') as fastafile:
	for line in fastafile:
		if line[0] == '>':
			curRead = line.split(' ')[0][1:]
		elif (mapped[curRead] != 1):
			outfile.write(">" + curRead + '\n' + line + '\n')
outfile.close()	
