from collections import defaultdict

mapped = defaultdict(int)
with open("magic_blast_reference_output.txt", 'r') as blastfile:
	for line in blastfile:
		if line[0] != '#':
			mapped[line.split('\t')[0][:-2]] = 1

seqs = defaultdict(str)
outfile = open("testout.fasta", 'w')

with open("SRR7445584.fasta", 'r') as fastafile:
	for line in fastafile:
		if line[0] == '>':
			curRead = line.split(' ')[0][1:]
		elif (mapped[curRead] != 1):
			outfile.write(">" + curRead + '\n' + line + '\n')
outfile.close()	
