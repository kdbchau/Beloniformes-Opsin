#consensus.py

import sys

#rev_comp = {'A':'T', 'T':'A', 'C':'G', 'G':'C'}
seqs = {}
pileup_file = sys.argv[1]
all_lines = []

for line in open(pileup_file):
	all_lines.append(line.split())
#print(all_lines)
for i in all_lines:
	ID,pos,ref,num_reads,bases,quals,check = i
	if ID not in seqs:
		seqs[ID] = ''
	if num_reads == 0:
		seqs[ID] += "-"
	else:
		counts = {'A':0, 'C':0, 'G':0, 'T':0, "N":0}
		for c in bases:
			if c in {'.',','}:
				counts[ref] += 1
			elif c in {'A', 'C', 'G', 'T', 'N'}:
				counts[c] += 1
			elif c in {'a', 'c', 'g', 't'}:
				counts[c.upper()] += 1
		if sum(counts.values()) == 0:
			seqs[ID] += "-"
		else:
			seqs[ID] += max(counts, key=lambda key: counts[key])


species = seqs.keys()[0]
sequence = seqs.values()[0]

with open(pileup_file.split(".mpileup")[0]+"_consensus.fa", "w") as outfile:
	outfile.write(">" + species + "\n")
	outfile.write(sequence)
