import sys

mapfile = open(sys.argv[1], 'r')
deetfile = open(sys.argv[2],'r')

outfile = open(sys.argv[3], "a")
outfile.write('seqname'+'\t'+'source'+'\t'+'feature'+'\t'+'start'+'\t'+'end'+'\t'+'score'+'\t'+'strand'+'\t'+'frame'+'\t'+'attributes'+'\n')

annDict = {}
for note in deetfile:

	note = note.strip("\n\r")
	npy = note.split('\t')
	annDict[npy[0]] = npy[2]

count = 0
for line in mapfile:
	count += 1
	line = line.strip("\n\r")
	lpy = line.split('\t')
	ARGH = (lpy[2]).split('#')
	if ARGH[-1] in annDict:
		outfile.write(line + '\t' + str(annDict[ARGH[-1]])+'\n')