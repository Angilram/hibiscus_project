import sys

exonfile = open(sys.argv[1], 'r')

outfile = open(sys.argv[2], "w")
outfile.write('seqname'+'\t'+'source'+'\t'+'feature'+'\t'+'start'+'\t'+'end'+'\t'+'score'+'\t'+'strand'+'\t'+'frame'+'\t'+'attributes'+'\n')

exonmap = {}

class Gene:
	def __init__(self):
		
		self.seqname = []
		self.source = []
		self.feature = [] #match length
		self.start = [] #match start
		self.end = [] #match end
		self.score = []
		self.strand = []
		self.frame = []
		self.attributes = []
		self.length = []

#test = Gene()
#test.start = 'x'
#print(test.S())

exonmap={}
#create summary for each gene
for line in exonfile:

	#print(line)

	line = line.strip("\n\r")
	line = line.split('\t')

	if line[0] == 'seqname':
		continue

	else:
		#print(line[0])
		expy = Gene()
		l = line[2].split('#')
		line[2] = l[1]
		#print(line[3])
		expy.seqname = line[0]
		expy.source = line[1]
		expy.feature = line[2]
		expy.start = int(line[3])
		expy.end = int(line[4])
		expy.score = float(line[5])
		expy.strand = line[6]
		expy.frame = line[7]
		try:
			expy.attributes = line[9]
		except:
			pass
		expy.length = abs(int(line[3])-int(line[4]))
		#print(line[9])

		#adding to dictionary
		try:
			exonmap[line[2]].seqname = expy.seqname
			exonmap[line[2]].source = expy.source
			exonmap[line[2]].feature = expy.feature
			exonmap[line[2]].start.append(expy.start)
			exonmap[line[2]].end.append(expy.end)
			exonmap[line[2]].score.append(expy.score)
			exonmap[line[2]].strand = expy.strand
			exonmap[line[2]].frame = expy.frame
			try:
				exonmap[line[2]].attributes = expy.attributes
			except:
				pass
			exonmap[line[2]].length.append(expy.length)

		except:
			gene = Gene()
			exonmap[line[2]] = gene
			exonmap[line[2]].seqname = expy.seqname
			exonmap[line[2]].source = expy.source
			exonmap[line[2]].feature = expy.feature
			exonmap[line[2]].start.append(expy.start)
			exonmap[line[2]].end.append(expy.end)
			exonmap[line[2]].score.append(expy.score)
			exonmap[line[2]].strand = expy.strand
			exonmap[line[2]].frame = expy.frame
			try:
				exonmap[line[2]].attributes = expy.attributes
			except:
				pass
			exonmap[line[2]].length.append(expy.length)

		continue

print(len(exonmap))

ls=[]

#determine gene starts and ends
for key in exonmap:
	
	#calculate "start" and "end"
	bases = exonmap[key].start+exonmap[key].end
	START = min(bases)
	END = max(bases)
	exonmap[key].start = START
	exonmap[key].end = END
	
	#weight scores by length
	d = sum(exonmap[key].length)
	SCORE = []
	for i in range(len(exonmap[key].score)):
		s = exonmap[key].score[i]*exonmap[key].length[i]
		SCORE.append(s)
	exonmap[key].score = sum(SCORE)/d

	ls.append(exonmap[key])

#reorder dictionary
messymap = {}

def S(g):
	return g.start

#split into sublists by contig
for item in ls:
	try:
		messymap[item.seqname].append(item)
	except:
		messymap[item.seqname] = [item]

#print(messymap)

#sort those sublists by start

sls = []

for entry in messymap:
	if len(messymap[entry])!=1:
		tls = messymap[entry]
		tls.sort(key=S)
		#print(tls)
		# then concatenate
		sls = sls + tls
	else:
		tls = messymap[entry]
		# then concatenate
		sls = sls + tls

#for i in sls:
#	print(i.seqname)
#	print(i.start)
#	print('\n')

#add to dictionary
try:
	for tidy in sls:
		outfile.write(tidy.seqname+'\t'+tidy.source+'\t'+tidy.feature+'\t'+str(tidy.start)+'\t'+str(tidy.end)+'\t'+str(tidy.score)+'\t'+tidy.strand+'\t'+tidy.frame+'\t'+tidy.attributes+'\n')
except:
	for tidy in sls:
		outfile.write(tidy.seqname+'\t'+tidy.source+'\t'+tidy.feature+'\t'+str(tidy.start)+'\t'+str(tidy.end)+'\t'+str(tidy.score)+'\t'+tidy.strand+'\t'+tidy.frame+'\n')