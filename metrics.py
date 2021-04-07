#python3 metrics.py Tc_genes.gff unannotated.txt Tc_rna.fa Tc_metrics.csv

import sys
import fasta_reader

#read in mapped genes
genes = open(sys.argv[1],'r')

#need seqname, feature, start, and end
class Gene:
	def __init__(self):
		
		self.seqname = []
		self.feature = []
		self.start = []
		self.end = []
		
	def _repr(self):

		return repr((self.seqname,self.feature,self.start,self.end))

#create dictionary of objects by key characteristic 'seqname'
mDict={}

for line in genes:
	line = line.strip("\n\r")
	line = line.split("\t")

	try:
		gene=Gene()
		gene.seqname = line[0]
		gene.feature = line[2]
		gene.start = float(line[3])
		gene.end = float(line[4])

	except:
		pass

#let each corresponding item be the full object
	if gene.seqname in mDict:
		(mDict[gene.seqname]).append(gene)
	else:
		mDict[gene.seqname] = [gene]

#generate dictionary with keys 'seqname'...
nDict = {}

for key in mDict:

#...and items from length of of items in mDict
	nDict[key] = len(mDict[key])

#calculate total number of genes
total = 0
record = []
for key in nDict:
	total += nDict[key]
	record.append(nDict[key])

#calculate distances between genes
alldist = []

#for each key in mDict, extract item
alllen = []

for entry in mDict:
	raw = mDict[entry]

#sort list by characteristic 'start'

	try:
		tidy = sorted(raw, key = lambda i: i.start)
	except:
		for thing in raw:
			if type(thing.start) != float:
				print(thing.start)
				print(type(thing.start))
				print(thing.end)
				print(thing.feature)

#flag overlapping genes
	count1 = 0
	newtidy = tidy

	for each in newtidy:

		if count1 != 0:
			if legacy <= each.start: #if end of previous gene comes before the start of the current gene
				legacy = each.end #record the end of the new gene

#combine those items to give new 'start' and 'end'
			else:
				combo = Gene()
				combo.seqname = entry
				combo.feature = ((newtidy[count1-1]).feature + "/" + (newtidy[count1]).feature)
				combo.start = (newtidy[count1-1]).start
				combo.end = (newtidy[count1]).end
				newtidy[count1-1] = "Oops"
				newtidy[count1] = combo
				legacy = combo.end #record end of combined gene

#nb: each == newtidy[count]
		else:
			legacy = each.end #record end of first gene
		
		count1 += 1

	newtidy = list(filter(("Oops").__ne__,newtidy))
	newlen = len(newtidy)
	alllen.append(newlen)

#iterate through new list
	count2 = 0
	dist = 0

	for each in newtidy:
		if count2 != 0:

#subtract 'end' from new 'start'
			temp = (newtidy[count2]).start-(newtidy[count2-1]).end
			dist += temp
			
			if temp < 0:
				print(temp)

		else:
			pass		
		
		count2 += 1

#divide sum by iteration count
	avdist = dist/count2
	alldist.append(avdist)

#calculate average across all contigs iyl
totavdist = sum(alldist)/len(mDict)

#read in length of unannotated genes
unann = open(sys.argv[2],'r')
unlen = -1

for line in unann:
	unlen += 1

#read in length of transcriptome fasta dictionary
trans = fasta_reader.readfa(sys.argv[3])
tralen = len(trans)

#simple division!
propun = (unlen/tralen)*100

#writing to file
totav = str(total/len(mDict))
total = str(total)
for i in range(len(record)):
	record[i] = str(record[i])

totgen = str(sum(alllen))
avgen = str(sum(alllen)/len(mDict))
alllen = str(alllen)

totavdist = str(totavdist)
for i in range(len(alldist)):
	alldist[i] = str(alldist[i])

propun = str(propun)

output = open(sys.argv[4],'w')
output.write("Mapped transcripts in genome"+"\n"+total+"\n"+"Average mapped transcripts per contig"+"\n"+totav+"\n"+"Transcripts mapped to each contig"+"\n"+",".join(record)+"\n")
output.write("Mapped genes in genome"+"\n"+totgen+"\n"+"Average mapped genes per contig"+"\n"+avgen+"\n"+"Genes mapped to each contig" + ",".join(record)+"\n")
output.write("Average intergenic distance in genome"+"\n"+totavdist+"\n"+"Average intergenic distance in each contig"+"\n"+",".join(alldist)+"\n")
output.write("Percentage of transcripts unmapped"+"\n"+propun)