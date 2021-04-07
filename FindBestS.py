import sys
import fasta_reader
import datetime

#a script to run post-BLAST

#1: makeblastdb -in C_rna.fa -out C_rna.fa -dbtype 'nucl'

#2: blastn -db C_rna.fa -query C_gen.fa -evalue 1e-3 -perc_identity 95 -num_threads 6 -max_target_seqs 20000 -out Cannabinus.rawblast -outfmt '6 qseqid qlen sseqid slen frames pident nident length mismatch gapopen qstart qend sstart send evalue bitscore'

#defining hit properties
class BlastHit:
	def __init__(self):
		
		self.genome_hits = [] #hit in genome
		self.percent_identity = [] #percent identity
		self.length = [] #match length
		self.start = [] #match start
		self.end = [] #match end
		self.score = []

	#defining perc_iden within class to sort matches for a given transcript
	def perc_iden(self):

		owzat = [i for i, x in enumerate(self.length) if int(x)<100] #selecting indices of short matches
		temp = []
		excluded = []

		temp = [x for i,x in enumerate(self.score) if i not in owzat]

		sorted_scores = sorted(temp,reverse=True) #sorting all matches by percent identity
		pos_of_interest = []
		x = []

		for i in range(5):
			try:
				x.append(sorted_scores[i])
				pos_of_interest.append(self.score.index(sorted_scores[i]))
			except:
				pass
		
		best_genome_hits = []
		for j in range(5):
			try:
				best_genome_hits.append(self.genome_hits[pos_of_interest[j]]) #record best hit in the genome
			except:
				pass

		pos_of_best_genome_hits = [i for i, x in enumerate(self.genome_hits) if x in best_genome_hits] #record its index
		#nb assuming "hit" refers to a physical location in the genome, handily groups transcripts

		#Step 1: make dictionary of contigs hit
		shortlist = {}
		for contig in pos_of_best_genome_hits:
			key = self.genome_hits[contig]

		#Step 2: add values corresponding to %id
			try:
				shortlist[key].append(float(self.percent_identity[contig]))
			except:
				shortlist[key] = [float(self.percent_identity[contig])]

		#Step 3: sum id% for each contigs
		try:
			best_score = 0
			for key in shortlist:
				shortlist[key] = sum(shortlist[key])
				if shortlist[key] > best_score:
					best_score = shortlist[key]
					winner = key
			#Step 4: extract top-scoring			
			pos_of_best = [i for i, x in enumerate(self.genome_hits) if x == winner]
		except:
			pos_of_best = 'ERROR'
		return pos_of_best


#opening file
#if len(sys.argv) != 2:
#	print(("python " + sys.argv[0] + " BlastFile"))
#	sys.exit()
File = open(sys.argv[1],"r")

#Step 1: filter by the transcriptome hits to the genome, use this to find what the most likely transcriptome is

#0 is transcriptome, 2 is genome, 5 percent identity, 7 length, 10 is qstart, 11 is qend, 12 is sstart, 13 send
#count = 0
HASH = {} #dictionary that will have the transcript and the BLAST hit object
#reading in the opened file
check = []
for line in File:
	line = line.strip("\r\n")
	array = line.split("\t")
	check.append(array[0])
	Blast = BlastHit()
	#recording progress
	#count += 1
	#sys.stderr.write("processing line " + str(count) + '\r')

	#Associate the BlastHit object with its corresponding transcript

	#first testing whether we can append results to array in dictionary
	try:
		HASH[array[0]].genome_hits.append(array[2])
		HASH[array[0]].percent_identity.append(array[5])
		HASH[array[0]].length.append(array[7])
		HASH[array[0]].start.append(array[12])
		HASH[array[0]].end.append(array[13])
		score = float(array[5])*int(array[7])
		HASH[array[0]].score.append(score)
	
	#if we can't - set Blast as a BlastHit object, then append!
	except:
		Blast = BlastHit()
		HASH[array[0]] = Blast
		HASH[array[0]].genome_hits.append(array[2])
		HASH[array[0]].percent_identity.append(array[5])
		HASH[array[0]].length.append(array[7])
		HASH[array[0]].start.append(array[12])
		HASH[array[0]].end.append(array[13])
		score = float(array[5])*int(array[7])
		HASH[array[0]].score.append(score)

	#print(count,datetime.datetime.now())

print(len(HASH))

#Step 2: Figure out which is the best genome hit for the transcriptome:

with open(sys.argv[3], "w") as outfile:
	outfile.write('seqname')
	outfile.write('\t')
	outfile.write('source')
	outfile.write('\t')
	outfile.write('feature')
	outfile.write('\t')
	outfile.write('start')
	outfile.write('\t')
	outfile.write('end')
	outfile.write('\t')
	outfile.write('score')
	outfile.write('\t')
	outfile.write('strand')
	outfile.write('\t')
	outfile.write('frame')
outfile = open(sys.argv[3],"a")

with open('unannotated.txt',"w") as oopsfile:
	oopsfile.write("Unannotated")
	oopsfile.write('\n')
oopsfile = open('unannotated.txt',"a")

count = 0
errors = 0

for transcript in HASH:
	try:
		array_points = []
		#by retrieving the indices of their best hits in the genome
		array_points = HASH[transcript].perc_iden()
		#count += 1
		#if 10 < count: 
		#	sys.exit()

		#use the positions to summarize location and hit info
		starts = []
		ends = []
		length = []
		perc_iden = []
		genome_seq = ""

		#going through each saved match index
		for hit in array_points:
			#saving data to lists
			genome_seq = HASH[transcript].genome_hits[hit] 
			starts.append(HASH[transcript].start[hit])
			ends.append(HASH[transcript].end[hit])
			length.append(HASH[transcript].length[hit])
			perc_iden.append(HASH[transcript].percent_identity[hit])

		#printing each match to log file
		for i in range(0,len(perc_iden)):
			yay = '#'.join(("exon",transcript))
			#print(yay)
			outfile.write('\n')
			outfile.write(genome_seq)
			outfile.write('\t')
			outfile.write('TGF')
			outfile.write('\t')
			outfile.write(yay)
			outfile.write('\t')
			outfile.write(starts[i])
			outfile.write('\t')
			outfile.write(ends[i])
			outfile.write('\t')
			outfile.write(perc_iden[i])
			outfile.write('\t')
			outfile.write('.')
			outfile.write('\t')
			outfile.write('.')
			outfile.write('\t')

		
			#print((genome_seq + "\t" + transcript + "\t" + starts[i] + "\t" + ends[i] + "\t" + length[i] + "\t" + perc_iden[i]))
		count += 1

	except:
		error = 'error'
		#print("ERROR")
		oopsfile.write(transcript)
		oopsfile.write('\n')
		errors +=1

extras = fasta_reader.readfa(sys.argv[2])
for key in extras:
	#print(key)
	K = key.split(' ')
	K = K[0]
	k = K.strip('>')
	#print(k)
	if k not in HASH:
		oopsfile.write(k)
		oopsfile.write('\n')
		errors += 1

print(count)
print(errors)