#import prerequisites
import sys
import csv
import fasta_reader
import kmer_script
import kmer_stream
import indexing
import reverse
import complement
import datetime

#define mapping function
def gene(reads, k, score, scaffold):

#set up new class
    class nsp:
        def __init__(self, ID, sequence,terminus):
            self.ID = ID
            self.sequence = sequence
            self.terminus = terminus

#establish dictionary of kmer transcripts
    reads = fasta_reader.readfa(reads)

    kmerDict = {}
    
    for key in reads:
        x = len(reads[key])
        y = (reads[key])[0:k]
        z = (reads[key])[x-k:x]
        #initalise keys for first and last kmers
        kmerDict[y] = []
        kmerDict[z] = []
        #calculate keys for new kmers
        values_y = nsp(key,reads[key],'front')
        values_z = nsp(key,reads[key],'back')
        #add values to keys
        kmerDict[y] += [values_y]
        kmerDict[z] += [values_z]
    
        #print(((kmerDict[y])[0]).ID)
    #print(kmerDict)

#set up scoring dictionary
    scoreDict = {}

#set up mapping list
    genemap = []
    
#go through kmers of each sequence in scaffold file
    indices = indexing.index(scaffold)
                
    for thingy in range(0,len(indices)-1):
                    
        entry = kmer_stream.stream(scaffold,indices,thingy)
        contig = entry[1]
        #print(entry[0])
        #print(entry[1])
        contig_r = reverse.reverse(contig)
        contig_c = complement.complement(contig)
        contig_r_c = complement.complement(contig_r)

#repeat calculations for reverse and for reverse complement

        for each in [contig,contig_r,contig_c,contig_r_c]:
            
            #print('NEW')

            if each == contig:
                name = 'contig'
            if each == contig_r:
                name = 'contig_r'
            if each == contig_c:
                name = 'contig_c'
            if each == contig_r_c:
                name = 'contig_r_c'

            #print(len(each))
            
            for base in range(0,len(each)-k+1):
                #print(base)
    #compare to terminal kmers from transcripts
                for key in kmerDict:
                    #print(key)
#calculate match score
                    boolean = []
                    ATCG = (each)[base:(base+k)]
                    #print(ATCG)
                    #print(key)

                    for char in range(0,len(ATCG)):
                        comp = [key[char]==ATCG[char]]
                        boolean += comp
                        #print(boolean)
                        #print(key)
                    #print(sum(boolean))    

                    if sum(boolean)>score:
                        #print(sum(boolean))
                        seqname = str(entry[0]+name)
                            
                        ambigs = len(kmerDict[key])
                        #print(ambigs)
                        for item in range(0,ambigs):
                            source = (kmerDict[key])[item].ID
                            
                        start = base+1
                        end = base+k
                        seqs = ATCG
                        boo = sum(boolean)
                        #print(boo)
                        #print(key)

#store best scores for each kmer
                    #is this better than the previous score?
                        if key in scoreDict:
                        #print(scoreDict[key])
                            if ((scoreDict[key])[4])<boo:
                                scoreDict[key] = [seqname,source,start,end,boo]
                        #is this the same as the previous best?        
                            if ((scoreDict[key])[4])==boo:
                                scoreDict[key] = [[seqname,(scoreDict[key])[0]],[source,(scoreDict[key])[1]],[start,(scoreDict[key])[2]],[end,(scoreDict[key])[3]],[boo,(scoreDict[key])[4]]]   
                        #if not, do nothing
                        #if key doesn't exist yet, put it in the dictionary                     
                        else:
                            #print(key)
                            scoreDict[key] = [seqname,source,start,end,boo]
    
    #print(scoreDict)
        
#pull out genes, put in new dictionary
        #searching scoreDict by source; let's assume single pair of kmers per tanscript for convenience
        #what to do about multiple matches???
    geneDict = {}
    for key in scoreDict:
        #pull source
        source = (scoreDict[key])[1]
    #    print(source)
        if source in geneDict:
            (geneDict[source])[1] = [(geneDict[source])[1],(scoreDict[key])[2]]
            (geneDict[source])[2] = [(geneDict[source])[2],(scoreDict[key])[3]]
            (geneDict[source])[3] = [(geneDict[source])[3],(scoreDict[key])[4]]
        else:
            geneDict[source] = [(scoreDict[key])[0],(scoreDict[key])[2],(scoreDict[key])[3],(scoreDict[key])[4]]
    #print(geneDict)
#we now have a dictionary hashed by transcript id, containing gff-relevant data for each id
#still need to separate values by whether they are start or end.. or do we?
    for key in geneDict:
        gene_seqname = (geneDict[key])[0]
        gene_source = key
        #want min and max values
        bases = (geneDict[key])[1]+(geneDict[key])[2]
        gene_start = min(bases)
        gene_end = max(bases)
        gene_score = max((geneDict[key])[3])
        #print(gene_start)
        #print(gene_end)
        #print(gene_score)

#write line describing gene
        gene_seqname = gene_seqname.replace('>','')
        gene_source = gene_source.replace('>','')
        gene = [str(gene_seqname),str(gene_source),'gene',str(gene_start),str(gene_end),str(boo),'.','.']
        genestring = '\t'.join(gene)

#write gene to map      
        genemap += [genestring]
        print([len(genemap),datetime.datetime.now()])
        
#write all strings to table
    with open('gene_map.gff', "w") as outfile:
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
        outfile.write("\n")
                
        for rows in genemap:
            outfile.write(rows)
            outfile.write("\n")
            
            
    
    return(genemap)

#allow parsing from command line
if __name__ == '__main__':
    gene(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])