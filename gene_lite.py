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

#calculate kmers
    reads = fasta_reader.readfa(reads)
    kmerDict = reads
    for each in kmerDict:
        x = len(kmerDict[each])
        kmerDict[each]=[(kmerDict[each])[0:k],(kmerDict[each])[x-k:x]]
    
    #print(kmerDict)
    
#set up mapping list(s)
    genemap = []
    
#take each key's values
    for key in kmerDict.keys():
        #print(key)
        #print([(kmerDict[key])[0]]+[(kmerDict[key])[len(kmerDict[key])-1]])
        
#go through first and last kmers
        for kmer in [(kmerDict[key])[0],(kmerDict[key])[len(kmerDict[key])-1]]:
            #print(kmer)
            seqname = []
            source = []
            ktest = []
            start = []
            end = []
            seqs = []
            boo = []
            bootlace = []
            
#compare to the kmers of each sequence in scaffold file
            indices = indexing.index(scaffold)
                
            for thingy in range(0,len(indices)-1):
                    
                entry = kmer_stream.stream(scaffold,indices,thingy)
                contig = entry[1]
                    
                for base in range(0,len(contig)-k):
                    
                    #calculate match score
                    boolean = []
                    ATCG = (contig)[base:(base+k)]
                    #print(ATCG)
                    for char in range(0,len(ATCG)):
                        #print(char)
                        #print(kmer[char])
                        #print(ATCG[char])
                        comp = [kmer[char]==ATCG[char]]
                        boolean += comp
                        #print(boolean)

                    if sum(boolean)>score:
                     
                        #print(kmer)
                        seqname += [entry[0]]
                        source += [key]
                        ktest += [kmer]
                        start += [str(base+1)]
                        end += [str(base+k+1)]
                        seqs += [ATCG]
                        boo += [sum(boolean)]
                        bootlace += [str(sum(boolean))]
                        #print([str(base)])
                
                #repeat calculations for reverse and for reverse complement
                contig = reverse.reverse(contig)
                
                for base in range(0,len(contig)-k):
                    
                    #calculate match score
                    boolean = []
                    ATCG = (contig)[base:(base+k)]
                    #print(ATCG)
                    for char in range(0,len(ATCG)):
                        #print(char)
                        #print(kmer[char])
                        #print(ATCG[char])
                        comp = [kmer[char]==ATCG[char]]
                        boolean += comp
                        #print(boolean)

                    if sum(boolean)>score:
                     
                        #print(kmer)
                        seqname += [entry[0]]
                        source += [key+'REVERSE']
                        ktest += [kmer]
                        start += [str(base+1)]
                        end += [str(base+k+1)]
                        seqs += [ATCG]
                        boo += [sum(boolean)]
                        bootlace += [str(sum(boolean))]
                
                contig = complement.complement(contig)
                
                for base in range(0,len(contig)-k):
                    
                    #calculate match score
                    boolean = []
                    ATCG = (contig)[base:(base+k)]
                    #print(ATCG)
                    for char in range(0,len(ATCG)):
                        #print(char)
                        #print(kmer[char])
                        #print(ATCG[char])
                        comp = [kmer[char]==ATCG[char]]
                        boolean += comp
                        #print(boolean)

                    if sum(boolean)>score:
                     
                        #print(kmer)
                        seqname += [entry[0]]
                        source += [key+'REVERSE_COMPLEMENT']
                        ktest += [kmer]
                        start += [str(base+1)]
                        end += [str(base+k+1)]
                        seqs += [ATCG]
                        boo += [sum(boolean)]
                        bootlace += [str(sum(boolean))]
                
                contig = reverse.reverse(contig)
                
                for base in range(0,len(contig)-k):
                    
                    #calculate match score
                    boolean = []
                    ATCG = (contig)[base:(base+k)]
                    #print(ATCG)
                    for char in range(0,len(ATCG)):
                        #print(char)
                        #print(kmer[char])
                        #print(ATCG[char])
                        comp = [kmer[char]==ATCG[char]]
                        boolean += comp
                        #print(boolean)

                    if sum(boolean)>score:
                     
                        #print(kmer)
                        seqname += [entry[0]]
                        source += [key+'COMPLEMENT']
                        ktest += [kmer]
                        start += [str(base+1)]
                        end += [str(base+k+1)]
                        seqs += [ATCG]
                        boo += [sum(boolean)]
                        bootlace += [str(sum(boolean))]

#write bases and seq for best match scores for given kmer
            if len(boo) != 0:
                for match in range(0,len(start)):
                
                    if boo[match] == max(boo):
                        if kmer == (kmerDict[key])[0]: #basically, if this is the start
                            gene_start = start[match]
                            boo_start = boo[match]
                        if len(boo_start) != 0:                    
                            if kmer == (kmerDict[key])[len(kmerDict[key])-1]: #basically, if this this is the end
                                gene_end = end[match]
                                gene_seqname = seqname[match]
                                gene_source = source[match]
                                avboo = str((boo[match]+boo_start)/2)
        
#write line describing gene
        if len(boo) != 0:
            gene_seqname = gene_seqname.replace('>','')
            gene_source = gene_source.replace('>','')
            gene = [''.join(gene_seqname),''.join(gene_source),'gene',''.join(gene_start),''.join(gene_end),''.join(avboo),'.','.']
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
        
        #for rows in maps:
        #    outfile.write(rows)
        #    outfile.write("\n")
        
        for rows in genemap:
            outfile.write(rows)
            outfile.write("\n")
            
            
    
    #return(maps)
    return(genemap)

#allow parsing from command line
if __name__ == '__main__':
    gene(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])