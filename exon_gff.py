#import prerequisites
import sys
import csv
import fasta_reader
import kmer_script
import gene_gff

#define function
def exons(gene_id, reads, scaffold, k, K, score):

#parse gene gff output
    genemap = gene_gff.gene(reads, K, score, scaffold)
    #print(genemap)
    
    for lines in genemap:
        lines = lines.split('\t')
        #print(lines)
        if lines[1] == gene_id:
            bounds = [int(lines[3])-1,int(lines[4])-1]
        #    print(bounds)
            contig = lines[0]
        #    print(contig)

#read in scaffold
    scaffold = fasta_reader.readfa(scaffold)

#extract gene from scaffold
    gene = (scaffold['>'+contig])[bounds[0]:bounds[1]+1]
    #print(gene)

#calculate kmers for transcript
    kmerDict = kmer_script.kmers(reads,k)
    transcripts = kmerDict['>'+gene_id]
    #print(transcripts)
    
    match_start = []
    match_end = []
    exon = []
    exonmap = []

#take each transcript kmer in turn
    for kmer in transcripts:
        #print(kmer)
    
        #take each gene kmer in turn
        for base in range(0,len(gene)-k):
            #print(base)
        
            #if there is a match...
            if gene[base:base+k]==kmer:
            
                #...extract start and end
                match_start += [base]
                #print(base)
                match_end += [base+k]
                
                #if start != previous match's start + 1...
                if len(match_start)>1:
                    #print(len(match_start))
                
                    if match_start[-1] != (match_start[-2]+1):            
                    
                        #define start and end of old exon
                        exon_end = str(match_end[-2])
                        exon_start = str(match_start[0])
                    
                        #write to line
                        exon = [''.join(contig),''.join(gene_id),'exon',''.join(exon_start),''.join(exon_end),'.','.','.']
                        exonstring = '\t'.join(exon)
                        exonmap += [exonstring]                    
                
                        #reset scores
                        match_start = [base]
                        match_end = [base+k]
                        
    #define start and end of old exon
    exon_end = str(match_end[-1])
    exon_start = str(match_start[0])
                    
    #write to line
    exon = [''.join(contig),''.join(gene_id),'exon',''.join(exon_start),''.join(exon_end),'.','.','.']
    exonstring = '\t'.join(exon)
    exonmap += [exonstring] 
                
    #write exon starts and ends to gff

    with open('exon_map.gff', "w") as outfile:
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
            
        for rows in exonmap:
            outfile.write(rows)
            outfile.write("\n")           
    
    #return(maps)
    return(exonmap)

#allow parsing from command line
if __name__ == '__main__':
    exons(sys.argv[1], sys.argv[2], sys.argv[3], int(sys.argv[4]), int(sys.argv[5]), int(sys.argv[6]))