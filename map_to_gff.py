#import prerequisites
import sys
import csv
import fasta_reader
import kmer_script

#define mapping function

def gff(reads, k, score, scaffold):

#read in scaffold
    scaffold = fasta_reader.readfa(scaffold)

#calculate kmers
    kmerDict = kmer_script.kmers(reads,k)
    #print(kmerDict)
    
#set up mapping list
    maps = []
    
#take each key's values
    for key in kmerDict.keys():
        #print(key)
        
#go through each kmer in list
        for kmer in kmerDict[key]:
        
#go through each scaffold contig            
            for contig in scaffold:
                line = [contig,key,'transcript kmer',]
                start = []
                end = []
                seqs = []
                boo = []
                #print(kmer)
            
#compare each kmer to scaffold
                for each in range(0,round(len(scaffold[contig])/k)):
                #print(scaffold[base+k])
                    base = each*k
                
                    #calculate match score
                    boolean = []
                    ATCG = (scaffold[contig])[base:(base+k)]
                    #print(ATCG)
                    for char in range(0,len(ATCG)):
                        #print(char)
                        comp = [kmer[char]==ATCG[char]]
                        boolean += comp
                #print(boolean)

                    if sum(boolean)>score:
                     
                        #print(kmer)
                        start += [str(base)]
                        end += [str(base+k)]
                        seqs += [ATCG]
                        boo += [str(sum(boolean))]
                        #print([str(base)])
                
#write bases and seq to kmer line
                for match in range(0,len(start)):
                    line += [''.join(start[match]),''.join(end[match]),''.join(boo[match]),'.','.']
                    stringline = '\t'.join(line)
                    del line[-5:-1]
                    del line[-1]
                    #print(line)

#write line to map
                    maps += [stringline]  

#write all strings to table
    with open('test.gff', "w") as outfile:
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
        
        for rows in maps:
            outfile.write(rows)
            outfile.write("\n")
    
    return(maps)

#allow parsing from command line
if __name__ == '__main__':
    gff(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4])