#import prerequisites
import sys
import csv
import fasta_reader
import kmer_script

#define mapping function

def gene(reads, k, score, scaffold, ol):

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
        
#go through first and last kmers
        for kmer in ((kmerDict[key])[0:ol]+(kmerDict[key])[(len(kmerDict[key])-ol):len(kmerDict[key])]):
            #print((kmerDict[key])[0:ol])
            #print((kmerDict[key])[(len(kmerDict[key])-ol):len(kmerDict[key])])
            seqname = []
            source = []
            ktest = []
            start = []
            end = []
            seqs = []
            boo = []
            bootlace = []
        
#go through each scaffold contig            
            for contig in scaffold.keys():
                #print(kmer)
            
#compare to scaffold
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
                        seqname += [contig]
                        source += [key]
                        ktest += [kmer]
                        start += [str(base)]
                        end += [str(base+k)]
                        seqs += [ATCG]
                        boo += [sum(boolean)]
                        bootlace += [str(sum(boolean))]
                        #print([str(base)])
                        
                      
                        #need to sum boolean across contigs
                
#write bases and seq for best match scores
            for match in range(0,len(start)):
                if boo[match]>=max(boo):
                    line = [''.join(seqname[match]),''.join(source[match]),''.join(ktest[match]),'transcript kmer',''.join(start[match]),''.join(end[match]),''.join(bootlace[match]),'.','.']
                    stringline = '\t'.join(line)
                    del line[-5:-1]
                    del line[-1]
                            #print(line)

#write line to map
                    maps += [stringline]  

#write all strings to table
    with open('gene_test.gff', "w") as outfile:
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
    gene(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]), sys.argv[4], int(sys.argv[5]))