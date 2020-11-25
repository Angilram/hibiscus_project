#import prerequisites
import sys
import csv
import fasta_reader
import kmer_script
import gene_gff

#define function for a single 
def synteny(genomes, K, score):

    genomes = genomes.split(',')
    temp = []
    for pair in genomes:
        pair = pair.split('#')
        temp += [pair]
    genomes = temp
    #print(genomes)
        
    synmap = []

#find genes in each genome
    for each in genomes:
        scaffold = each[0]
        reads = each[1]
        genemap = gene_gff.gene(reads, K, score, scaffold)
        #print(genemap)
        
        #map transcript names to locus start
        synDict = {}
        for lines in genemap:
            #print(lines)
            lines = lines.split('\t')
            synDict[lines[1]]=lines[3]
            
        #gene function already sorts by locus start (within a contig)
        #print(synDict)
    
#write values to file
        line = list(synDict.keys())
        #print(line)
        stringline = ','.join(line)
        #print(stringline)
        synmap += [stringline]
        
        with open('synteny_map.csv', "w") as outfile:
            for rows in synmap:
                outfile.write(rows)
                outfile.write("\n")
            
    return(synmap)

if __name__ == '__main__':
    synteny(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))