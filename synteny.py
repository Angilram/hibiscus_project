#import prerequisites
import sys
import csv
import fasta_reader
import kmer_script
import gene_sonic
import length

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
        genemap = gene_sonic.gene(reads, K, score, scaffold)
        #print(genemap)
        
        #map transcript names to locus start
        synDict = {}
        for lines in genemap:
            #matching bases in reverse to normal contigs
            if '_r' in lines[0]: 
                temp=lines[0].split("contig_r")
                temp=temp[0]
                length = length.length(scaffold,temp)
                lines[3] = length - int(lines[3])
                lines[4] = length - int(lines[4])
            #print(lines)
            lines = lines.split('\t')
            synDict[lines[3]]=lines[1]
            
        #need to sort genes
        synDict_items = synDict.items()
        sorted_items = sorted(synDict.items)
    
#write values to file
        line = list(sorted_items())
        #print(line)
        stringline = '\t'.join(line)
        #print(stringline)
        synmap += [stringline]
        
        with open('synteny_map.csv', "w") as outfile:
            for rows in synmap:
                outfile.write(rows)
                outfile.write("\n")
            
    return(synmap)

if __name__ == '__main__':
    synteny(sys.argv[1], int(sys.argv[2]), int(sys.argv[3]))