import sys
import csv
import fasta_reader
import kmer_script
import kmer_stream
import indexing
import reverse
import complement
import datetime

def matchmake(kmerDict, k, scaffold, indices, thingy):

    scoreDict = {}
                    
    entry = kmer_stream.stream(scaffold,indices,thingy)
    contig = entry[1]
    #print(entry[0])
    #print(entry[1])
    contig_r = reverse.reverse(contig)
    contig_c = complement.complement(contig)
    contig_r_c = complement.complement(contig_r)
    print(contig_r_c)

    #repeat calculations for reverse and for reverse complement

    for each in [contig,contig_r,contig_c,contig_r_c]:

        print('NEW')
        print(each)
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
                ATCG = (each)[base:(base+k)]
                #print(ATCG)
                #print(key)

                if key == ATCG:

                    seqname = str(entry[0]+name)
                        
                    ambigs = len(kmerDict[key])
                    #print(ambigs)
                    for item in range(0,ambigs):
                        source = (kmerDict[key])[item].ID
                        
                    start = base+1
                    end = base+k
                    seqs = ATCG
                    boo = k
                    #print(key)

        #store best scores for each kmer
                    if key in scoreDict:
                    #print(scoreDict[key])
                        scoreDict[key] = [[seqname,(scoreDict[key])[0]],[source,(scoreDict[key])[1]],[start,(scoreDict[key])[2]],[end,(scoreDict[key])[3]],[boo,(scoreDict[key])[4]]]   
                    #if key doesn't exist yet, put it in the dictionary                     
                    else:
                        #print(key)
                        scoreDict[key] = [seqname,source,start,end,boo]
    
    print(scoreDict)
    return(scoreDict)