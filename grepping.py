import sys
import os
import complement
import reverse
import fasta_reader
import regex as re

test = 'dog'
#print(re.search("(dog){s<=3}", "cat and dog"))
print(re.search("(%s){s<=3}" % (test), "cat and dog"))

#take a given transcript
def krep(transcript,scaffold):
 
    #read in transcript file
    d = fasta_reader.readfa(transcript)
    for key in d:
        read = d[key]
    print(read)

    #reformat scaffold file
    FastaDict = fasta_reader.readfa(scaffold)
    with open('formatted.fa', "w") as outfile:
        for key in FastaDict:
            outfile.write(key)
            outfile.write("\n")
            outfile.write(FastaDict[key])
            outfile.write("\n")

    print("DONE")

    #set up list for match frequency
    matches = []

    for i in range(1,50):
        print("Thank u, NEXT")
        matches.append(0)
        #apply grep to kmers of increasing length
        kmer = read[0:i]
        kmer_r = reverse.reverse(kmer)
        kmer_c = complement.complement(kmer)
        kmer_r_c = reverse.reverse(kmer_c)
        file = open("formatted.fa", "r")
        for line in file:
            #print('NEW')
            if re.search('(%s){e<=7}' % (kmer), line):
                matches[(i-1)] = matches[(i-1)]+1
            if re.search('(%s){e<=7}' % (kmer_r), line):
                matches[(i-1)] = matches[(i-1)]+1
            if re.search('(%s){e<=7}' % (kmer_c), line):
                matches[(i-1)] = matches[(i-1)]+1
            if re.search('(%s){e<=7}' % (kmer_r_c), line):
                matches[(i-1)] = matches[(i-1)]+1
    #ditto their reverse/complement/etc
    print(matches)
    k = matches.index(1)
    k = k+1
    print(k)
        
#return kmer length(s) that give exactly one match
    return(k)

if __name__ == '__main__':
    krep(sys.argv[1],sys.argv[2])