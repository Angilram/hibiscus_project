import sys
import fasta_reader

def kmers(file,k):
    #read in fasta file
    FastaDict = fasta_reader.readfa(file)

    for fasta in FastaDict:
        code=FastaDict[fasta]
            
    #set up empty list and dictionary
        mers = list()
        kmerDict=FastaDict

    #loop through sequence
        for x in range(0,len(code)-k):
            mers.append(code[x:(x+k)])
        
    #save results to dictionary
        kmerDict[fasta] = mers

    #print results
    for key in kmerDict.keys():
        print(key)
        print(kmerDict[key])
        
    return(kmerDict)
    
if __name__ == '__main__':
    kmers(sys.argv[1], int(sys.argv[2]))