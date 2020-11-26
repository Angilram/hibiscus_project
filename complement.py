import sys

def complement(sequence):

    flipped = []
    
    for base in sequence:
        if base == 'A':
            baze = 'T'
        if base == 'T':
            baze = 'A'
        if base == 'G':
            baze = 'C'
        if base == 'C':
            baze = 'G'
        if base == 'N':
            baze == 'N'
        flipped += baze
    
    flipped = ''.join(flipped)
    
    #print(flipped)
    return(flipped)

if __name__ == '__main__':
    complement(sys.argv[1])
            
        