import sys

def length(file,contig):

    OpenedFile = open(file, 'r')

    seq = ''
    #print('>'+contig)
    for line in OpenedFile:

        line = line.strip("\n\r")
        #print(line)
        
#only recognise
        if line == '>'+contig:
            #print('found it')
            continue

        if line[0] != '>' and line[0] != [' ']:
            #print(line)
            seq = seq+line

    length = len(seq)
    print(seq)
    print(length)

    return(length)
    
if __name__ == '__main__':
    length(sys.argv[1],sys.argv[2])