import sys
import indexing

def stream(file,indices,x):

#pull specifc sequence using indices
    if x != (len(indices)-1):
        #print(indices)
        name = list(open(file, 'r'))[indices[x]]
        name = name.strip("\n\r")
        a = indices[x]+1
        b = indices[x+1]
        SeqX = list(open(file, 'r'))[a:b]
        SeqY = []
        for each in SeqX:
            each = each.strip("\n\r")
            SeqY += each
            SeqY = ''.join(SeqY)
        #print(SeqY)
        SeqZ = [name,SeqY]
        #print(SeqZ)
    else:
        SeqZ = 0
    
    return(SeqZ)
        
if __name__ == '__main__':
    stream(sys.argv[1],sys.argv[2],int(sys.argv[3]))