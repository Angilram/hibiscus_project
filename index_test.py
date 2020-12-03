import sys

def index(file,x):

    OpenedFile = open(file, 'r')
    lino=-1
    indices = []
#go through every line
    for line in OpenedFile:
        lino = lino+1
        line = line.strip("\n\r")
        #print(line)
    
#stop at sequence name
        if line [0] == ">":
        
            indices += [lino]
        
    indices += [lino]    
#we now have a record of where every sequence starts and stops  
    del(OpenedFile)
    print(indices)
#pull specific sequence using indices
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
    print(SeqZ)
    return(SeqZ)
    
if __name__ == '__main__':
    index(sys.argv[1],int(sys.argv[2]))