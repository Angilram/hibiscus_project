import sys

def index(file):

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
    return(indices)

if __name__ == '__main__':
    index(sys.argv[1])