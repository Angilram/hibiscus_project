import sys

def readfa(file):

    OpenedFile = open(file, 'r')

    FastaDict = {}
    for line in OpenedFile:

        line = line.strip("\n\r")
    
#only recognise
        if line [0] == ">":
        
            FastaDict[line] = ""
            name = line  
            
        else:
            FastaDict[name] += line
    return FastaDict
    print(FastaDict)
    
if __name__ == '__main__':
    readfa(sys.argv[1])