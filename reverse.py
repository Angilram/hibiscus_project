import sys

def reverse(sequence):

     flipped = sequence[::-1]
     print(flipped)
     return flipped

if __name__ == '__main__':
    reverse(sys.argv[1])