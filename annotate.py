import sys

print(sys.argv[1])
print(sys.argv[2])
print(sys.argv[3])
print(sys.argv[4])


transfile = open(sys.argv[1], 'r')
deetfile = open(sys.argv[2],'r')
notfile = open(sys.argv[3],'r')

modfile = open(sys.argv[4], 'w')

annDict = {}
for note in deetfile:
	#print(note)
	note = note.strip("\n\r")
	npy = note.split('\t')
	annDict[npy[0]] = npy[2]

#print(annDict)

arghDict = {}

for line in transfile:
	line = line.strip("\n\r")
	lpy = line.split('\t')
	ortho = lpy[2]
	#print(ortho)
	#print(lpy[0])
	if lpy[0] in annDict:
		#print('yes')
		arghDict[ortho] = str(annDict[lpy[0]])

print(arghDict)

for each in notfile:
	each = each.strip("\n\r")
	epy = each.split('\t')
	if epy[2] in arghDict:
		annote = arghDict[epy[2]]
		modfile.write(each + '\t' + str(annote)+'\n')
	elif epy[0]=='seqname':
		annote = "attributes"
		modfile.write(each + '\t' + str(annote)+'\n')
	else:
		annote = "Not mapped to H. trionum"
		modfile.write(each + '\t' + str(annote)+'\n')