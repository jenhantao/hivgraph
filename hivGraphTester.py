# arguments: matrix file, sequence file, threshold for edges, threshold for discarding edges

#null hypothesis is the collapsed tree
#have to test once for each version of the collapsed tree -> three tests for each triple
#write a new fasta file each tree variant? or load all sequences into memory

# IMPORTS
import sys
import re

# CLASSES
class Edge:
    def __init__(self,first,second):
	self.first = str(first)
        self.second = str(second)
    def toString(self):
        return self.first+"--"+self.second

# read in arguments
matrixFileName = sys.argv[1]
sequenceFileName = sys.argv[2]
edgeThreshold = sys.argv[3]
discardThreshold = sys.argv[4]

# read in matrix
with open(matrixFileName) as f:
    matrixFile = f.readlines()
f.close()
# store the order of sequences
sequenceIndex = []
for name in matrixFile[0].split("\t"):
    sequenceIndex.append(name.strip())

# make the edges of the graph
firstVertexHash = {} # key: the first vertex listed in an edge, value: an array list of edges with given first vertex
secondVertexHash ={} # key: the second vertex listed in an edge, value: an array list of edges with given second vertex
numberPossible = 0
for i in range(1,len(matrixFile)): # skip the first line
    tokens = matrixFile[i].split("\t")
    for j in range(i+1,len(tokens)):
        score = float(tokens[j].strip())
        numberPossible = numberPossible + 1
        if score < edgeThreshold: # if the distance between two sequences is less than the threshold, create an edge
            firstVertex = sequenceIndex[i-1]
            secondVertex = sequenceIndex[j-1]
            edge = Edge(firstVertex, secondVertex);
            if not firstVertex in firstVertexHash:
                firstVertexHash[firstVertex] = [];
            firstVertexHash[firstVertex].append(edge)
            if not secondVertex in secondVertexHash:
                secondVertexHash[secondVertex] = [];
            secondVertexHash[secondVertex].append(edge)            
      
# read in sequences
sequenceHash = {} #key: sequence name, value: sequence
with open(sequenceFileName) as f:
    sequenceFile = f.readlines()
f.close()
name = ""
sequence = ""
for line in sequenceFile:
    if not re.match(line,"[\n\r]+"):
        if line.startswith(">"):
            if len(name) > 0 and len(sequence) > 0:
                sequenceHash[name] = sequence
                name = ""
                sequence = ""
            name = line.strip()
        else:
            sequence = sequence + line.strip()

# determine which sequences to test (find all triplets in graph)





# print summary
print("Number of possible edges: " + str(numberPossible))
print("There were originally: " + str(len(firstVertexHash.keys())) + " edges")
