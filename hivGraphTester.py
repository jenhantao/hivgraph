# arguments: matrix file, sequence file, threshold for edges, threshold for discarding edges

#null hypothesis is the collapsed tree
#have to test once for each version of the collapsed tree -> three tests for each triple
#write a new fasta file each tree variant? or load all sequences into memory

# IMPORTS
import sys

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

# read in matrix to make edges
with open(matrixFileName) as f:
    matrixFile = f.readlines()
f.close()
# store the order of sequences
sequenceIndex = []
for name in matrixFile[0].split("\t"):
    sequenceIndex.append(name.strip())

firstVertexHash = {} # key: the first vertex listed in an edge, value: an array list of edges with given first vertex
secondVertexHash ={} # key: the second vertex listed in an edge, value: an array list of edges with given second vertex

for i in range(1,len(matrixFile)):
    tokens = matrixFile[i].split("\t")
    for j in range(i+1,len(tokens)-1):
	score = float(tokens[j].strip())
        if score < edgeThreshold: # if the distance between two sequences is less than the threshold, create an edge
            firstVertex = sequenceIndex[i]
            secondVertex = sequenceIndex[j]
            edge = Edge(firstVertex, secondVertex);
            if not firstVertex in firstVertexHash:
                firstVertexHash[firstVertex] = [];
            firstVertexHash[firstVertex].append(edge)
            if not secondVertex in secondVertexHash:
                secondVertexHash[secondVertex] = [];
            secondVertexHash[secondVertex].append(edge)            
      
# read in sequences
sequenceDict = {} #key: sequence name, value: sequence
with open(sequenceFileName) as f:
    sequenceFile = f.readlines()


