# arguments: matrix file, sequence file, threshold for edges, threshold for discarding edges, output file name

#null hypothesis is the collapsed tree
#have to test once for each version of the collapsed tree -> three tests for each triple
#write a new fasta file each tree variant? or load all sequences into memory

# IMPORTS
import sys
import re
import os
import subprocess

# CLASSES
class Edge:
    def __init__(self,first,second):
	self.first = str(first)
        self.second = str(second)
    def toString(self):
        return self.first+"--"+self.second
    def __eq__(self, other):
        if (self.first == other.first and self.second == other.second):
            return True
        elif (self.first == other.second and self.second == other.first):
            return True
        else:
            return False
    def __cmp__(self, other):
        if self.__eq__(self,other):
            return 0
        else:
            return -1
# read in arguments
matrixFileName = sys.argv[1]
sequenceFileName = sys.argv[2]
edgeThreshold = float(sys.argv[3])
discardThreshold = float(sys.argv[4])
outputFileName = sys.argv[5]

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
numberPossible = 0
numEdges = 0;
allEdges = []
for i in range(1,len(matrixFile)): # skip the first line
    tokens = matrixFile[i].split("\t")
    firstVertex = sequenceIndex[i-1]
    if not firstVertex in firstVertexHash:
        firstVertexHash[firstVertex] = [];
    for j in range(i+1,len(tokens)):
        score = float(tokens[j].strip())
        numberPossible = numberPossible + 1
        if score < edgeThreshold: # if the distance between two sequences is less than the threshold, create an edge
            secondVertex = sequenceIndex[j-1]
            edge = Edge(firstVertex, secondVertex);
            allEdges.append(edge)
            firstVertexHash[firstVertex].append(edge)
            numEdges = numEdges + 1


# read in sequences
sequenceHash = {} #key: sequence name, value: sequence
with open(sequenceFileName) as f:
    sequenceFile = f.readlines()
f.close()
name = ""
sequence = ""
for line in sequenceFile:
    if not re.match(line,"^[\n\r]+"):
        if line.startswith(">"):
            if len(name) > 0 and len(sequence) > 0:
                sequenceHash[name[1:]] = sequence
                name = ""
                sequence = ""
            name = line.strip()
        else:
            sequence = sequence + line.strip()
sequenceHash[name[1:]] = sequence # write the last entry

# determine which sequences to test (find all triplets in graph)
triples = set()
for first in firstVertexHash.keys():
    for firstEdge in firstVertexHash[first]:
        for secondEdge in firstVertexHash[firstEdge.second]:
            for thirdEdge in firstVertexHash[firstEdge.first]:
                if secondEdge.second == thirdEdge.second:
                    triples.add(frozenset([firstEdge.first, firstEdge.second, secondEdge.second]))
# write temp file and then run test
vertexNames = ['A', 'B', 'C']
for triple in triples:
    # rotate which sequence is the second
    tripleList = list(triple)
    for i in range(3): # there are always three sequences
        currentFirstVertex = tripleList[0]
        # write the fasta file for this order
        f = open('temp.fasta','w')
        for j in range(3):
            sequenceName = tripleList[j]
            vertexName = vertexNames[j]
            f.write(">"+vertexName+"\n"+sequenceHash[sequenceName]+"\n")
        f.close()
        # call HyPhy to perform test, if we can reject collapsing A, then remove edge B-C
        pvalue = subprocess.check_output(["HYPHYMP","testTriple.bf"])
        if pvalue > discardThreshold:
            # try removing the edge in both orientations
            toRemove = Edge(str(tripleList[1]),str(tripleList[2]))
            try:
                allEdges.remove(toRemove)
            except:
                toRemove = Edge(str(tripleList[2]),str(tripleList[1]))
                try:
                    allEdges.remove(toRemove)
                except:
                    pass #edge must have been removed
        # reorder the triple 
        tripleList.remove(currentFirstVertex)
        tripleList.append(currentFirstVertex)

# construct adjacency matrix
adjacencyMatrix = []
for i in range(len(sequenceHash)):
    adjacencyMatrix.append(["0"]*len(sequenceHash))
for edge in allEdges:
    i = sequenceIndex.index(edge.first)    
    j = sequenceIndex.index(edge.second)
    adjacencyMatrix[i][j] = "1"
    adjacencyMatrix[j][i] = "1"

# write output file
f = open(outputFileName,'w')
f.write("\t".join(sequenceIndex))
for i in range(len(adjacencyMatrix)):
    f.write(sequenceIndex[i]+"\t"+"\t".join(adjacencyMatrix[i]))
f.close()
# print summary
print("Number of possible edges: " + str(numberPossible))
print("There were originally: " + str(numEdges) + " edges")
print("Now there are: " + str(len(allEdges)) + " edges")




