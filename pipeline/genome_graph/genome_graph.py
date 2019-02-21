'''
Implementation of genome graphs.

Data structure : 
g.nodes : dict to store node attributes (seq, name)

g.edges : adjacency list of the nodes. 
       A negative nodeId indicates that the contig is reversed. In that case, 
       the source is the start of the contig, and the target is the end.

'''

import re


class GenomeNode:

       def __init__(self,nodeSeq,nodeName):
              # self.nodeId = nodeId
              self.nodeSeq = nodeSeq
              self.nodeName = nodeName
              # self.active = True  # Alternative way to remove nodes ?
       
       def __eq__(self, other):
              if isinstance(other, GenomeNode):
                     return self.nodeSeq == other.nodeSeq and self.nodeName == other.nodeName 
              return False

       def __ne__(self, other):
              return not self.__eq__(other)

       def __hash__(self):
              return hash((self.nodeSeq, self.nodeName))

class GenomeGraph:

       def __init__(self):
              self.nodes = {}
              self.edges = {}
              self.overlap = 0

       def nNodes(self):
              return(len(self.nodes))
       
       def nEdges(self):
              return(len([j for i in self.edges.values() for j in i])/2)

       def add_node(self,nodeName,nodeSeq):
              newNode = GenomeNode(nodeSeq,nodeName)
              try:
                     assert newNode not in self.nodes.values()
              except AssertionError:
                     print("Node already in graph")
                     
              if len(self.nodes)>0:
                     nodeId = max(self.nodes.keys())+1
              else: 
                     nodeId = 1
              self.nodes[nodeId] = newNode
              self.edges[nodeId] = set()
              self.edges[-nodeId] = set()

       def add_edge(self, src,dst):
              self.edges[src].add(dst)
              self.edges[-dst].add(-src)

       def rem_edge(self,src,dst):
              try:
                     self.edges[src].remove(dst)
                     self.edges[-dst].remove(-src)
              except KeyError:
                     print("Edge not in graph")
              
       def rem_node(self,nodeId):
              try:
                     self.nodes.pop(nodeId)
              except KeyError:
                     print("Node not in graph")

              for n in self.get_neighbors(nodeId).copy():
                     self.rem_edge(-n,-nodeId)

              for n in self.get_neighbors(-nodeId).copy():
                     self.rem_edge(-n,nodeId)

              self.edges.pop(nodeId)
              self.edges.pop(-nodeId)

       def get_neighbors(self,nodeId):
              return(self.edges[nodeId])

           
       @classmethod
       def read_gfa(self,file):
              
              g = GenomeGraph()
              nodeIds = {}  # Used to retrieve a node Id from a name. Is it useful though?
              nlines = 0
              with open(file) as f:
                     for line in f:
                            nlines += 1
                            if re.match(r"S.*", line):
                                   nodeName = line.split("\t")[1]
                                   nodeSeq = line.split("\t")[2]
                                   g.add_node(nodeName,nodeSeq)
                                   nodeIds[nodeName] = g.nNodes()
                            elif re.match(r"L.*",line):
                                   startName,startDir,endName,endDir,overlap = line.split("\t")[1:]

                                   if g.overlap == 0:
                                          g.overlap = overlap
                                   
                                   startNode = nodeIds[startName]
                                   endNode = nodeIds[endName]

                                   if startDir == "-":
                                          startNode = - startNode
                                   if endDir == "-":
                                          endNode = - endNode

                                   g.add_edge(startNode,endNode)                                  
              return(g)


