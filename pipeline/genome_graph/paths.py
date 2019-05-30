from utils import reverse_complement
from copy import deepcopy
#from pipeline.genome_graph.genome_graph import *

class Path:
    def __init__(self,g,nodeId):
        self.nodes = [g.nodes[abs(nodeId)]]
        self.nodeIds = [nodeId]
        self.nNodes = 1
    
    def extend_right(self,g):
        lastNode = self.nodeIds[-1]
        neighbors = g.edges[lastNode].copy()
        extendedPaths = set()

        if len(neighbors)>0:
            for neighbor in neighbors:
                if abs(neighbor) not in [abs(n) for n in self.nodeIds]: # Node has not been visited by path
                    rev_neighbors = g.edges[-neighbor]
                    assert -lastNode in rev_neighbors
                    newP = deepcopy(self)
                    newP.nodes.append(g.nodes[abs(neighbor)]) 
                    newP.nodeIds.append(neighbor)
                    newP.nNodes += 1
                    extendedPaths.add(newP)
        if len(extendedPaths) > 0:
            return(extendedPaths)
        else:
            return(False)
            

    def extend_linear_right(self,g):
        lastNode = self.nodeIds[-1]
        neighbors = g.edges[lastNode].copy()
        if len(neighbors)==1:
            neighbor = neighbors.pop()
            if abs(neighbor) not in [abs(n) for n in self.nodeIds]: # Node has not been visited by path
                rev_neighbors = g.edges[-neighbor]
                if rev_neighbors == {-lastNode}:   
                    self.nodes.append(g.nodes[abs(neighbor)]) 
                    self.nodeIds.append(neighbor)
                    self.nNodes += 1
                    return(True)        
        return(False)

    def extend_linear_left(self,g):
        firstNode = self.nodeIds[0]
        neighbors = g.edges[-firstNode].copy()
        if len(neighbors)==1:
            neighbor = neighbors.pop()
            if abs(neighbor) not in [abs(n) for n in self.nodeIds]: # Node has not been visited by path
                rev_neighbors = g.edges[-neighbor] 
                if rev_neighbors == {firstNode}: # In the case of linear Paths, we also check there is no branching in the opposite direction
                    self.nodes.insert(0,g.nodes[abs(neighbor)]) 
                    self.nodeIds.insert(0,-neighbor)
                    self.nNodes += 1
                    return(True)
        return(False)

    def extend_left(self,g):
        # If the path is extendable,  returns a set containing all the extended paths
        firstNode = self.nodeIds[0]
        neighbors = g.edges[-firstNode].copy()
        extendedPaths = set()
        if len(neighbors)>0:
            for neighbor in neighbors:
                if abs(neighbor) not in [abs(n) for n in self.nodeIds]: # Node has not been visited by path
                    rev_neighbors = g.edges[-neighbor]
                    assert firstNode in rev_neighbors
                    newP = deepcopy(self)
                    newP.nodes.insert(0,g.nodes[abs(neighbor)]) 
                    newP.nodeIds.insert(0,-neighbor)
                    newP.nNodes += 1
                    extendedPaths.add(newP)
        if len(extendedPaths) > 0:
            return(extendedPaths)
        else:
            return(False)

    def getSeq(self,g):
        seq = ''
        for node in self.nodeIds:
            nodeSeq = g.nodes[abs(node)].nodeSeq
            if node < 0:
                nodeSeq = reverse_complement(nodeSeq)
            if seq != '':
                nodeSeq = nodeSeq[g.overlap:]    
            seq = seq + nodeSeq
        return(seq)

    def getName(self,g):
        name = ''
        rev = False
        for node in self.nodeIds:
            nodeName = g.nodes[abs(node)].nodeName
            if node < 0:
                rev = not rev
            if rev:
                nodeName = nodeName + "_Rc"
            if name != '':
                nodeName = "_" + nodeName   
            name = name + nodeName
        return(name)

    def merge(self,g):
        newName = self.getName(g)
        newSeq = self.getSeq(g)
        g.add_node(newName,newSeq)
        newId = g.maxId

        neighbors_left = g.get_neighbors(-self.nodeIds[0])
        for n in neighbors_left:
            g.add_edge(-newId,n)

        neighbors_right = g.get_neighbors(self.nodeIds[-1])
        for n in neighbors_right:
            g.add_edge(newId,n)
        
        for node in self.nodeIds:
            g.rem_node(abs(node))

            
def setExtend(paths,g):
# extends a set of paths in both directions
    extendedPaths = set()
    for p in paths:
        extension = p.extend_right(g)
        if extension != False:  # path has been extended
            extendedPaths.update(extension)
        else: 
            extendedPaths.add(p) # We keep the unextended path

    # We do the same on the left side
    paths = deepcopy(extendedPaths)
    extendedPaths = set()
    for p in paths:
        extension = p.extend_left(g)
        if extension != False:  # path has been extended
            extendedPaths.update(extension)
        else:
            extendedPaths.add(p) # We keep the unextended path
    
    return(extendedPaths)

