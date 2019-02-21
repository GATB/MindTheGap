import unittest

from pipeline.genome_graph.genome_graph import *


class GenomeGraphTests(unittest.TestCase):

    def test_read_gfa(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        assert(g.nNodes()==28)
        assert(g.nEdges()==36)
        assert(g.edges[11]=={2})
        assert(g.edges[13]=={-1})
    def test_add_node(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")

        g.add_node("newNode","newSeq")
        assert(g.nNodes()==29)
        g.add_edge(1,29)
        assert(29 in g.get_neighbors(1))

    def test_rem_node(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        g.rem_node(1)

        assert(1 not in g.nodes.keys())
        assert(1 not in g.get_neighbors(11))

    def test_add_edge(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        g.add_edge(1,2)
        assert(-1 in g.get_neighbors(-2))
        assert(2 in g.get_neighbors(1))

    def test_rem_edge(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        g.rem_edge(1,11)
        assert(1 not in g.get_neighbors(11))
        assert(11 not in g.get_neighbors(1))     

    def test_in_nodes(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        assert(g.nodes[1] in g.nodes.values())
        assert(GenomeNode('ATCG','node') not in g.nodes.values() )

    def test_simple_bubles_removal(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        g.pop_all_bubbles()
        assert(len(g.nodes)==19)
        for node in g.nodes.keys():
            assert len(g.get_neighbors(node)) <= 1
            assert len(g.get_neighbors(-node)) <= 1


if __name__ == '__main__':
    unittest.main()
