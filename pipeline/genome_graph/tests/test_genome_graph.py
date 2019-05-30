import unittest

from pipeline.genome_graph.genome_graph import *
from pipeline.genome_graph.paths import *
from pipeline.genome_graph.utils import *


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


    def test_branching_path(self):
        # Simple graph only contains branching nodes
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        p = LinearPath(g,1)
        for node in g.nodes.keys():
            p1 = LinearPath(g,node)
            p2 = LinearPath(g,-node)
            assert p1.extend_right(g)==False
            assert p2.extend_right(g)==False

    def test_one_linear_path(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        g.pop_all_bubbles()
        p = LinearPath(g,1)
        assert p.extend_right(g) == True

    def test_path_extend(self):
        # Assert that left and right extend have the same output
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        g.pop_all_bubbles()
        p1 = LinearPath(g,1)
        p1.extend_right(g)
        p2 = LinearPath(g,11)
        p2.extend_left(g)
        p1.nodeIds ==  p2.nodeIds == [1,11]

    def test_simplify_graph(self):
        g = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple.gfa")
        g.pop_all_bubbles()
        g.merge_all_linear_paths()
        assert(len(g.nodes)==1)

        g2 = GenomeGraph.read_gfa("pipeline/genome_graph/data/simple4.gfa")
        g2.pop_all_bubbles()
        g2.merge_all_linear_paths()


if __name__ == '__main__':
    unittest.main()
