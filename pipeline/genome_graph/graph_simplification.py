#!/usr/bin/env python3

"""
Remove redundant gapfillings in a MindTheGap assembly graph
"""

import argparse

import genome_graph

op = argparse.ArgumentParser()
op.add_argument("infile")
op.add_argument("outfile")
opts = op.parse_args()

g = genome_graph.GenomeGraph.read_gfa(opts.infile)

g.pop_all_bubbles()

g.write_gfa(opts.outfile)