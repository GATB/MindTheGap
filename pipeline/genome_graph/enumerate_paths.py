#!/usr/bin/env python3

"""
Enumerate all paths going thorugh the longest node of a gfa graph
"""

import argparse

import genome_graph

op = argparse.ArgumentParser()
op.add_argument("infile")
op.add_argument("outfile")
opts = op.parse_args()

print("Loading graph")
g = genome_graph.GenomeGraph.read_gfa(opts.infile) 

maxLen = 0
for n in g.nodes:
    seq = g.nodes[n].nodeSeq
    name = g.nodes[n].nodeName
    if len(seq) > maxLen:
        longestNode = n
        maxLen = len(seq)

paths = g.find_all_paths(n)

ofile = open(opts.outfile, "w")
npath = 1
for p in paths:
    seq = p.getSeq(g)
    ofile.write(">" + "path" + str(npath) + "\n" + seq + "\n")
    npath += 1
ofile.close()

