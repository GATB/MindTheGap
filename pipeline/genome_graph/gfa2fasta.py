#!/usr/bin/env python3

"""
Output all sequences of a GFA graph longuer than minlength to a fasta file
"""

import argparse

import genome_graph

op = argparse.ArgumentParser()
op.add_argument("infile")
op.add_argument("outfile")
op.add_argument("minlength")
opts = op.parse_args()

print("Loading graph")
g = genome_graph.GenomeGraph.read_gfa(opts.infile) 

ofile = open(opts.outfile, "w")
for n in g.nodes.values():
    seq = n.nodeSeq
    name = n.nodeName
    if len(seq) > int(opts.minlength):
        ofile.write(">" + name + "\n" + seq + "\n")
ofile.close()

