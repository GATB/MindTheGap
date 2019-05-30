''' ***********************************************
    
    Pairwise alignment using BioPython
    
    Usage:
    identity = PairAlign(sequence1, sequence2, 10, -5, -5)
    
    *********************************************** '''

from Bio import pairwise2

def nb_match(aln):
    i = 0
    for a in range(0,len(aln[0])):
        #s = aln[:,a]
        if aln[0][a] == aln[1][a]:
            i += 1
    return i

def PairAlign(seq1,seq2,match,mismatch,gap):
    alignments = pairwise2.align.globalms(seq1,seq2,match,mismatch,gap,gap)
    seq_length = max(len(seq1), len(seq2))
    matches = nb_match(alignments[0])
    percent_match = (matches / seq_length)
    return(percent_match)


