#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''
Program to generate reads from a fasta file (for tests). P. Marijon based one
get_subsequence.py write by C. Lemaitre.
'''

import getopt, sys
import random
        

def usage():
    '''Usage'''
    print "-----------------------------------------------------------------------------"
    print sys.argv[0]," : sub-sequence"
    print "-----------------------------------------------------------------------------"
    print "usage: ",sys.argv[0]," -f fasta_file -n numbre -l length"
    print "  -f: input fasta file"
    print "  -n: numbre of read you want(def : 1)"
    print "  -l: length of sub-sequence (def : 1)"
    print "  -h: help"
    print "-----------------------------------------------------------------------------"
    sys.exit(2)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hf:n:l:", ["help", "fasta=", "num=", "len="])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    # Default parameters
    read_len=1
    fasta_file=0
    num_loop=1
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-f", "--fasta"):
            fasta_file = arg
        elif opt in ("-l", "--len"):
            read_len = int(arg)
        elif opt in ("-n", "--num"):
            num_loop = int(arg)
        else:
            assert False, "unhandled option"

    if fasta_file == 0 :
        print "Missing arguments"
        usage()
        return 2
    
    else:
        header = "base_header"
        sequence = ""
        filin = open(fasta_file,"r")

        for line in filin:
            if line[0] ==">":
                # header
                header = (line.lstrip(">")).rstrip("\n").rstrip(" ")
            else:
                sequence+=(line.rstrip("\n")).upper()

        filin.close()

        sequence_len = len(sequence)
        
        if sequence_len == 0 :
            print "warning we didn't find fasta sequence in file."
            return 1

        if sequence_len < read_len :
            print "warning read length is upper than sequence length we can't generate read."
            return 1
        
        for i in range(num_loop) :
            pos = sequence_len
            while pos + read_len > sequence_len :
                pos = random.randint(0, sequence_len)

            print ">"+header+"_read"+str(i)+"_pos_"+str(pos)+":"+str(pos+read_len)
            print sequence[pos:pos+read_len]

    return 0
            
if __name__ == "__main__":
    exit(main())
    
# exemple :
# ./get_subsequence -f random_1Mb.fa -b 44444

