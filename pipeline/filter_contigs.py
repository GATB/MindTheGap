#!/usr/bin/env python3

import fileinput
import sys
import os

minLength = int(sys.argv[1])
length=0
for line_in in sys.stdin:
    line = line_in.strip()
    if line.startswith('>'):
        if length > minLength:
            print(seq)
        seq = line+"\n"
        length=0
    else:
        length+=len(line)
        seq+=line
        seq+="\n"
