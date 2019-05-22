import subprocess
import os
import sys

alt_map = {'ins':'0'}
complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'} 

def reverse_complement(dna):
    complement = {'A': 'T', 'C': 'G', 'G': 'C', 'T': 'A'}
    return ''.join([complement[base] for base in dna[::-1]])

def compare_strings(str1,str2):
    # returns the position of last common letter between two strings
    i = 0
    cont = True
    imax = min(len(str1),len(str2))

    while cont:
        eq = str1[i]==str2[i]
        if eq:
            i += 1
            if i == imax:
                return(i)
        else:
            return(i)

def locate_nw_binary():
    scriptDir = os.path.abspath(os.path.dirname(sys.argv[0]))
    nwCommand = os.path.join(scriptDir,"../../build/bin/nwalign")
    if os.path.isfile(nwCommand)==False:
        print("No nwAlign binary found in " + nwCommand)
    return(nwCommand)

def nw_align(seq1,seq2,nwCommand):
    seq = seq1+"\n"+seq2
    p = subprocess.Popen(nwCommand,stdin=subprocess.PIPE,stdout=subprocess.PIPE)# ,stdout=out,stderr=out)
    
    try:
        out = p.communicate(input=str.encode(seq))[0].decode()
        id = float(out.strip())
        return(id)
    except:
        print("Error during nwAlign")

            