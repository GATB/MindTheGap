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

            