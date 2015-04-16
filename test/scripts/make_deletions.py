#!/usr/bin/env python
# -*- coding: utf-8 -*-
'''Program to generate random deletions in a given genome. C. Lemaitre'''


# 25/07/2013 modif pour avoir un format bed pour les positions de la délétion
import getopt, sys
import random
        

def usage():
    '''Usage'''
    print "-----------------------------------------------------------------------------"
    print sys.argv[0]," : deletion simulator"
    print "-----------------------------------------------------------------------------"
    print "usage: ",sys.argv[0]," -g genome_file -o output [-n nb_del] [-m min_length] [-M max_length] [-s separator] [-N] [-b]"
    print "  -g: input fasta file containing the genome sequence(s)"
    print "  -o: file preffix for output files (.fasta : the new genome, .del.fasta : the deleted sequences, .del.txt : the positions of deletions)"
    print "  -n: number of deletions to generate, default=1"
    print "  -m: min size of the deletions (in bp), default = 100"
    print "  -M: max size of the deletions (in bp), default = 500"
    print "  -s: min distance between two consecutive deletions (in bp), default = 1, only positive separator are taken into account. For now it can not generate overlapping deletions."
    print "  -N: autorize N inside the deletion (but still not in the borders)"
    print "  -b: bed format for output of the positions of deletions (instead of .del.txt will be .del.bed)"
    print "  -h: help"
    print "-----------------------------------------------------------------------------"
    sys.exit(2)


def main():
    try:
        opts, args = getopt.getopt(sys.argv[1:], "hg:n:o:m:M:s:Nb", ["help", "genome=", "num=", "output=", "min=", "max=", "sep=", "enableN", "bed"])
    except getopt.GetoptError, err:
        # print help information and exit:
        print str(err)  # will print something like "option -a not recognized"
        usage()
        sys.exit(2)

    # Default parameters
    sep=1
    min_length=100
    max_length=500
    genome_file=0
    nb_del=1
    output=0
    enableN=0
    bed_format=0
    for opt, arg in opts:
        if opt in ("-h", "--help"):
            usage()
            sys.exit()
        elif opt in ("-o", "--output"):
            output = arg
        elif opt in ("-g", "--genome"):
            genome_file = arg
        elif opt in ("-n", "--num"):
            nb_del = int(arg)
        elif opt in ("-m", "--min"):
            min_length = int(arg)
        elif opt in ("-M", "--max"):
            max_length = int(arg)
        elif opt in ("-s", "--sep"):
            sep = int(arg)
            if sep<0:
                sep=0
        elif opt in ("-N", "--enableN"):
            enableN=1
        elif opt in ("-b", "--bed"):
            bed_format=1
        else:
            assert False, "unhandled option"

    if genome_file == 0 or output==0:
        print "Missing arguments"
        usage()
        sys.exit(2)
    elif min_length<=0 or max_length<min_length:
        print "Error in parameters : deletion length must respect the following condition : 0 < min_length <= max_length"
        sys.exit(2)
    elif nb_del<=0:
        print "Error in parameters : number of deletions should be greater than 0"
        sys.exit(2)
    else:
        
        # defining special characters :
        char_seq="$"
        char_begin="b"
        char_end="e"
        char_N="N"
        char_set="1234567890"+char_seq+char_begin+char_end;
        if enableN==0:
            char_set=char_set+char_N
                
        
        # Getting the genome sequence
        cumul_seq=""
        total_length=0
        names=[]
        nb_seq=0
        nb_col_fasta=0
        filin=open(genome_file,"r")
        for line in filin:
            if line[0] ==">":
                # header
                current_name=(line.lstrip(">")).rstrip("\n")
                names.append(current_name)
                if nb_seq>0:
                    cumul_seq+="$"
                nb_seq+=1
            else:
                current_seq=(line.rstrip("\n")).upper()
                le=len(current_seq)
                cumul_seq+=current_seq
                total_length+=le
                # we record the fasta format (width) to write output in the same format
                if le > nb_col_fasta:
                    nb_col_fasta=le
        filin.close()
        cumul_seq+="$"
        total_length+=nb_seq #  to account for the $


#        # To test :
#        test_pos=[25,5,26]
#        test_length=[5,12,20]
#        remaining_length=total_length
#        new_seq=cumul_seq
#        for i in range(len(test_pos)):
#            left_pos=test_pos[i]
#            del_length=test_length[i]


        # Placing the deletions
        nb_ok=0  # nb of placed deletions
        nb_boucle=0  # to limit the running time if deletions are too difficult to place
        remaining_length=total_length
        new_seq=cumul_seq
        while nb_ok < nb_del and nb_boucle < (20*nb_del):
            left_pos=random.randint(1,remaining_length)
            del_length=random.randint(min_length,max_length)            
            right_pos=left_pos+del_length
            # checking if deletion is possible
            # extend the region from both sides
            if right_pos+sep<remaining_length and left_pos-sep>=0:
                del_seq_ext=new_seq[(left_pos-sep):(right_pos+sep)]
                borders=new_seq[(left_pos-sep):(left_pos+sep)]+new_seq[(right_pos-sep):(right_pos+sep)]
                if not containsAny(del_seq_ext,char_set) and not char_N in borders:
                    # good deletion : does not contain any end of sequence (chr) or begin or end of previous deletions
                    # and near the borders does not have any N
                    nb_ok+=1
                    # inserting the deletion : replacing the delerted sequence by b125e (for an del of length 125)
                    tmp_seq=new_seq[0:left_pos]+char_begin+str(del_length)+char_end+new_seq[right_pos:]
                    new_seq=tmp_seq
                    remaining_length+=2+len(str(del_length))-del_length# new_seq length
            nb_boucle+=1
                
        if nb_ok < nb_del:
            print "Warning: too difficult to place ",str(nb_del)," deletions, only ",str(nb_ok)," placed"

#        print "placed deletions : "+str(nb_ok)


        # Printing the output
        filout_fasta=open(output+".fasta","w")
        filout_del=open(output+".del.fasta","w")
        if bed_format==1:
            filout_pos=open(output+".del.bed","w")
        else:
            filout_pos=open(output+".del.txt","w")

        # Reading the cumul sequence to get the positions and write the 2 fasta files
        # Note : positions are 0-based and begin and end positions are included in the deletion : [5-10] length = 11
        index_seq=0
        chr_seq=""
        old_pos=0 # position in initial sequence (without deletion)
        new_pos=0 # position in new sequence (with the deletions)
        del_begin=0
        inside_del=False
        compt_del=0
        current_name=names[index_seq]
        header=current_name
        del_length_str=""

        # position file header for txt format
        if bed_format==0:
            filout_pos.write("id\tname\tpos\tlength\tinit.inf\tinit.sup\n")
        for c in new_seq:
            if c==char_seq:
                # writing previous sequence :
                fasta_lines=writeFasta(chr_seq,header,nb_col_fasta)
                filout_fasta.write(fasta_lines)
                index_seq+=1
                if index_seq<nb_seq:  # if not the last $
                    current_name=names[index_seq]
                    header=current_name
                    cumul_seq=cumul_seq[old_pos+1:]
                    old_pos=0
                    new_pos=0
                    chr_seq=""
            elif c==char_begin:
                del_begin=old_pos
                del_length_str=""
            elif c==char_end:
                compt_del+=1
                del_length=int(del_length_str)
                del_end=del_begin+del_length
                # Getting the sequence to delete using the coordinates on the initial sequence
                del_seq=cumul_seq[del_begin:del_end]
                old_pos+=del_length
                # Writing the deleted sequence :
                del_header="deletion_"+str(compt_del)+" : "+current_name+"_"+str(new_pos)
                fasta_lines=writeFasta(del_seq,del_header,nb_col_fasta)
                filout_del.write(fasta_lines)
                # Writing the deletion positions in text files
                if bed_format==1:
                    ## ordre des champs : ch_name, pos_begin, pos_end, id, length, init.begin, init.end
                    filout_pos.write(current_name+"\t"+str(new_pos)+"\t"+str(new_pos+1)+"\t"+str(compt_del)+"\t"+str(del_end-del_begin)+"\t"+str(del_begin)+"\t"+str(del_end)+"\n")
                else:
                    filout_pos.write(str(compt_del)+"\t"+current_name+"\t"+str(new_pos)+"\t"+str(del_end-del_begin)+"\t"+str(del_begin)+"\t"+str(del_end)+"\n")
            elif c in "0123456789":
                del_length_str+=c
            else:
                chr_seq+=c
                old_pos+=1
                new_pos+=1

        filout_fasta.close()
        filout_del.close()
        filout_pos.close()

#        print "recovered and printed deletions : "+str(compt_del)

def containsAny(str, set):
    """Check whether 'str' contains ANY of the chars in 'set'"""
    return 1 in [c in str for c in set]

def writeFasta(sequence,name,ncol):
    fasta_lines=">"+name+"\n"
    if ncol>0:
        le=len(sequence)
        ind=0
        while ind<le:
            fasta_lines+=sequence[ind:(ind+ncol)]+"\n"
            ind+=ncol
    else:
        fasta_lines+=sequence+"\n"
    return fasta_lines
    

if __name__ == "__main__":
    main()

# exemple :
# ./make_deletions.py -g essai.fasta -n 1 -m 10 -M 12 -o output

