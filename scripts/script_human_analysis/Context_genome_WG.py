from gatb import Graph
import csv
import sys
import numpy as np
import pandas as pd
from Bio.Seq import Seq
from Bio import SeqIO
import statistics
import pandas as pd
import getopt


def main():
    print(sys.argv[1:])
    try:
        opts, args = getopt.getopt(sys.argv[1:], "g:p:c:b:s:t:v:o:m:", ["graph=", "genome_parser=", "branching_outp=", "bkpt_file=", "truth_vcf=", "outp_context=","vcf_fill=","bkpt_outp=","threshold="])
    except getopt.GetoptError:
        # print help information and exit:
		#print ('error')  # will print something like "option -a not recognized"
        sys.exit(2)

    # Default parameters
	#print(opts)
    graph = ""
    genome_parser = ""
    branching_outp = ""
    bkpt_file = ""
    truth_vcf = ""
    outp_context = ""
    vcf_fill=""
    dic_parse = {}
    bkpt_outp=""
    threshold=0.80
    for opt, arg in opts:
        print(opt, arg)
        if opt in ('-g', "--graph"):
            graph = arg
	#print(i)
        elif opt in ('-p', "--genome_parser"):
            genome_parser = arg
	#print(r)
        elif opt in ('-c', "--branching_outp"):
            branching_outp = arg
	#print(i)
        elif opt in ('-b', "--bkpt_file"):
            bkpt_file = arg
        elif opt in ('-s', "--truth_vcf"):
            truth_vcf = arg
        elif opt in ('-t', "--outp_context"):
            outp_context = arg
        elif opt in ('-v', "--vcf_fill"):
            vcf_fill = arg
        elif opt in ('-o', "--bkpt_outp"):
            bkpt_outp = arg
        elif opt in ('-m', "--threshold"):
            threshold =int(arg)
        else:
	        assert False, "unhandled option"
    #dic_parse=parsing_genome_branching(graph, genome_parser)
    #parsing_genome_branching2(graph,genome_parser,branching_outp)
    #analyze_genomic_context(bkpt_file, branching_outp,bkpt_outp)
    #analyze_genomic_context2(bkpt_file,branching_outp,truth_vcf,vcf_fill)
    #write_context_genomic(dico_TP, dico_FP,outp_context)

    analyze_genomic_context_direct(bkpt_file, graph,genome_parser,bkpt_outp,threshold)
    
    
def analyze_genomic_context_direct(bkpt, graph_h5,genome,outp_bkpt,threshold):
    forma = "-in "+graph_h5
    graph = Graph(forma)
    graph
    genomes_parser = SeqIO.parse(open(genome), "fasta")
    genome_parser = SeqIO.parse(open(bkpt), "fasta")
    outp_find = csv.writer(open(outp_bkpt, 'w'), delimiter='\n')
    liste_good=[]
    liste_chrom_seen=[]
    liste_position = []
    dico_first={}
    dico_second={}
    total_bkpt=0
    count=0
    for element in genome_parser :
        if count%2==0 :
            dico_first.setdefault(element.description.split('_')[1], []).append(int(element.description.split('_')[3]))
            total_bkpt+=1
        count+=1
	
    for elt in genomes_parser :
        str_chromosome = str(elt.seq)
        id_chrom=str(elt.description)
        if id_chrom in dico_first :
            for value in dico_first[id_chrom] :
                 #print(value)
                 sum_degree=[]
                 for i in range (50) :
                    kmer = str_chromosome[value-i-31:value-i]
                    node = graph[kmer]
                    bytes(node)
                    assert node.reversed == node
                    sum_degree.append(node.out_degree)
                    sum_degree.append(node.in_degree)
                 percentage_concatenate = (sum_degree.count(1)+sum_degree.count(2))/(len(sum_degree))
                    #print(percentage_concatenate)
                 if percentage_concatenate > threshold :  
                    dico_second.setdefault(id_chrom, []).append(int(value))
    a=0
    for i in dico_second:
        a+=len(dico_second[i])
    print ("total breakpoints kept : ", a, " on ", total_bkpt)
    genome_parser = SeqIO.parse(open(bkpt), "fasta")
    for element in genome_parser :
        #print(element)
        if int(element.description.split('_')[3]) in dico_second[element.description.split('_')[1]] :
            outp_find.writerow([">"+element.description,element.seq])


def parsing_genome_branching(graph_h5, genome):
    forma = "-in "+graph_h5
    graph = Graph(forma)
    graph
    genome_parser = SeqIO.parse(open(genome), "fasta")
    dico_parse={}
    for chromosome in genome_parser:
        str_chromosome = str(chromosome.seq)
        #print(str_chromosome[0:31])
        for i in range(len(str_chromosome)-31):
            if "N" not in str_chromosome[i:i+31]:
                #print (str_chromosome[i:i+31])
                kmer = str_chromosome[i:i+31]
                node = graph[kmer]
                bytes(node)
                assert node.reversed == node
                position = i
                out_deg = node.out_degree
                in_deg = node.in_degree
                #if out_deg>1 or in_deg>1
                dico_parse.setdefault(chromosome.description, []).append((position, in_deg, out_deg))
            else:
                dico_parse.setdefault(chromosome.description, []).append((i, 0, 0))
    for a in dico_parse :
        print(a)
    return (dico_parse)


def parsing_genome_branching2(graph_h5, genome, outpt):
    forma = "-in "+graph_h5
    graph = Graph(forma)
    graph
    genome_parser = SeqIO.parse(open(genome), "fasta")
    output_bed = csv.writer(open(outpt, "w"), delimiter="\t")
    output_bed.writerow(["chr","position", "in_degree", "out_degree"])

    

    for chromosome in genome_parser:
        str_chromosome = str(chromosome.seq)
        #print(str_chromosome[0:31])
        for i in range(len(str_chromosome)-31):
            if "N" not in str_chromosome[i:i+31]:
                #print (str_chromosome[i:i+31])
                kmer = str_chromosome[i:i+31]
                node = graph[kmer]
                bytes(node)
                assert node.reversed == node
                position = i
                out_deg = node.out_degree
                in_deg = node.in_degree
                #if out_deg>1 or in_deg>1:
                output_bed.writerow([chromosome.description,position, in_deg, out_deg])
            else:
                output_bed.writerow([chromosome.description,i, 0, 0])


def analyze_genomic_context2(bkpt, branching_bed,truth_file,vcf_file):
    genome_parser = SeqIO.parse(open(bkpt), "fasta")
    input_bed = pd.read_csv(branching_bed, sep='\t')
    sum_FP = 0
    sum_TP = 0
    good_tp=0
    bad_tp=0
    good_fp=0
    bad_fp=0
    fp_remove=0
    fp_kept=0
    tp_kept=0
    outp_fill = csv.reader(open(vcf_file, 'r'), delimiter='\t')
    liste_fill=[]
    for elt in outp_fill :
        if "#" not in elt[0] and "@" not in elt[0] :
            liste_fill.append(int(elt[1]))
    #print(liste_fill)



    truth_parser = csv.reader(open(truth_file, 'r'), delimiter='\t')
    next(truth_parser, None)  # skip header
    liste_position = []
    liste_truth = []
    i=0
    for elements in truth_parser:
        if i % 2 == 0:
            liste_position.append(int(elements[2]))
            liste_truth.append(int(elements[2]))
        i += 1
    for element in genome_parser:
        #Avoid repetition from couple breakpoint
        if 'left' in str(element.description):
            #positon changes if it is a back up insertion
            if element.description.split('_')[2] == 'backup':
                pos = int(element.description.split('_')[4])
            else:
                pos = int(element.description.split('_')[3])
            test_in = input_bed[(input_bed['chr'] == element.description.split('_')[1]) & (input_bed['position'] > pos-100) & (input_bed['position'] <= pos-31)]
            liste_in=test_in.in_degree.tolist()
            liste_out=test_in.out_degree.tolist()
            percentage_concatenate = (liste_in.count(1)+liste_in.count(2)+liste_out.count(1)+liste_out.count(2))/(len(liste_in)+len(liste_out))
            if percentage_concatenate > 0.80:
                    sum_TP += 1
                   
                    if pos in liste_truth or pos-1 in liste_truth or pos+1 in liste_truth:
                        good_tp+=1
                        if pos in liste_fill or pos-1 in liste_fill or pos+1 in liste_fill:
                            #print('FP failed filled', element.description)
                            tp_kept += 1
                    else : 
                        bad_tp+=1
                        #print('FP failed', element.description,percentage_concatenate)
                        if pos in liste_fill or pos-1 in liste_fill or pos+1 in liste_fill:
                            #print('FP failed filled', element.description)
                            fp_kept += 1

                        

            else : 
                sum_FP+=1
                #print('FP',element.description, percentage_concatenate)
                if pos in liste_truth or pos-1 in liste_truth or pos+1 in liste_truth:
                    bad_fp += 1
                    #print('TP failed', element.description, percentage_concatenate)
                else : 
                    good_fp+=1
                    if pos in liste_fill or pos-1 in liste_fill or pos+1 in liste_fill:
                            #print('FP failed filled', element.description)
                            fp_remove+=1
    print("sum sup 0.5 ", sum_TP,"sum below 0.5 ", sum_FP)
    print('TP predicted ', good_tp,'TP_failed ',bad_tp,'FP predicted ', good_fp,'FP failed ',bad_fp, " fp removed ", fp_remove, "fp_kept",fp_kept,"tp_kept",tp_kept)

def analyze_genomic_context(bkpt, branching_bed,outp_bkpt):
    genome_parser = SeqIO.parse(open(bkpt), "fasta")
    input_bed = pd.read_csv(branching_bed, sep='\t')
    outp_find = csv.writer(open(outp_bkpt, 'w'), delimiter='\n')
    liste_good=[]

    liste_position = []

    for element in genome_parser:
        #Avoid repetition from couple breakpoint
        if 'left' in str(element.description):
            #positon changes if it is a back up insertion
            if element.description.split('_')[2] == 'backup':
                pos = int(element.description.split('_')[4])
            else:
                pos = int(element.description.split('_')[3])
            test_in = input_bed[(input_bed['chr'] == element.description.split('_')[1]) & (input_bed['position'] > pos-100) & (input_bed['position'] <= pos-31)]
            liste_in=test_in.in_degree.tolist()
            liste_out=test_in.out_degree.tolist()
            percentage_concatenate = (liste_in.count(1)+liste_in.count(2)+liste_out.count(1)+liste_out.count(2))/(len(liste_in)+len(liste_out))
            if percentage_concatenate > 0.50:
                liste_good.append(element.description.split('_')[0])

    print(liste_good)
    genome_parser = SeqIO.parse(open(bkpt), "fasta")
    for element in genome_parser :
        print(element)
        if element.description.split('_')[0] in liste_good :
            outp_find.writerow([">"+element.description,element.seq])

    #break
#print ("HEY",dico_FP)


if __name__ == "__main__":
    main()
