import sys
import csv
from Bio.Seq import Seq
from Bio import SeqIO
from collections import defaultdict
import getopt
csv.field_size_limit(100000000)


def main():
	print(sys.argv[1:])
	try:
		opts, args = getopt.getopt(sys.argv[1:], "s:g:o:", ["snp=", "genome=", "genome_altered="])
	except getopt.GetoptError:
	# print help information and exit:
	#print ('error')  # will print something like "option -a not recognized"
		sys.exit(2)
	
	genome_parser=""
	vcf_reader=""
	out_m=""
	dic_snp={}
	
	for opt, arg in opts:
		print(opt, arg)
		if opt in ('-s', "--snp"):
			vcf_reader = arg
		#print(i)
		elif opt in ('-g', "--genome"):
			genome_parser = arg
		#print(r)
		elif opt in ('-o', "--genome_altered"):
			out_m = arg
			dic_snp=extract_snp(vcf_reader)	
			alter_genome(genome_parser,dic_snp,out_m)
	
def is_valid(inser) :
	allowed="ATCGatcg"
	if all(c in allowed for c in inser ) :
		return True
	else : 
		return False
	
def insert_snp(genome,chromosome,dic_snp) :
	nb_error=0
	#print (chromosome)
	#print ("len before", len(genome))
	if chromosome in dic_snp :
		for listing in dic_snp[chromosome] :
			if listing[1]==genome[listing[0]-1] :
				genome[listing[0]-1]=listing[2]
			else :
				print( "Error SNP: in genome ",genome[listing[0]-2:listing[0]+2], "in vcf ", listing[1],"at position", listing[0] )
				nb_error+=1
	print (" nb SNP substitution failed",nb_error)
	return genome

def extract_snp(vcf_reader):
	vcf_readers=csv.reader(open(vcf_reader,'r'),delimiter='\t')
	dic_snp={}
	for elits in vcf_readers :
		if '#' not in elits[0] and '@' not in elits[0] :
			if len(elits[3])==1 and len(elits[4])==1  and is_valid(elits[3])==True and is_valid(elits[4])==True:
				dic_snp.setdefault(elits[0],[]).append((int(elits[1]),elits[3],elits[4]))
	return dic_snp
			
def alter_genome(genome_parser,dic_snp,out_m):
	genome_parsers=SeqIO.parse(genome_parser, "fasta")
	out_ms=csv.writer(open(out_m,"w"),delimiter="\n")
	for record in genome_parsers :
			elts=str(record.description)
			head=">" + str(record.description)
			old_sequence=list(str(record.seq).upper())
			old_length=len(old_sequence)
			work_seq=old_sequence
			sequence=""
			unmasked_seq=""
			num_exon=0
			total_length=0
			work_seq=insert_snp(work_seq,elts,dic_snp)
			finals=''.join(work_seq)
			print( 'final', len(finals))
			out_m_lists=[head,finals]
			out_ms.writerow(out_m_lists)

if __name__ == "__main__":
    main()
