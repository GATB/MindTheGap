
# Wrapper for the MindTheGap reference guided assembly pipeline

# 1) Map reads on a remote reference genome using BWA MEM
# 2) Assemble reads using minia
# 3) Use MindTheGap to fill the gaps between contigs
# 4) (not included in the wrapper) Clean the GFA graph

import os, math
from os import listdir
from os.path import isfile, join, splitext
import sys, argparse
import subprocess
import shutil
import logging
import time

#from MtgMin_utils import MtgParser, ArgumentFormatterMtg, read_mapping_header, ProgressBar
from pipeline_utils import MtgParser, ArgumentFormatterMtg, contig_stats


#os.chdir(os.path.split(os.path.realpath(__file__))[0])

#-------------------------------------------------------------------------------------------------------------
# Arg parser
#-------------------------------------------------------------------------------------------------------------
parser = MtgParser(formatter_class=ArgumentFormatterMtg)

parserMain = parser.add_argument_group("[main options]")
parserMapping = parser.add_argument_group("[mapping options]")
parserAssembly = parser.add_argument_group("[assembly options]")
parserGapfilling = parser.add_argument_group("[gapfilling options]")
parserContinue = parser.add_argument_group("[continue options]")
parserCore = parser.add_argument_group("[core options]")

parserMain.add_argument('-in', action="store", dest="input_file", help="input reads file", required=False)
parserMain.add_argument('-1', action="store", dest="input_file1", help="input reads first file", required=False)
parserMain.add_argument('-2', action="store", dest="input_file2", help="input reads second file", required=False)
parserMain.add_argument('-fof', action="store", dest="input_fof", help="input file of read files (if paired files, 2 columns tab-separated)", required=False)
parserMain.add_argument('-out', action="store", dest="out", default="./mtg_results", help="output directory for result files")

parserMapping.add_argument('-ref', action="store", dest="ref_genome", help="bwa index", required=False)

parserAssembly.add_argument('-minia-bin', action="store", dest="minia_bin", help="path to Minia binary")
parserAssembly.add_argument('-assembly-kmer-size', action="store", dest="minia_kmer_size", help="kmer size used for Minia assembly (should be given even if bypassing minia assembly step, usefull knowledge for gap-filling)", default="31")
parserAssembly.add_argument('-assembly-abundance-min', action="store", dest="minia_abundance", help="Minimal abundance of kmers used for assembly", default="auto")
parserAssembly.add_argument('-min-contig-size', action="store", dest="min_contig_size", default="0", help="minimal size for a contig to be used in gapfilling")

parserGapfilling.add_argument('-mtg-dir', action="store", dest="mtg_dir", help="path to MindTheGap build directory", required=False)
parserGapfilling.add_argument('-gapfilling-kmer-size', action="store", dest="mtg_kmer_size", help="kmer size used for gapfilling", default="31")
parserGapfilling.add_argument('-gapfilling-abundance-min', action="store", dest="mtg_abundance", help="Minimal abundance of kmers used for gapfilling", default="auto")
parserGapfilling.add_argument('-max-nodes', action="store", dest="max_nodes", help="Maximum number of nodes in contig graph", default="100")
parserGapfilling.add_argument('-max-length', action="store", dest="max_length", help="Maximum length of gapfilling (nt)", default="10000")

parserContinue.add_argument('-contigs',action="store",dest="continue_contigs",help="Contigs in fasta format - override mapping and assembly")
parserContinue.add_argument('-graph',action="store",dest="continue_h5",help="Graph in h5 format - override graph creation")

parserCore.add_argument('-nb-cores', action="store", dest="nb_cores", help="number of cores", default="0")


args =  parser.parse_args()

if args.input_file is None and (args.input_file1 is None or args.input_file2 is None) and (args.input_fof is None) and (args.continue_h5 is None or args.continue_contigs is None):
    parser.error("Please supply reads as -in or -1/-2 or -fof or -graph")

if args.continue_contigs is None and args.ref_genome is None:
    parser.error("Either -ref or -contigs is required")

if args.continue_contigs is None and args.minia_bin is None and shutil.which("minia") is None :
    parser.error("minia is not in PATH, please supply -minia-bin argument")

if (shutil.which("MindTheGap") is None or shutil.which("dbgh5") is None) and args.mtg_dir is None:
    parser.error("MindTheGap and dbgh5 are not in PATH, please supply -mtg-dir argument")

#if (args.continue_h5 is None or args.continue_contigs is None) and args.input_file is None:
    #parser.error("-in is required if neither -graph and -contigs are supplied")


#-------------------------------------------------------------------------------------------------------------
# MindTheGap pipeline
#-------------------------------------------------------------------------------------------------------------

#Create some dirs and filenamesF
if not os.path.exists(args.out): os.makedirs(args.out)
outDir = args.out

logsDir = os.path.join(outDir, "logs")
if not os.path.exists(logsDir): os.makedirs(logsDir)
logging.basicConfig(
    format="%(asctime)s [%(threadName)-12.12s] [%(levelname)-5.5s]  %(message)s",
    handlers=[logging.FileHandler(os.path.join(logsDir,"pipeline.log")),logging.StreamHandler()],
    level=logging.INFO)
logger = logging.getLogger(__name__)
#-------------------------------------------------------------------------------------------------------------
# Mapping
#-------------------------------------------------------------------------------------------------------------

startTime = time.time()
if args.continue_contigs is None:

    mappingDir = os.path.join(outDir, "mapping")
    if not os.path.exists(mappingDir): os.makedirs(mappingDir)

    input1 = []
    input2 = []
    #Getting input files
    if args.input_file is not None:
        input1.append(args.input_file)
        input2.append("")
    if args.input_file1 is not None:
        input1.append(args.input_file1)
        input2.append(args.input_file2)
    if args.input_fof is not None:
        with open(args.input_fof,"r") as fofFile:
            for line in fofFile:
                line = line.rstrip()
                sp = line.split("\t")
                input1.append(sp[0])
                if len(sp)>1:
                    input2.append(sp[1])
                else:
                    input2.append("")

    mappingLog = os.path.join(logsDir,"mapping.log")


    #Check if ref is indexed
    if not os.path.exists(args.ref_genome+".bwt"):
        logger.info("Indexing reference genome")
        indexCommand = "bwa index " + args.ref_genome
        with open(mappingLog,"wb") as log:
            p = subprocess.Popen(indexCommand.split(),stdout=log,stderr=log)
        p.wait()
        if p.returncode != 0: logger.error("Indexing ref genome failed"); exit(1)

    logger.info("Starting mapping")


    fqFile = os.path.join(mappingDir,"mapped_reads.fastq")
    catCommand = "cat "
    nbTotReads = 0

    for i in range(len(input1)):
        name="file" + str(i)
        #Create commands
        mappingCommand = "bwa mem"
        mappingCommand += " -t " + args.nb_cores
        mappingCommand += " " + args.ref_genome
        mappingCommand += " " + input1[i] + " " + input2[i]

        sam2bamCommand = "samtools view -b -F 4 -"
        bamFile = os.path.join(mappingDir, name + ".bam")

        tmpFqFile = os.path.join(mappingDir, name + "_mapped_reads.fastq")
        bam2fqCommand = "samtools bam2fq " + bamFile

        # Executing
        logger.info("\tCall : "+ mappingCommand)

        with open(mappingLog,"wb") as out:
            p1 = subprocess.Popen(mappingCommand.split(),stdout=subprocess.PIPE,stderr=out)
            p2 = subprocess.Popen(sam2bamCommand.split(),stdin=p1.stdout,stdout=open(bamFile,"w"),stderr=out)
        p2.wait()
        if p2.returncode != 0: logger.error("Mapping failed"); exit(1)

        with open(tmpFqFile,"wb") as out,open(mappingLog,"wb") as log:
            p = subprocess.Popen(bam2fqCommand.split(),stdout=out,stderr=subprocess.PIPE)
            out.close()
        p.wait()
        if p.returncode != 0: logger.error("Conversion to fastq failed"); exit(1)

        logger.info(name + " : Mapping done")
        logLine = p.stderr.readlines()[-1]
        nbReads = [int(s) for s in logLine.split() if s.isdigit()]
        logger.info(name + " : " + str(nbReads)+" reads mapped")

        catCommand += " " + tmpFqFile

        nbTotReads += sum(nbReads)

    #Concatenation of fastq files
    with open(fqFile,"wb") as out, open(mappingLog,"wb") as log:
        #print(catCommand)
        p = subprocess.Popen(catCommand.split(),stdout=out,stderr=log)
    p.wait()
    if p.returncode != 0: logger.error("Cat fastq files failed"); exit(1)

    logger.info("All mapping done")
    logger.info(str(nbTotReads)+" reads mapped")


mappingTime = time.time()
mappingDuration = round(mappingTime - startTime,1)


#-------------------------------------------------------------------------------------------------------------
# Assembly
#-------------------------------------------------------------------------------------------------------------
if args.continue_contigs is None:

    assemblyDir = os.path.join(outDir, "assembly")
    if not os.path.exists(assemblyDir): os.makedirs(assemblyDir)
    assemblyName = "minia_k"+args.minia_kmer_size+"_abundancemin_"+args.minia_abundance
    assemblyPrefix = os.path.join(assemblyDir,assemblyName)

    if args.minia_bin is None:
        assemblyCommand = "minia"
    else:
        assemblyCommand = args.minia_bin
    assemblyCommand += " -nb-cores " + args.nb_cores
    assemblyCommand += " -in " + fqFile
    assemblyCommand += " -kmer-size " + args.minia_kmer_size
    assemblyCommand += " -abundance-min " + args.minia_abundance
    assemblyCommand += " -out " + assemblyPrefix
    assemblyCommand += " -out-tmp " + assemblyDir
    assemblyLog = os.path.join(logsDir,"assembly.log")

    logger.info("Starting assembly")
    logger.info("\tCall : "+ assemblyCommand)
    logger.info("\tLog file : "+assemblyLog)

    with open(assemblyLog,"wb") as out:
        p = subprocess.Popen(assemblyCommand.split(),stdout=out,stderr=out)
    p.wait()
    if p.returncode != 0: logger.error("Assembly failed"); exit(1)

    contigFile = assemblyPrefix + ".contigs.fa"

    #Check if assembly succeeded :
    if not os.path.exists(contigFile):
        logger.error("Assembly failed: no contig output (try with other parameters)"); exit(1)

    logger.info("Assembly done")


else:
    contigFile = args.continue_contigs
    assemblyPrefix =  args.continue_contigs
    assemblyName = os.path.basename(contigFile)

contig_stats(contigFile)

logger.info("Filtering contigs")
scriptPath = sys.path[0]
filteringCommand = os.path.join(scriptPath,"filter_contigs.py")
filteringCommand += " " + args.min_contig_size

filteredFile = assemblyPrefix + "_filtered_"+args.min_contig_size+".fa"
with open(contigFile,"r") as input, open(filteredFile,"w") as output:
    p = subprocess.Popen(filteringCommand.split(),stdin=input,stdout=output)
p.wait()
logger.info("Contigs filtered")
contig_stats(filteredFile)

assemblyTime = time.time()
assemblyDuration = round(assemblyTime - mappingTime,1)


#-------------------------------------------------------------------------------------------------------------
# Graph construction
#-------------------------------------------------------------------------------------------------------------

gapfillingDir = os.path.join(outDir, "gapfilling")
if not os.path.exists(gapfillingDir): os.makedirs(gapfillingDir)
h5Name = "graph_k"+args.mtg_kmer_size+"_abundancemin_"+args.mtg_abundance+".h5"
h5File = os.path.join(gapfillingDir,h5Name)
gapfillingName = assemblyName+"_filtered_"+args.min_contig_size+"_gapfilling_k"+args.mtg_kmer_size+"_abundancemin_"+args.mtg_abundance
gapfillingPrefix = os.path.join(gapfillingDir,gapfillingName)


if args.continue_h5 is None:
    if args.mtg_dir is None:
        h5Command = ["dbgh5"]
    else:
        h5Command = [os.path.join(args.mtg_dir,"ext/gatb-core/bin/dbgh5")]
    if args.input_file1 is not None:
        h5Command.extend(["-in",args.input_file1+","+args.input_file2])
    if args.input_file is not None:
        h5Command.extend(["-in",args.input_file])
    if args.input_fof is not None:
        input_file_list=[]
        with open(args.input_fof,"r") as fofFile:
            for line in fofFile:
                line = line.rstrip()
                sp = line.split("\t")
                input_file_list.append(sp[0])
                if len(sp)>1:
                    input_file_list.append(sp[1])
        h5Command.extend(["-in",",".join(input_file_list)])

    h5Command.extend(["-kmer-size",args.mtg_kmer_size])
    h5Command.extend(["-abundance-min",args.mtg_abundance])
    h5Command.extend(["-out",h5File])
    h5Command.extend(["-nb-cores",args.nb_cores])
    h5Log = os.path.join(logsDir,"dbgh5.log")


    logger.info("Building graph")
    logger.info("\tCall : "+ " ".join(h5Command))
    logger.info("\tLog file : "+ h5Log)

    with open(h5Log,"wb") as out:
        p = subprocess.Popen(h5Command,stdout=out,stderr=out)
    p.wait()

graphTime = time.time()
graphDuration = round(graphTime - assemblyTime,1)

#-------------------------------------------------------------------------------------------------------------
# Gapfilling
#-------------------------------------------------------------------------------------------------------------

if args.mtg_dir is None:
    mtgCommand = ["MindTheGap"]
else:
    mtgCommand = [os.path.join(args.mtg_dir,"bin/MindTheGap")]
mtgCommand.append("fill" )
mtgCommand.extend(["-contig",filteredFile])
if args.continue_h5 is None:
    mtgCommand.extend(["-graph",h5File])
else:
    mtgCommand.extend(["-graph",args.continue_h5])
#mtgCommand.extend(["-abundance-min",args.mtg_abundance])
mtgCommand.extend(["-overlap",args.minia_kmer_size])
mtgCommand.extend(["-out",gapfillingPrefix])
mtgCommand.extend(["-nb-cores",args.nb_cores])
mtgCommand.extend(["-max-length",args.max_length])
mtgCommand.extend(["-max-nodes",args.max_nodes])
mtgLog = os.path.join(logsDir,"gapfilling.log")

logger.info("Gapfilling")
logger.info("\tCall : "+ " ".join(mtgCommand))

with open(mtgLog,"wb") as out:
    p = subprocess.Popen(mtgCommand,stdout=out,stderr=out)
p.wait()
gapfillingTime = time.time()
gapfillingDuration = round(gapfillingTime - graphTime,1)



# Output Mtg results
logFile=open(mtgLog,'r')
switch = False
for line in logFile:
    if line.startswith("Results"):
        switch = True
    if switch:
        logger.info(line.rstrip())

logger.info("Runtime :")
logger.info("\tMapping : " + str(mappingDuration))
logger.info("\tAssembly : " + str(assemblyDuration))
logger.info("\tGraph creation : " + str(graphDuration))
logger.info("\tGapfilling : " + str(gapfillingDuration))
