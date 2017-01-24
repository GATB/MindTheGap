/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G.Rizk
 *
 *  This program is free software: you can redistribute it and/or modify
 *  it under the terms of the GNU Affero General Public License as
 *  published by the Free Software Foundation, either version 3 of the
 *  License, or (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU Affero General Public License for more details.
 *
 *  You should have received a copy of the GNU Affero General Public License
 *  along with this program.  If not, see <http://www.gnu.org/licenses/>.
 *****************************************************************************/

#include "Finder.hpp"
#include <FindBreakpoints.hpp>
#include <IFindObserver.hpp>
#include <FindBackup.hpp>
#include <FindDeletion.hpp>
#include <FindHeteroInsertion.hpp>
#include <FindInsertion.hpp>
#include <FindSNP.hpp>
#include <limits> //for std::numeric_limits

//#define PRINT_DEBUG
/********************************************************************************/


Finder::~Finder()
{
    if(_refBank != 0) {_refBank->forget();}
    // delete _graph ?
}



void HelpFinder(void* target)
{
	if(target!=NULL)
	{
		Finder * obj = (Finder *) target;
		obj->FinderHelp();
	}
}
/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
Finder::Finder ()  : Tool ("MindTheGap find")
{
    _refBank = 0;
    _kmerSize = 31;
    _max_repeat = 0;
    _het_max_occ = 1;
    _snp_min_val = 5;
    _nbCores = 0;
    _breakpoint_file_name = "";
    _vcf_file_name = "";
    _nb_homo_clean = 0;
    _nb_homo_fuzzy = 0;
    _nb_hetero_clean = 0;
    _nb_hetero_fuzzy = 0;
    _nb_fuzzy_deletion = 0;
    _nb_clean_deletion = 0;
    _nb_solo_snp = 0;
    _nb_multi_snp = 0;
    _nb_backup = 0;
    
    _homo_only = false;
    _homo_insert = true;
    _hete_insert = true;
    _snp = true;
    _backup = false;
    _deletion = true;
	
	setHelp(&HelpFinder);
	setHelpTarget(this);
	
    // Option parser, with several sub-parsers
    setParser (new OptionsParser ("MindTheGap find"));

    IOptionsParser* generalParser = new OptionsParser("General");
	
    // generalParser->push_back (new OptionNoParam (STR_VERSION, "version", false)); // move this option in the main.cpp
    generalParser->push_front (new OptionOneParam (STR_VERBOSE,     "verbosity level",      false, "1"  ));
    generalParser->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory for graph building (in MBytes)", false, "2000"));
    generalParser->push_front (new OptionOneParam (STR_MAX_DISK, "max disk for graph building (in MBytes)", false, "0"));
    generalParser->push_front (new OptionOneParam (STR_NB_CORES,    "number of cores",      false, "0"  ));

    IOptionsParser* inputParser = new OptionsParser("Input / output");
    inputParser->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
	inputParser->push_front (new OptionOneParam (STR_URI_OUTPUT_TMP, "prefix for output temporary files", false, "."));
	
	
    inputParser->push_front (new OptionOneParam (STR_URI_REF, "reference genome file", true,""));
    inputParser->push_front (new OptionOneParam (STR_URI_GRAPH, "input graph file (likely a hdf5 file)",  false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_INPUT, "input read file(s)",  false, ""));

    IOptionsParser* finderParser = new OptionsParser("Detection");
    finderParser->push_front (new OptionNoParam (STR_INSERT_ONLY, "search only insertion breakpoints (do not report other variants)", false));
    //finderParser->getParser(STR_INSERT_ONLY)->setVisible(false);
    finderParser->push_front (new OptionOneParam (STR_HET_MAX_OCC, "maximal number of occurrences of a kmer in the reference genome allowed for heterozyguous breakpoints", false,"1"));
    //allow to find heterozyguous breakpoints in n-repeated regions of the reference genome
    finderParser->push_front (new OptionOneParam (STR_MAX_REPEAT, "maximal repeat size detected for fuzzy sites", false, "5"));
    finderParser->push_front (new OptionNoParam (STR_HOMO_ONLY, "search only homozygous breakpoints", false));

    //Options not for the common user
    finderParser->push_front (new OptionOneParam (STR_SNP_MIN_VAL, "minimal number of kmers to validate a SNP", false, "5"));
    finderParser->getParser(STR_SNP_MIN_VAL)->setVisible(false);


    //Options usefull only for debugging
    finderParser->push_front (new OptionNoParam (STR_SNP_ONLY, "search only SNPs", false));
    finderParser->getParser(STR_SNP_ONLY)->setVisible(false);
    finderParser->push_front (new OptionNoParam (STR_DELETION_ONLY, "search only deletion variants", false));
    finderParser->getParser(STR_DELETION_ONLY)->setVisible(false);
    finderParser->push_front (new OptionNoParam (STR_HETERO_ONLY, "search only heterozygous insertion breakpoints", false));
    finderParser->getParser(STR_HETERO_ONLY)->setVisible(false);
    finderParser->push_front (new OptionNoParam (STR_NO_SNP, "do not search SNPs", false));
    finderParser->getParser(STR_NO_SNP)->setVisible(false);
    finderParser->push_front (new OptionNoParam (STR_NO_INSERT, "do not search insertion breakpoints", false));
    finderParser->getParser(STR_NO_INSERT)->setVisible(false);
    finderParser->push_front (new OptionNoParam (STR_NO_DELETION, "do not search deletions", false));
    finderParser->getParser(STR_NO_DELETION)->setVisible(false);
    finderParser->push_front (new OptionNoParam (STR_NO_HETERO, "do not search heterozygous insertion breakpoints", false));
    finderParser->getParser(STR_NO_HETERO)->setVisible(false);
    finderParser->push_front (new OptionNoParam (STR_WITH_BACKUP, "report also unusual breakpoints (gap size is larger than kmer-size/2 and does not validate a common variant)", false));
    finderParser->getParser(STR_WITH_BACKUP)->setVisible(false);


    IOptionsParser* graphParser = new OptionsParser("Graph building");
    string abundanceMax = Stringify::format("%ld", std::numeric_limits<CountNumber>::max()); //to be sure in case CountNumber definition changes
    graphParser->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "maximal abundance threshold for solid kmers", false, abundanceMax));
    //TODO for release : change 3 -> auto
    graphParser->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "minimal abundance threshold for solid kmers", false, "auto"));
    graphParser->push_front (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));
	//IOptionsParser* graphParser = SortingCountAlgorithm<>::getOptionsParser(false);
    //graphParser->setName ("Graph building");

    /** We hide some options. */
//    const char* optionNames[] = {
//        STR_URI_INPUT, STR_KMER_ABUNDANCE_MIN_THRESHOLD, STR_HISTOGRAM_MAX, STR_SOLIDITY_KIND,
//        STR_URI_SOLID_KMERS, STR_URI_OUTPUT, STR_URI_OUTPUT_DIR, STR_MINIMIZER_TYPE, STR_MINIMIZER_SIZE, STR_REPARTITION_TYPE
//    };
//    for (size_t i=0; i<sizeof(optionNames)/sizeof(optionNames[0]); i++)
//    {
//        if (IOptionsParser* p = graphParser->getParser(optionNames[i]))  { p->setVisible(false);  }
//    }

    getParser()->push_front(generalParser);
    getParser()->push_front(finderParser);
    getParser()->push_front(graphParser);
    getParser()->push_front(inputParser);
    
}


void Finder::FinderHelp()
{
	cout << endl << "Usage:  MindTheGap find (-in <reads.fq> | -graph <graph.h5>) -ref <reference.fa> [options]" << endl;
	OptionsHelpVisitor v(cout);
	getParser()->accept(v);
	throw Exception(); // to get out with EXIT_FAILURE
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
void Finder::execute ()
{
    


    // Checks mandatory options
    if ((getInput()->get(STR_URI_GRAPH) != 0 && getInput()->get(STR_URI_INPUT) != 0) || (getInput()->get(STR_URI_GRAPH) == 0 && getInput()->get(STR_URI_INPUT) == 0))
    {
        throw OptionFailure(getParser(), "ERROR: options -graph and -in are incompatible, but at least one of these is mandatory");
    }
    
    if (getInput()->get(STR_URI_REF) == 0){
    	throw OptionFailure(getParser(), "ERROR: option -ref is mandatory");
    }


    // If outputPrefix is not provided we create one using the current date-time
    if (getInput()->get(STR_URI_OUTPUT) == 0)
    {
        time_t     now = time(0);
        struct tm  tstruct;
        char       buf[80];
        tstruct = *localtime(&now);
        strftime(buf, sizeof(buf), "%Y-%m-%d.%I:%M", &tstruct);
        string outputPrefix="MindTheGap_Expe-"+string(buf);
        
        getInput()->add (0, STR_URI_OUTPUT, outputPrefix);
    }
    

    // Getting the graph
	

    // Case 1 : -in option, we create the graph from read files
    if (getInput()->get(STR_URI_INPUT) != 0)
    {
        //fprintf(log,"Creating the graph from file(s) %s\n",getInput()->getStr(STR_URI_INPUT).c_str());
        
        // We need to add the options of dbgh5/Graph that were masked to the user
    	getInput()->add(0,STR_SOLIDITY_KIND, "sum"); //way to consider a solid kmer with several datasets (sum, min or max)

        getInput()->add(0,STR_BANK_CONVERT_TYPE,"tmp");
        getInput()->add(0,STR_URI_OUTPUT_DIR, ".");
		//getInput()->add(0,STR_URI_OUTPUT_TMP, ".");
        getInput()->add(0,STR_BLOOM_TYPE, "basic"); //neighbor basic cache
        getInput()->add(0,STR_DEBLOOM_TYPE, "original"); //DO NOT use cascading : generates too many FP inside  pas bien car bcp plus de FP non critique au milieur trou
        getInput()->add(0,STR_DEBLOOM_IMPL, "basic"); //minimizer => STR_BLOOM_TYPE = neighbor
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_INTEGER_PRECISION, "0");
        getInput()->add(0,STR_MPHF_TYPE, "none");
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_MINIMIZER_SIZE, "8");
        getInput()->add(0,STR_REPARTITION_TYPE, "0");
        getInput()->add(0,STR_MINIMIZER_TYPE, "0");
        getInput()->add(0,STR_HISTOGRAM_MAX, "10000");
        getInput()->add(0,STR_KMER_ABUNDANCE_MIN_THRESHOLD,"3");
        //getInput()->add(0,STR_URI_SOLID_KMERS, ""); //DONOT uncomment this line, otherwise solid kmers are stored in file ./.h5 outside from the considered graph h5 file.
        
        //Warning if kmer size >128 cascading debloom does not work
        if(getInput()->getInt(STR_KMER_SIZE)>128){
            getInput()->get(STR_DEBLOOM_TYPE)->value="original";
        }
        
        //de Bruijn graph building
        _graph = Graph::create (getInput());

		
        _kmerSize = getInput()->getInt(STR_KMER_SIZE);
        
    }
    
    // Case 2 : -graph option, we load the graph from a .h5 file
    if (getInput()->get(STR_URI_GRAPH) != 0)
    {
        //fprintf(log,"Loading the graph from file %s\n",getInput()->getStr(STR_URI_GRAPH).c_str());
        _graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));
        _kmerSize = _graph.getKmerSize();
    }


    // Preparing the output files
    _breakpoint_file_name = getInput()->getStr(STR_URI_OUTPUT)+".breakpoints";
    _breakpoint_file = fopen(_breakpoint_file_name.c_str(), "w");
    if(_breakpoint_file == NULL){
        //cerr <<" Cannot open file "<< _output_file <<" for writting" << endl;
        string message = "Cannot open file "+ _breakpoint_file_name + " for writting";
        throw Exception(message.c_str());
    }

    _vcf_file_name = getInput()->getStr(STR_URI_OUTPUT)+".othervariants.vcf";
    _vcf_file = fopen(_vcf_file_name.c_str(), "w");
    if(_vcf_file == NULL){
    	//cerr <<" Cannot open file "<< _output_file <<" for writting" << endl;
    	string message = "Cannot open file "+ _vcf_file_name + " for writting";
    	throw Exception(message.c_str());
    }
    writeVcfHeader();

    // Getting the reference genome
    //_refBank = new BankFasta(getInput()->getStr(STR_URI_REF));
    _refBank = Bank::open(getInput()->getStr(STR_URI_REF)); // more general can be a list or a file of files
    _refBank->use(); //to be able to use the bank several times (do not forget at the end to do _refBank->forget() = delete)
    
    //Getting other parameters
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_repeat = getInput()->getInt(STR_MAX_REPEAT);
    _het_max_occ=getInput()->getInt(STR_HET_MAX_OCC);
    _snp_min_val=getInput()->getInt(STR_SNP_MIN_VAL);

    if(_het_max_occ<1){
    	_het_max_occ=1;
    }

    if(getInput()->get(STR_HOMO_ONLY) != 0)
    {
	_homo_only = true;
	_homo_insert = true;
	_hete_insert = false;
	_snp = true;
	_backup = false;
	_deletion = true;
    }
    
    if(getInput()->get(STR_INSERT_ONLY) != 0)
    {
	_homo_only = false;
	_homo_insert = true;
	_hete_insert = true;
	_snp = false;
	_backup = false;
	_deletion = false;
    }

    if(getInput()->get(STR_SNP_ONLY) != 0)
    {
	_homo_only = true;
	_homo_insert = false;
	_hete_insert = false;
	_snp = true;
	_backup = false;
	_deletion = false;
    }

    if(getInput()->get(STR_DELETION_ONLY) != 0)
    {
	_homo_only = true;
	_homo_insert = false;
	_hete_insert = false;
	_snp = false;
	_backup = false;
	_deletion = true;
    }

    if(getInput()->get(STR_HETERO_ONLY) != 0)
    {
	_homo_only = false;
	_homo_insert = false;
	_hete_insert = true;
	_snp = false;
	_backup = false;
	_deletion = false;
    }

    if(getInput()->get(STR_WITH_BACKUP) != 0)
    {
	_backup = true;
    }

    if(getInput()->get(STR_NO_SNP) != 0)
    {
	_snp = false;
    }

    if(getInput()->get(STR_NO_INSERT) != 0)
    {
	_homo_insert = false;
    }

    if(getInput()->get(STR_NO_DELETION) != 0)
    {
	_deletion = false;
    }

    if(getInput()->get(STR_NO_HETERO) != 0)
    {
	_hete_insert = false;
    }

    // Now do the job
    time_t start_time = time(0);
    // According to the kmer size,  we call one fillBreakpoints method.
    Integer::apply<runFindBreakpoints,Finder*> (_kmerSize, this);
    time_t end_time = time(0);
    double seconds=difftime(end_time,start_time);

    //cout << "in MTG" <<endl;
    // We gather some statistics.
    fclose(_breakpoint_file);
    fclose(_vcf_file);

    // Printing result informations (ie. add info to getInfo(), in Tool Info is printed automatically after end of execute() method
    getInfo()->add(1,"version",_mtg_version);
    getInfo()->add(1,"gatb-core-library",STR_LIBRARY_VERSION);
    getInfo()->add(1,"supported_kmer_sizes","%s", KSIZE_STRING);
    //getInfo()->add (1, &LibraryInfo::getInfo());
    resumeParameters();
    resumeResults(seconds);
}

void Finder::resumeParameters(){
    
    getInfo()->add(0,"Parameters");
    getInfo()->add(1,"Input data");
    if (getInput()->get(STR_URI_INPUT) != 0){
        getInfo()->add(2,"Reads",getInput()->getStr(STR_URI_INPUT).c_str());
    }
    if (getInput()->get(STR_URI_GRAPH) != 0){
        getInfo()->add(2,"Graph",getInput()->getStr(STR_URI_GRAPH).c_str());
    }
    getInfo()->add(2,"Reference",getInput()->getStr(STR_URI_REF).c_str());
    getInfo()->add(1,"Graph");
    getInfo()->add(2,"kmer-size","%i", _kmerSize);

    //In MindTheGap, solidity-kind always at "sum" (not tunable)
    //getInfo()->add(2,"solidity_kind",_graph.getInfo().getStr("solidity_kind").c_str());
    try { // use try/catch because this key is present only if auto asked
    	getInfo()->add(2,"abundance_min (auto inferred)",_graph.getInfo().getStr("cutoffs_auto.values").c_str());
    } catch (Exception e) {
    	// doing nothing
    }
    string min_abundance;
    //if(_graph.getInfo().getStr("solidity_kind")=="sum"){
    int thre = _graph.getInfo().getInt("thresholds"); //with getInt obtains the first number (if sum : threshold = 4 4 if -in had 2 input files, less confusing for the user if only one value shown)
    stringstream ss;
    ss << thre;
    min_abundance = ss.str();
    //    }
    //    else{
    //    	min_abundance = _graph.getInfo().getStr("thresholds").c_str();
    //    }
    getInfo()->add(2,"abundance_min (used)",min_abundance);

    try { // version actuelle info manquante si -graph
    	getInfo()->add(2,"abundance_max",_graph.getInfo().getStr("abundance_max").c_str());
    } catch (Exception e) {
    	// doing nothing
    }
    try { // entour try/catch ici au cas ou le nom de la cle change dans gatb-core
    	getInfo()->add(2,"nb_solid_kmers",_graph.getInfo().getStr("kmers_nb_solid").c_str());
    	getInfo()->add(2,"nb_branching_nodes",_graph.getInfo().getStr("nb_branching").c_str());
    } catch (Exception e) {
    	// doing nothing
    }

    getInfo()->add(1,"Breakpoint detection options");
    getInfo()->add(2,"max_repeat","%i", _max_repeat);
    getInfo()->add(2,"hetero_max_occ","%i", _het_max_occ);
    getInfo()->add(2,"homo_insertions","%s", _homo_insert ? "yes" : "no");
    getInfo()->add(2,"hete_insertions","%s", _hete_insert ? "yes" : "no");
    getInfo()->add(2,"snp","%s", _snp ? "yes" : "no");
    //getInfo()->add(2,"backup","%s", _backup ? "yes" : "no");
    getInfo()->add(2,"deletion","%s", _deletion ? "yes" : "no");    
}

void Finder::resumeResults(double seconds){
    getInfo()->add(0,"Results");
    getInfo()->add(1,"Insertion breakpoints");
    getInfo()->add(2,"homozygous","%i", _nb_homo_clean+_nb_homo_fuzzy);
    getInfo()->add(3,"clean","%i", _nb_homo_clean);
    getInfo()->add(3,"fuzzy","%i", _nb_homo_fuzzy);
    getInfo()->add(2,"heterozygous","%i", _nb_hetero_clean+_nb_hetero_fuzzy);
    getInfo()->add(3,"clean","%i", _nb_hetero_clean);
    getInfo()->add(3,"fuzzy","%i", _nb_hetero_fuzzy);
    getInfo()->add(1,"Other variants");
    getInfo()->add(2,"deletions","%i", _nb_clean_deletion+_nb_fuzzy_deletion);
    //getInfo()->add(3,"clean", "%i", _nb_clean_deletion);
    //getInfo()->add(3,"fuzzy", "%i", _nb_fuzzy_deletion);
    getInfo()->add(2,"SNPs","%i", _nb_solo_snp+_nb_multi_snp);
    //getInfo()->add(3,"isolated", "%i", _nb_solo_snp);
    //getInfo()->add(3,"close", "%i", _nb_multi_snp);
    //getInfo()->add(2,"backup","%i", _nb_backup);
    getInfo()->add(1,"Time", "%.1f s",seconds);
    getInfo()->add(1,"Output files");
    if(getInput()->get(STR_URI_INPUT) != 0){
            getInfo()->add(2,"graph_file", "%s.h5",getInput()->getStr(STR_URI_OUTPUT).c_str());
        }
    getInfo()->add(2,"breakpoint_file","%s",_breakpoint_file_name.c_str());
    getInfo()->add(2,"othervariants_file","%s",_vcf_file_name.c_str());


}

void Finder::writeVcfHeader(){

	string sample="";
	if (getInput()->get(STR_URI_INPUT) != 0){
		sample = getInput()->getStr(STR_URI_INPUT);
	}
	if (getInput()->get(STR_URI_GRAPH) != 0){
		sample=getInput()->getStr(STR_URI_GRAPH);
	}

	//getting the date
	time_t current_time;
	char* c_time_string;
	current_time = time(NULL);
	c_time_string = ctime(&current_time);

	fprintf(_vcf_file,
			"##fileformat=VCFv4.1\n\
##filedate=%s\
##source=MindTheGap find version %s\n\
##SAMPLE=file:%s\n\
##REF=file:%s\n\
##INFO=<ID=TYPE,Number=1,Type=String,Description=\"SNP, INS, DEL or .\">\n\
##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"variant size\">\n\
##INFO=<ID=FUZZY,Number=1,Type=Integer,Description=\"repeat size at the breakpoint, only for INS and DEL\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tG1\n",
c_time_string, _mtg_version, sample.c_str(),getInput()->getStr(STR_URI_REF).c_str());
}

template<size_t span>
void Finder::runFindBreakpoints<span>::operator ()  (Finder* object)
{
	FindBreakpoints<span> findBreakpoints(object);
	
	/* Add Gap observer */
	if(object->_snp)
	{
		findBreakpoints.addGapObserver(new FindSoloSNP<span>(&findBreakpoints));
		//findBreakpoints.addGapObserver(new FindFuzzySNP<span>(&findBreakpoints));
		findBreakpoints.addGapObserver(new FindMultiSNP<span>(&findBreakpoints));
		findBreakpoints.addGapObserver(new FindMultiSNPrev<span>(&findBreakpoints));

	}
	
	if(object->_deletion)
	{
		findBreakpoints.addGapObserver(new FindDeletion<span>(&findBreakpoints));
	}
	
	if(object->_homo_insert)
	{
		findBreakpoints.addGapObserver(new FindCleanInsertion<span>(&findBreakpoints));
		findBreakpoints.addGapObserver(new FindFuzzyInsertion<span>(&findBreakpoints));
	}
	
	if(object->_backup)
	{
		findBreakpoints.addGapObserver(new FindBackup<span>(&findBreakpoints));
	}
	
	/* Add kmer observer*/
	if(object->_hete_insert)
	{
		findBreakpoints.addKmerObserver(new FindHeteroInsertion<span>(&findBreakpoints));
	}
	
	/* Run */
	findBreakpoints();
}
