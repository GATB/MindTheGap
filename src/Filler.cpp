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

#include <Filler.hpp>
#include <Utils.hpp>
#include <GraphAnalysis.hpp>
#include <limits> // for numeric_limits
#include <unordered_map>

#define PRINT_DEBUG
/********************************************************************************/

// We define some constant strings for names of command line parameters
//static const char* STR_FOO = "-foo";


void HelpFiller(void* target)
{
    if(target!=NULL)
    {
        Filler * obj = (Filler *) target;
        obj->FillerHelp();
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
Filler::Filler ()  : Tool ("MindTheGap fill") , _progress(0)
{

    //TODO rajouter les parametres
    _nb_mis_allowed = 2;
    _nb_gap_allowed = 0;
    _nb_breakpoints = 0;
    _nb_filled_breakpoints = 0;
    _nb_multiple_fill = 0;


    setHelp(&HelpFiller);
    setHelpTarget(this);


    setParser (new OptionsParser ("MindTheGap fill"));

    IOptionsParser* generalParser = new OptionsParser("General");
    //generalParser->push_front (new OptionNoParam (STR_HELP, "help", false));
    // generalParser->push_back (new OptionNoParam (STR_VERSION, "version", false)); // move this option in the main.cpp
    generalParser->push_front (new OptionOneParam (STR_VERBOSE,     "verbosity level",      false, "1"  ));
    generalParser->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory for graph building (in MBytes)", false, "2000"));
    generalParser->push_front (new OptionOneParam (STR_MAX_DISK, "max disk for graph building   (in MBytes)", false, "0"));
    generalParser->push_front (new OptionOneParam (STR_NB_CORES,    "number of cores",      false, "0"  ));

    IOptionsParser* inputParser = new OptionsParser("Input / output");
    inputParser->push_front (new OptionOneParam (STR_CONTIG_OVERLAP, "Overlap between input contigs",  false, "0"));
    inputParser->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_BKPT, "breakpoint file", false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_CONTIG, "contig file", false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_GRAPH, "input graph file (likely a hdf5 file)",  false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_INPUT, "input read file(s)",  false, ""));

    IOptionsParser* fillerParser = new OptionsParser("Assembly");
    //TODO HERE PUT THE FILL OPTIONS
    fillerParser->push_front (new OptionOneParam (STR_MAX_DEPTH, "maximum length of insertions (nt)", false, "10000"));
    fillerParser->push_front (new OptionOneParam (STR_MAX_NODES, "maximum number of nodes in contig graph (nt)", false, "100"));

    IOptionsParser* graphParser = new OptionsParser("Graph building");
    string abundanceMax = Stringify::format("%ld", std::numeric_limits<CountNumber>::max()); //to be sure in case CountNumber definition changes
    graphParser->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "maximal abundance threshold for solid kmers", false, abundanceMax));
    graphParser->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "minimal abundance threshold for solid kmers", false, "auto"));
    graphParser->push_front (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));


    getParser()->push_front(generalParser);
    getParser()->push_front(fillerParser);
    getParser()->push_front(graphParser);
    getParser()->push_front(inputParser);

}



void Filler::FillerHelp()
{
    cout << endl << "Usage:  MindTheGap fill (-in <reads.fq> | -graph <graph.h5>) -bkpt <breakpoints.fa or -contig <contig.fa> [options]" << endl;
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
void Filler::execute ()
{


    if ((getInput()->get(STR_URI_GRAPH) != 0 && getInput()->get(STR_URI_INPUT) != 0) || (getInput()->get(STR_URI_GRAPH) == 0 && getInput()->get(STR_URI_INPUT) == 0))
    {

        throw OptionFailure(getParser(), "options -graph and -in are incompatible, but at least one of these is mandatory");

    }

    if (getInput()->get(STR_URI_BKPT) == 0 && getInput()->get(STR_URI_CONTIG) == 0)
        {
                throw OptionFailure(getParser(), "option -bkpt or -contig is mandatory");
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

        // We need to add the options of dbgh5/Graph that were masked to the user (or we could create a new Properties object)
        getInput()->add(0,STR_SOLIDITY_KIND, "sum"); //way to consider a solid kmer with several datasets (sum, min or max)

        getInput()->add(0,STR_BANK_CONVERT_TYPE,"tmp");
        getInput()->add(0,STR_URI_OUTPUT_DIR, ".");
        getInput()->add(0,STR_URI_OUTPUT_TMP, ".");
        getInput()->add(0,STR_BLOOM_TYPE, "basic"); //neighbor basic cache
        getInput()->add(0,STR_DEBLOOM_TYPE, "original"); //cascading  pas bien car bcp plus de FP non critique au milieur trou
        getInput()->add(0,STR_DEBLOOM_IMPL, "basic"); //minimizer => STR_BLOOM_TYPE = neighbor
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_INTEGER_PRECISION, "0");
       // getInput()->add(0,STR_MPHF_TYPE, "none");
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_MINIMIZER_SIZE, "8");
        getInput()->add(0,STR_REPARTITION_TYPE, "0");
        getInput()->add(0,STR_MINIMIZER_TYPE, "0");
        getInput()->add(0,STR_HISTOGRAM_MAX, "10000");
        getInput()->add(0,STR_KMER_ABUNDANCE_MIN_THRESHOLD,"3");
        getInput()->add(0,STR_STORAGE_TYPE,"hdf5");
        //getInput()->add(0,STR_URI_SOLID_KMERS, ""); //surtout ne pas decommenter cette ligne, sinon les kmers solids sont stockes dans le fichier ./.h5 et les infos ne sont plus dans le output.h5

        //Warning if kmer size >128 cascading debloom does not work
        if(getInput()->getInt(STR_KMER_SIZE)>128){
            getInput()->get(STR_DEBLOOM_TYPE)->value="original";
        }

        _graph = Graph::create (getInput());
        _kmerSize = getInput()->getInt(STR_KMER_SIZE);

    }

    // Case 2 : -graph option, we load the graph from a .h5 file
    if (getInput()->get(STR_URI_GRAPH) != 0)
    {
        //fprintf(log,"Loading the graph from file %s\n",getInput()->getStr(STR_URI_GRAPH).c_str());

    	fprintf(stderr,"Loading the graph..."); //TODO better a progress bar
    	fflush(stderr);
        _graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));
        _kmerSize = _graph.getKmerSize();
        fprintf(stderr,"done\n");
        fflush(stderr);
    }


    // Output file names

    _insert_file_name = getInput()->getStr(STR_URI_OUTPUT)+".insertions.fasta";
    _insert_file = fopen(_insert_file_name.c_str(), "w");
    if (getInput()->get(STR_URI_CONTIG) != nullptr)
    {
        _gfa_file_name = getInput()->getStr(STR_URI_OUTPUT)+".gfa";
        _gfa_file = fopen(_gfa_file_name.c_str(),"w");
        if(_gfa_file == NULL){
            string message = "Cannot open file "+ _gfa_file_name + " for writting";
            throw Exception(message.c_str());
        }
    }
    if(_insert_file == NULL){
        string message = "Cannot open file "+ _insert_file_name + " for writting";
        throw Exception(message.c_str());
    }

    _insert_info_file_name = getInput()->getStr(STR_URI_OUTPUT)+".info.txt";
    _insert_info_file = fopen(_insert_info_file_name.c_str(), "w");
    if(_insert_info_file == NULL){
        string message = "Cannot open file "+ _insert_info_file_name + " for writting";
        throw Exception(message.c_str());
    }
	
    if (getInput()->get(STR_URI_BKPT) != nullptr)
    {
        _vcf_file_name = getInput()->getStr(STR_URI_OUTPUT)+".insertions.vcf";
        _vcf_file = fopen(_vcf_file_name.c_str(), "w");
        if(_vcf_file == NULL){
            //cerr <<" Cannot open file "<< _output_file <<" for writting" << endl;
            string message = "Cannot open file "+ _vcf_file_name + " for writting";
            throw Exception(message.c_str());
        }
        writeVcfHeader();
    }
    
    //Getting the breakpoint sequences
    if (getInput()->get(STR_URI_BKPT) != nullptr)
       {
           _breakpointBank = new BankFasta(getInput()->getStr(STR_URI_BKPT));

       }
    if (getInput()->get(STR_URI_CONTIG) != nullptr)
       {
           _breakpointBank = new BankFasta(getInput()->getStr(STR_URI_CONTIG));

        }
    
    //Getting other parameters
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_depth = getInput()->getInt(STR_MAX_DEPTH);
    _max_nodes = getInput()->getInt(STR_MAX_NODES);
    // Now do the job
    time_t start_time = time(0);

    Integer::apply<fillAny,Filler*> (_kmerSize, this);
    time_t end_time = time(0);

    //job is done, closing the output files
    fclose(_insert_file);
    fclose(_insert_info_file);

    if (getInput()->get(STR_URI_BKPT) != nullptr)
    {
        fclose(_vcf_file);
    }
    if (getInput()->get(STR_URI_CONTIG) != nullptr)
    {
        fclose(_gfa_file);
    }
    
    // We gather some info/statistics to print in stdout
    resumeParameters();

    double seconds=difftime(end_time,start_time);
    resumeResults(seconds);
}

void Filler::writeVcfHeader(){

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
##source=MindTheGap fill version %s\n\
##SAMPLE=file:%s\n\
##REF=file:%s\n\
##INFO=<ID=TYPE,Number=1,Type=String,Description=\"INS\">\n\
##INFO=<ID=LEN,Number=1,Type=Integer,Description=\"variant size\">\n\
##INFO=<=AVK,Number=.,Type=Float,Description=\"Average k-mer coverage along the insertion\">\n\
##INFO=<=MDK,Number=.,Type=Float,Description=\"Median k-mer coverage along the insertion\">\n\
##INFO=<ID=GTT,Number=1,Type=String,Description=\"Unphased  genotypes\">\n\
##INFO=<ID=NPOS,Number=1,Type=Integer,Description=\"number of alternative positions for the insertion site (= size of repeat (fuzzy) +1)\">\n\
##FORMAT=<ID=GT,Number=1,Type=String,Description=\"Genotype\">\n\
#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tG1\n",
c_time_string, _mtg_version, sample.c_str(),getInput()->getStr(STR_URI_OUTPUT).c_str());
    
    //TODO later :
    //##INFO=<=QUA,Number=.,Type=Integer,Description=\"Quality of the insertion\">\n\
    //##INFO=<=NS,Number=1,Type=String,Description=\"Solution index\">\n\

}

void Filler::resumeParameters(){

    getInfo()->add(0,"MindTheGap fill");

    //getInfo()->add(1,"version",getVersion());
    getInfo()->add(1,"version",_mtg_version);
    getInfo()->add(1,"gatb-core-library",System::info().getVersion().c_str());
    getInfo()->add(1,"supported_kmer_sizes","%s", KSIZE_STRING);
    //getInfo()->add (1, &LibraryInfo::getInfo());
    
    getInfo()->add(0,"Parameters");
    getInfo()->add(1,"Input data");
    if (getInput()->get(STR_URI_INPUT) != 0){
        getInfo()->add(2,"Reads",getInput()->getStr(STR_URI_INPUT).c_str());
    }
    if (getInput()->get(STR_URI_GRAPH) != 0){
        getInfo()->add(2,"Graph",getInput()->getStr(STR_URI_GRAPH).c_str());
    }
    
    if (getInput()->get(STR_URI_BKPT) != nullptr)
    {
        getInfo()->add(2,"Breakpoints",getInput()->getStr(STR_URI_BKPT).c_str());
    }
    if (getInput()->get(STR_URI_CONTIG) != nullptr)
    {
        getInfo()->add(2,"Contigs",getInput()->getStr(STR_URI_CONTIG).c_str());
    }
    getInfo()->add(1,"Graph");
    getInfo()->add(2,"kmer-size","%i", _kmerSize);
    try { // use try/catch because this key is present only if auto asked
        getInfo()->add(2,"abundance_min (auto inferred)",_graph.getInfo().getStr("cutoffs_auto.values").c_str());
    } catch (Exception e) {
        // doing nothing
    }
    string min_abundance;
    int thre = _graph.getInfo().getInt("thresholds"); //with getInt obtains the first number (if sum : threshold = 4 4 if -in had 2 input files, less confusing for the user if only one value shown)
    stringstream ss;
    ss << thre;
    min_abundance = ss.str();
    getInfo()->add(2,"abundance_min (used)",min_abundance);

    try {     // entour try/catch ici au cas ou le nom de la cle change dans gatb-core
            getInfo()->add(2,"nb_solid_kmers",_graph.getInfo().getStr("kmers_nb_solid").c_str());
            getInfo()->add(2,"nb_branching_nodes",_graph.getInfo().getStr("nb_branching").c_str());
        } catch (Exception e) {
            // doing nothing
        }

    getInfo()->add(1,"Assembly options");
    getParser()->push_front (new OptionOneParam (STR_MAX_DEPTH, "maximum length of insertions (nt)", false, "10000"));
    getInfo()->add(2,"max_depth","%i", _max_depth);
    getInfo()->add(2,"max_nodes","%i", _max_nodes); //todo

}

void Filler::resumeResults(double seconds){

    getInfo()->add(0,"Results");
    getInfo()->add(1,"Breakpoints");
    getInfo()->add(2,"nb_input","%i", _nb_breakpoints);
    getInfo()->add(2,"nb_filled","%i", _nb_filled_breakpoints);
    getInfo()->add(3,"unique_sequence","%i", _nb_filled_breakpoints-_nb_multiple_fill);
    getInfo()->add(3,"multiple_sequence","%i", _nb_multiple_fill);
    getInfo()->add(1,"Time", "%.1f s",seconds);
    getInfo()->add(1,"Output files");
    getInfo()->add(2,"assembled sequence file","%s",_insert_file_name.c_str());
    if (getInput()->get(STR_URI_CONTIG) != nullptr)
    {
        getInfo()->add(2,"assembly graph file","%s",_gfa_file_name.c_str());
    }
    if (getInput()->get(STR_URI_BKPT) != nullptr)
    {
        getInfo()->add(2,"insertion variant vcf file","%s",_vcf_file_name.c_str());
    }
    getInfo()->add(2,"assembly statistics file","%s",_insert_info_file_name.c_str());

}


template<size_t span>
class contigFunctor
{

    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type  Type;

public:
    void operator() (Sequence& sequence)
    {

    string sourceSequence = string(sequence.getDataBuffer(),sequence.getDataSize());
    string seedk = sourceSequence;

    string seedName = sequence.getComment();
    string infostring;
    bool isRc = !seedName.compare (seedName.length() - 3, 3, "_Rc");

    //cout << seedName << endl;

    //bool begin_kmer_repeated = false;
    //bool end_kmer_repeated = false;
    bool is_anchor_repeated = false;
    bool reverse = false;
    string conc_targetSequence;
    bkpt_dict_t targetDictionary;
    int nb_mis_allowed = 2; // To fix

    for (auto its=_all_targetDictionary.begin(); its !=_all_targetDictionary.end(); ++its)
    {
        string& targetName= its->second.first;
        string tempName = targetName;
        if(its->second.second) tempName += "_Rc";
        if (tempName.compare(seedName) != 0 )
        {
            conc_targetSequence.append(its->first);
            targetDictionary.insert({its->first,its->second});
        }
     }
        
    std::vector<filled_insertion_t> filledSequences;
     _object->contigGapFill<span>(infostring,_tid, sourceSequence, conc_targetSequence,filledSequences, nb_mis_allowed, targetDictionary, false );
        
     infostring +=   Stringify::format ("\t%d", filledSequences.size()) ;


     // Write insertions to file
     _object->writeFilledBreakpoint(filledSequences,seedName,infostring,is_anchor_repeated,false);
     _object->writeToGFA(filledSequences,sourceSequence,seedName,isRc,is_anchor_repeated);
        

     _nb_breakpoints++;

     //progress bar
     _nbBreakpointsProgressDone++;
     if (_nbBreakpointsProgressDone > 0)   {  _object->_progress->inc (_nbBreakpointsProgressDone);  _nbBreakpointsProgressDone = 0;  }
    }


    //constructor
    contigFunctor(Filler* object, int * nb_living, int * global_nb_breakpoints, bkpt_dict_t all_targetDictionary) : _object(object),_global_nb_breakpoints(global_nb_breakpoints),_all_targetDictionary(all_targetDictionary)
    {
        _nb_living =nb_living;
        _tid =  __sync_fetch_and_add (_nb_living, 1);
        _nb_breakpoints = 0;
        //printf("creating thread id %i \n",_tid);

    }

    contigFunctor(contigFunctor const &r)
    {
        _all_targetDictionary = r._all_targetDictionary;
        _nb_living = r._nb_living;
        _object= r._object;
        _global_nb_breakpoints = r._global_nb_breakpoints;
        _nb_breakpoints = 0;
        _tid =  __sync_fetch_and_add (_nb_living, 1);
        
        //printf("CC creating thread id %i \n",_tid);

    }

    ~contigFunctor()
    {
         __sync_fetch_and_add (_global_nb_breakpoints, _nb_breakpoints);
    }
private:
    Filler* _object;
    int _tid;
    int _nb_breakpoints;
    int * _global_nb_breakpoints;
    int * _nb_living;

    Sequence _previousSeq;
    u_int64_t _nbBreakpointsProgressDone = 0;
    bkpt_dict_t _all_targetDictionary;
};


template<size_t span>
class gapfillerFunctor
{
    typedef typename gatb::core::kmer::impl::Kmer<span>::ModelCanonical ModelCanonical;
    typedef typename Kmer<span>::Count Count;
    typedef typename Kmer<span>::Type  Type;

public:
    void operator() (Sequence& sequence)
    {
        if( (sequence.getIndex() & 1) == 0)
        {
            _previousSeq = sequence; //first sequence (source) of the breakpoint
        }
        else //second sequence (target) of the breakpoint
        {
            //compute pair
            //	printf("%s  id %i  %s  id %i   tid %i  \n",_previousSeq.getCommentShort().c_str(), _previousSeq.getIndex()
            //		   ,sequence.getCommentShort().c_str(), sequence.getIndex(), _tid);

            string sourceSequence =  string(_previousSeq.getDataBuffer(),_previousSeq.getDataSize());//previously L
            string seedk = sourceSequence; //used for coverage computation

            string breakpointName = string(_previousSeq.getCommentShort());

            string infostring; //to store some statitistics about the gap-filling process

            bool begin_kmer_repeated = _previousSeq.getComment().find("REPEATED") !=  std::string::npos;


            string targetSequence =  string(sequence.getDataBuffer(),sequence.getDataSize());//previously R

            string breakpointName_R = string(sequence.getCommentShort());
            bool end_kmer_repeated = sequence.getComment().find("REPEATED") !=  std::string::npos;
            bool is_anchor_repeated = begin_kmer_repeated || end_kmer_repeated;

            int nb_mis_allowed = _nb_mis_allowed;

            if(is_anchor_repeated){
                nb_mis_allowed=0;
            }
            //printf("nb_mis_allowed %i \n",nb_mis_allowed);

            //Initialize set of filled sequences
            std::vector<filled_insertion_t> filledSequences;
            bkpt_dict_t targetDictionary;
            // If Source and Target sequences are larger than kmer-size, resize to kmer-size :
            if(sourceSequence.size()> _object->_kmerSize ){
                sourceSequence.substr(sourceSequence.size()- _object->_kmerSize,_object->_kmerSize); //suffix of size _kmerSize
            }

            if(targetSequence.size()>_object->_kmerSize){
                targetSequence.substr(0,_object->_kmerSize); //prefix of size _kmerSize
            }
            targetDictionary.insert ({targetSequence, std::make_pair(breakpointName_R, false)});

            //_object->gapFill<span>(infostring,_tid,sourceSequence,targetSequence,filledSequences,begin_kmer_repeated,end_kmer_repeated);
            _object->contigGapFill<span>(infostring,_tid, sourceSequence, targetSequence,filledSequences, nb_mis_allowed, targetDictionary, false);

            //If gap-filling failed in one direction, try the other direction (from target to source in revcomp)
            if(filledSequences.size()==0){
                string targetSequence2 = revcomp_sequence(sourceSequence);
                targetDictionary.clear();
                targetDictionary.insert({targetSequence2, std::make_pair(breakpointName, false)});
                string sourceSequence2 = revcomp_sequence(targetSequence);
                breakpointName= breakpointName_R;


                //_object->GapFill<span>(infostring,_tid,sourceSequence2,targetSequence2,filledSequences,begin_kmer_repeated,end_kmer_repeated,true);
                _object->contigGapFill<span>(infostring,_tid, sourceSequence2, targetSequence2,filledSequences, nb_mis_allowed, targetDictionary, true);

            }

            infostring +=   Stringify::format ("\t%d", filledSequences.size()) ;

            _object->writeFilledBreakpoint(filledSequences,breakpointName,infostring,is_anchor_repeated,true);
            _object->writeVcf(filledSequences,breakpointName,seedk);

            
            // We increase the breakpoint counter.
            _nb_breakpoints++;

            //progress bar
            _nbBreakpointsProgressDone++;
            if (_nbBreakpointsProgressDone > 50)   {  _object->_progress->inc (_nbBreakpointsProgressDone);  _nbBreakpointsProgressDone = 0;  }
        }

    }

    //constructor
    gapfillerFunctor(Filler* object, int * nb_living, int * global_nb_breakpoints) : _object(object),_global_nb_breakpoints(global_nb_breakpoints)
    {
        _nb_living =nb_living;
        _tid =  __sync_fetch_and_add (_nb_living, 1);
        _nb_breakpoints = 0;

        //printf("creating thread id %i \n",_tid);

    }

    gapfillerFunctor(gapfillerFunctor const &r)
    {
        _nb_living = r._nb_living;
        _object= r._object;
        _global_nb_breakpoints = r._global_nb_breakpoints;
        _nb_breakpoints = 0;
        _tid =  __sync_fetch_and_add (_nb_living, 1);

    //	printf("CC creating thread id %i \n",_tid);

    }

    ~gapfillerFunctor()
    {
         __sync_fetch_and_add (_global_nb_breakpoints, _nb_breakpoints);
    }
private:
    Filler* _object;

    int _tid;
    int _nb_breakpoints;
    int * _global_nb_breakpoints;
    int * _nb_living;
    int _nb_mis_allowed = 2; // To fix, should be read from global parameters
    Sequence _previousSeq;
    u_int64_t _nbBreakpointsProgressDone = 0;

};


template<size_t span>
struct Count2TypeAdaptor  {  typename Kmer<span>::Type& operator() (typename Kmer<span>::Count& c)  { return c.value; }  };
//template method : enabling to deal with all sizes of kmer <KSIZE_4

// fill use for both -contig and -bkpt
template<size_t span>
void Filler::fillAny<span>::operator () (Filler* object)
{
    size_t kmerSize = object->_kmerSize;

    BankFasta::Iterator itSeq (*object->_breakpointBank);
    int nbBreakpoints=0;

    if (object->getInput()->get(STR_URI_CONTIG) != nullptr){

        bkpt_dict_t seedDictionary;
        bkpt_dict_t all_targetDictionary;

        // seed sequences will be written on disk to create a seed Bank
        // Original contigs written as nodes of the GFA file
        ofstream seedFile;
        string seedFileName = object->getInput()->getStr(STR_URI_OUTPUT)+"_seed_dictionary.fasta";
        seedFile.open (seedFileName);
        int overlap = object->getInput()->getInt(STR_CONTIG_OVERLAP);
        for (itSeq.first(); !itSeq.isDone(); itSeq.next())
        {
            std::string seedSequence = string(itSeq->getDataBuffer(),itSeq->getDataSize());

            // Write the contigs to GFA
            // Contig has not been trimmed of overlap
            fprintf(object->_gfa_file,"S\t%s\t%s\n", itSeq->getComment().c_str(),seedSequence.c_str());

            // Remove overlap supplied as parameter
            // TODO : Remove small contigs
            seedSequence = seedSequence.substr(overlap, itSeq->getDataSize() - 2*overlap);


             // create seed and targetDictionary

            if(seedSequence.size() > (2*kmerSize)){
                std::string seedSequence_f = seedSequence.substr(seedSequence.size()-(2*kmerSize), kmerSize);
                std::string name = string(itSeq->getCommentShort());
                std::string targetSequence_f =seedSequence.substr(kmerSize,kmerSize);
                std::string seedSequence_Rc = revcomp_sequence(seedSequence);
                std::string seedSequence_Rc_f= seedSequence_Rc.substr(seedSequence_Rc.size() - (2 *kmerSize),kmerSize);
                std::string targetSequence_Rc_f = seedSequence_Rc.substr(kmerSize,kmerSize);
                seedDictionary.insert ({{seedSequence_f,std::make_pair(name, false)},
                                        {seedSequence_Rc_f,std::make_pair(name, true)}
                                       });
                all_targetDictionary.insert ({{targetSequence_f, std::make_pair(name, false)},
                                          {targetSequence_Rc_f, std::make_pair(name, true)}
                                         });

                // Output seed dictionary to a file
                seedFile << ">"+string(itSeq->getCommentShort())+"\n";
                seedFile << seedSequence_f << endl;
                seedFile << ">"+string(itSeq->getCommentShort())+"_Rc\n";
                seedFile << seedSequence_Rc_f << endl;


            }
            else{
                cerr << "Warning contig not used (too short): " << string(itSeq->getCommentShort()) << " of size "<< seedSequence.size() << " nt" << endl;
            }
        }
        seedFile.close();

        u_int64_t nbBreakpointsEstimated = seedDictionary.size() ;  // Number of seeds to iterate
        u_int64_t nbBreakpointsProgressDone = 0;

        object->setProgress (new ProgressSynchro (
                                              object->createIteratorListener (nbBreakpointsEstimated, "Filling the breakpoints"),
                                              System::thread().newSynchronizer())
                         );
        object->_progress->init ();


        BankFasta inbank (seedFileName);
        BankFasta::Iterator it (inbank);

        int nb_living=0;
        

        Dispatcher(object->getInput()->getInt(STR_NB_CORES)).iterate(it, contigFunctor<span>(object,&nb_living,&object->_nb_breakpoints,all_targetDictionary),30);

        object->_nb_breakpoints = object->_nb_breakpoints ;
        object->_progress->finish ();

    }


    if (object->getInput()->get(STR_URI_BKPT) != nullptr){
        u_int64_t nbBreakpointsEstimated = object->_breakpointBank->estimateNbItems()  / 2 ;  // 2 seq per breakpoint
        u_int64_t nbBreakpointsProgressDone = 0;

        object->setProgress (new ProgressSynchro (
                                          object->createIteratorListener (nbBreakpointsEstimated, "Filling breakpoints"),
                                          System::thread().newSynchronizer())
                     );
        object->_progress->init ();


        int nb_living=0;

        Dispatcher(object->getInput()->getInt(STR_NB_CORES)).iterate(itSeq, gapfillerFunctor<span>(object,&nb_living,&object->_nb_breakpoints),30);
        //WARNING : 30 = number of sequences sent to each thread, keep it *even* (2 sequences for one breakpoint)

        object->_nb_breakpoints = object->_nb_breakpoints ;
        object->_progress->finish ();


    }
}

template<size_t span>
void Filler::contigGapFill(std::string & infostring, int tid, string sourceSequence, string targetSequence, std::vector<filled_insertion_t>& filledSequences,int nb_mis_allowed, bkpt_dict_t targetDictionary, bool reverse ){
    typedef typename gatb::core::kmer::impl::Kmer<span>::ModelCanonical ModelCanonical;


    //object used to mark the traversed nodes of the graph (note : it is reset at the beginning of construct_linear_seq)
    BranchingTerminator terminator (_graph);
    IterativeExtensions<span> extension (_graph, terminator, TRAVERSAL_CONTIG, ExtendStopMode_until_max_depth, SearchMode_Breadth, false, _max_depth, _max_nodes);
    //todo check param dontOutputFirstNucl=false ??
    //todo put these two above lines in fillBreakpoints and pass object extension in param

    //prepare temporary file names
    string rev_str=""; //used to have distinct contig file names for forward and reverse extension (usefull if investigating one breakpoint and code lines with remove command are commented)
    if (reverse){
        rev_str="_rev";
    }

    string file_name_suffix=rev_str+ Stringify::format("%i",getpid()) + "_t" + Stringify::format("%i",tid); //to have unique file names when breakpoints are treated in parallel
    string contig_file_name = "contigs.fasta" + file_name_suffix;
    string contig_graph_file_prefix="contig_graph" + file_name_suffix;
    //std::cout << contig_file_name << std::endl;
    //std::cout << contig_graph_file_prefix << std::endl;

    //Build contigs and output them in a file in fasta format
    extension.construct_linear_seqs(sourceSequence,targetSequence,contig_file_name,true); //last param : swf will be true

    // connect the contigs into a graph and write the resulting graph in a temporary file
    GraphOutputDot<span> graph_output(_kmerSize,contig_graph_file_prefix);
    graph_output.load_nodes_extremities(contig_file_name,infostring);
    graph_output.first_id_els = graph_output.construct_graph(contig_file_name,"LEFT");
    graph_output.close();

    set< info_node_t > terminal_nodes_with_endpos = find_nodes_containing_multiple_R(targetDictionary, contig_file_name, nb_mis_allowed, _nb_gap_allowed);
    // printf("nb contig with target %zu \n",terminal_nodes_with_endpos.size());

    //also convert it to list of node id for traditional use
    set<int> terminal_nodes;
    for (set< info_node_t >::iterator it = terminal_nodes_with_endpos.begin(); it != terminal_nodes_with_endpos.end(); it++)
    {
        terminal_nodes.insert((*it).node_id);
    }

    // analyze the graph to find a satisfying gap sequence between L and R

    GraphAnalysis graph = GraphAnalysis(graph_output.get_dot_file_name(),_kmerSize);
    graph.debug = true;

    infostring +=   Stringify::format ("\t%d", terminal_nodes.size()) ;
    /*if(terminal_nodes.size()==0)
    {
        //if(verb)
            //printf("Right anchor not found.. gapfillling failed... \n");
        return ;
    }*/

    //if(verb)    fprintf(stderr," analysis..");
    bool success;
    set<pair<unlabeled_path,bkpt_t>> paths = graph.find_all_paths(terminal_nodes_with_endpos, success);
    

    // We build a map to sort paths leading to the same target
    unordered_map<string,set<unlabeled_path>> paths_to_compare;
    for (set<pair<unlabeled_path,bkpt_t>>::iterator it = paths.begin(); it!=paths.end();it++)
    {
        string key = it->second.first;
        if (it->second.second){
            key += "_Rc";
        }
        paths_to_compare[key].insert(it->first);
    }

    int nbTotal_filled_insertions = 0;
    int nbTotal_reported_insertions = 0; //after reducing redundancy of highly similar sequences
    for (unordered_map<string,set<unlabeled_path>>::iterator it = paths_to_compare.begin(); it!=paths_to_compare.end();++it)
    {
        std::vector<filled_insertion_t> tmpSequences;
        set<unlabeled_path> current_paths = it->second;
        tmpSequences = graph.paths_to_sequences(current_paths,terminal_nodes_with_endpos);
        
        int nb_filled_insertions = tmpSequences.size();
        nbTotal_filled_insertions += nb_filled_insertions;
        
        // remove almost identical sequences
        if (tmpSequences.size() > 1)
        {
            //if(verb)     printf(" [SUCCESS]\n");
            remove_almost_identical_solutions(tmpSequences,90);
        }
        
        int nb_reported_insertions = tmpSequences.size();
        nbTotal_reported_insertions += nb_reported_insertions;

        //Here add information for each filled insertion : coverage, quality, revcomp if reverse
        
        
        
        filledSequences.insert(filledSequences.end(),tmpSequences.begin(),tmpSequences.end());
     
    }
    
    if((nbTotal_filled_insertions>0) | reverse)
    {
        infostring +=   Stringify::format ("\t%d", nbTotal_filled_insertions) ;
    }
    
    
    /////////compute coverage of filled sequences // make sure to compute before reverse !!!
    ModelCanonical model (_kmerSize);
    for (std::vector<filled_insertion_t>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
    {
        typename ModelCanonical::Iterator itk (model);
        std::string cseq = sourceSequence + it->seq;
        Data data ((char*)cseq.c_str());
        
        itk.setData (data);
        std::vector<unsigned int> vec_abundances;
        
        u_int64_t sum = 0;
        int nbkmers =0;
        
        // We iterate the kmers of this seq
        for (itk.first(); !itk.isDone(); itk.next())
        {
            Node node(Node::Value(itk->value()));
            unsigned int cov = _graph.queryAbundance(node);
            //TODO add warning if cov==0
            if(cov==0){
                cerr << "WARNING Unknown kmer : " << model.toString(itk->value()) << endl;
            }
            sum+= cov; nbkmers++;
            vec_abundances.push_back(cov);
        }
        it->median_coverage = median(vec_abundances);
        it->avg_coverage  = sum /(float) nbkmers;
        
        if(reverse)
        {
            it->reverse();
        }
    }
    

    remove(contig_file_name.c_str());
    remove((contig_graph_file_prefix+".graph").c_str());

}


void Filler::writeFilledBreakpoint(std::vector<filled_insertion_t>& filledSequences,  string seedName, std::string info,bool is_anchor_repeated, bool breakpointMode){
    
    //printf("-- writeFilledBreakpoint --\n");
    
    //printf("found %zu seq \n",filledSequences.size());
    
    flockfile(_insert_file);
    
    int nbInsertions = 0;
    int nbTotalInsertions = 0;
    //cout << "enter filled insertion t" << endl;
    for (std::vector<filled_insertion_t>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
    {
        
        string insertion = it->seq;
        //cout << "Insert" << it->seq << endl;
        
        //std::cout <<"\n insertion" << insertion << std::endl;
        int llen = insertion.length() ;
        if(llen > 0) nbTotalInsertions++;
    }
    
    
    for (std::vector<filled_insertion_t>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
    {
        string insertion = it->seq;
        int llen = insertion.length() ;// - (int) R.length() - (int) L.length() - 2*hetmode;
        bkpt_t targetId = it->targetId_anchor;
        //printf("Insertion %i  %s \n",nbContig,insertion.c_str() );
        //discards insert too long (can happen when last node is very large)
        //		if(llen > max_insertions_size)
        //			continue;
        
        // save sequences to results file
        if(llen > 0)
        {
            
            std::ostringstream osolu_i;
            osolu_i <<   "solution " <<    nbInsertions+1 << "/" << nbTotalInsertions ;
            string solu_i = nbTotalInsertions >1 ?  osolu_i.str() : "" ;
            
            int qual = it->compute_qual(is_anchor_repeated);
            
            if(filledSequences.size()>1 && qual>10) qual = 15; // if multiple solutions and 0 err or repeated in ref
            
            //else if(filledSequences.size()>1 && qual ==50) qual = 2;
            
            
            //int bkptid;
            //parse bkpt header
            //get bkpt id
            //sscanf(breakpointName.c_str(),"bkpt%i*",&bkptid );
            //const char * end_header = strstr(breakpointName.c_str(), "kmer_");
            
            // bkpt%i insertion_len_%d_%s
            
            //Writting sequence header
            if(breakpointMode) //-bkpt mode, to keep the same header name as before
            {
                fprintf(_insert_file,">%s_len_%d_qual_%i_avg_cov_%.2f_median_cov_%.2f   %s\n",
                        seedName.c_str(),llen,qual,solu_i.c_str()
                        ,it->avg_coverage,it->median_coverage);
            }
            else{
                string targetName = targetId.first;
                if(targetId.second) {
                    targetName.append("_Rc");
                }
                int cov = it->median_coverage + 0.5;
                
                string insertionName = ">"+seedName+";"+targetName+";len_"+to_string(llen)+"_qual_"+to_string(qual)+"_median_cov_"+to_string(cov)+"\t"+solu_i+"\n";
                fprintf(_insert_file,"%s",insertionName.c_str());
            }
            
            //Writting DNA sequence
            fprintf(_insert_file,"%.*s\n",(int)llen,insertion.c_str() );
        }
        nbInsertions++;
    }
    
    if(nbInsertions>0){
        __sync_fetch_and_add(& _nb_filled_breakpoints,1);
        if(nbInsertions>1){
            __sync_fetch_and_add(& _nb_multiple_fill,1);
        }
    }
    
    funlockfile(_insert_file);
    
    
    //breakpoint name  info
    flockfile(_insert_info_file);
    
    fprintf(_insert_info_file,"%s\t%s\n",seedName.c_str(),info.c_str());
    
    funlockfile(_insert_info_file);
}

void Filler::writeVcf(std::vector<filled_insertion_t>& filledSequences, string breakpointName, string seedk){
    //warning : may not work if breakpoint file is not formated by MindTheGap find (to find info in the sequence headers (such as CHR, pos...)
    
    flockfile(_vcf_file);
    
    for (std::vector<filled_insertion_t>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
    {
        string insertion = it->seq;
        int llen = insertion.length() ;// - (int) R.length() - (int) L.length() - 2*hetmode;
        if(llen > 0)
        {
            vector<char> left(seedk.begin(), seedk.end());
            vector<char> filled(it->seq.begin(), it->seq.end());
            int fuzzy=1;
            bool loop=true;
            //std::string common;
            string fille=it->seq;
            //std::cout <<breakpointName<<std::endl;
            //std::cout <<filled.size()<<std::endl;
            //std::cout<<left.size()<<std::endl;
            //loop to identify alternative position and left normalized position
            while (loop)
            {
                
                if ((filled.size()> 1) & (left[seedk.size()-fuzzy] == filled[fille.size()-fuzzy]))
                {
                    fuzzy++;
                    continue;
                }
                if ((filled.size()==1) & (left[seedk.size()-fuzzy] == filled[fille.size()-1]))
                {
                    
                    fuzzy++;
                    continue;
                }
                else
                {
                    loop=false;
                }
                
                
            }
            
            //left normalization: if fuzzy>0 alt field is the concatenation of the char before insertion, the repeated sequence (of size=fuzzy) and the assembled sequence, the ref field is the char before insertion + the repeated sequence (of size=fuzzy)
            insertion=seedk.substr(seedk.size()-fuzzy,seedk.size()-1)+fille;
            string ref = seedk.substr(seedk.size()-fuzzy,seedk.size()-1);
            
            
            string token;
            istringstream iss(breakpointName);
            std::vector<string> tokens;
            // we split the header and put it in a vector tokens
            // split header breakpoint to extract information
            while(getline(iss,token,'_')) tokens.push_back(token);
            
            //cerr << tokens.size() << endl;
            
            string bkpt=breakpointName;
            string position = ".";
            string chromosome = ".";
            string GT = "./.";
            if (tokens.size()==7){ //MindTheGap find expected format
                bkpt = tokens[0].c_str();
                int pos = atoi(tokens[3].c_str())-ref.size()+1;
                position = to_string(pos);
                chromosome = tokens[1].c_str();
                string genotype = tokens[6].c_str();
                GT = genotype.compare("HOM")==0 ?  "1/1" : "0/1" ;
            }
            
            
            int size =insertion.size()-ref.size();

            // write in vcf format
            fprintf(_vcf_file,"%s\t%s\t%s\t%s\t%s\t.\tPASS\tTYPE=INS;LEN=%i;NPOS=%i;AVK=%.2f;MDK=%.2f\tGT\t%s\n",chromosome.c_str(),position.c_str(),bkpt.c_str(),ref.c_str(),insertion.c_str(),size,fuzzy,it->avg_coverage,it->median_coverage,GT.c_str());
            
        }
        
    }
    
    funlockfile(_vcf_file);
    
}

void Filler::writeToGFA(std::vector<filled_insertion_t>& filledSequences, string sourceSequence, string seedName, bool isRc, bool is_anchor_repeated){

    int overlap = getInput()->getInt(STR_CONTIG_OVERLAP);
    string seedDirection = "+";
    string targetDirection;
    string seedNameNode = seedName;
    string targetNameNode;

    if (isRc) // Seed contig has been read Rc
    {
        // Remove "_Rc" tag of contig name
        seedName = seedName.substr(0, seedName.size()-3);
        seedDirection = "-";
    }

    flockfile(_gfa_file);
    int nbInsertions = 0;
    int nbTotalInsertions = 0;
    //cout << "enter filled insertion t" << endl;
    for (std::vector<filled_insertion_t>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
    { // Duplicate of WriteFilledBreakpoint code
        string insertion = it->seq;
        int llen = insertion.length() ;
        if(llen > 0) nbTotalInsertions++;
    }

    // Write gapfilling as GFA node + 2 edges
    for (std::vector<filled_insertion_t>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
    {
        int qual = it->compute_qual(is_anchor_repeated);
        string insertion = it->seq;
        int llen = insertion.length() ;// - (int) R.length() - (int) L.length() - 2*hetmode;

        if(llen > 0)
        {

            std::ostringstream osolu_i;
            osolu_i <<   "solution " <<    nbInsertions+1 << "/" << nbTotalInsertions ;
            string solu_i = nbTotalInsertions >1 ?  osolu_i.str() : "" ;

            // Add seed and target kmers to insertion seq
            //insertion = sourceSequence + insertion + it->targetId
            bkpt_t targetId = it->targetId_anchor;
            string targetName = targetId.first;

            if(targetId.second) {
               targetDirection = "-";
               targetNameNode = targetName + "_Rc";
            } else {
               targetDirection = "+";
               targetNameNode = targetName;
            }

            // Write node
            int cov = it->median_coverage + 0.5;
            string nodeName = seedNameNode+";"+targetNameNode+";len_"+to_string(llen)+"_qual_"+to_string(qual)+"_median_cov_"+to_string(cov)+" "+solu_i; // Name could be computed once for gfa and fasta
            fprintf(_gfa_file,"S\t%s\t%s\n",nodeName.c_str(),insertion.c_str());

            // Write link between nodes

            // From seed to gapfilling
            fprintf(_gfa_file,"L\t%s\t%s\t%s\t+\t%iM\n",seedName.c_str(),seedDirection.c_str(),nodeName.c_str(),_kmerSize+overlap);

            // From gapfilling to seed
            fprintf(_gfa_file,"L\t%s\t+\t%s\t%s\t%iM\n",nodeName.c_str(),targetName.c_str(),targetDirection.c_str(),_kmerSize+overlap);

            nbInsertions++;
        }
    }
    funlockfile(_gfa_file);
}


set< info_node_t >  Filler::find_nodes_containing_multiple_R(bkpt_dict_t targetDictionary, string linear_seqs_name, int nb_mis_allowed, int nb_gaps_allowed)
{
    bool debug = false;
    //set<int> terminal_nodes;
    set< info_node_t >  terminal_nodes;
    BankFasta* Nodes = new BankFasta((char *)linear_seqs_name.c_str());

    long nodeNb = 0;
    char * nodeseq;
    size_t nodelen;

    // heuristics: R has to be seen entirely in the node up to nb_mis_allowed errors, in the forward strand
    BankFasta::Iterator  * itSeq  =  new BankFasta::Iterator  (*Nodes);

    // We loop over sequences.
    for (itSeq->first(); !itSeq->isDone(); itSeq->next())
    {
        int anchor_size = _kmerSize;

        nodelen = (*itSeq)->getDataSize();
        if (nodelen < _kmerSize )
        {
            cout << "Too short" << endl;
            nodeNb++;
            continue;
        }
        nodeseq =  (*itSeq)->getDataBuffer();

        int best_match=0;
        bkpt_t best_id;
        int position=0;
        bool arret=false;
        for (unsigned int j = 0; j < nodelen-_kmerSize+1 && !arret; j++)
        {
            //best_match = 0; //commented by CL 06/06/18 bug fix (this line prevented to find target kmers with errors in the middle of the contig)
            
            for (auto it=targetDictionary.begin(); it !=targetDictionary.end() && !arret; ++it)
                {
                    string ide = (it->second).first;
                    const char * anchor =(it->first).c_str();
                    int nbmatch=0;

                    for (int i = 0; i < anchor_size; i++)
                    {
                        nbmatch += identNT(nodeseq[j+i], anchor[i]);
                    }


                    if (nbmatch > best_match && nbmatch >=(anchor_size - nb_mis_allowed))
                    {
                        //cout << "nb match" << nbmatch << endl;
                        //cout << " target seq" << anchor << endl;
                        best_id = it->second;
                        position=j;
                        best_match = nbmatch;
                        if (nbmatch == anchor_size) {
                            arret = true;
                            break;
                        }
                    }
                }
        }
        if (best_match != 0)
        {
            //cout << endl;
            //cout << "cible" << best_id.first << " position " << position << "nodeId" << nodeNb << "  len " << nodelen << endl;
            //cout << " nodeseq" << nodeseq << endl;


            terminal_nodes.insert((info_node_t) {(int)nodeNb,(int)position, anchor_size - best_match, best_id}); // nodeNb,  j pos of beginning of right anchor

        }


        nodeNb++;
    }
    /*if (debug)
    {
        for(set< info_node_t >::iterator it = terminal_nodes.begin(); it != terminal_nodes.end() ; ++it)
            fprintf(stderr," (node %d pos %d) ",(*it).node_id,(*it).pos);
    }*/

    delete itSeq;
    delete Nodes;
    return terminal_nodes;
}
