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

#include <Finder.hpp>
#include <FindBreakpoints.hpp>
#include <FindObserver.hpp>
#include <IFindObserver.hpp>

//#define PRINT_DEBUG
/********************************************************************************/

// We define some constant strings for names of command line parameters
static const char* STR_FOO = "-foo";

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
    
    _nb_homo_clean=0;
    _nb_homo_fuzzy=0;
    _nb_hetero_clean=0;
    _nb_hetero_fuzzy=0;
    

    // Option parser, with several sub-parsers
    setParser (new OptionsParser ("MindTheGap find"));

	IOptionsParser* generalParser = new OptionsParser("General");
	generalParser->push_front (new OptionNoParam (STR_HELP, "help", false));
	// generalParser->push_back (new OptionNoParam (STR_VERSION, "version", false)); // move this option in the main.cpp
	generalParser->push_front (new OptionOneParam (STR_VERBOSE,     "verbosity level",      false, "1"  ));
	generalParser->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory for graph building (in MBytes)", false, "2000"));
	generalParser->push_front (new OptionOneParam (STR_MAX_DISK, "max disk for graph building (in MBytes)", false, "0"));
	generalParser->push_front (new OptionOneParam (STR_NB_CORES,    "number of cores",      false, "0"  ));

	IOptionsParser* inputParser = new OptionsParser("Input / output");
    inputParser->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_REF, "reference genome file", false,""));
    inputParser->push_front (new OptionOneParam (STR_URI_GRAPH, "input graph file (likely a hdf5 file)",  false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_INPUT, "input read file(s)",  false, ""));

	IOptionsParser* finderParser = new OptionsParser("Detection");
    finderParser->push_front (new OptionOneParam (STR_MAX_REPEAT, "maximal repeat size detected for fuzzy site", false, "5"));
    finderParser->push_front (new OptionNoParam (STR_HOMO_ONLY, "only search for homozygous breakpoints", false));

    IOptionsParser* graphParser = new OptionsParser("Graph building");
    graphParser->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "minimal abundance threshold for solid kmers", false, "3"));
  	graphParser->push_front (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));

	getParser()->push_front(generalParser);
	getParser()->push_front(finderParser);
	getParser()->push_front(graphParser);
	getParser()->push_front(inputParser);
    
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
    
	if (getInput()->get(STR_HELP) != 0){
		cout << endl << "Usage:  MindTheGap find (-in <reads.fq> | -graph <graph.h5>) -ref <reference.fa> [options]" << endl;
		OptionsHelpVisitor v(cout);
		getParser()->accept(v);
		throw Exception(); // to get out with EXIT_FAILURE
	}

	// Checks mandatory options
    if ((getInput()->get(STR_URI_GRAPH) != 0 && getInput()->get(STR_URI_INPUT) != 0) || (getInput()->get(STR_URI_GRAPH) == 0 && getInput()->get(STR_URI_INPUT) == 0))
    {
        throw OptionFailure(getParser(), "options -graph and -in are incompatible, but at least one of these is mandatory");
    }
    
    if (getInput()->get(STR_URI_REF) == 0){
    	throw OptionFailure(getParser(), "option -ref is mandatory");
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
    	getInput()->add(0,STR_KMER_ABUNDANCE_MAX, "4294967295"); //maximal abundance threshold for solid kmers

        getInput()->add(0,STR_BANK_CONVERT_TYPE,"tmp");
        getInput()->add(0,STR_URI_OUTPUT_DIR, ".");
        getInput()->add(0,STR_BLOOM_TYPE, "basic"); //neighbor basic cache
        getInput()->add(0,STR_DEBLOOM_TYPE, "original"); //DO NOT use cascading : generates too many FP inside  pas bien car bcp plus de FP non critique au milieur trou
        getInput()->add(0,STR_DEBLOOM_IMPL, "basic"); //minimizer => STR_BLOOM_TYPE = neighbor
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_INTEGER_PRECISION, "0");
        getInput()->add(0,STR_MPHF_TYPE, "none");
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


    // Preparing the output file
    _breakpoint_file_name = getInput()->getStr(STR_URI_OUTPUT)+".breakpoints";
    _breakpoint_file = fopen(_breakpoint_file_name.c_str(), "w");
    if(_breakpoint_file == NULL){
        //cerr <<" Cannot open file "<< _output_file <<" for writting" << endl;
        string message = "Cannot open file "+ _breakpoint_file_name + " for writting";
        throw Exception(message.c_str());

    }

    
    //Getting other parameters
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_repeat = getInput()->getInt(STR_MAX_REPEAT);
    _homo_only=getInput()->get(STR_HOMO_ONLY) !=0;

    // Getting the reference genome
    _refBank = new BankFasta(getInput()->getStr(STR_URI_REF));
    
    // Now do the job

    // According to the kmer size,  we call one fillBreakpoints method.
    if (_kmerSize < KSIZE_1) { runFindBreakpoints<KSIZE_1>(); }
    else if (_kmerSize < KSIZE_2) { runFindBreakpoints<KSIZE_2>(); }
    else if (_kmerSize < KSIZE_3) { runFindBreakpoints<KSIZE_3>(); }
    else if (_kmerSize < KSIZE_4) { runFindBreakpoints<KSIZE_4>(); }
    else  { throw Exception ("unsupported kmer size %d", _kmerSize);  }


    //cout << "in MTG" <<endl;
    // We gather some statistics.
    fclose(_breakpoint_file);

    // Printing result informations (ie. add info to getInfo(), in Tool Info is printed automatically after end of execute() method
    //getInfo()->add(1,"version",getVersion());
    getInfo()->add (1, &LibraryInfo::getInfo());
    resumeParameters();
    resumeResults();
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
    try { // entour try/catch ici au cas ou le nom de la cle change dans gatb-core
            getInfo()->add(2,"abundance_min",_graph.getInfo().getStr("abundance_min").c_str());
            getInfo()->add(2,"nb_solid_kmers",_graph.getInfo().getStr("kmers_nb_solid").c_str());
            getInfo()->add(2,"nb_branching_nodes",_graph.getInfo().getStr("nb_branching").c_str());
        } catch (Exception e) {
            // doing nothing
        }

    getInfo()->add(1,"Breakpoint detection options");
    getInfo()->add(2,"max_repeat","%i", _max_repeat);
    getInfo()->add(2,"homo_only","%i", _homo_only); //todo
    
}

void Finder::resumeResults(){
	getInfo()->add(0,"Results");
	getInfo()->add(1,"Breakpoints");
	getInfo()->add(2,"homozygous","%i", _nb_homo_clean+_nb_homo_fuzzy);
	getInfo()->add(3,"clean","%i", _nb_homo_clean);
	getInfo()->add(3,"fuzzy","%i", _nb_homo_fuzzy);

	getInfo()->add(1,"output file","%s",_breakpoint_file_name.c_str());

}

void Finder::writeBreakpoint(int bkt_id, string& chrom_name, uint64_t position, string& kmer_begin, string& kmer_end, int repeat_size){
	fprintf(_breakpoint_file,">left_contig_%i_%s_pos_%lli_repeat_%i\n%s\n>right_contig_%i_%s_pos_%lli_repeat_%i\n%s\n",
			bkt_id,
			chrom_name.c_str(),
			position,
			repeat_size,
			kmer_begin.c_str(),
			bkt_id,
			chrom_name.c_str(),
			position,
			repeat_size,
			kmer_end.c_str()
	);
}

template<size_t span>
void Finder::runFindBreakpoints()
{
    FindBreakpoints<span> findBreakpoints(this);

    /* Add observer */
    findBreakpoints.addObserver(new FindCleanInsert<span>(&findBreakpoints));
    findBreakpoints.addObserver(new FindFuzzyInsert<span>(&findBreakpoints));

    // Add this Observer after all other 
    findBreakpoints.addObserver(new FindEndSolid<span>(&findBreakpoints));
    findBreakpoints.addObserver(new FindEndGap<span>(&findBreakpoints));

    /* Run */
    findBreakpoints();
}
