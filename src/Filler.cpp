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

#define PRINT_DEBUG
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
Filler::Filler ()  : Tool ("MindTheGap find")
{
    
    getParser()->push_back (new OptionNoParam (STR_VERSION, "version", false));
    
    // from Graph.cpp
    getParser()->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory (in MBytes)", false, "2000"));
    getParser()->push_front (new OptionOneParam (STR_MAX_DISK, "max disk   (in MBytes)", false, "0"));
    
    //TODO HERE PUT THE FILL OPTIONS
    getParser()->push_front (new OptionOneParam (STR_MAX_DEPTH, "maximum length of insertions (nt)", false, "10000"));
    getParser()->push_front (new OptionOneParam (STR_MAX_NODES, "maximum number of nodes in contig graph (nt)", false, "100"));
    
    getParser()->push_front (new OptionOneParam (STR_SOLIDITY_KIND, "way to consider a solid kmer with several datasets (sum, min or max)", false, "sum"));
    getParser()->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "maximal abundance threshold for solid kmers", false, "4294967295"));
    getParser()->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "minimal abundance threshold for solid kmers", false, "3"));

    getParser()->push_front (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
    getParser()->push_front (new OptionOneParam (STR_URI_BKPT, "breakpoint file", true));
    getParser()->push_front (new OptionOneParam (STR_URI_GRAPH, "input graph file (likely a hdf5 file)",  false, ""));
    getParser()->push_front (new OptionOneParam (STR_URI_INPUT, "input read file(s)",  false, ""));
    
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
        getInput()->add(0,STR_BANK_CONVERT_TYPE,"tmp");
        getInput()->add(0,STR_URI_OUTPUT_DIR, ".");
        getInput()->add(0,STR_BLOOM_TYPE, "basic"); //neighbor basic cache
        getInput()->add(0,STR_DEBLOOM_TYPE, "original"); //cascading  pas bien car bcp plus de FP non critique au milieur trou
        getInput()->add(0,STR_DEBLOOM_IMPL, "basic"); //minimizer => STR_BLOOM_TYPE = neighbor
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_INTEGER_PRECISION, "0");
        getInput()->add(0,STR_MPHF_TYPE, "none");
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
        
        _graph = Graph::load (getInput()->getStr(STR_URI_GRAPH));
        _kmerSize = _graph.getKmerSize();
    }

    //Getting the breakpoint sequences
    _breakpointBank = new BankFasta(getInput()->getStr(STR_URI_BKPT));
    
    //Getting other parameters
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_depth = getInput()->getInt(STR_MAX_DEPTH);
    _max_nodes = getInput()->getInt(STR_MAX_NODES);
    
    // Now do the job
    // According to the kmer size, we call one fillBreakpoints method.
    if (_kmerSize < KSIZE_1)  { fillBreakpoints<KSIZE_1>  ();  }
    else if (_kmerSize < KSIZE_2)  { fillBreakpoints<KSIZE_2>  ();  }
    else if (_kmerSize < KSIZE_3)  { fillBreakpoints<KSIZE_3>  ();  }
    else if (_kmerSize < KSIZE_4)  { fillBreakpoints<KSIZE_4> ();  }
    else  { throw Exception ("unsupported kmer size %d", _kmerSize);  }

    //cout << "in MTG Fill" <<endl;
    // We gather some statistics.

    //getInfo()->add(1,"version",getVersion());
    getInfo()->add (1, &LibraryInfo::getInfo());
    resumeParameters();
    resumeResults();
}

void Filler::resumeParameters(){
    
    //Properties resumeParams;
    getInfo()->add(0,"Parameters");
    getInfo()->add(1,"Input data");
    if (getInput()->get(STR_URI_INPUT) != 0){
        getInfo()->add(2,"Reads",getInput()->getStr(STR_URI_INPUT).c_str());
    }
    if (getInput()->get(STR_URI_GRAPH) != 0){
        getInfo()->add(2,"Graph",getInput()->getStr(STR_URI_GRAPH).c_str());
    }
    getInfo()->add(2,"Breakpoints",getInput()->getStr(STR_URI_BKPT).c_str());
    getInfo()->add(1,"Graph");
    getInfo()->add(2,"kmer-size","%i", _kmerSize);
    try { // entour try/catch ici au cas ou le nom de la cle change dans gatb-core
            getInfo()->add(2,"abundance_min",_graph.getInfo().getStr("abundance_min").c_str());
            getInfo()->add(2,"abundance_max",_graph.getInfo().getStr("abundance_max").c_str());
            getInfo()->add(2,"solidity_kind",_graph.getInfo().getStr("solidity_kind").c_str());
            getInfo()->add(2,"nb_solid_kmers",_graph.getInfo().getStr("kmers_nb_solid").c_str());
            getInfo()->add(2,"nb_branching_nodes",_graph.getInfo().getStr("nb_branching").c_str());
        } catch (Exception e) {
            // doing nothing
        }

    getInfo()->add(1,"Breakpoint Filler options");
    getParser()->push_front (new OptionOneParam (STR_MAX_DEPTH, "maximum length of insertions (nt)", false, "10000"));
    getInfo()->add(2,"max_depth","%i", _max_depth);
    getInfo()->add(2,"max_nodes","%i", _max_nodes); //todo
    
}

void Filler::resumeResults(){
	getInfo()->add(0,"Results");
//	getInfo()->add(1,"Breakpoints");
//	getInfo()->add(2,"homozygous","%i", _nb_homo_clean+_nb_homo_fuzzy);
//	getInfo()->add(3,"clean","%i", _nb_homo_clean);
//	getInfo()->add(3,"fuzzy","%i", _nb_homo_fuzzy);

}

//template method : enabling to deal with all sizes of kmer <KSIZE_4
template<size_t span>
void Filler::fillBreakpoints(){

	// We create an iterator over the breakpoint bank.
	BankFasta::Iterator itSeq (*_breakpointBank);

	int nbBreakpoints=0;

	// We loop over sequences.
	for (itSeq.first(); !itSeq.isDone(); itSeq.next())
	{
		//iterate by pair of sequences (WARNING : no verification same breakpoint id)
		string sourceSequence =  string(itSeq->getDataBuffer());//previously L
		itSeq.next();
		if(itSeq.isDone()){
			throw Exception("Wrong breakpoint file: odd number of sequences...");
		}
		string targetSequence =  string(itSeq->getDataBuffer());//previously R

		//Initialize set of filled sequences
		set<string> filledSequences;

		// Resize to kmer-size :
		if(sourceSequence.size()>_kmerSize){
			sourceSequence.substr(sourceSequence.size()-_kmerSize,_kmerSize); //suffix of size _kmerSize
		}
		if(targetSequence.size()>_kmerSize){
			targetSequence.substr(0,_kmerSize); //prefix of size _kmerSize
		}

		gapFill<span>(sourceSequence,targetSequence,filledSequences);

		// TODO if set vide et reverse : reverse les 2 sequences + gapFill

		// TODO ecrire les resultats dans le fichier (method)
		writeFilledBreakpoint();

		// We increase the breakpoint counter.
		nbBreakpoints++;
	}

	cout << "nb breakpoints=" << nbBreakpoints <<endl;

}

//template method : enabling to deal with all sizes of kmer <KSIZE_4
template<size_t span>
void Filler::gapFill(string sourceSequence, string targetSequence, set<string>& filledSequences){

	//object used to mark the traversed nodes of the graph (note : it is reset at the beginning of construct_linear_seq)
	BranchingTerminator terminator (_graph);
	IterativeExtensions<span> extension (_graph, terminator, TRAVERSAL_CONTIG, ExtendStopMode_until_max_depth, SearchMode_Breadth, false, _max_depth, _max_nodes);
	//todo check param dontOutputFirstNucl=false ??
	// todo put these two above lines in fillBreakpoints and pass object extension in param

	//Build contigs and output them in a file in fasta format
	string contig_file_name = "contigs.fasta";
	extension.construct_linear_seqs(sourceSequence,targetSequence,contig_file_name,false); //last param : swf=stopWhenFound

    // connect the contigs into a graph
	string contig_graph_file_prefix="contig_graph";
	GraphOutputDot<span> graph_output(_kmerSize,contig_graph_file_prefix);
	graph_output.load_nodes_extremities(contig_file_name);
	graph_output.first_id_els = graph_output.construct_graph(contig_file_name,"LEFT");
	graph_output.close();
}

void Filler::writeFilledBreakpoint(){
//	fprintf(_breakpoint_file,">left_contig_%i_%s_pos_%lli_repeat_%i\n%s\n>right_contig_%i_%s_pos_%lli_repeat_%i\n%s\n",
//			bkt_id,
//			chrom_name.c_str(),
//			position,
//			repeat_size,
//			kmer_begin.c_str(),
//			bkt_id,
//			chrom_name.c_str(),
//			position,
//			repeat_size,
//			kmer_end.c_str()
//	);
}
