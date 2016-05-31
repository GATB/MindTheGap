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
Filler::Filler ()  : Tool ("MindTheGap fill") , _progress(0)
{
    
	//TODO rajouter les parametres
	_nb_mis_allowed = 0;
	_nb_gap_allowed = 0;
	_nb_breakpoints = 0;
	_nb_filled_breakpoints = 0;
	_nb_multiple_fill = 0;


    setParser (new OptionsParser ("MindTheGap fill"));

	IOptionsParser* generalParser = new OptionsParser("General");
	generalParser->push_front (new OptionNoParam (STR_HELP, "help", false));
	// generalParser->push_back (new OptionNoParam (STR_VERSION, "version", false)); // move this option in the main.cpp
	generalParser->push_front (new OptionOneParam (STR_VERBOSE,     "verbosity level",      false, "1"  ));
	generalParser->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory for graph building (in MBytes)", false, "2000"));
	generalParser->push_front (new OptionOneParam (STR_MAX_DISK, "max disk for graph building   (in MBytes)", false, "0"));
	generalParser->push_front (new OptionOneParam (STR_NB_CORES,    "number of cores",      false, "0"  ));

	IOptionsParser* inputParser = new OptionsParser("Input / output");
    inputParser->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
    inputParser->push_front (new OptionOneParam (STR_URI_BKPT, "breakpoint file", false, ""));
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
    
	if (getInput()->get(STR_HELP) != 0){
		cout << endl << "Usage:  MindTheGap fill (-in <reads.fq> | -graph <graph.h5>) -bkpt <breakpoints.fa> [options]" << endl;
		OptionsHelpVisitor v(cout);
		getParser()->accept(v);
		throw Exception(); // to get out with EXIT_FAILURE
	}

    if ((getInput()->get(STR_URI_GRAPH) != 0 && getInput()->get(STR_URI_INPUT) != 0) || (getInput()->get(STR_URI_GRAPH) == 0 && getInput()->get(STR_URI_INPUT) == 0))
    {

        throw OptionFailure(getParser(), "options -graph and -in are incompatible, but at least one of these is mandatory");

    }
    
    if (getInput()->get(STR_URI_BKPT) == 0){
    	throw OptionFailure(getParser(), "option -bkpt is mandatory");
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
        getInput()->add(0,STR_MPHF_TYPE, "none");
        getInput()->add(0,STR_BRANCHING_TYPE, "stored");
        getInput()->add(0,STR_MINIMIZER_SIZE, "8");
        getInput()->add(0,STR_REPARTITION_TYPE, "0");
        getInput()->add(0,STR_MINIMIZER_TYPE, "0");
        getInput()->add(0,STR_HISTOGRAM_MAX, "10000");
        getInput()->add(0,STR_KMER_ABUNDANCE_MIN_THRESHOLD,"3");
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

	
	_insert_file_name = getInput()->getStr(STR_URI_OUTPUT)+".insertions.fasta";
	_insert_file = fopen(_insert_file_name.c_str(), "w");
	if(_insert_file == NULL){
		string message = "Cannot open file "+ _insert_file_name + " for writting";
		throw Exception(message.c_str());
	}
	
	
    //Getting the breakpoint sequences
    _breakpointBank = new BankFasta(getInput()->getStr(STR_URI_BKPT));
    
    //Getting other parameters
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_depth = getInput()->getInt(STR_MAX_DEPTH);
    _max_nodes = getInput()->getInt(STR_MAX_NODES);
    
    // Now do the job
    time_t start_time = time(0);
    // According to the kmer size,  we call one fillBreakpoints method.
    Integer::apply<fillBreakpoints,Filler*> (_kmerSize, this);
    time_t end_time = time(0);
    double seconds=difftime(end_time,start_time);

    //cout << "in MTG Fill" <<endl;
    // We gather some statistics.

	fclose(_insert_file);
	
    //getInfo()->add(1,"version",getVersion());
	getInfo()->add(1,"version",_mtg_version);
	getInfo()->add(1,"gatb-core-library",STR_LIBRARY_VERSION);
	getInfo()->add(1,"supported_kmer_sizes","%s", KSIZE_STRING);
	
    //getInfo()->add (1, &LibraryInfo::getInfo());
    resumeParameters();
    resumeResults(seconds);
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
	getInfo()->add(1,"Output file","%s",_insert_file_name.c_str());

}

//template method : enabling to deal with all sizes of kmer <KSIZE_4
template<size_t span>
void Filler::fillBreakpoints<span>::operator ()  (Filler* object)
{
	//TODO count the number of filled insertions (in an attribute of class Filler), to output results in resumeResults()

	// We create an iterator over the breakpoint bank.
	BankFasta::Iterator itSeq (*object->_breakpointBank);

	int nbBreakpoints=0;

	
	
	u_int64_t nbBreakpointsEstimated = object->_breakpointBank->estimateNbItems()  / 2 ;  // 2 seq per breakpoint
	u_int64_t nbBreakpointsProgressDone = 0;
	
	object->setProgress (new ProgressSynchro (
									  object->createIteratorListener (nbBreakpointsEstimated, "Filling breakpoints"),
									  System::thread().newSynchronizer())
				 );
	object->_progress->init ();
	
	
	
	
	// We loop over sequences.
	for (itSeq.first(); !itSeq.isDone(); itSeq.next())
	{
		

		//iterate by pair of sequences (WARNING : no verification same breakpoint id)
		string sourceSequence =  string(itSeq->getDataBuffer(),itSeq->getDataSize());//previously L

		string breakpointName = string(itSeq->getComment());

		itSeq.next();
		if(itSeq.isDone()){
			throw Exception("Wrong breakpoint file: odd number of sequences...");
		}

		string targetSequence =  string(itSeq->getDataBuffer(),itSeq->getDataSize());//previously R
		

		//printf("break %i L : %s  R: %s \n",nbBreakpoints,sourceSequence.c_str(),targetSequence.c_str());
		
		//Initialize set of filled sequences
		set<string> filledSequences;

		// Resize to kmer-size :
		if(sourceSequence.size()>object->_kmerSize){
			sourceSequence.substr(sourceSequence.size()-object->_kmerSize,object->_kmerSize); //suffix of size _kmerSize
		}

		if(targetSequence.size()>object->_kmerSize){
			targetSequence.substr(0,object->_kmerSize); //prefix of size _kmerSize
		}

		object->gapFill<span>(sourceSequence,targetSequence,filledSequences);

		//Can be modified : could do in reverse mode even if filledSequences is not empty (new filled sequences are inserted into the set : to verify)
		if(filledSequences.size()==0){
			string sourceSequence2 = revcomp_sequence(targetSequence);
			string targetSequence2 = revcomp_sequence(sourceSequence);
			object->gapFill<span>(sourceSequence2,targetSequence2,filledSequences,true);
		}
		
		//Checks if all sequences are roughly the same :
		if (all_consensuses_almost_identical(filledSequences,90))
		{
			//if(verb)     printf(" [SUCCESS]\n");
			if (filledSequences.size() > 1) {
				//				stringstream ss;
				//				ss << "cons" <<filledSequences.size();
				//				breakpointName=breakpointName+ss.str();
				filledSequences.erase(++(filledSequences.begin()),filledSequences.end()); // keep only one consensus sequence
			}
		}
		else
			;
		//if(verb)   printf(" [MULTIPLE SOLUTIONS]\n");
		

		// TODO ecrire les resultats dans le fichier (method) : attention checker si mode Une ou Multiple Solutions
		object->writeFilledBreakpoint(filledSequences,breakpointName);

		// We increase the breakpoint counter.
		nbBreakpoints++;

		//progress bar
		nbBreakpointsProgressDone++;
		if (nbBreakpointsProgressDone > 50)   {  object->_progress->inc (nbBreakpointsProgressDone);  nbBreakpointsProgressDone = 0;  }
	}

	object->_progress->finish ();

	object->_nb_breakpoints = nbBreakpoints;

	cout << "nb breakpoints=" << nbBreakpoints <<endl;
}

//template method : enabling to deal with all sizes of kmer <KSIZE_4
template<size_t span>
void Filler::gapFill(string sourceSequence, string targetSequence, set<string>& filledSequences, bool reversed){


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

	//find targetSequence in nodes of the contig graph
	set< std::pair<int,int> > terminal_nodes_with_endpos = find_nodes_containing_R(targetSequence, contig_file_name, _nb_mis_allowed, _nb_gap_allowed);
	

	//printf("nb contig with target %zu \n",terminal_nodes_with_endpos.size());
	
	//also convert it to list of node id for traditional use
	set<int> terminal_nodes;
	for (set< std::pair<int,int> >::iterator it = terminal_nodes_with_endpos.begin(); it != terminal_nodes_with_endpos.end(); it++)
	{
		terminal_nodes.insert((*it).first);
	}


	// analyze the graph to find a satisfying gap sequence between L and R

	GraphAnalysis graph = GraphAnalysis(graph_output.get_dot_file_name(),_kmerSize);
	graph.debug = false;

	if(terminal_nodes.size()==0)
	{
		//if(verb)
		//	printf("Right anchor not found.. gapfillling failed... \n");
		return ;
	}

	//if(verb)    fprintf(stderr," analysis..");
	bool success;
	set<unlabeled_path> paths = graph.find_all_paths(terminal_nodes, success);

	//now this func also cuts the last node just before the beginning of the right anchor
	set<string> tmpSequences = graph.paths_to_sequences(paths,terminal_nodes_with_endpos);
	
	if(reversed)
	{
		set<string>::iterator its;
		for (its = tmpSequences.begin(); its != tmpSequences.end(); ++its)
		{
			filledSequences.insert (revcomp_sequence(*its));
		}
		return;
	}
		
	

	filledSequences.insert(tmpSequences.begin(),tmpSequences.end());


	remove(contig_file_name.c_str());
	remove((contig_graph_file_prefix+".graph").c_str());

}

void Filler::writeFilledBreakpoint(set<string>& filledSequences, string breakpointName){
	
	//printf("-- writeFilledBreakpoint --\n");
	
	//printf("found %zu seq \n",filledSequences.size());
	
	int nbInsertions = 0;
	int nbTotalInsertions = 0;

	for (set<string>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
	{
		string insertion = *it;
		int llen = insertion.length() ;
		if(llen > 0) nbTotalInsertions++;
	}
	
	
	for (set<string>::iterator it = filledSequences.begin(); it != filledSequences.end() ; ++it)
	{
		string insertion = *it;
		int llen = insertion.length() ;// - (int) R.length() - (int) L.length() - 2*hetmode;
		
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
			
			
			int bkptid;
			//parse bkpt header
			//get bkpt id
			sscanf(breakpointName.c_str(),"bkpt%i*",&bkptid );
			const char * end_header = strstr(breakpointName.c_str(), "kmer_");

			fprintf(_insert_file,">bkpt%i insertion_len_%d_%s  %s\n",bkptid,llen,end_header+5,solu_i.c_str());

			//fprintf(_insert_file,"> insertion ( len= %d ) for breakpoint \"%s\"  %s  \n",llen, breakpointName.c_str(),solu_i.c_str());
			//todo check  revcomp here
			fprintf(_insert_file,"%.*s\n",(int)llen,insertion.c_str() );
		}
		nbInsertions++;
	}
	
	if(nbInsertions>0){
		_nb_filled_breakpoints++;
		if(nbInsertions>1){
			_nb_multiple_fill++;
		}
	}


}

set< std::pair<int,int> >  Filler::find_nodes_containing_R(string targetSequence, string linear_seqs_name, int nb_mis_allowed, int nb_gaps_allowed)
{
    bool debug = false;
    //set<int> terminal_nodes;
    set< std::pair<int,int> >  terminal_nodes;
    BankFasta* Nodes = new BankFasta((char *)linear_seqs_name.c_str());

    long nodeNb = 0;
    char * nodeseq;
    size_t nodelen;
    int nbmatch =0;

    const char * anchor = targetSequence.c_str();
    int anchor_size=  targetSequence.size();

    // heuristics: R has to be seen entirely in the node up to nb_mis_allowed errors, in the forward strand

    BankFasta::Iterator  * itSeq  =  new BankFasta::Iterator  (*Nodes);

    // We loop over sequences.
    for (itSeq->first(); !itSeq->isDone(); itSeq->next())
    {
    	nodelen = (*itSeq)->getDataSize();
        if (nodelen < targetSequence.size())
        {
            nodeNb++;
            continue;
        }

		nodeseq =  (*itSeq)->getDataBuffer();

		
        if (debug)
            printf("searching %s (size=%zu nt) in %s (size=%zu nt)\n",anchor,targetSequence.size(),nodeseq,nodelen);

		int best_err = 100000000;
		int curr_err = 0;
		int   best_j = -1;
		// int bg, bm; debug only
        for (unsigned int j = 0; j < nodelen-targetSequence.size()+1; j++)
        {
            nbmatch=0;

			int nw_mis   = 0 ;
			int nw_gaps  = 0 ;
			int nw_match = 0 ;

			if(nb_gaps_allowed>0)
			{


			    std::string nodestring(nodeseq + j , min(anchor_size + nb_gaps_allowed, (int)(nodelen-j)));
				int deb = j==128;
				needleman_wunsch(targetSequence,nodestring,&nw_match,&nw_mis, &nw_gaps);
				curr_err = nw_mis + nw_gaps;

//				printf("-----------   ---------\n");
//				printf("testing pos %i:   mis %i gap %i  curr_err %i\n",j,nw_mis,nw_gaps,curr_err);
//				printf("%s  \n",R.c_str());
//				printf("%s  \n",nodestring.c_str());
//

				if( (nw_mis <= nb_mis_allowed && nw_gaps <= nb_gaps_allowed )  && (curr_err < best_err))
				{
					best_err = curr_err;
					best_j  = j;
					// bg = nw_gaps; bm = nw_mis;
				//	printf("found correct pos %i nb mis %i nb gaps %i \n",j,nw_mis,nw_gaps);

				}
			}
			else
			{
				for (int i = 0; i < anchor_size; i++)
				{
					nbmatch += identNT(nodeseq[j+i], anchor[i]);
				}
				if(nbmatch >= (anchor_size - nb_mis_allowed))
				{
					terminal_nodes.insert( std::make_pair (nodeNb,j) ); // nodeNb,  j pos of beginning of right anchor
				//	printf("found pos %i nbmatch %i  \n",j,nbmatch);

					break;
				}

			}
        }

		if(nb_gaps_allowed>0 && best_j != -1)
		{
			terminal_nodes.insert( std::make_pair (nodeNb,best_j) ); // nodeNb,  j pos of beginning of right anchor
		//	printf("found pos %i nb mis %i nb gaps %i \n",best_j,bm,bg);
			break;
		}


        nodeNb++;
    }


    if (debug)
    {
        for(set< std::pair<int,int> >::iterator it = terminal_nodes.begin(); it != terminal_nodes.end() ; ++it)
            fprintf(stderr," (node %d pos %d) ",(*it).first,(*it).second);
    }

	delete itSeq;
	delete Nodes;

    return terminal_nodes;
}
