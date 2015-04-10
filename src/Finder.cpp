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
    finderParser->push_front (new OptionOneParam (STR_HET_MAX_OCC, "maximal number of occurrences of a kmer in the reference genome allowed for heterozyguous breakpoints", false,"1"));
    //allow to find heterozyguous breakpoints in n-repeated regions of the reference genome

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

    // Getting the reference genome
    _refBank = new BankFasta(getInput()->getStr(STR_URI_REF));

    
    //Getting other parameters
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_repeat = getInput()->getInt(STR_MAX_REPEAT);
    _homo_only=getInput()->get(STR_HOMO_ONLY) !=0;
    _het_max_occ=getInput()->getInt(STR_HET_MAX_OCC);

    if(!_homo_only){
    	// Building the index of reference (k-1)-mers occurring more than _het_max_occ + 1 times
    	string tempFileName="trashme";
    	stringstream commandLine;
    	commandLine << STR_URI_INPUT << " " << getInput()->getStr(STR_URI_REF) << " " <<  //start from fasta file (and not from Bank : can not be used several times)
    			STR_KMER_ABUNDANCE_MIN << " " << _het_max_occ + 1 << " " <<
    			STR_KMER_SIZE << " " << _kmerSize-1 << " " <<
    			STR_DEBLOOM_TYPE << " none " <<
    			STR_URI_OUTPUT << " " << tempFileName << " ";

    	cout << commandLine.str() << endl;
    	_ref_graph = Graph::create(commandLine.str().c_str());
    	cout << "ok" <<endl;
    	System::file().remove(tempFileName+".h5");
    }

    // Now do the job

    // According to the kmer size,  we call one fillBreakpoints method.
    if (_kmerSize < KSIZE_1)  { findBreakpoints<KSIZE_1>  ();  }
    else if (_kmerSize < KSIZE_2)  { findBreakpoints<KSIZE_2>  ();  }
    else if (_kmerSize < KSIZE_3)  { findBreakpoints<KSIZE_3>  ();  }
    else if (_kmerSize < KSIZE_4)  { findBreakpoints<KSIZE_4> ();  }
    else  { throw Exception ("unsupported kmer size %d", _kmerSize);  }

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

//template method : enabling to deal with all sizes of kmer <KSIZE_4
template<size_t span>
void Finder::findBreakpoints(){

	uint64_t bkt_id=0; // id of detected breakpoint

	int nbKmers=0;
	int nbSequences=0; // not used

	uint64_t nb_ref_solid = 0; // not used
	uint64_t nb_ref_notsolid = 0; //not used
	uint64_t solid_stretch_size = 0; //size of current stretch of 1 (ie kmer indexed)
	uint64_t gap_stretch_size = 0; //size of current stretch of 0 (ie kmer not indexed)
	uint64_t previous_gap_stretch_size = 0;

	typedef typename gatb::core::kmer::impl::Kmer<span>::ModelCanonical KmerModel; // ModelCanonical ModelDirect
	typedef typename gatb::core::kmer::impl::Kmer<span>::ModelCanonical::Iterator KmerIterator;
	typedef typename gatb::core::kmer::impl::Kmer<span>::Type        kmer_type;
	typedef typename gatb::core::kmer::impl::Kmer<span>::Count       kmer_count; // not used

#ifdef PRINT_DEBUG
	string deb01;
#endif

	kmer_type kmer_begin;
	kmer_type kmer_end;
	kmer_type previous_kmer;

	// We declare a kmer model with a given span size.
	KmerModel model (_kmerSize);
	//std::cout << "span: " << model.getSpan() << std::endl;
	// We create an iterator over this bank.
	BankFasta::Iterator itSeq (*_refBank);
	// We declare an iterator on a given sequence.
	KmerIterator itKmer (model);
	// We loop over sequences.
	for (itSeq.first(); !itSeq.isDone(); itSeq.next())
	{
		solid_stretch_size = 0;
		gap_stretch_size = 0;

		//Method : an homozyguous breakpoint is detected as a gap_stretch (ie. consecutive kmers on the sequence, that are not indexed in the graph) of particular sizes.
		//BUT there are some False Positives when we query the graph : when we ask the graph if a kmer is indexed (when starting from a previous not indexed kmer) it may wrongly answer yes
		//(the gatb dbg is exact only if the kmer we query is the neighbor of a truly solid kmer)
		//FP are likely to be isolated, ie. surrounded by not indexed kmers, therefore they can be detected as solid_stretches of size 1.
		
#ifdef PRINT_DEBUG
		deb01.clear();
#endif

		
		// We set the data from which we want to extract kmers, to the kmer iterator
		itKmer.setData (itSeq->getData());
		char* chrom_sequence = itSeq->getDataBuffer();
		string chrom_name = itSeq->getComment();
		uint64_t position=0;

		// We iterate the kmers.
		for (itKmer.first(); !itKmer.isDone(); itKmer.next(), position++)
		{
			nbKmers++;

			//we need to convert the kmer in a node to query the graph.
			Node node(Node::Value(itKmer->value()));

			if (_graph.contains(node)) //the kmer is indexed
			{
#ifdef PRINT_DEBUG
				deb01+= "1";
#endif
				nb_ref_solid++;
				solid_stretch_size++;

				if(solid_stretch_size > 1){ //to be sure the gap is finished, it must be followed by a solid stretch of size >1, otherwise it could be a false positive of the graph (likely isolated)
					if(gap_stretch_size == (_kmerSize-1)){
						// clean insert site

//						cout << "gap size k-1" << endl;
//						cout << "position " << position -1 << endl;
//						cout << "kmer begin " << model.toString (kmer_begin) << endl;
//						cout << "kmer end " << model.toString (kmer_end) << endl;

						string kmer_begin_str = model.toString (kmer_begin);
						string kmer_end_str = model.toString (kmer_end);
						writeBreakpoint(bkt_id,chrom_name,position-1,kmer_begin_str, kmer_end_str,0);
						bkt_id++;
						_nb_homo_clean++;
					}
					else if(gap_stretch_size < _kmerSize - 1 && gap_stretch_size >= _kmerSize -1 -_max_repeat){
						// Fuzzy site, position and kmer_end are impacted by the repeat

						int repeat_size = _kmerSize - 1 - gap_stretch_size;
						string kmer_begin_str = model.toString (kmer_begin);
						string kmer_end_str = string(&chrom_sequence[position-1+repeat_size], _kmerSize);
						writeBreakpoint(bkt_id,chrom_name,position -1 + repeat_size,kmer_begin_str,kmer_end_str, repeat_size);
						bkt_id++;
						_nb_homo_fuzzy++;
					}
					else if(gap_stretch_size>0) {
						//for debug
						//cout << "gap_stretch_size = " << gap_stretch_size << " in sequence " << chrom_name << " position " << position -1 << endl;
					}

					gap_stretch_size = 0;// gap stretch size is re-set to 0 only when we are sure that the end of the gap is not due to an isolated solid kmer (likely FP)
				}
				if (solid_stretch_size==1) kmer_end = itKmer->forward(); // kmer_end should be the first kmer indexed after a gap (the first kmer of a solid_stretch is when solid_stretch_size=1)
			}
			else //kmer is not indexed, measure the size of the zone not covered by kmers of the reads (= a gap)
			{
#ifdef PRINT_DEBUG
				deb01+= "0";
#endif
				nb_ref_notsolid++;
				if(solid_stretch_size==1) // it means that the previous position was an isolated solid kmer
				{
					gap_stretch_size = gap_stretch_size + solid_stretch_size ; //if previous position was an isolated solid kmer, we need to add 1 to the gap_stretch_size (as if replacing the FP by a non indexed kmer)
				}
				if(solid_stretch_size > 1) // previous position was a solid kmer, but not an isolated one  => we are at the beginning of a gap
				{
					kmer_begin = previous_kmer ;
				}
				gap_stretch_size ++;
				solid_stretch_size =0;

			}
			previous_kmer = itKmer->forward();

		}
		
#ifdef PRINT_DEBUG
		cout << deb01 << endl;
#endif
		
		// We increase the sequence counter.
		nbSequences++;
	}

//	cout << "nb sequences=" << nbSequences <<endl;
//	cout << "nb kmers=" << nbKmers <<endl;

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
