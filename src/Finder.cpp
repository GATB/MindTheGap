/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G. Rizk
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


    /* Commented this until there is a function hide() for parser
     because we want some of the options of Graph to be masked to the user
     plus it enables to master the order of appearance of options in the Help message
     //Getting the options of graph (ie. dbgh5)
     //OptionsParser parser = Graph::getOptionsParser(false);
     //getParser()->add (parser);
     */
    
    // Here getParser() has already the options inherited from Tool : nb-cores, verbose et help
    // push_front options, in the order of more and more importance
    getParser()->push_back (new OptionNoParam (STR_VERSION, "version", false));
    
    // from Graph.cpp
    getParser()->push_front (new OptionOneParam (STR_MAX_MEMORY, "max memory (in MBytes)", false, "2000"));
    getParser()->push_front (new OptionOneParam (STR_MAX_DISK, "max disk   (in MBytes)", false, "0"));
    
    getParser()->push_front (new OptionOneParam (STR_MAX_REPEAT, "maximal repeat size detected for fuzzy site", false, "5"));
    getParser()->push_front (new OptionNoParam (STR_HOMO_ONLY, "only search for homozygous breakpoints", false));
    
    getParser()->push_front (new OptionOneParam (STR_SOLIDITY_KIND, "way to consider a solid kmer with several datasets (sum, min or max)", false, "sum"));
    getParser()->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MAX, "maximal abundance threshold for solid kmers", false, "4294967295"));
    getParser()->push_front (new OptionOneParam (STR_KMER_ABUNDANCE_MIN, "minimal abundance threshold for solid kmers", false, "3"));

    getParser()->push_front (new OptionOneParam (STR_KMER_SIZE, "size of a kmer", false, "31"));
    getParser()->push_front (new OptionOneParam (STR_URI_OUTPUT, "prefix for output files", false, ""));
    getParser()->push_front (new OptionOneParam (STR_URI_REF, "reference genome file", true));
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
void Finder::execute ()
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
        getInput()->add(0,STR_BLOOM_TYPE, "neighbor"); //neighbor
        getInput()->add(0,STR_DEBLOOM_TYPE, "cascading");
        getInput()->add(0,STR_DEBLOOM_IMPL, "minimizer"); //minimizer => STR_BLOOM_TYPE = neighbor
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

    string breakpoint_file_name = getInput()->getStr(STR_URI_OUTPUT)+".breakpoints";
    _breakpoint_file = fopen(breakpoint_file_name.c_str(), "w");
    if(_breakpoint_file == NULL){
        //cerr <<" Cannot open file "<< _output_file <<" for writting" << endl;
        string message = "Cannot open file "+ breakpoint_file_name + " for writting";
        throw Exception(message.c_str());

    }

    
    //Getting other parameters
    _nbCores = getInput()->getInt(STR_NB_CORES);
    _max_repeat = getInput()->getInt(STR_MAX_REPEAT);
    _homo_only=getInput()->get(STR_HOMO_ONLY) !=0;
    cout << _homo_only <<endl;

    // Getting the reference genome
    _refBank = new BankFasta(getInput()->getStr(STR_URI_REF));
    
    // Now do the job

    /** According to the kmer size, we instantiate one DSKAlgorithm class and delegate the actual job to it. */
    if (_kmerSize < KSIZE_1)  { findBreakpoints<KSIZE_1>  ();  }
    else if (_kmerSize < KSIZE_2)  { findBreakpoints<KSIZE_2>  ();  }
    else if (_kmerSize < KSIZE_3)  { findBreakpoints<KSIZE_3>  ();  }
    else if (_kmerSize < KSIZE_4)  { findBreakpoints<KSIZE_4> ();  }
    else  { throw Exception ("unsupported kmer size %d", _kmerSize);  }


    //cout << "in MTG" <<endl;
    // We gather some statistics.

    fclose(_breakpoint_file);
    //getInfo()->add(1,"version",getVersion());
    getInfo()->add (1, &LibraryInfo::getInfo());
    resumeParameters();

}

void Finder::resumeParameters(){
    
    //Properties resumeParams;
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
            getInfo()->add(2,"abundance_max",_graph.getInfo().getStr("abundance_max").c_str());
            getInfo()->add(2,"solidity_kind",_graph.getInfo().getStr("solidity_kind").c_str());
            getInfo()->add(2,"nb_solid_kmers",_graph.getInfo().getStr("kmers_nb_solid").c_str());
            getInfo()->add(2,"nb_branching_nodes",_graph.getInfo().getStr("nb_branching").c_str());
        } catch (Exception e) {
            // doing nothing
        }

    getInfo()->add(1,"Breakpoint Finder options");
    getInfo()->add(2,"max_repeat","%i", _max_repeat);
    getInfo()->add(2,"homo_only","%i", "T"); //todo
    
}

template<size_t span>
void Finder::findBreakpoints(){

	uint64_t bkt_id=0;
	int nbKmers=0;
	int nbSequences=0;
	//TODO gerer le span automatique (cf. template)

	uint64_t nb_ref_solid = 0;
	uint64_t nb_ref_notsolid = 0;
	uint64_t solid_stretch_size = 0; //size of current stretch of 1 (ie kmer indexed)
	uint64_t gap_stretch_size = 0; //size of current stretch of 0 (ie kmer not indexed)


	typedef typename gatb::core::kmer::impl::Kmer<span>::ModelDirect KmerModel;
	typedef typename gatb::core::kmer::impl::Kmer<span>::ModelDirect::Iterator KmerIterator;
	typedef typename gatb::core::kmer::impl::Kmer<span>::Type        kmer_type;
	typedef typename gatb::core::kmer::impl::Kmer<span>::Count       kmer_count;

	kmer_type kmer_begin;
	kmer_type kmer_end;
	kmer_type previous_kmer;

	// We declare a kmer model with a given span size.
	KmerModel model (_kmerSize);
	std::cout << "span: " << model.getSpan() << std::endl;
	// We create an iterator over this bank.
	BankFasta::Iterator itSeq (*_refBank);
	// We declare an iterator on a given sequence.
	KmerIterator itKmer (model);
	// We loop over sequences.
	for (itSeq.first(); !itSeq.isDone(); itSeq.next())
	{
		solid_stretch_size = 0;
		gap_stretch_size = 0;

		// We set the data from which we want to extract kmers.
		itKmer.setData (itSeq->getData());
		char* chrom_sequence = itSeq->getDataBuffer();
		string chrom_name = itSeq->getComment();
		uint64_t position=0;
		// We iterate the kmers.
		for (itKmer.first(); !itKmer.isDone(); itKmer.next(), position++)
		{
			nbKmers++;
			Node node(Node::Value(itKmer->value()));
			bool ok=_graph.contains(node);
			if( !ok){
				cout << "kmer absent " << model.toString (itKmer->value()) << endl;
			}

			if (_graph.contains(node)) //kmer is indexed
			{
				nb_ref_solid++;
				solid_stretch_size++;

				if(solid_stretch_size > 1){
					if(gap_stretch_size == (_kmerSize-1)){
						// clean insert site
						cout << "trou size k-1" << endl;
						cout << "position " << position -1 << endl;
						cout << "kmer begin " << model.toString (kmer_begin) << endl;
						cout << "kmer end " << model.toString (kmer_end) << endl;
						//TODO : ecrire breakpoint dans fichier
						string kmer_begin_str = model.toString (kmer_begin);
						string kmer_end_str = model.toString (kmer_end);
						writeBreakpoint(bkt_id,chrom_name,position-1,kmer_begin_str, kmer_end_str);
						bkt_id++;
						_nb_homo_clean++;
					}
					else if(gap_stretch_size < _kmerSize - 1 && gap_stretch_size >= _kmerSize -1 -_max_repeat){
						// Fuzzy site, position and kmer_end are impacted by the repeat
						int repeat_size = _kmerSize - 1 - gap_stretch_size;
						string kmer_begin_str = model.toString (kmer_begin);
						string kmer_end_str = string(chrom_sequence[position-1+repeat_size], _kmerSize);
						writeBreakpoint(bkt_id,chrom_name,position -1 + repeat_size,kmer_begin_str,kmer_end_str);
						bkt_id++;
						_nb_homo_fuzzy++;
					}
				}
				if (solid_stretch_size > 1) gap_stretch_size = 0; // du coup on sort le trou a tai indexed ==2, gap_stretch_size pas remis a 0 par solide isole (FP)
				if (solid_stretch_size==1) kmer_end = itKmer->value(); // kmer_end should be first kmer indexed after a hole

			}
			else //kmer is not indexed, measure size of the zone not covered by kmers of the reads
			{
				nb_ref_notsolid++;
				// if(gap_stretch_size == 0 && solid_stretch_size==1) //si zone indexe prec est ==1, probable FP, merge size, keep old kmer_begin
				if(solid_stretch_size > 1) // begin of not indexed zone
				{
					kmer_begin = previous_kmer ;
				}
				gap_stretch_size ++;
				solid_stretch_size =0;

			}
			previous_kmer = itKmer->value();

		}
		// We increase the sequences counter.
		nbSequences++;
	}

	cout << "nb sequences=" << nbSequences <<endl;
	cout << "nb kmers=" << nbKmers <<endl;

}

void Finder::writeBreakpoint(int bkt_id, string& chrom_name, uint64_t position, string& kmer_begin, string& kmer_end){
	fprintf(_breakpoint_file,">left_contig_%i_%s_pos_%i\n%s\n>right_contig_%i_%s_pos_%i\n%s\n",
			bkt_id,
			chrom_name.c_str(),
			position,
			kmer_begin.c_str(),
			bkt_id,
			chrom_name.c_str(),
			position,
			kmer_end.c_str()
	);
}
