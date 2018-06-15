/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G.Rizk, P.Marijon
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

/**
 * \file FindBreakpoins.hpp
 * \date 09/04/2015
 * \author pmarijon
 * \brief FindBreakpoint definition class
 */
#ifndef _TOOL_FindBreakpoints_HPP_
#define _TOOL_FindBreakpoints_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <Finder.hpp>
#include <FindSNP.hpp>
#include "CircularBuffer.hpp"
/********************************************************************************/

template<size_t type>
class IFindObserver;

/**
 * \brief An observable functor for find gaps in reference genome
 *
 * This class associated with IFindObserver inherit, to find gaps in reference genome 
 */
template<size_t span>
class FindBreakpoints
{
public :

    typedef typename gatb::core::kmer::impl::Kmer<span> Kmer;

    typedef typename Kmer::ModelCanonical KmerModel;
    typedef typename Kmer::Type KmerType;
    typedef typename Kmer::Count KmerCount;
    typedef typename Kmer::KmerCanonical KmerCanonical;

    typedef typename KmerModel::Iterator KmerIterator;

    /** Variables for the heterozyguous mode */
    // structure to store information about a kmer : kmer that will be treated as kmer_begin of a breakpoint
    typedef struct info_type
    {
	KmerType kmer;
	int nb_in;
	int nb_out;
	bool is_repeated; // is the k-1 suffix of this kmer is repeated in the reference genome
    } info_type;

public :

    /** Constructor
     * \param[in] find : A pointeur one Finder instance
     */
    FindBreakpoints(Finder * find);

    /** Destructor. */
    virtual ~FindBreakpoints();
    
    //Functor
    /** overloading operator ()
     * Read reference genome, and find gaps
     */
    void operator()();

    // Observable
    /** Notify gap observer
     * \param[in] If kmer is in graph in_graph is true else is false
     */
    void notify(Node node, bool is_valid);

    /** Add observer call after a gap detection
     */
    void addGapObserver(IFindObserver<span>* new_obs);

    /** Add observer call after a gap detection
     */
    void addKmerObserver(IFindObserver<span>* new_obs);
    
    /** writes a given breakpoint in the output file
     */
    void writeBreakpoint(int bkt_id, string& chrom_name, uint64_t position, string& kmer_begin, string& kmer_end, int repeat_size, string type,bool repeat_in_genome_kmer_begin =false, bool repeat_in_genome_kmer_end = false);

    /** writes a given variant in the output vcf file
     */
    void writeVcfVariant(int bkt_id, string& chrom_name, uint64_t position, char* ref_char, char* alt_char, int repeat_size, string type);


    /*Getter*/
    /** Return the number of found breakpoints
     */
    uint64_t breakpoint_id();

    /** Return the position of first pb of actual read kmer
     */
    uint64_t position();

    /** Return the reference chromosome sequence
     */
    char* chrom_seq();

    /** Return the comment of sequence
     */
    string& chrom_name();

    /** Return the model of Kmer
     */
    KmerModel& model();

    /** Return the size of kmer used for gap search 
     */
    size_t kmer_size();

    /** Return the numbre of max repeat 
     */
    int max_repeat();

    /** Return the number of minimal kmer required to validate snp
     */
    int snp_min_val();

    
    /** The last solid kmer before gap
     */
    KmerCanonical& kmer_begin();

    /** The first solid kmer after gap
     */
    KmerCanonical& kmer_end();

    /** Size of current solid stretch
     */
    uint64_t solid_stretch_size();
    
    /** Size of current gap stretch
     */
    uint64_t gap_stretch_size();

    /** MindTheGap run with homo-only flag 
     */
    bool homo_only();

    /** Mask to get easily the k-1 prefix/suffix of a kmer (in KmerType unit)
     */
    KmerType kminus1_mask();

    /** Information one the next kmer
     */
    info_type& current_info();

    /** if >0, the precedent hetero site was too close, to avoid very close hetero sites
     */
    int recent_hetero();

    /**
     */
    bool kmer_end_is_repeated();
	
	bool kmer_begin_is_repeated();

    /** Get info_type at the index in circular buffer, rq : limits the kmerSize <256
     */
    info_type& het_kmer_history(unsigned char index);

    /**
     */
    unsigned char het_kmer_begin_index();
	unsigned char het_kmer_end_index();

    /**
     */
    bool graph_contains(Node& kmer_node);
    
    /*Iterater*/
    /** Incremente the value of breakpoint_id counter
     */
    uint64_t breakpoint_id_iterate();

    /** Incremente the value of homo_fuzzy_iterate
     */
    int homo_fuzzy_iterate();

    /** Incremente the value of homo_clean_iterate
     */
    int homo_clean_iterate();

    /** Incremente the value of hetero_fuzzy_iterate
     */
    int hetero_fuzzy_iterate();

    /** Incremente the value of hetero_clean_iterate
     */
    int hetero_clean_iterate();

    /** Incremente the value of fuzzy_deletion_iterate
     */
    int fuzzy_deletion_iterate();

    /** Incremente the value of fuzzy_deletion_iterate
     */
    int clean_deletion_iterate();

    /** Incremente the value of solo_snp_iterate
     */
    int solo_snp_iterate();

    /** Incremente the value of multi_snp_iterate
     */
    int multi_snp_iterate();

    /** Incremente the value of backup_iterate
     */
    int backup_iterate();
    
    /*Setter*/
    /** Set value of recent_hetero
     */
    void recent_hetero(int value);
    
private :

    IBloom<KmerType>* fillRefBloom();
    void store_kmer_info(Node node);

private :

    /*Observable membre*/
    std::vector<IFindObserver<span>* > gap_obs;
    std::vector<IFindObserver<span>* > kmer_obs;

    /*Find breakpoint membre*/
    /*Write breakpoint*/
    uint64_t m_breakpoint_id;
    uint64_t m_position;
    char* m_chrom_sequence;
    string m_chrom_name;

    /*Kmer related object*/
    KmerModel m_model;
    KmerCanonical m_previous_kmer;
    KmerIterator m_it_kmer;

    /*Kmer related object*/
    KmerCanonical m_kmer_begin;
    KmerCanonical m_kmer_end;

    /*Gap type detection*/
    uint64_t m_solid_stretch_size;
    uint64_t m_gap_stretch_size;
    
    /*Finder access*/
    Finder* finder;

    /*Hetero mode*/
    info_type m_het_kmer_history[256]; 
    unsigned char m_het_kmer_end_index; // index in history, must remain an unsigned char = same limit as the history array
    unsigned char m_het_kmer_begin_index;
    info_type m_current_info;
    int m_recent_hetero;
    bool m_kmer_end_is_repeated;
	bool m_kmer_begin_is_repeated;

	
//	CircularBuffer<info_type> m_het_kmer_history_CB;
//	typedef typename CircularBuffer<info_type>::itCB  iterCB;
//	iterCB* m_het_kmer_end_index_CB;
//	iterCB* m_het_kmer_begin_index_CB;
	
	
    /** Bloom of the repeated kmers of the reference genome 
     */
    IBloom<KmerType>* m_ref_bloom;

    //Please didn't add other friends please
    friend bool FindMultiSNP<span>::update();
	
	friend bool FindMultiSNPrev<span>::update();

	
	/** Handle on the progress information. */
	gatb::core::tools::dp::IteratorListener* _progress;
	void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }
	
};

template<size_t span>
FindBreakpoints<span>::FindBreakpoints(Finder * find) : gap_obs(), m_model(find->_kmerSize), m_it_kmer(m_model), _progress (0)
{
	this->m_breakpoint_id = 1;
	this->m_position = 0;
	this->m_chrom_sequence = NULL;
	this->m_chrom_name = "";
	this->m_kmer_begin = KmerCanonical(); // init kmerbegin and kmerend otherwise not init when checking this->_find->kmer_begin().isValid() in update
	this->m_kmer_end = KmerCanonical();
	
	//m_het_kmer_end_index_CB = new iterCB (&m_het_kmer_history_CB);
	//m_het_kmer_begin_index_CB = new iterCB (&m_het_kmer_history_CB);
	
	/*Homozygote usage*/
	this->m_solid_stretch_size = 0;
	this->m_gap_stretch_size = 0;
	
	this->finder = find;
	
	/*Heterozygote usage*/ //always fill repeat ref bloom
	//if(this->finder->_hete_insert)
	{
		this->m_ref_bloom = this->fillRefBloom();
		this->m_ref_bloom->use();
	}
}

template<size_t span>
FindBreakpoints<span>::~FindBreakpoints()
{
	for(typename std::vector<IFindObserver<span>* >::iterator it = this->kmer_obs.begin(); it != this->kmer_obs.end(); it++)
	{
		(*it)->forget();
	}
	
	for(typename std::vector<IFindObserver<span>* >::iterator it = this->gap_obs.begin(); it != this->gap_obs.end(); it++)
	{
		(*it)->forget();
	}
	
	//if(this->finder->_hete_insert)  //always fill repeat ref bloom
		this->m_ref_bloom->forget();
	
	setProgress             (0);
	
}

template<size_t span>
void FindBreakpoints<span>::operator()()
{
	// We create an iterator over this bank
	Iterator<Sequence>* it_seq = this->finder->_refBank->iterator();
	LOCAL(it_seq);
	
	
	u_int64_t  totalsize = this->finder->_refBank->estimateSequencesSize();
	u_int64_t nbkmersdone = 0;
	
	//printf("bank size %lli \n",totalsize);
	setProgress (new ProgressSynchro (
									  finder->createIteratorListener (totalsize, "Finding breakpoints"), //bon sang le createIteratorListener est dans le tool finder
									  System::thread().newSynchronizer())
				 );
	_progress->init ();
	
	
	
	// We loop over sequences
	for (it_seq->first(); !it_seq->isDone(); it_seq->next())
	{

		this->m_kmer_begin = KmerCanonical();
		this->m_kmer_end = KmerCanonical();
		//DEBUG
		//cout<<"sequence "<< (*it_seq)->getCommentShort() << endl;

		//Reintialize stretch_size for each sequence
		this->m_solid_stretch_size = 0;
		this->m_gap_stretch_size = 0;
		
		// for hetero mode:
		memset(this->m_het_kmer_history, 0, sizeof(info_type)*256);
		//m_het_kmer_history_CB.clear();
		
		this->m_het_kmer_end_index = this->finder->_kmerSize +1;
		this->m_het_kmer_begin_index = 1;
		
		//m_het_kmer_end_index_CB->set(this->finder->_kmerSize +1);
		//m_het_kmer_end_index_CB->set(1);
		
		this->m_recent_hetero = 0;
		
		// We set the data from which we want to extract kmers.
		m_it_kmer.setData ((*it_seq)->getData());
		this->m_chrom_sequence = (*it_seq)->getDataBuffer();
		this->m_chrom_name = (*it_seq)->getCommentShort();
		this->m_position = 0;
		
		// We iterate the kmers.
		for (m_it_kmer.first(); !m_it_kmer.isDone(); m_it_kmer.next(), m_position++, m_het_kmer_begin_index++, m_het_kmer_end_index++
			 ) //,m_het_kmer_begin_index_CB++, m_het_kmer_end_index++
		{
			if(!(*m_it_kmer).isValid())
			{
				this->m_solid_stretch_size = 0;
                this->m_gap_stretch_size = 0;
                this->m_kmer_begin = KmerCanonical();
                this->m_kmer_end = KmerCanonical();
                //DEBUG
				//cout<<"n";
			}
			else 
			{

			//we need to convert the kmer in a node to query the graph.
			Node node(Node::Value(m_it_kmer->value()), m_it_kmer->strand());// strand is necessary for hetero mode (in/out degree depends on the strand
			
			uint64_t save_position = m_position; // m_position can be modified by observer (multisnp rev)
			
			//we notify all observer
			this->notify(node, (*m_it_kmer).isValid());
			
			m_position = save_position;
			
			//save actual kmer for potential False Positive
			m_previous_kmer = *m_it_kmer;
			
			//if(!graph_contains(node) & (*m_it_kmer).isValid()) {cout << m_position << endl;}
			nbkmersdone++;
			if (nbkmersdone > 1000)   {  _progress->inc (nbkmersdone);  nbkmersdone = 0;  }
			}
		}

		//DEBUG
		//cout<<endl;
	}
	
	  _progress->finish ();
}

template<size_t span>
void FindBreakpoints<span>::notify(Node node, bool is_valid)
{
	bool in_graph = this->graph_contains(node);
	this->store_kmer_info(node);
	
	for(typename std::vector<IFindObserver<span>* >::iterator it = this->kmer_obs.begin(); it != this->kmer_obs.end(); it++)
	{
		(*it)->update();
	}
	
	// Kmer is in graph incremente scretch size
	if(in_graph && is_valid)
	{
		//DEBUG
		//cout<<"1";

		m_solid_stretch_size++;
		
		if(m_solid_stretch_size > 1 && m_gap_stretch_size > 0)
		{
			// Call each readonly observer
			for(typename std::vector<IFindObserver<span>* >::iterator it = this->gap_obs.begin(); it != this->gap_obs.end(); it++)
			{
				//DEBUG
				//cout << m_gap_stretch_size << endl;
				if((*it)->update())
				{
					break;
				}
			}
			
			// gap stretch size is re-set to 0 only when we are sure that the end of the gap is not due to an isolated solid kmer (likely FP)
			this->m_gap_stretch_size = 0;
		}
		
		if (this->m_solid_stretch_size==1)
		{
			// kmer_end should be the first kmer indexed after a gap (the first kmer of a solid_stretch is when m_solid_stretch_size=1)
			this->m_kmer_end = *this->m_it_kmer;
		}
	}
		
	// Kmer isn't in graph incremente gap size and reset solid size
	if(!in_graph && is_valid)
	{
		//DEBUG
		//cout<<"0";

		if(this->m_solid_stretch_size==1)
		{
			this->m_gap_stretch_size = this->m_gap_stretch_size + this->m_solid_stretch_size; //if previous position was an isolated solid kmer, we need to add 1 to the m_gap_stretch_size (as if replacing the FP by a non indexed kmer)
		}
		if(this->m_solid_stretch_size > 1 && this->m_previous_kmer.isValid()) // begin of not indexed zone
		{
			this->m_kmer_begin = this->m_previous_kmer;
			this->m_kmer_begin_is_repeated =  this->m_current_info.is_repeated ;
		}
		
		m_gap_stretch_size++;
		m_solid_stretch_size = 0;
	}
}

template<size_t span>
void FindBreakpoints<span>::addGapObserver(IFindObserver<span>* new_obs)
{
    new_obs->use();
    // Add observer in tables use unique_ptr for safety destruction
    this->gap_obs.push_back(new_obs);
}

template<size_t span>
void FindBreakpoints<span>::addKmerObserver(IFindObserver<span>* new_obs)
{
    new_obs->use();
    // Add observer in tables use unique_ptr for safety destruction
    this->kmer_obs.push_back(new_obs);
}

template<size_t span>
void FindBreakpoints<span>::writeBreakpoint(int bkt_id, string& chrom_name, uint64_t position, string& kmer_begin, string& kmer_end, int repeat_size, string type, bool repeat_in_genome_kmer_begin, bool repeat_in_genome_kmer_end  ){
    fprintf(this->finder->_breakpoint_file,">bkpt%i_%s_pos_%lli_fuzzy_%i_%s %s left_kmer\n%s\n>bkpt%i_%s_pos_%lli_fuzzy_%i_%s %s right_kmer\n%s\n",
	    bkt_id,
	    chrom_name.c_str(),
	    position+1, //switch to 1-based
	    repeat_size,
	    type.c_str(),
		repeat_in_genome_kmer_begin ? "REPEATED" : "",
	    kmer_begin.c_str(),
	    bkt_id,
	    chrom_name.c_str(),
	    position+1, //switch to 1-based
	    repeat_size,
	    type.c_str(),
		repeat_in_genome_kmer_end ? "REPEATED" : "",
	    kmer_end.c_str()
	);
}

template<size_t span>
void FindBreakpoints<span>::writeVcfVariant(int bkt_id, string& chrom_name, uint64_t position, char* ref_char, char* alt_char, int repeat_size, string type){
	//cout << ref_char << alt_char << endl;
	// NOTE : currently all positions coming from FindObservers are 0-based, VCF is supposed to be 1-based, so we add +1
	int variant_size=1;
	if (strcmp(type.c_str(),STR_DEL_TYPE)==0){
		variant_size = strlen(ref_char) - 1;
	}
	fprintf(this->finder->_vcf_file,"%s\t%lli\tbkpt%i\t%s\t%s\t.\tPASS\tTYPE=%s;LEN=%i;FUZZY=%i\tGT\t1/1\n",
			chrom_name.c_str(),
			position+1,  //switch to 1-based
			bkt_id,
			ref_char,
			alt_char,
			type.c_str(),
			variant_size,
			repeat_size
	);
}

/*Getter*/
template<size_t span>
uint64_t FindBreakpoints<span>::breakpoint_id()
{
    return this->m_breakpoint_id;
}

template<size_t span>
uint64_t FindBreakpoints<span>::position()
{
    return this->m_position;
}

template<size_t span>
char * FindBreakpoints<span>::chrom_seq()
{
    return this->m_chrom_sequence;
}

template<size_t span>
string& FindBreakpoints<span>::chrom_name()
{
    return this->m_chrom_name;
}

template<size_t span>
typename FindBreakpoints<span>::KmerModel& FindBreakpoints<span>::model()
{
    return this->m_model;
}

template<size_t span>
size_t FindBreakpoints<span>::kmer_size()
{
    return this->finder->_kmerSize;
}

template<size_t span>
int FindBreakpoints<span>::max_repeat()
{
    return this->finder->_max_repeat;
}

template<size_t span>
int FindBreakpoints<span>::snp_min_val()
{
    return this->finder->_snp_min_val;
}

/*Kmer related object*/
template<size_t span>
typename FindBreakpoints<span>::KmerCanonical& FindBreakpoints<span>::kmer_begin()
{
    return this->m_kmer_begin;
}

template<size_t span>
typename FindBreakpoints<span>::KmerCanonical& FindBreakpoints<span>::kmer_end()
{
    return this->m_kmer_end;
}

template<size_t span>
uint64_t FindBreakpoints<span>::solid_stretch_size()
{
    return this->m_solid_stretch_size();
}

template<size_t span>
uint64_t FindBreakpoints<span>::gap_stretch_size()
{
    return this->m_gap_stretch_size;
}

template<size_t span>
bool FindBreakpoints<span>::homo_only()
{
    return this->finder->_homo_only;
}

template<size_t span>
typename FindBreakpoints<span>::info_type& FindBreakpoints<span>::current_info()
{
    return this->m_current_info;
}

template<size_t span>
int FindBreakpoints<span>::recent_hetero()
{
    return this->m_recent_hetero;
}

template<size_t span>
bool FindBreakpoints<span>::kmer_end_is_repeated()
{
    return this->m_kmer_end_is_repeated;
}

template<size_t span>
bool FindBreakpoints<span>::kmer_begin_is_repeated()
{
	return this->m_kmer_begin_is_repeated;
}

template<size_t span>
typename FindBreakpoints<span>::info_type& FindBreakpoints<span>::het_kmer_history(unsigned char index)
{
    return this->m_het_kmer_history[index];
	
//	return  index(); //with index of type iterCB

}

template<size_t span>
unsigned char FindBreakpoints<span>::het_kmer_begin_index()
{
    return this->m_het_kmer_begin_index;
}


template<size_t span>
unsigned char FindBreakpoints<span>::het_kmer_end_index()
{
	return this->m_het_kmer_end_index;
}




template<size_t span>
bool FindBreakpoints<span>::graph_contains(Node& kmer_node)
{
    return this->finder->_graph.contains(kmer_node);
	//keep tips and internal node sonly
	//return	( this->finder->_graph.contains(kmer_node) && (this->finder->_graph.indegree(kmer_node)>=1 || this->finder->_graph.outdegree(kmer_node)>=1 ));
	
	//keep internal nodes only
	//	 return	( this->finder->_graph.contains(kmer_node) && (this->finder->_graph.indegree(kmer_node)>=1 && this->finder->_graph.outdegree(kmer_node)>=1 ));
	

}

/*Iterater*/
template<size_t span>
uint64_t  FindBreakpoints<span>::breakpoint_id_iterate()
{
    return this->m_breakpoint_id++;
}

template<size_t span>
int FindBreakpoints<span>::homo_fuzzy_iterate()
{
    return this->finder->_nb_homo_fuzzy++;
}

template<size_t span>
int FindBreakpoints<span>::homo_clean_iterate()
{
    return this->finder->_nb_homo_clean++;
}

template<size_t span>
int FindBreakpoints<span>::hetero_fuzzy_iterate()
{
    return this->finder->_nb_hetero_fuzzy++;
}

template<size_t span>
int FindBreakpoints<span>::hetero_clean_iterate()
{
    return this->finder->_nb_hetero_clean++;
}

template<size_t span>
int FindBreakpoints<span>::fuzzy_deletion_iterate()
{
    return this->finder->_nb_fuzzy_deletion++;
}

template<size_t span>
int FindBreakpoints<span>::clean_deletion_iterate()
{
    return this->finder->_nb_clean_deletion++;
}

template<size_t span>
int FindBreakpoints<span>::solo_snp_iterate()
{
    return this->finder->_nb_solo_snp++;
}

template<size_t span>
int FindBreakpoints<span>::multi_snp_iterate()
{
    return this->finder->_nb_multi_snp++;
}

template<size_t span>
int FindBreakpoints<span>::backup_iterate()
{
    return this->finder->_nb_backup++;
}

/*Setter*/
template<size_t span>
void FindBreakpoints<span>::recent_hetero(int value)
{
    this->m_recent_hetero = value;
}


//todo later replace this by mphf+ abundance per kmer
template<size_t span>
IBloom<typename FindBreakpoints<span>::KmerType>* FindBreakpoints<span>::fillRefBloom(){
	
	//Bloom of the repeated (k-1)mers of the reference genome
	IBloom<KmerType>* ref_bloom = 0;
	
	//solid kmers must be stored in a file
	string tempFileName = this->finder->getInput()->getStr(STR_URI_OUTPUT)+"_trashme.h5";
	
	// Parameters for SortingCountAlgorithm // all defaults
	IProperties* props = SortingCountAlgorithm<>::getDefaultProperties();
	props->setInt (STR_KMER_ABUNDANCE_MIN, this->finder->_het_max_occ+1);
	props->setInt (STR_KMER_SIZE,          this->finder->_kmerSize-1);
	props->setStr (STR_URI_OUTPUT,         tempFileName);
	//Remark : could re-use MAX_DISK or others from Finder options ? not necessary here, small counting in theory
	//props->setStr (STR_MAX_DISK, this->finder->getInput()->getStr(STR_MAX_DISK));
	
	/** We create a DSK (kmer counting) instance and execute it. */
	SortingCountAlgorithm<span> sortingCount (this->finder->_refBank,props);
	
	sortingCount.getInput()->add (0, STR_VERBOSE, 0);//do not show progress bar
	sortingCount.execute();
	
	// OLD WAY : Partition<KmerCount> & solidCollection = storage->root().getGroup("dsk").getPartition<KmerCount> ("solid");
	Partition<KmerCount> & solidCollection = * sortingCount.getSolidCounts();
	
	/** We get the number of solid kmers. */
	u_int64_t nb_solid = solidCollection.getNbItems();
	
	/** parameters of the Bloom filter */
	float NBITS_PER_KMER = 12;
	u_int64_t estimatedBloomSize = (u_int64_t) ((double)nb_solid * NBITS_PER_KMER * 2); //TODO *3 ?
	if (estimatedBloomSize ==0 )
	{
		estimatedBloomSize = 1000;
	}
	
	size_t nbHash = (int)floorf (0.7*NBITS_PER_KMER);
	
	//iterator of KmerCount
	Iterator<KmerCount>* itKmers = this->finder->createIterator(
																solidCollection.iterator(),
																nb_solid
																);
	LOCAL (itKmers);
	
	// building the bloom
	BloomBuilder<span> builder (estimatedBloomSize, nbHash, this->finder->_kmerSize-1, BLOOM_CACHE, this->finder->getDispatcher()->getExecutionUnitsNumber(), this->finder->_het_max_occ+1);
	ref_bloom = builder.build (itKmers);
	//cout << typeid(*ref_bloom).name() << endl;  // to verify the type of bloom
	
	System::file().remove(tempFileName);
	
	return ref_bloom;
}

template<size_t span>
void FindBreakpoints<span>::store_kmer_info(Node node)
{
	KmerType one; one.setVal(1);
	KmerType kminus1_mask = (one << ((this->finder->_kmerSize-1)*2)) - one;
	
	this->m_current_info.kmer = this->m_it_kmer->forward();
	if (this->finder->_graph.contains(node))
	{
		this->m_current_info.nb_in = this->finder->_graph.indegree (node);
		this->m_current_info.nb_out = this->finder->_graph.outdegree (node);
	}
	else
	{
		this->m_current_info.nb_in = 0;
		this->m_current_info.nb_out = 0;
	}
	
	//checking if the k-1 suffix is repeated
	KmerType suffix = this->m_it_kmer->forward() & kminus1_mask ; // getting the k-1 suffix (because putative kmer_begin)
	KmerType suffix_rev = revcomp(suffix,this->finder->_kmerSize-1); // we get its reverse complement to compute the canonical value of this k-1-mer
	
	//if(this->finder->_hete_insert) //alwayss fill repeat info
		this->m_current_info.is_repeated = this->m_ref_bloom->contains(min(suffix,suffix_rev));
	
	//filling the history array with the current kmer information
	this->m_het_kmer_history[m_het_kmer_end_index] = m_current_info;
	//m_het_kmer_end_index_CB->item() = m_current_info ;
	
	//checking if the k-1 prefix is repeated
	KmerType prefix = (this->m_it_kmer->forward() >> 2) & kminus1_mask; // getting the k-1 prefix (applying kminus1_mask after shifting of 2 bits to get the prefix)
	KmerType prefix_rev = revcomp(prefix,this->finder->_kmerSize-1); // we get its reverse complement to compute the canonical value of this k-1-mer
	
//	if(this->finder->_hete_insert) //alwayss fill repeat info
		this->m_kmer_end_is_repeated = this->m_ref_bloom->contains(min(prefix,prefix_rev));
}

#endif /* _TOOL_FindBreakpoints_HPP_ */

