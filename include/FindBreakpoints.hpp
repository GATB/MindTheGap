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

/**
 * \file FindBreakpoins.hpp
 * \date 09/04/2015
 * \author pmarijon
 * \brief FindBreakpoint definition class
 */
#ifndef _TOOL_FindBreakpoints_HPP_
#define _TOOL_FindBreakpoints_HPP_

/********************************************************************************/
#include <memory>
#include <gatb/gatb_core.hpp>
#include <Finder.hpp>

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

    typedef typename gatb::core::kmer::impl::Kmer<span>::ModelCanonical KmerModel;
    typedef typename KmerModel::Iterator KmerIterator;
    typedef typename gatb::core::kmer::impl::Kmer<span>::Type KmerType;

public :

    /** Constructor
     * \param[in] find : A pointeur one Finder instance
     */
    FindBreakpoints(Finder * find, IFindObserver<span>* backup);

    /** Destructor. */
    virtual ~FindBreakpoints() {}
    
    //Functor
    /** overloading operator ()
     * Read reference genome, and find gaps
     */
    void operator()();

    // Observable
    /** Notify all observer
     * \param[in] If kmer is in graph in_graph is true else is false
     */
    void notify(bool in_graph);

    /** Add an observer in the observer list
     */
    void addObserver(IFindObserver<span>* new_obs);

    /** writes a given breakpoint in the output file
     */
    void writeBreakpoint(int bkt_id, string& chrom_name, uint64_t position, string& kmer_begin, string& kmer_end, int repeat_size);

    /*Getter*/
    /** Return the number of found breakpoints
     */
    uint64_t breakpoint_id();

    /** Return the position of first pb of actual read kmer
     */
    uint64_t position();

    /** Return the reference sequence of this kmer 
     */
    char * chrom_seq();

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

    /** The last solid kmer before gap
     */
    KmerType& kmer_begin();

    /** The first solid kmer after gap
     */
    KmerType& kmer_end();

    /** Size of current solid stretch
     */
    uint64_t solid_stretch_size();
    
    /** Size of current gap stretch
     */
    uint64_t gap_stretch_size();
    
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

private :

    /*Observable membre*/
    std::vector<std::unique_ptr<IFindObserver<span> > > list_obs;
    std::unique_ptr<IFindObserver<span> > m_backup;

    /*Find breakpoint membre*/
    /*Write breakpoint*/
    uint64_t m_breakpoint_id;
    uint64_t m_position;
    char * m_chrom_sequence;
    string m_chrom_name;

    /*Kmer related object*/
    KmerModel m_model;
    KmerType m_previous_kmer;
    KmerIterator m_it_kmer;

    /*Kmer related object*/
    KmerType m_kmer_begin;
    KmerType m_kmer_end;

    /*Gap type detection*/
    uint64_t m_solid_stretch_size;
    uint64_t m_gap_stretch_size;
    
    /*Finder access*/
    Finder * finder;
};

template<size_t span>
FindBreakpoints<span>::FindBreakpoints(Finder * find, IFindObserver<span>* backup) : list_obs(), m_model(find->_kmerSize), m_it_kmer(m_model), m_backup(backup)
{
    this->m_breakpoint_id = 0;
    this->m_position = 0;
    this->m_chrom_sequence = NULL;
    this->m_chrom_name = "";

    this->m_solid_stretch_size = 0;
    this->m_gap_stretch_size = 0;

    this->finder = find;
}

template<size_t span>
void FindBreakpoints<span>::operator()()
{
    // We create an iterator over this bank
    Iterator<Sequence>* it_seq = this->finder->_refBank->iterator();
    LOCAL(it_seq);
    
    // We loop over sequences
    for (it_seq->first(); !it_seq->isDone(); it_seq->next())
    {
	//Reintialize stretch_size for each sequence
	this->m_solid_stretch_size = 0;
	this->m_gap_stretch_size = 0;

        //Method : an homozyguous breakpoint is detected as a gap_stretch (ie. consecutive kmers on the sequence, that are not indexed in the graph) of particular sizes.
	//BUT there are some False Positives when we query the graph : when we ask the graph if a kmer is indexed (when starting from a previous not indexed kmer) it may wrongly answer yes
	//(the gatb dbg is exact only if the kmer we query is the neighbor of a truly solid kmer)
	//FP are likely to be isolated, ie. surrounded by not indexed kmers, therefore they can be detected as solid_stretches of size 1.
	// See FindObserver for code
	
	
	// We set the data from which we want to extract kmers.
	m_it_kmer.setData ((*it_seq)->getData());
	this->m_chrom_sequence = (*it_seq)->getDataBuffer();
	this->m_chrom_name = (*it_seq)->getComment();
	this->m_position = 0;
	
	// We iterate the kmers.
	for (m_it_kmer.first(); !m_it_kmer.isDone(); m_it_kmer.next(), m_position++)
	{
	    //we need to convert the kmer in a node to query the graph.
	    Node node(Node::Value(m_it_kmer->value()));

	    //we notify all observer
	    this->notify(this->finder->_graph.contains(node));

	    //save actual kmer for potential False Positive
	    m_previous_kmer = m_it_kmer->forward();
	}
    }
}

template<size_t span>
void FindBreakpoints<span>::notify(bool in_graph)
{
    // Kmer is in graph incremente scretch size
    if(in_graph)
    {
	m_solid_stretch_size++;

	if(m_solid_stretch_size > 1)
	{
	    bool one_observer_ret_true = false;
	    // Call each readonly observer
	    for(auto it = this->list_obs.begin(); it != this->list_obs.end(); it++)
	    {
		bool current_observer_ret = (*it)->update();
		if(!one_observer_ret_true && current_observer_ret)
		{
		    one_observer_ret_true = true;
		}
	    }

	    if(!one_observer_ret_true)
	    {
		this->m_backup->update();
	    }

	    // gap stretch size is re-set to 0 only when we are sure that the end of the gap is not due to an isolated solid kmer (likely FP)
	    this->m_gap_stretch_size = 0; 
	}

	if (this->m_solid_stretch_size==1)
	{
	    // kmer_end should be the first kmer indexed after a gap (the first kmer of a solid_stretch is when m_solid_stretch_size=1)
	    this->m_kmer_end = this->m_it_kmer->forward();
	}
    }
    
    // Kmer isn't in graph incremente gap size and reset solid size
    if(!in_graph)
    {
	if(this->m_solid_stretch_size==1)
	{
	    this->m_gap_stretch_size = this->m_gap_stretch_size + this->m_solid_stretch_size; //if previous position was an isolated solid kmer, we need to add 1 to the m_gap_stretch_size (as if replacing the FP by a non indexed kmer)
	}
	if(this->m_solid_stretch_size > 1) // begin of not indexed zone
	{
	    this->m_kmer_begin = this->m_previous_kmer;
	}
	
        m_gap_stretch_size++;
	m_solid_stretch_size = 0;
    }
}

template<size_t span>
void FindBreakpoints<span>::addObserver(IFindObserver<span>* new_obs)
{
    // Add observer in tables use unique_ptr for safety destruction
    this->list_obs.push_back(std::unique_ptr<IFindObserver<span> >(new_obs));
}

template<size_t span>
void FindBreakpoints<span>::writeBreakpoint(int bkt_id, string& chrom_name, uint64_t position, string& kmer_begin, string& kmer_end, int repeat_size){
	fprintf(this->finder->_breakpoint_file,">left_contig_%i_%s_pos_%lli_repeat_%i\n%s\n>right_contig_%i_%s_pos_%lli_repeat_%i\n%s\n",
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

/*Kmer related object*/
template<size_t span>
typename FindBreakpoints<span>::KmerType& FindBreakpoints<span>::kmer_begin()
{
    return this->m_kmer_begin;
}

template<size_t span>
typename FindBreakpoints<span>::KmerType& FindBreakpoints<span>::kmer_end()
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

#endif /* _TOOL_FindBreakpoints_HPP_ */
