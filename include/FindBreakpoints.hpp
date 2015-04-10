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

#ifndef _TOOL_FindBreakpoints_HPP_
#define _TOOL_FindBreakpoints_HPP_

/********************************************************************************/
#include <memory>
#include <gatb/gatb_core.hpp>
#include <Finder.hpp>

/********************************************************************************/

template<size_t type>
class IFindObserver;

template<size_t span>
class FindBreakpoints
{
public :

    // Constructor
    FindBreakpoints(Finder * find);

    //Functor
    void operator()();

    // Observable
    void notify(bool in_graph);
    void addObserver(IFindObserver<span>* new_obs);

public :

    typedef typename gatb::core::kmer::impl::Kmer<span>::ModelCanonical KmerModel;
    typedef typename KmerModel::Iterator KmerIterator;
    typedef typename gatb::core::kmer::impl::Kmer<span>::Type KmerType;

public :

    /*Write breakpoint*/
    uint64_t breakpoint_id;
    uint64_t position;
    char * chrom_sequence;
    string chrom_name;

    /*Kmer related object*/
    KmerModel model;
    KmerType kmer_begin;
    KmerType kmer_end;
    KmerType previous_kmer;
    KmerIterator it_kmer;

    /*Gap type detection*/
    uint64_t solid_stretch_size;
    uint64_t gap_stretch_size;
    uint64_t previous_gap_stretch_size;

    Finder * finder;

private :

    std::vector<std::unique_ptr<IFindObserver<span> > > list_obs;
};

template<size_t span>
FindBreakpoints<span>::FindBreakpoints(Finder * find) : list_obs(), model(find->_kmerSize), it_kmer(model)
{
    this->breakpoint_id = 0;
    this->position = 0;
    this->chrom_sequence = NULL;
    this->chrom_name = "";

    this->solid_stretch_size = 0;
    this->gap_stretch_size = 0;
    this->previous_gap_stretch_size = 0;

    this->finder = find;
}

template<size_t span>
void FindBreakpoints<span>::operator()()
{
    BankFasta::Iterator it_seq(*(finder->_refBank));

    for (it_seq.first(); !it_seq.isDone(); it_seq.next())
    {
	this->solid_stretch_size = 0;
	this->gap_stretch_size = 0;
	this->previous_gap_stretch_size = 0;
		
	// We set the data from which we want to extract kmers.
	it_kmer.setData (it_seq->getData());
	this->chrom_sequence = it_seq->getDataBuffer();
	this->chrom_name = it_seq->getComment();
	this->position = 0;
	
	// We iterate the kmers.
	for (it_kmer.first(); !it_kmer.isDone(); it_kmer.next(), position++)
	{
	    //we need to convert the kmer in a node to query the graph.
	    Node node(Node::Value(it_kmer->value()));
	    this->notify(this->finder->_graph.contains(node));
	    previous_kmer = it_kmer->forward();
	}
    }
}

template<size_t span>
void FindBreakpoints<span>::notify(bool in_graph)
{
    if(in_graph)
    {
	solid_stretch_size++;
    }

    for(auto it = this->list_obs.begin(); it != this->list_obs.end(); it++)
    {
	(*it)->update(in_graph);
    }

    if(!in_graph)
    {
        gap_stretch_size++;
	solid_stretch_size = 0;
    }
}

template<size_t span>
void FindBreakpoints<span>::addObserver(IFindObserver<span>* new_obs)
{
    this->list_obs.push_back(std::unique_ptr<IFindObserver<span> >(new_obs));
}

#endif /* _TOOL_FindBreakpoints_HPP_ */
