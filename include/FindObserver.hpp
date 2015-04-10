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

#ifndef _TOOL_FindObserver_HPP_
#define _TOOL_FindObserver_HPP_

/*******************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>

template<size_t span>
class FindCleanInsert : public IFindObserver<span>
{
public :
    
    FindCleanInsert(FindBreakpoints<span> * find);

    void update(bool in_graph);
};

template<size_t span>
FindCleanInsert<span>::FindCleanInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
void FindCleanInsert<span>::update(bool in_graph)
{
    if(in_graph)
    {
	if(this->_find->solid_stretch_size > 1){
	    if(this->_find->gap_stretch_size == (this->_find->finder->_kmerSize-1)){
		// clean insert site
		string kmer_begin_str = this->_find->model.toString(this->_find->kmer_begin);
		string kmer_end_str = this->_find->model.toString(this->_find->kmer_end);
		this->_find->finder->writeBreakpoint(this->_find->breakpoint_id, this->_find->chrom_name, this->_find->position - 1, kmer_begin_str, kmer_end_str, 0);
		this->_find->breakpoint_id++;
	    }
	}
    }
}

template<size_t span>
class FindFuzzyInsert : public IFindObserver<span>
{
public :

    FindFuzzyInsert(FindBreakpoints<span> * find);

    void update(bool in_graph);
};

template<size_t span>
FindFuzzyInsert<span>::FindFuzzyInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
void FindFuzzyInsert<span>::update(bool in_graph)
{
    if(in_graph)
    {
	if(this->_find->solid_stretch_size > 1){
	    if(this->_find->gap_stretch_size < this->_find->finder->_kmerSize - 1 && this->_find->gap_stretch_size >= this->_find->finder->_kmerSize - 1 - this->_find->finder->_max_repeat){
		// Fuzzy site, position and kmer_end are impacted by the repeat

		int repeat_size = this->_find->finder->_kmerSize - 1 - this->_find->gap_stretch_size;
		string kmer_begin_str = this->_find->model.toString(this->_find->kmer_begin);
		string kmer_end_str = string(&(this->_find->chrom_sequence[this->_find->position - 1 + repeat_size]), this->_find->finder->_kmerSize);
		this->_find->finder->writeBreakpoint(this->_find->breakpoint_id, this->_find->chrom_name, this->_find->position - 1 + repeat_size, kmer_begin_str, kmer_end_str, repeat_size);
		this->_find->breakpoint_id++;
		this->_find->finder->_nb_homo_fuzzy++;
	    }
	}
    }
}

template<size_t span>
class FindEndSolid : public IFindObserver<span>
{
public :

    FindEndSolid(FindBreakpoints<span> * find);

    void update(bool in_graph);
};

template<size_t span>
FindEndSolid<span>::FindEndSolid(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
void FindEndSolid<span>::update(bool in_graph)
{
    if(in_graph)
    {
	if (this->_find->solid_stretch_size > 1) this->_find->gap_stretch_size = 0; // du coup on sort le trou a tai indexed ==2, gap_stretch_size pas remis a 0 par solide isole (FP)
	if (this->_find->solid_stretch_size==1) this->_find->kmer_end = this->_find->it_kmer->forward(); // kmer_end should be first kmer indexed after a hole
	if(this->_find->gap_stretch_size) this->_find->previous_gap_stretch_size = this->_find->gap_stretch_size;
    }
}

template<size_t span>
class FindEndGap : public IFindObserver<span>
{
public :

    FindEndGap(FindBreakpoints<span> * find);

    void update(bool in_graph);
};

template<size_t span>
FindEndGap<span>::FindEndGap(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
void FindEndGap<span>::update(bool in_graph)
{
    if(!in_graph)
    {
	if(this->_find->solid_stretch_size==1)
	{
	    this->_find->gap_stretch_size = this->_find->previous_gap_stretch_size + this->_find->solid_stretch_size - 1; //inutile maintenant il me semble, car tai_not_indexed non reset par FP
	}
	if(this->_find->solid_stretch_size > 1) // begin of not indexed zone
	{
	    this->_find->kmer_begin = this->_find->previous_kmer ;
	}
    }
}
#endif /* _TOOL_FindObserver_HPP_ */
