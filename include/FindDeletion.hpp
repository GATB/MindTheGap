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

#ifndef _TOOL_FindDeletion_HPP_
#define _TOOL_FindDeletion_HPP_

/*****************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>


template<size_t span>
class FindDeletion : public IFindObserver<span>
{
public :
    typedef typename gatb::core::kmer::impl::Kmer<span> Kmer;

    typedef typename Kmer::ModelCanonical KmerModel;
    typedef typename KmerModel::Iterator KmerIterator;
    
public :

    /** \copydoc IFindObserver<span>
     */
    FindDeletion(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::IFindObserver
     */
    bool update();
};

template<size_t span>
FindDeletion<span>::FindDeletion(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindDeletion<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }

    std::string begin = this->_find->model().toString(this->_find->kmer_begin().forward());
    std::string end = this->_find->model().toString(this->_find->kmer_end().forward());

    // Detect if deletion is a fuzzy deletion
    unsigned int repeat_size;
    for(repeat_size = 0; begin.substr(begin.length() - repeat_size, 1) == end.substr(repeat_size, 1); repeat_size++);
    if(repeat_size != 0)
    {
	unsigned char hist_end_pos = this->_find->position() - 2 - repeat_size % 256;
	unsigned char hist_begin_pos = hist_end_pos - this->_find->gap_stretch_size() + repeat_size;
	begin = this->_find->model().toString(this->_find->het_kmer_history(hist_begin_pos).kmer);
	end = this->_find->model().toString(this->_find->het_kmer_history(hist_end_pos).kmer);
    }
    
    unsigned int del_size = this->_find->gap_stretch_size() - this->_find->kmer_size();
    
    std::string seq = begin + end;
    
    KmerModel local_m(this->_find->kmer_size());
    KmerIterator local_it(local_m);
    Data local_d(const_cast<char*>(seq.c_str()));
    local_d.setRef(const_cast<char*>(seq.c_str()), (size_t)seq.length());
    local_it.setData(local_d);

    for(local_it.first(); !local_it.isDone(); local_it.next())
    {
        if(!this->contains(local_it->forward()))
    	{
    	    return false;
    	}
    }
    
    // obtains the kmer sequence
    string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
    string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());

    this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - del_size - 2 + repeat_size, kmer_begin_str, kmer_end_str, 0, STR_DEL_TYPE);
    this->_find->breakpoint_id_iterate();

    if(repeat_size != 0)
	this->_find->fuzzy_deletion_iterate();
    else
	this->_find->clean_deletion_iterate();

    return true;
}

#endif /* _TOOL_FindDeletion_HPP_ */
