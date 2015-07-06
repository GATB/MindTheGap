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

#ifndef _TOOL_FindHetero_HPP_
#define _TOOL_FindHetero_HPP_

/*******************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>

template<size_t span>
class FindHeteroInsert : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindHeteroInsert(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindHeteroInsert<span>::FindHeteroInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindHeteroInsert<span>::update()
{
    if(!this->_find->homo_only())
    {
        // hetero site detection
	if(!this->_find->kmer_end_is_repeated() && this->_find->current_info().nb_in == 2 && !this->_find->recent_hetero())
	{
	    //loop over putative repeat size (0=clean, >0 fuzzy), reports only the smallest repeat size found.
	    for(int i = 0; i <= this->_find->max_repeat(); i++)
	    {
		if(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 2 && !this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).is_repeated)
		{
		    //hetero breakpoint found
		    string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).kmer);
		    string kmer_end_str = this->_find->model().toString(this->_find->current_info().kmer);
		    this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position()-1+i, kmer_begin_str, kmer_end_str,i, STR_HET_TYPE);

		    this->_find->breakpoint_id_iterate();

		    if(i==0)
		    {
			this->_find->hetero_clean_iterate();
		    }
		    else
		    {
			this->_find->hetero_fuzzy_iterate();
		    }
		    
		    this->_find->recent_hetero(this->_find->max_repeat()); // we found a breakpoint, the next hetero one mus be at least _max_repeat apart from this one.
		    return true; //reports only the smallest repeat size found.
		}
	    }
	}

	this->_find->recent_hetero(max(0, this->_find->recent_hetero() - 1));  // when recent_hetero=0 : we are sufficiently far from the previous hetero-site
    }

    return false;
}

#endif /* _TOOL_FindHetero_HPP_ */
