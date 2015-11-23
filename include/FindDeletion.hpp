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
    
public:

    /** \copydoc IFindObserver<span>
     */
    FindDeletion(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::IFindObserver
     */
    bool update();

private:

    /** Detect if the end of a kmer is equal to the begin of other
     * \param[in] begin first kmer
     * \param[in] end the other kmer
     * \return The size of repetition						
     */
    unsigned int fuzzy_site(std::string begin, std::string end);
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

    // Test if deletion is a fuzzy deletion
    std::string begin = this->_find->model().toString(this->_find->kmer_begin().forward());
    std::string end = this->_find->model().toString(this->_find->kmer_end().forward());

    unsigned int repeat_size = this->fuzzy_site(begin, end);

    if(repeat_size > (unsigned)this->_find->max_repeat())
    {
	return false;
    }
    
    if(repeat_size != 0)
    {	
	begin = begin.substr(0, begin.length() - repeat_size);
    }

    // Compute del_size
    unsigned int del_size = this->_find->gap_stretch_size() - this->_find->kmer_size() + repeat_size + 1;

    // Create a sequence maybe is in graphe
    std::string seq = begin + end;

    // Create variable required for iterate on kmer
    KmerModel local_m(this->_find->kmer_size());
    KmerIterator local_it(local_m);
    Data local_d(const_cast<char*>(seq.c_str()));

    // Init this variable
    local_d.setRef(const_cast<char*>(seq.c_str()), (size_t)seq.length());
    local_it.setData(local_d);

    bool is_deletion = true;
    for(local_it.first(); !local_it.isDone(); local_it.next())
    {
        if(!this->contains(local_it->forward()))
    	{
    	    is_deletion = false;
	    break;
	}
    }

    if(is_deletion == false)
    {
	if(repeat_size == 0)
	{
	    return false;
	}
	else // Maybee isn't a fuzzy deletion
	{
	    seq = this->_find->model().toString(this->_find->kmer_begin().forward()) + end;
	    local_d.setRef(const_cast<char*>(seq.c_str()), (size_t)seq.length());
	    local_it.setData(local_d);

	    for(local_it.first(); !local_it.isDone(); local_it.next())
	    {
		if(!this->contains(local_it->forward()))
		{
		    return false;
		}
	    }
 
	    del_size -= repeat_size;
	    repeat_size = 0;
	}
    }
    
    // Write the breakpoint
    //this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - del_size - 1, begin, end, repeat_size, STR_DEL_TYPE);
    //NOTE : position will always be the left-most when repeat_size>0.
    size_t del_start_pos = this->_find->position() - 2 - del_size; //begining position of the deletion -1 (0-based): because in VCF we need to put the letter just before the deleted sequence
    //cout << "start pos = " << del_start_pos << "size = " << del_size << endl;
    char *del_sequence = new char[del_size+2];
    sprintf(del_sequence,"%.*s", del_size+1, this->_find->chrom_seq()+del_start_pos);
    char *alt_char = new char[2];
    sprintf(alt_char,"%.*s", 1, del_sequence);
    //cout << del_sequence << endl;
    //cout << alt_char << endl;
    // here position is 0-based
    this->_find->writeVcfVariant(this->_find->breakpoint_id(), this->_find->chrom_name(), del_start_pos, del_sequence, alt_char, repeat_size, STR_DEL_TYPE);

    delete(del_sequence);
    delete(alt_char);

    this->_find->breakpoint_id_iterate();

    if(repeat_size != 0)
	this->_find->fuzzy_deletion_iterate();
    else
	this->_find->clean_deletion_iterate();

    return true;
}

/*
  with max_repeat = 5
  good case 1 + 5 + 1 = 6 operation exemple AAAAATTCGG TTCGGCCCCC
*/
template<size_t span>
unsigned int FindDeletion<span>::fuzzy_site(std::string begin, std::string end)
{
    for(unsigned int i = this->_find->max_repeat(); i != 0; i--)
	for(unsigned int j = 1; begin.substr(begin.length() - i, j) == end.substr(0, j); j++)
	    if(i == j)
		return j;
	    
    return 0;
}

#endif /* _TOOL_FindDeletion_HPP_ */














