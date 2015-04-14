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
    
    /** \copydoc IFindObserver::IFindObserver
     */
    FindCleanInsert(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindCleanInsert<span>::FindCleanInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindCleanInsert<span>::update()
{
    if(this->_find->gap_stretch_size() == (this->_find->kmer_size()-1)) //Check size of gap 
    {
	// obtains the kmer sequence
	string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin());
	string kmer_end_str = this->_find->model().toString(this->_find->kmer_end());

	this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1, kmer_begin_str, kmer_end_str, 0);

	// iterate counter
	this->_find->breakpoint_id_iterate();
	this->_find->homo_clean_iterate();

	return true;
    }
    
    return false;
}

template<size_t span>
class FindFuzzyInsert : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindobserver
     */
    FindFuzzyInsert(FindBreakpoints<span> * find);
    
    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindFuzzyInsert<span>::FindFuzzyInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindFuzzyInsert<span>::update()
{
    if(this->_find->gap_stretch_size() < this->_find->kmer_size() - 1 && this->_find->gap_stretch_size() >= this->_find->kmer_size() - 1 - this->_find->max_repeat())
    {
	// Fuzzy site, position and kmer_end are impacted by the repeat
	int repeat_size = this->_find->kmer_size() - 1 - this->_find->gap_stretch_size();
	    
	// obtains the kmer sequence
	string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin());
	string kmer_end_str = string(&(this->_find->chrom_seq()[this->_find->position() - 1 + repeat_size]), this->_find->kmer_size());
	    
	this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1 + repeat_size, kmer_begin_str, kmer_end_str, repeat_size);

	//iterate counter
	this->_find->breakpoint_id_iterate();
	this->_find->homo_fuzzy_iterate();

	return true;
    }

    return false;
}

template<size_t span>
class FindSoloSNP : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindSoloSNP(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindSoloSNP<span>::FindSoloSNP(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindSoloSNP<span>::update()
{
    if(this->_find->gap_stretch_size() == this->_find->kmer_size())
    {
	//for each kmer :
	//    if cor_nuc != undef && in_graph(corect(kmer, cor_nuc))
	//         kmer_corrected++
	//     else
	//         
	std::cout<<"Potential solo SNP"<<std::endl;
	return true;
    }

    return false;
}

template<size_t span>
class FindBackup : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindBackup(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindBackup<span>::FindBackup(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindBackup<span>::update()
{

}

#endif /* _TOOL_FindObserver_HPP_ */
