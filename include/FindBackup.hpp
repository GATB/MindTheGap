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

#ifndef _TOOL_FindBackup_HPP_
#define _TOOL_FindBackup_HPP_

/*******************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>

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
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }

    if(this->_find->gap_stretch_size() > (this->_find->kmer_size() / 2)) {
	string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
	string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());
	string chrom_name_bak = this->_find->chrom_name()+"_backup";

	this->_find->writeBreakpoint(this->_find->breakpoint_id(), chrom_name_bak, this->_find->position() - 1, kmer_begin_str, kmer_end_str, 0, STR_BKP_TYPE);

	this->_find->breakpoint_id_iterate();
	this->_find->backup_iterate();
	
	return true;
    }

    return false;
}

#endif /* _TOOL_FindBackup_HPP_ */
