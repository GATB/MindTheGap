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

#ifndef _TOOL_IFindObserver_HPP_
#define _TOOL_IFindObserver_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <FindBreakpoints.hpp>

/********************************************************************************/

template<size_t span>
class FindBreakpoints;

template<size_t span>
class IFindObserver
{
public:

IFindObserver(FindBreakpoints<span>* find);

    virtual void update(bool in_graph) = 0;

protected :

    FindBreakpoints<span>* _find;
};

template<size_t span>
IFindObserver<span>::IFindObserver(FindBreakpoints<span>* find)
{
    this->_find = find;
}

#endif /* _TOOL_IFindObserver_HPP_ */
