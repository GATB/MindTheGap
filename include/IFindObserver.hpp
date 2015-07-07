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
 * \file IFindObserver.hpp
 * \date 09/04/2015
 * \author pmarijon
 * \brief Interface definition for FindBreakpoints observer
 */
#ifndef _TOOL_IFindObserver_HPP_
#define _TOOL_IFindObserver_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <FindBreakpoints.hpp>

/********************************************************************************/

template<size_t span>
class FindBreakpoints;

/** \brief Interface for FindBreakpoints observer
 * Implementation can be add in FindBreakpoints observer list and call by update.
 */
template<size_t span>
class IFindObserver : public SmartPointer
{
public:
    typedef typename gatb::core::kmer::impl::Kmer<span> Kmer;

    typedef typename Kmer::Type KmerType;
    
public:

    /** Constructor.
     * \param[in,out] The FindBreakpoints instance, implementation can read and 
     * write information one this instance.
     */
    IFindObserver(FindBreakpoints<span>* find);

    /** Destructor.
     */
    virtual ~IFindObserver() {}

    /** Called when FindBreakpoints::notify is called
     * \param[in] kmer is in graph or not
     */
    virtual bool update() = 0;

protected :

    /** Pointer one FindBreakpoints instance
     */
    FindBreakpoints<span>* _find;
    bool contains(KmerType kmer);
};

template<size_t span>
IFindObserver<span>::IFindObserver(FindBreakpoints<span>* find)
{
    this->_find = find;
}

template<size_t span>
bool IFindObserver<span>::contains(KmerType kmer)
{
    kmer = std::min(kmer, revcomp(kmer, this->_find->kmer_size()));
    Node node = Node(Node::Value(kmer));
    return this->_find->graph_contains(node);
}

#endif /* _TOOL_IFindObserver_HPP_ */

