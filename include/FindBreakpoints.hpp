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
#include <gatb/gatb_core.hpp>
#include <Finder.hpp>

/********************************************************************************/

class IFindBreakpoints
{
public :

    virtual void operator()() = 0;
};

template<size_t span>
class FindBreakpoints : public IFindBreakpoints
{
public :

    // Constructor
    FindBreakpoints(Finder * find);

    // Observable

    //Functor
    void operator()();

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

    /*Gap type detection*/
    uint64_t solid_stretch_size;
    uint64_t gap_stretch_size;
    uint64_t previous_gap_stretch_size;

private :

    Finder * finder;
};

template class FindBreakpoints<KSIZE_1>;
template class FindBreakpoints<KSIZE_2>;
template class FindBreakpoints<KSIZE_3>;
template class FindBreakpoints<KSIZE_4>;
#endif /* _TOOL_FindBreakpoints_HPP_ */
