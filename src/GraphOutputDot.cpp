/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen
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

#include <GraphOutputDot.hpp>
#include <iostream>
#include <fstream>

#include <cstdio>
#define DEBUG(a)    //printf a

using namespace std;



/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
GraphOutputDot<span>::GraphOutputDot (size_t kmerSize, const string& prefix)
    : IGraphOutput<span> (kmerSize,prefix)
{
    init (true);
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE :
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
void GraphOutputDot<span>::init(bool erase)
{

	_dot_file_suffix = ".graph"; //should be static (but pb with span)
    _dot_file_name  = (this->_prefix + _dot_file_suffix);
    _graph_file    = fopen (_dot_file_name.c_str(),          erase ? "w":"a");
    fprintf(_graph_file,"digraph dedebruijn {\n");
}

/*********************************************************************
 ** METHOD  :
 ** PURPOSE : write graph file or sequence file
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
void GraphOutputDot<span>::close()
{
    fprintf(_graph_file,"}\n");
    fclose(_graph_file);
}

/*********************************************************************
** METHOD  :
** PURPOSE : print head for a starter
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void GraphOutputDot<span>::print_starter_head (int index, char* sequence, size_t sequenceLen)
{

}

/*********************************************************************
** METHOD  :
** PURPOSE : write mark for end of nodes list
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void GraphOutputDot<span>::print_starter_end() // output a single node to a file
{

}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void GraphOutputDot<span>::print_sequence_head (const string& linear_seqs_name, const string& direction)
{

}

/*********************************************************************
** METHOD  :
** PURPOSE :
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void GraphOutputDot<span>::print_sequence_end ()
{

}

/*********************************************************************
** METHOD  :
** PURPOSE : output a single node to a file
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void GraphOutputDot<span>::print_node (long index, const string& seq) // output a single node to a file
{

	fprintf(_graph_file,"%ld [label=\"%s\"];\n",index,seq.c_str());

}

/*********************************************************************
** METHOD  :
** PURPOSE : output a single edges to a file
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void GraphOutputDot<span>::print_edge (long index, long id, long id2, const string& label, const string& comment)
{
	fprintf(_graph_file,"%ld -> %ld [label=\"%s\"];\n",id,id2,label.c_str());

}

// WARNING !!! The following code is not generic !!!
// It is designed to cope with 4 values of supported kmer size.

template class GraphOutputDot <KMER_SPAN(0)>;
template class GraphOutputDot <KMER_SPAN(1)>;
template class GraphOutputDot <KMER_SPAN(2)>;
template class GraphOutputDot <KMER_SPAN(3)>;


