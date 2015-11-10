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

#include <IGraphOutput.hpp>
#include <iostream>
#include <fstream>

#include <cstdio>
#define DEBUG(a)    //printf a

using namespace std;

/*********************************************************************
 ** METHOD  :
 ** PURPOSE : Initialize first elements and files  (files are erasing)
 ** INPUT   :
 ** OUTPUT  :
 ** RETURN  :
 ** REMARKS :
 *********************************************************************/
template<size_t span>
IGraphOutput<span>::IGraphOutput (size_t kmerSize, const string& prefix)
    : original(true), _modelKmer(kmerSize), _modelKmerMinusOne(kmerSize-1), _prefix(prefix)
{
    first_id_els.node = 0;
    first_id_els.edge = 0;
}

/*********************************************************************
** METHOD  :
** PURPOSE : load nodes extremities
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void IGraphOutput<span>::load_nodes_extremities (const string& linear_seqs_name)
{
    kmer_links.clear();

    /** We open the bank that holds the extensions. */
    IBank* Nodes = Bank::open (linear_seqs_name);    LOCAL (Nodes);

    long nb_nodes = first_id_els.node;

    DEBUG (("[GraphOutput::load_nodes_extremities]  kmerSize=%ld  name=%s\n", _modelKmerMinusOne.getKmerSize(), linear_seqs_name.c_str()));

    Iterator<Sequence>* itSeq = Nodes->iterator();   LOCAL (itSeq);
    for (itSeq->first(); !itSeq->isDone(); itSeq->next())
    {
        char* rseq    = itSeq->item().getDataBuffer();
        int   readlen = itSeq->item().getDataSize();

        DEBUG (("[GraphOutput::load_nodes_extremities]  seq.size=%ld\n", readlen));

        ModelKmer leftKmer  = _modelKmerMinusOne.codeSeed (rseq, Data::ASCII, 0);
        ModelKmer rightKmer = _modelKmerMinusOne.codeSeed (rseq, Data::ASCII, readlen-_modelKmerMinusOne.getKmerSize());

        kmer_links[leftKmer. value()].insert (node_strand(nb_nodes, leftKmer.strand(),  LEFT));
        kmer_links[rightKmer.value()].insert (node_strand(nb_nodes, rightKmer.strand(), RIGHT));

        nb_nodes++;
    }

    DEBUG (("[GraphOutput::load_nodes_extremities]  nbNodes=%ld\n", nb_nodes));
}

/*********************************************************************
** METHOD  :
** PURPOSE : construct node file and edge file for graph file
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
id_els IGraphOutput<span>::construct_graph (const string& linear_seqs_name, const string& direction)
{
    DEBUG (("[GraphOutput::construct_graph]  linear_seqs_name=%s   direction=%s\n", linear_seqs_name.c_str(), direction.c_str() ));

    /** We open the bank that holds the extensions. */
    IBank* Nodes = Bank::open (linear_seqs_name);    LOCAL (Nodes);

    id_els nb_els = first_id_els;
    bool found = false;

    print_sequence_head (linear_seqs_name, direction);

    Iterator<Sequence>* itSeq = Nodes->iterator();   LOCAL (itSeq);
    for (itSeq->first(); !itSeq->isDone(); itSeq->next())
    {
        /** We get the current sequence as a string object.
         * NOTE: we can't rely on strlen() on itSeq->item().getDataBuffer() because the buffer may be not
         * 0 terminated (it's just a buffer, not a C-like string). */
        string seq;  seq.assign (itSeq->item().getDataBuffer(), itSeq->item().getDataSize());

        ModelKmer leftKmer  = _modelKmerMinusOne.codeSeed (seq.c_str(), Data::ASCII, 0);
        ModelKmer rightKmer = _modelKmerMinusOne.codeSeed (seq.c_str(), Data::ASCII, seq.size()-_modelKmerMinusOne.getKmerSize());

        print_edges (leftKmer,  seq.size(), LEFT,  nb_els);  // left edges (revcomp extensions)
        print_edges (rightKmer, seq.size(), RIGHT, nb_els);  // right edges

        print_node (nb_els.node, seq);
        nb_els.node++;

        found = true ;

    } /* end of for (itSeq->first(); !itSeq->isDone(); itSeq->next()) */

    if (found)  {  print_sequence_end ();  }

    return nb_els;
}

/*********************************************************************
** METHOD  :
** PURPOSE : construct node file and edge file for graph file
** INPUT   :
** OUTPUT  :
** RETURN  :
** REMARKS :
*********************************************************************/
template<size_t span>
void IGraphOutput<span>::print_edges (const ModelKmer& kmer, size_t seqLen, LeftOrRight direction, id_els& nb_els)
{
    /** Shortcut. */
    std::set<node_strand>& nodes = kmer_links[kmer.value()];

    size_t sizeKmer = _modelKmerMinusOne.getKmerSize();

    static const char* table0[] = { "R", "F" };
    static const char* table1[] = { "F", "R" };

    for (typename set<node_strand>::iterator it = nodes.begin(); it != nodes.end(); it++)
    {
        long        cur_node          = it->node;
        Strand      cur_strand        = it->strand;
        LeftOrRight cur_left_or_right = it->left_or_right;

        // prevent self loops on same kmer
        if (cur_node == nb_els.node)  {  if (seqLen == sizeKmer)  { continue; } }

        string label = table0[direction];

        if (cur_left_or_right == direction)
        {
            if (cur_strand != kmer.strand())    {  label += table1[direction];  }
            else                                {  continue;                    }
        }
        else
        {
            if (cur_strand == kmer.strand())    {  label += table0[direction];  }
            else                                {  continue;                    }
        }

        print_edge (nb_els.edge, nb_els.node, cur_node, label, "");
        nb_els.edge++;
    }
}

// WARNING !!! The following code is not generic !!!
// It is designed to cope with 4 values of supported kmer size.

template class IGraphOutput <KMER_SPAN(0)>;
template class IGraphOutput <KMER_SPAN(1)>;
template class IGraphOutput <KMER_SPAN(2)>;
template class IGraphOutput <KMER_SPAN(3)>;


