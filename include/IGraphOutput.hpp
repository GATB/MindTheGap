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

#ifndef _IGRAPHOUTPUT_H
#define _IGRAPHOUTPUT_H

/********************************************************************************/
#include <gatb/gatb_core.hpp>
/********************************************************************************/

#ifdef DONT_USE_TR1
    #include <unordered_map>
    #include <functional>

    #define NS_TR1_BEGIN
    #define NS_TR1_END

    #define NS_TR1_PREFIX std
#else
    #include <tr1/unordered_map>
    #include <tr1/functional>

    #define NS_TR1_BEGIN  namespace tr1  {
    #define NS_TR1_END    }

    #define NS_TR1_PREFIX std::tr1
#endif

/********************************************************************************/

//structure for print id nodes and edges in graph output
struct id_els
{
    long node;
    long edge;
};

/********************************************************************************/

// hash functions for unordered_map with various kmer_type's
namespace std
{
NS_TR1_BEGIN
template <>  struct hash<Kmer<KSIZE_1>::Type> : public unary_function<Kmer<KSIZE_1>::Type, size_t>
{
    size_t operator()(const Kmer<KSIZE_1>::Type& elem) const  {  return hash1(elem);  }
};
template <>  struct hash<Kmer<KSIZE_2>::Type> : public unary_function<Kmer<KSIZE_2>::Type, size_t>
{
    size_t operator()(const Kmer<KSIZE_2>::Type& elem) const  {  return hash1(elem);  }
};
template <>  struct hash<Kmer<KSIZE_3>::Type> : public unary_function<Kmer<KSIZE_3>::Type, size_t>
{
    size_t operator()(const Kmer<KSIZE_3>::Type& elem) const  {  return hash1(elem);  }
};
template <>  struct hash<Kmer<KSIZE_4>::Type> : public unary_function<Kmer<KSIZE_4>::Type, size_t>
{
    size_t operator()(const Kmer<KSIZE_4>::Type& elem) const  {  return hash1(elem);  }
};
NS_TR1_END
}

/********************************************************************************/

template<size_t span>
class IGraphOutput
{
public:
    id_els first_id_els;

    // The extended kmer comes originally from the starter (true), or (false) if is it a degenerated kmer (one substitution or one indel).
    bool original;

public:

    /** Constructor.
     * \param[in] kmerSize : size of the kmer
     * \param[in] prefix : prefix of the file name
     * */
    IGraphOutput (size_t kmerSize, const std::string& prefix);

    /** Destructor. */
    virtual ~IGraphOutput() {}

    /** */
    void load_nodes_extremities (const std::string& linear_seqs_name);

    /** */
    id_els construct_graph (const std::string& linear_seqs_name, const std::string& direction);

    /** Finish the output. */
    virtual void close() = 0;

    virtual void print_starter_head (int index, char* sequence, size_t sequenceLen)  = 0;
    virtual void print_starter_end  () = 0;

    virtual void print_sequence_head (const std::string& filename, const std::string& direction) = 0;
    virtual void print_sequence_end  () = 0;

    virtual void print_node (long index, const std::string& seq) = 0;
    virtual void print_edge (long index, long id, long id2, const std::string& label, const std::string& comment) = 0;

    /** */
    void reset () {  first_id_els.node = first_id_els.edge = 0; }

protected:

    typedef typename Kmer<span>::Type  kmer_type;
    typedef typename Kmer<span>::ModelCanonical Model;
    typedef typename Model::Kmer  ModelKmer;

    enum LeftOrRight { LEFT=0, RIGHT=1 };

    Model _modelKmer;
    Model _modelKmerMinusOne;

    virtual void print_edges (const ModelKmer& kmer, size_t seqLen, LeftOrRight direction, id_els& nb_els);

    struct node_strand {
        long node;
        Strand strand;
        LeftOrRight left_or_right;
        node_strand(long node, Strand strand, LeftOrRight left_or_right) : node(node), strand(strand), left_or_right(left_or_right) {}
        bool operator<(const node_strand &other) const {
            if (node != other.node)
                return (node < other.node);
            if (left_or_right != other.left_or_right)
                return left_or_right < other.left_or_right;
            return (strand < other.strand);
        }
    };

    NS_TR1_PREFIX::unordered_map < kmer_type, std::set<node_strand> > kmer_links;

    std::string _prefix;
};

#endif //_IGRAPHOUTPUT_H
