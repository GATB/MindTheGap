/*****************************************************************************
 *   GATB : Genome Assembly Tool Box
 *   Copyright (C) 2014  INRIA
 *   Authors: R.Chikhi, G.Rizk, E.Drezen, C.Lemaitre
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

#ifndef _GRAPHOUTPUTDOT_H
#define _GRAPHOUTPUTDOT_H

/********************************************************************************/
#include <IGraphOutput.hpp>
/********************************************************************************/




template<size_t span>
class GraphOutputDot : public IGraphOutput<span>
{
public:

    /** Constructor.
     * \param[in] kmerSize : size of the kmer
     * \param[in] prefix : prefix of the file name
     * */
    GraphOutputDot (size_t kmerSize, const std::string& prefix);

    /** Finish the output. */
    virtual void close();

    virtual void print_starter_head (int index, char* sequence, size_t sequenceLen);
    virtual void print_starter_end  ();

    virtual void print_sequence_head (const std::string& filename, const std::string& direction);
    virtual void print_sequence_end  ();

    virtual void print_node (long index, const std::string& seq);
    virtual void print_edge (long index, long id, long id2, const std::string& label, const std::string& comment);

    std::string get_dot_file_name() {return _dot_file_name;};

private:

    void init (bool erase);

    FILE* _graph_file;

    std::string _dot_file_name;

    std::string _dot_file_suffix;
};

/********************************************************************************/

#endif //_GRAPHOUTPUTDOT_H
