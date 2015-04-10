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

#ifndef _TOOL_Finder_HPP_
#define _TOOL_Finder_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
using namespace std;

/********************************************************************************/

static const char* STR_URI_REF = "-ref";
static const char* STR_MAX_REPEAT = "-max-rep";
static const char* STR_HOMO_ONLY = "-homo-only";


class Finder : public Tool
{
public:

    // Constructor
    Finder ();
    
    size_t _kmerSize;
    Graph _graph;
    int _max_repeat;
    int _nbCores;
    bool _homo_only;
    BankFasta* _refBank;
    string _breakpoint_file_name;
    FILE * _breakpoint_file;

    int _nb_homo_clean;
    int _nb_homo_fuzzy;
    int _nb_hetero_clean;
    int _nb_hetero_fuzzy;


    // Actual job done by the tool is here
    void execute ();

    /** writes a given breakpoint in the output file
         */
    void writeBreakpoint(int bkt_id, string& chrom_name, uint64_t position, string& kmer_begin, string& kmer_end, int repeat_size);
    
private:
    
    /** fills getInfo() with parameters informations
     */
    void resumeParameters();

    /** fills getInfo() with results informations
         */
    void resumeResults();

    template<size_t span>
    void runFindBreakpoints();
};

/********************************************************************************/

#endif /* _TOOL_Finder_HPP_ */

