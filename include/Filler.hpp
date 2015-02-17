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

#ifndef _TOOL_Filler_HPP_
#define _TOOL_Filler_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
#include <GraphOutputDot.hpp>
using namespace std;

/********************************************************************************/

static const char* STR_URI_BKPT = "-bkpt";
static const char* STR_MAX_DEPTH = "-max-length";
static const char* STR_MAX_NODES = "-max-nodes";

class Filler : public Tool
{
public:

    // Constructor
    Filler ();
    
    size_t _kmerSize;
    Graph _graph;
    int _nbCores;
    BankFasta* _breakpointBank;
    int _max_depth;
    int _max_nodes;


    // Actual job done by the tool is here
    void execute ();

private:
    
    /** fills getInfo() with parameters informations
     */
    void resumeParameters();

    /** fills getInfo() with results informations
         */
    void resumeResults();

    /** Main function to fill the gaps of all breakpoints
         */
    template<size_t span>
    void fillBreakpoints();

    /** Fill one gap
            */
    template<size_t span>
    void gapFill(string sourceSequence, string targetSequence, set<string>& filledSequences);

    /** writes a given breakpoint in the output file
         */
    void writeFilledBreakpoint();

};

/********************************************************************************/

#endif /* _TOOL_Filler_HPP_ */

