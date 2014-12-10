/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G. Rizk
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
/********************************************************************************/

static const char* STR_URI_GRAPH = "-graph";
static const char* STR_URI_REF = "-ref";
static const char* STR_MAX_REPEAT = "-max-rep";
static const char* STR_ONLY_HOMO = "-hom-only";


class Finder : public Tool
{
public:

    // Constructor
    Finder ();
    
    size_t _kmerSize;
    Graph _graph;
    int _max_repeat;
    int _nbCores;

    // Actual job done by the tool is here
    void execute ();
    
    /** fills getInfo() with parameters informations
     */
    void resumeParameters();

};

/********************************************************************************/

#endif /* _TOOL_Finder_HPP_ */

