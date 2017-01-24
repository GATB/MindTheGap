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
#include <Utils.hpp>

using namespace std;

/********************************************************************************/

static const char* STR_URI_BKPT = "-bkpt";
static const char* STR_MAX_DEPTH = "-max-length";
static const char* STR_MAX_NODES = "-max-nodes";


 class info_node_t
{
public:
	int node_id;
	int pos;  // pos of beginning of right anchor
	int nb_errors;
	bool anchor_is_repeated;
	
	
	//required to be inserted in std::set
	bool operator< (const info_node_t & other) const
	{
		return ( (node_id < other.node_id) || ((node_id == other.node_id)  && (pos < other.pos))  );
	}
	
	bool operator> (const info_node_t & other) const
	{
		return ( (node_id > other.node_id) || ((node_id == other.node_id)  && (pos > other.pos))  );
	}
};





class Filler : public Tool
{

public:

    // Constructor
    Filler ();
	void FillerHelp();

	const char* _mtg_version;

	
    size_t _kmerSize;
    Graph _graph;
    int _nbCores;
    BankFasta* _breakpointBank;
    int _nb_breakpoints;
	int _nb_filled_breakpoints;
	int _nb_multiple_fill;

	
	string _insert_file_name;
	FILE * _insert_file;
	
	string _insert_info_file_name;
	FILE * _insert_info_file;

    int _max_depth;
    int _max_nodes;

    int _nb_mis_allowed;
    int _nb_gap_allowed;


    // Actual job done by the tool is here
    void execute ();


	
	//these two func moved to public because need access from functor gapfillerFunctor
	/** writes a given breakpoint in the output file
	 */
	void writeFilledBreakpoint(set<filled_insertion_t>& filledSequences, string breakpointName);

	/** Fill one gap
	 */
	template<size_t span>
	void gapFill(int tid,string sourceSequence, string targetSequence, set<filled_insertion_t>& filledSequences, bool begin_kmer_repeated, bool end_kmer_repeated
				 ,bool reversed =false);

	gatb::core::tools::dp::IteratorListener* _progress;

	
private:
	
    /** fills getInfo() with parameters informations
     */
    void resumeParameters();

    /** fills getInfo() with results informations
         */
    void resumeResults(double seconds);

    /** Main function to fill the gaps of all breakpoints
         */
    template<size_t span>
    struct fillBreakpoints {  void operator ()  (Filler* object); };



    /**
     * returns the nodes containing the targetSequence
     */
    set< info_node_t >  find_nodes_containing_R(string targetSequence, string linear_seqs_name, int nb_mis_allowed, int nb_gaps_allowed, bool anchor_is_repeated);

	/** Handle on the progress information. */
	void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }
	
	//tiens, on a pas encore de destructeur pour cette classe?
};

/********************************************************************************/

#endif /* _TOOL_Filler_HPP_ */

