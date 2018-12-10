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


static const char* STR_URI_CONTIG = "-contig";
static const char* STR_CONTIG_OVERLAP = "-overlap";
static const char* STR_URI_BKPT = "-bkpt";
static const char* STR_MAX_DEPTH = "-max-length";
static const char* STR_MAX_NODES = "-max-nodes";


 class info_node_t
{
public:
    int node_id;
    int pos;  // pos of beginning of right anchor
    int nb_errors;
    //bool anchor_is_repeated;
    bkpt_t targetId;

    //required to be inserted in std::set
    bool operator< (const info_node_t & other) const
    {
        return ( (node_id < other.node_id) || ((node_id == other.node_id)  && (pos < other.pos))  );
    }

    bool operator> (const info_node_t & other) const
    {
        return ( (node_id > other.node_id) || ((node_id == other.node_id)  && (pos > other.pos))  );
    }
    //required to be used in unordered_map
    bool operator==(const info_node_t & other) const
    {
        return(node_id == other.node_id
               && pos == other.pos
               && nb_errors == other.nb_errors
               && targetId == other.targetId);
    }
};

 struct nodeHasher
 {
   std::size_t operator()(const info_node_t& k) const
   {
     using std::size_t;
     using std::hash;
     using std::string;

     return ((std::hash<int>()(k.node_id)
              ^ (hash<int>()(k.pos) << 1)) >> 1)
              ^ (hash<int>()(k.nb_errors) << 1);
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

    //to print some statistics at the end
    int _nb_breakpoints;
    int _nb_filled_breakpoints;
    int _nb_multiple_fill;

    //useful to compute average abundance of each filled sequence
    Storage* _storage;

    //output file with filled sequences
    string _insert_file_name;
    FILE * _insert_file;

    //output file in GFA format (fot -contig)
    string _gfa_file_name;
    FILE * _gfa_file;

    //output file with statistics about each attempt of gap-filling
    string _insert_info_file_name;
    FILE * _insert_info_file;

    //parameters for dbg traversal (stop criteria)
    int _max_depth;
    int _max_nodes;

    //parameters for looking for the target sequence in the contig graph, with some mismatches and/or gaps
    int _nb_mis_allowed;
    int _nb_gap_allowed;

    int _overlap_length;

    string _vcf_file_name;
    FILE * _vcf_file;

    // Actual job done by the tool is here
    void execute ();



    //these two func moved to public because need access from functors breakpointFunctor and contigFunctor
    /** writes a given breakpoint in the output file
     */
    void writeFilledBreakpoint(std::vector<filled_insertion_t>& filledSequences, string breakpointName, std::string infostring, bool breakpointMode);
    void writeToGFA(std::vector<filled_insertion_t>& filledSequences, string sourceSequence, string SeedName, bool isRc);
    /** writes a given variant in the output vcf file
     */
    void writeVcf(std::vector<filled_insertion_t>& filledSequences, string breakpointName, string seedk);

    /** Fill one gap
     */
    /*template<size_t span>
    void gapFill(std::string & infostring,int tid,string sourceSequence, string targetSequence, set<filled_insertion_t>& filledSequences, bool begin_kmer_repeated, bool end_kmer_repeated
                 ,bool reversed =false);*/

    template<size_t span>
    void gapFillFromSource(std::string & infostring, int tid, string sourceSequence, string targetSequence, std::vector<filled_insertion_t>& filledSequences, bkpt_dict_t targetDictionary,bool is_anchor_repeated, bool reverse, string SeedName );

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

    template<size_t span>
    struct fillContig { void operator () (Filler* object); };
    template<size_t span>
    struct fillAny { void operator () (Filler* object); };

    /** writes the header of the vcf file
         */
    void writeVcfHeader();

    /**
     * returns the nodes containing the targetSequence (can be an approximate match)
     */
    set< info_node_t >  find_nodes_containing_R(string targetSequence, string linear_seqs_name, int nb_mis_allowed, int nb_gaps_allowed, bool anchor_is_repeated);
    set< info_node_t> find_nodes_containing_multiple_R(bkpt_dict_t targetDictionary, string linear_seqs_name, int nb_mis_allowed, int nb_gaps_allowed);

    /** Handle on the progress information. */
    void setProgress (gatb::core::tools::dp::IteratorListener* progress)  { SP_SETATTR(progress); }

    //tiens, on a pas encore de destructeur pour cette classe?
};

/********************************************************************************/

#endif /* _TOOL_Filler_HPP_ */

