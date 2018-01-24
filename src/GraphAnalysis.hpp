/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G.Rizk, R. Chikhi
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
#include <string>


//USE_NEW_CXX variable defined in CMakeList.txt of gatb-core : depending on the compil version unordered_map is not in the same location...
#include <unordered_map>
#define NS_TR1_PREFIX std


#include <set>
#include <vector>
#include <fstream>
#include <stdlib.h>
#include <Filler.hpp>

using namespace std;

// path type
typedef vector<int> unlabeled_path;


class GraphAnalysis {

public:
    static const int max_breadth = 20; //changed from 10 to 20

    string prefix;
    FILE *graph_file;

    int nb_nodes, nb_edges;

    string node_identifier(int node);
    int revcomp_node(int node);

    NS_TR1_PREFIX::unordered_map<int,string > node_sequences;
    NS_TR1_PREFIX::unordered_map<int,set<int> > out_edges;

	size_t _sizeKmer;


    GraphAnalysis(string graph_file_name,size_t kmerSize);

    

    set<pair<unlabeled_path,bkpt_t>> find_all_paths(set<info_node_t> terminal_nodes_with_endpos, bool &success);
    set<pair<unlabeled_path,bkpt_t>> find_all_paths(int start_node, set<info_node_t> terminal_nodes_with_endpos, unlabeled_path current_path, int &nb_calls, bool &success);
    
    static int debug; // 0: no debug, 1: node id debug, 2: ful sequence debug; useful to see the sequences of the traversed paths
    set<filled_insertion_t> paths_to_sequences(set<unlabeled_path> paths, set< info_node_t > terminal_nodes_with_endpos); // is it also used in mapsembler  ?
};

