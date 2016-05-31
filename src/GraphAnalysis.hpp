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
#ifdef USE_NEW_CXX 
    #include <unordered_map>
    #define NS_TR1_PREFIX std
#else
    #include <tr1/unordered_map>
	#define NS_TR1_PREFIX std::tr1
#endif

#include <set>
#include <vector>
#include <fstream>
#include <stdlib.h>

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

    

    set<unlabeled_path> find_all_paths(set<int> terminal_nodes, bool &success);
    set<unlabeled_path> find_all_paths(int start_node, set<int> terminal_nodes, unlabeled_path current_path, int &nb_calls, bool &success);
    
    static int debug; // 0: no debug, 1: node id debug, 2: ful sequence debug; useful to see the sequences of the traversed paths
    set<string> paths_to_sequences(set<unlabeled_path> paths, set< std::pair<int,int> > terminal_nodes_with_endpos); // is it also used in mapsembler  ?
};
