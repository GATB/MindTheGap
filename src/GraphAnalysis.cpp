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
#include <GraphAnalysis.hpp>
#include <Utils.hpp>

#include <algorithm> // for find()

int GraphAnalysis::debug = 0;

/*
 * the graph (produced by GraphOutput) is loaded as a directed graph where
 * each node and its revcomp are separated. For instance,
 *        node 0
 * becomes:
 *        node "0f" and node "0r"
 * enabling edges between, for instance
 *        0 -> 1 [label="fr"]
 * becomes:
 *        0f -> 1r
 *
 * internal representation of "0f" is 0, and "0r" is 0+nb_nodes
 *
 * the following functions performs the conversion
 */

string GraphAnalysis::node_identifier(int node)
{
    char node_id[10];
    sprintf(node_id,"%d%s",(node>nb_nodes)?(node-nb_nodes):node, (node<nb_nodes)?"f":"r");
    string res (node_id);
    return res;
}

int GraphAnalysis::revcomp_node(int node)
{
    return (node<nb_nodes)?(node+nb_nodes):(node-nb_nodes);
}

GraphAnalysis::GraphAnalysis(string graph_file_name,size_t kmerSize)
{
    _sizeKmer =kmerSize;
    ifstream graph_file (graph_file_name.c_str());
    string line;

    if (!graph_file.is_open() || graph_file.eof())
    {
        printf("error opening file %s\n",graph_file_name.c_str());
        exit(1);
    }

    bool parsing_nodes = true;
    nb_nodes = 0;
    nb_edges = 0;
    char *node_sequence = (char*)malloc(1000000);

    getline(graph_file, line);

    while (graph_file.good())
    {
        getline(graph_file, line);

        // hypothesis: nodes are numbered in GAP-LESS, STRICTLY INCREMENTAL order, i.e. 0,1,2,3...

        int node_a, node_b, nb_numbers_seen;
        nb_numbers_seen = sscanf(line.c_str(), "%d%*s%d",&node_a,&node_b); // sketchy but works on my system

        if (nb_numbers_seen == 1)
        {
            sscanf(line.c_str(), "%*d %*[^\"]%*[\"]%[A-Z]%*[\"]",node_sequence); // ugly regexp to get the node sequence
            string str_node_sequence(node_sequence);
            node_sequences[nb_nodes] = node_sequence;
            nb_nodes++;
        }
        if (nb_numbers_seen == 2)
        {
            char label[100]; //needs to be large enough for the regexp below
            sscanf(line.c_str(), "%*d %*s %*d %*[^\"]%*[\"]%s%*[\"]",label); // ugly regexp to get the label of the edge
            label[2]='\0';

            if (label[0] == 'R')
                node_a = revcomp_node(node_a);
            if (label[1] == 'R')
                node_b = revcomp_node(node_b);

            if (out_edges[node_a].find(node_b) == out_edges[node_a].end())
            {
                out_edges[node_a].insert(node_b);
                nb_edges++;
            }
        }
    }

    free(node_sequence);
    //printf("finished parsing %s: %d nodes and %d edges found\n", graph_file_name.c_str(), nb_nodes, nb_edges);
}

// wrapper
set<unlabeled_path> GraphAnalysis::find_all_paths(set<int> terminal_nodes, bool &success)
{
    success = true;
    unlabeled_path start_path;
    start_path.push_back(0);
    int nb_calls = 0;
    set<unlabeled_path> paths = find_all_paths(0, terminal_nodes, start_path, nb_calls, success);
    //std::cout << "PATHS0 \n"  << endl;

    return paths;
}

// precondition: terminal_nodes is non-empty
set<unlabeled_path> GraphAnalysis::find_all_paths(int start_node, set<int> terminal_nodes, unlabeled_path current_path, int &nb_calls, bool &success)
{
    set<unlabeled_path> paths;


    // don't explore for too long
    if (nb_calls++ > 10000000)
    {
       // printf("fail, max nb_calls reached \n");

        success = false;

        return paths;
    }

    if (terminal_nodes.find(start_node) != terminal_nodes.end()) //stops when reaches one of the terminal nodes.
    {
        paths.insert(current_path);
        return paths;

    }
    // visit all neighbors
    for(set<int>::iterator it_edge = out_edges[start_node].begin(); it_edge != out_edges[start_node].end(); it_edge++)
    {
        int next_node = *it_edge;

        // each node of the graph is used in a path at most once (as a consequence, we won't gapfill some tandem repeats)
        if (find(current_path.begin(), current_path.end(), next_node) == current_path.end())
        {
            unlabeled_path extended_path(current_path);
            extended_path.push_back(next_node);

            // some debug
            /*printf("calling with %d, curent path is:",next_node);
            for(vector<int>::iterator it_path = current_path.begin(); it_path != current_path.end(); it_path++)
                 printf(" %d",*it_path);
            printf("\n");*/


            // recursive call
            set<unlabeled_path> new_paths = find_all_paths(next_node, terminal_nodes, extended_path, nb_calls, success);
            paths.insert(new_paths.begin(), new_paths.end());


            // mark to stop we end up with too large breadth
            if (paths.size() >= max_breadth)
            {
                //printf("fail, max breadth reached \n");
                success = false;
            }
        }

        // propagate the stop if too many consensuses reached
        if (success == false)
            return paths;
    }
    return paths;
}

set<filled_insertion_t> GraphAnalysis::paths_to_sequences(set<unlabeled_path> paths , set< info_node_t > terminal_nodes_with_endpos )
{
    //debug =2;
    set<filled_insertion_t> sequences;
    int errs_in_anchor;
    bool anchor_repeated_in_ref;
    bkpt_t targetId_anchor;

    //printf("paths set size %i \n",paths.size());
    for (set<unlabeled_path>::iterator it = paths.begin(); it != paths.end(); it++)
    {
//        if (debug)
//            printf("processing path: \n");

        unlabeled_path p = *it;
        string sequence;
        //cout << "new path, length: " << p.size() << endl;
        for (unlabeled_path::iterator it_path = p.begin(); it_path != p.end(); it_path++)
        {
            int node = *it_path;
            int pos_anchor = 0;


//            if (debug)
//            {
//                string node = node_identifier(*it_path);
//                printf("\n %s \n",node.c_str());
//            }

            bool revcomp = node > nb_nodes;
            if (revcomp)
                node -= nb_nodes;

            string node_sequence = node_sequences[node];

            if (revcomp)
            {
                char node_str[node_sequence.length()+1];
                sprintf(node_str,"%s",node_sequence.c_str());
                revcomp_sequence(node_str,node_sequence.length());
                node_sequence = node_str;
            }
            if (debug > 1)
                printf("%s ",node_sequence.c_str());




            //trim everything after right anchor for last node (if at beginning, must also erase seq from previous overlap)
            if (it_path == (p.end() -1 ))
            {
                set<int> terminal_nodes;
                //retreive anchor pos
                for (set< info_node_t >::iterator it = terminal_nodes_with_endpos.begin(); it != terminal_nodes_with_endpos.end(); it++)
                {
                    if((*it).node_id == node)
                    {
                        pos_anchor = (*it).pos;
                        errs_in_anchor = it->nb_errors;
                        anchor_repeated_in_ref= it->anchor_is_repeated;
                        targetId_anchor = it ->targetId;
                        break;
                    }
                }
                //sprintf("last node of the path, pos anchor %i \n",pos_anchor);
                node_sequence = node_sequence.substr(0,pos_anchor);
                //cout << endl << node_sequence << " "<<node << endl;

                if(pos_anchor <= (_sizeKmer-1))
                {
                    sequence = sequence.substr(0, sequence.length() - ((_sizeKmer-1) - pos_anchor)); //nothing else to add
                }
                else
                {
                    //must still erase k-1 overlap

                    if (it_path != p.begin())
                    {
                        node_sequence = node_sequence.substr(_sizeKmer-1,node_sequence.npos);
                    }
                    else
                    {
                        node_sequence = node_sequence.substr(_sizeKmer,node_sequence.npos);
                    }

                    sequence += node_sequence;
                }
                break;

            }


            // trim (k-1) overlaps at beginning, except for the first node
            if (it_path != p.begin())
            {
                node_sequence = node_sequence.substr(_sizeKmer-1,node_sequence.npos); // check if  seq is now larger than sizeKmer-1
            }
            else // for first node, remove sizeKmer, the left anchor
            {

                node_sequence = node_sequence.substr(_sizeKmer,node_sequence.npos);

            }


            if (debug > 1)
                printf("(%s) ",node_sequence.c_str());

            // append to sequence
            sequence += node_sequence;
        }
        if (debug > 1)
            printf(" sequence: %s\n",sequence.c_str());

        if(sequence.length()>0) // test filtrage ici ?
           { sequences.insert(filled_insertion_t(sequence,errs_in_anchor,anchor_repeated_in_ref,targetId_anchor));
            //cout << targetId_anchor.first << endl;
            //cout  << sequence << endl;

           }
    }

    return sequences;
}
