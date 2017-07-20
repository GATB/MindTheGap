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

#ifndef _Utils_HPP_
#define _Utils_HPP_

#include <string>
#include <set>
#include <stdlib.h>
#include <vector>
#include <algorithm> 

using namespace std;

#include <unordered_map>
typedef pair<string, bool> bkpt_t;
typedef unordered_map<string, bkpt_t> bkpt_dict_t;

class filled_insertion_t
{
public:
	
    filled_insertion_t(string insert, int nb_errors, bool anchor_repeated_in_ref, bkpt_t targetId) : nb_errors_in_anchor(nb_errors),is_anchor_repeated(anchor_repeated_in_ref), targetId_anchor(targetId)	{
		seq = insert;
	}
    filled_insertion_t(string insert, int nb_errors, bool anchor_repeated_in_ref) : nb_errors_in_anchor(nb_errors),is_anchor_repeated(anchor_repeated_in_ref) {
        seq = insert;
    }
	
	string seq;
	int nb_errors_in_anchor;
	bool is_anchor_repeated;
	
	float avg_coverage;
	float median_coverage;
    bkpt_t targetId_anchor;

    //required to be inserted in set
    bool operator< (const filled_insertion_t & other) const
    {
        if (this->targetId_anchor != other.targetId_anchor)
            return this->targetId_anchor < other.targetId_anchor;
        else
            return this->seq < other.seq;
    }
	
	int compute_qual() const
	{
		
		if(is_anchor_repeated)
			return 25;
		
		if(nb_errors_in_anchor==2)
			return 5;
		
		if(nb_errors_in_anchor==1)
			return 10;
		
		return 50;
		
	}
};



// TODO factoriser
// this one is used in GraphAnalysis (modifies s)
void revcomp_sequence(char s[], int len);
// this one is used to reverse source and target Sequence (copies the sequence)
string revcomp_sequence(const string& dna);

/**
 * verifies if a and b are identical (tolerant to case), if one equals N returns false (even if both N)
 */
int identNT(char a, char b);

/**
 * gapped alignment
 * used by find_nodes_containingR : need to get the details of differences
 */
float needleman_wunsch(string a, string b, int * nbmatch,int * nbmis,int * nbgaps);




/**
 * returns true if all pairs of sequences have identity percent > threshold
 */
bool all_consensuses_almost_identical(set<filled_insertion_t> consensuses, int identity_threshold);

double median(std::vector<unsigned int> &v);

#endif /* _Utils_HPP_ */
