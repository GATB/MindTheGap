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

using namespace std;

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
bool all_consensuses_almost_identical(set<string> consensuses, int identity_threshold);

#endif /* _Utils_HPP_ */
