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

#include <Utils.hpp>

void revcomp_sequence(char s[], int len)
{
#define CHAR_REVCOMP(a,b) {switch(a){\
	case 'A': b='T';break;case 'C': b='G';break;case 'G': b='C';break;case 'T': b='A';break;default: b=a;break;}}
		  int i;
		  unsigned char t;
		  for (i=0;i<len/2;i++)
		  {
			  t=s[i];
			  CHAR_REVCOMP(s[len-i-1],s[i]);
			  CHAR_REVCOMP(t,s[len-i-1]);
		  }
		  if (len%2==1)
			  CHAR_REVCOMP(s[len/2],s[len/2]);

}


string revcomp_sequence(const string& dna) {

    //cout<<"dna="<<dna<<endl;
    string revComp= "";

    for (string::const_reverse_iterator it = dna.rbegin(); it != dna.rend(); it++) {
        switch(*it) { //each element in the temp vector to it's complement element.
            case 'a' :
                revComp += "t";
                break;
            case 't' :
                revComp += "a";
                break;
            case 'c' :
                revComp += "g";
                break;
            case 'g' :
                revComp += "c";
                break;
            case 'A' :
                revComp += "T";
                break;
            case 'T' :
                revComp += "A";
                break;
            case 'C' :
                revComp += "G";
                break;
            case 'G' :
                revComp += "C";
                break;
        }
    }

    //cout<<"revComp="<<revComp<<endl;
    return revComp;
}



int identNT(char a, char b)
{
    return   ( (a==b || a-b==32 || a-b ==-32) && a!='N');
}


float needleman_wunsch(string a, string b, int * nbmatch,int * nbmis,int * nbgaps)
{
    float gap_score = -5;
    float mismatch_score = -5;
    float match_score = 10;
#define nw_score(x,y) ( (x == y) ? match_score : mismatch_score )

    int n_a = a.length(), n_b = b.length();
    float ** score =  (float **) malloc (sizeof(float*) * (n_a+1));
    for (int ii=0; ii<(n_a+1); ii++)
    {
        score [ii] = (float *) malloc (sizeof(float) * (n_b+1));
    }

	// float score[n_a+1][n_b+1];  //stack is too small
    // float pointer[n_a+1][n_b+1];

    for (int i = 0; i <= n_a; i++)
        score[i][0] = gap_score * i;
    for (int j = 0; j <= n_b; j++)
        score[0][j] = gap_score * j;

    // compute dp
    for (int i = 1; i <= n_a; i++)
    {
        for (int j = 1; j <= n_b; j++)
        {
            float match = score[i - 1][j - 1] + nw_score(a[i-1],b[j-1]);
            float del =  score[i - 1][j] + gap_score;
            float insert = score[i][j - 1] + gap_score;
            score[i][j] = max( max(match, del), insert);
        }
    }

    // traceback
    int i=n_a, j=n_b;
    float identity = 0;
	int nb_mis = 0;
	int nb_gaps = 0;

	int nb_end_gaps = 0 ;
	bool end_gap = true;

    while (i > 0 && j > 0)
    {

        float score_current = score[i][j], score_diagonal = score[i-1][j-1], score_up = score[i][j-1], score_left = score[i-1][j];
        if (score_current == score_diagonal + nw_score(a[i-1], b[j-1]))
        {
            if (a[i-1]== b[j-1])
			{

                identity++;
			}
			else
			{
				nb_mis++;
			}
            i -= 1;
            j -= 1;


			end_gap = false;
        }
        else
        {
            if (score_current == score_left + gap_score)
			{
				i -= 1;
			}
            else if (score_current == score_up + gap_score)
			{
				j -= 1;
			}


			if(!end_gap) //pour ne pas compter gap terminal
			{

				nb_gaps++;
			}
        }
    }

	//pour compter gaps au debut  :
	nb_gaps += i+j;
	//if(deb==1)printf("add gaps i j %i %i \n",i,j);

    identity /= max( n_a, n_b); // modif GR 27/09/2013    max of two sizes, otherwise free gaps

    if(nbmatch!=NULL) *nbmatch = identity;
    if(nbmis!=NULL)  *nbmis = nb_mis;
    if(nbgaps!=NULL) *nbgaps = nb_gaps;

    for (int ii=0; ii<(n_a+1); ii++)
    {
        free (score [ii]);
    }
    free(score);

    //printf("---nw----\n%s\n%s -> %.2f\n--------\n",a.c_str(),b.c_str(),identity);
    return identity;
}

bool all_consensuses_almost_identical(set<filled_insertion_t> consensuses, int identity_threshold)
{
	for (set<filled_insertion_t>::iterator it_a = consensuses.begin(); it_a != consensuses.end(); it_a++)
    {
		set<filled_insertion_t>::iterator it_b = it_a;
        advance(it_b,1);
        while (it_b != consensuses.end())
        {
            if (needleman_wunsch(it_a->seq,it_b->seq, NULL, NULL, NULL) * 100 < identity_threshold)
                return false;
            advance(it_b,1);
        }
    }
    return true;
}


double median(std::vector<unsigned int> &v)
{
	size_t n = v.size() / 2;
	std::nth_element(v.begin(), v.begin()+n, v.end());
	unsigned int vn = v[n];
	if(v.size()%2 == 1)
	{
		return vn;
	}else
	{
		std::nth_element(v.begin(), v.begin()+n-1, v.end());
		return 0.5*(vn+v[n-1]);
	}
}
