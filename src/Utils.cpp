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
// Remove insertions that correspond to reference sequences that contain only SNP
void remove_false_Het_insertion (std::vector<filled_insertion_t>& consensuses, string info_hetero,string sourceSequence,string targetSequence)
{
    std::vector<filled_insertion_t> final_set;
    int smallest_seq_pos = 0;
    int counter=0;
    vector<char> smallest_seq_vector;
    int smallest_seq_size;
    string smallest_seq;
    string current_seq;
    int current_size_seq;
     //cout << "SEQ " <<  endl;
     //Smallest sequence filled is supposed to be the one corresponding to the reference
    for (std::vector<filled_insertion_t>::iterator it_a=consensuses.begin(); it_a!=consensuses.end(); ++it_a) // could be improved : no need to compare first seq
    {

        vector<char> filled_f(it_a->seq.begin(), it_a->seq.end());
        //cout << "\n SEQ " << it_a->seq << endl;
        current_seq=it_a->seq;
        current_size_seq = filled_f.size();
        if (counter==0)
        {
           smallest_seq_vector=filled_f;
           smallest_seq_size=current_size_seq;
           smallest_seq=current_seq;
           smallest_seq_pos=counter;

        }
        if (counter!=0 && current_size_seq<smallest_seq_size)
        {
            smallest_seq_vector=filled_f;
            smallest_seq_size=current_size_seq;
            smallest_seq=current_seq;
            smallest_seq_pos=counter;
        }
        counter++;

    }
    int it_pos=0;
    //cout << "Smallest seq" << smallest_seq << endl;
    bool Check_push_smallest_seq=false;
    int size_second;
    //Remove sequence corresponding to the SNP + sequence between SNP and start of insertion
     for (std::vector<filled_insertion_t>::iterator it_b=consensuses.begin(); it_b!=consensuses.end(); ++it_b)
     {
		 int old_temp_size = consensuses.size();
		it_b->solution_count_het_snp=old_temp_size;
		 //cout << "\n compare " << (it_b->seq).compare(smallest_seq) << endl;
        if ((it_b->seq).compare(smallest_seq)!=0)
         {
         counter=0;
         int smallest_seq_size_2=smallest_seq_size-1;
         // TODO : Chan
		 //cout << info_hetero << endl;
         if (info_hetero=="down")
         {
             vector<char> filled_second(it_b->seq.begin(), it_b->seq.end());
             //string Seq_second=it_b->seq;
             size_second = filled_second.size()-1;
             //cout << it_b->seq << "    " << smallest_seq << endl;
			 //cout << smallest_seq_size_2 << size_second << endl;
			 //cout << "END " << smallest_seq_vector[smallest_seq_size_2] << filled_second[size_second] << endl;


             while (smallest_seq_size_2>0 && size_second>=0) // Check if the shortest sequence is found at one of the end. It is done 1 nt by 1 nt
             {
				 //cout << "END " << smallest_seq_vector[smallest_seq_size_2] << filled_second[size_second] << endl;
             if (smallest_seq_vector[smallest_seq_size_2]==filled_second[size_second])
             {
                 //cout << smallest_seq_vector[smallest_seq_size_2] << filled_second[size_second] << endl;
                 smallest_seq_size_2--;
                 size_second--;
                 counter++;
                 continue;
             }
             else
             {
                 it_b->altered_seq=0;
                 it_b->targetSeq=targetSequence;
                 it_b->sourceSeq=sourceSequence;
                 final_set.push_back(*it_b);
                 if (!Check_push_smallest_seq)
                 {
                     it_b->seq=smallest_seq;
                     final_set.push_back(*it_b);
                     Check_push_smallest_seq=true;
                 }
                 //cout << " ALTE : " << it_b->altered_seq << " SEQ : " << it_b->seq << endl;
                 break;
             }
            }
             if (smallest_seq_size_2!=0) continue;
             else
             {
			
             it_b->seq=it_b->seq.substr(0,filled_second.size()-smallest_seq_vector.size());
             it_b->sourceSeq=sourceSequence;
			 //cout << " SEQA0 : " << it_b->seq << endl;
             // Have to change insertion position
             it_b->altered_seq=counter;

             it_b->targetSeq=smallest_seq+targetSequence.substr(0,it_b->sourceSeq.size()-smallest_seq_vector.size());
             final_set.push_back(*it_b);


             }
          }

         //cout << "SMall 2 : " << smallest_seq << " Longest " << it_b->seq << endl;

         if (info_hetero=="up")
         {

             vector<char> filled_second(it_b->seq.begin(), it_b->seq.end());
             size_second = filled_second.size();
             //cout << it_b->seq << endl;
             while (smallest_seq_size_2>0 && size_second>=0)
             {
                 //cout << smallest_seq_size_2 << "    " << size_second << endl;
                 if (smallest_seq_vector[counter]==filled_second[counter])
                 {
                     smallest_seq_size_2--;
                     size_second--;
                     counter++;
                 }
                 else
                 {
                     it_b->altered_seq=0;
                     it_b->targetSeq=targetSequence;
                     it_b->sourceSeq=sourceSequence;
                     final_set.push_back(*it_b);
                     if (!Check_push_smallest_seq)
                     {
                         it_b->seq=smallest_seq;
                         final_set.push_back(*it_b);
                         Check_push_smallest_seq=true;
                     }
                     //cout << " ALTE : " << it_b->altered_seq << " SEQ : " << it_b->seq << endl;

                     break;
                 }
             }
             if (smallest_seq_size_2!=0) continue;

             else
             {
             //cout << "TEST " << it_b->seq <<  endl;
             //cout << "First " << size_second << "  " << smallest_seq_vector.size() << endl;
             it_b->seq=it_b->seq.substr(smallest_seq_vector.size());
             //cout << "SMALLEST " << smallest_seq.size() <<"  SEQ : " << smallest_seq << endl;
             if (smallest_seq_vector.size()>sourceSequence.size())
             {
                 it_b->sourceSeq=smallest_seq.substr(smallest_seq_vector.size()-sourceSequence.size());
				 it_b->altered_seq=0;
             	it_b->targetSeq=targetSequence;
             	final_set.push_back(*it_b);
                 //cout << " original seq" << smallest_seq << " new ref "<< it_b->sourceSeq << endl;
             }
             else
             {
                 it_b->sourceSeq=sourceSequence.substr(smallest_seq_vector.size())+smallest_seq;
				 it_b->altered_seq=0;
             	it_b->targetSeq=targetSequence;
             	final_set.push_back(*it_b);
             }
             //cout << " SEQ old " << sourceSequence << " CUt " << sourceSequence.substr(smallest_seq_vector.size(),sourceSequence.size()) << " NEW " << it_b->sourceSeq << " SEQ filled " << it_b->seq << endl;
             //it_b->altered_seq=counter;
             //cout << "Second" << endl;

             
             //cout << it_a->seq << endl;
             //cout << " ALTE : " << it_b->altered_seq << " SEQ : " << it_b->seq << endl;

             }
         }
        it_pos++;
     }
     }
     consensuses=final_set;
}
void remove_almost_identical_solutions(std::vector<filled_insertion_t>& consensuses, int identity_threshold,string sourceSequence,string targetSequence)
{
    // heuristic : add first seq to final set. Compare every seq to final_seq and add to final_set if different

    std::vector<filled_insertion_t> final_set;
    final_set.push_back(*consensuses.begin()  );
    //std::cerr << "remove_almost... after first push_back" << std::endl;

    for (std::vector<filled_insertion_t>::iterator it_a=consensuses.begin(); it_a!=consensuses.end(); ++it_a) // could be improved : no need to compare first seq
    {
		if (it_a->sourceSeq=="" && it_a->targetSeq=="")
		{
        it_a->targetSeq=targetSequence;
        it_a->sourceSeq=sourceSequence;
        it_a->altered_seq=0;
		}
        bool found_a_similar_seq = false;
        for (std::vector<filled_insertion_t>::iterator it_b=final_set.begin(); it_b!=final_set.end(); ++it_b){
            if (it_a->seq.compare(it_b->seq) == 0 || needleman_wunsch(it_a->seq,it_b->seq, NULL, NULL, NULL) * 100 >= identity_threshold){ // time optimisation ? if identical sequences, will not run needleman

                //This insertion is removed, but we select the one with nb_errors_in_anchor minimal
                if(it_a->nb_errors_in_anchor < it_b->nb_errors_in_anchor){
                    it_b->seq = it_a->seq;
                    it_b->nb_errors_in_anchor = it_a->nb_errors_in_anchor;
                    it_b->targetSeq=targetSequence;
                    it_b->sourceSeq=sourceSequence;
                    it_b->altered_seq=0;
                }
                found_a_similar_seq = true;
                break;
            }
        }
        // if the sequence has a %id always < threashold for all sequences in final_set, we add it to the final set.
        if(!found_a_similar_seq){
            final_set.push_back(*it_a);
        }
    }

    consensuses = final_set;
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
