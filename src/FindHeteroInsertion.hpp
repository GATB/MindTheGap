/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G.Rizk, P.Marijon
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

#ifndef _TOOL_FindHetero_HPP_
#define _TOOL_FindHetero_HPP_

/*******************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>
template<size_t span>
class FindHeteroInsertion : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindHeteroInsertion(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindHeteroInsertion<span>::FindHeteroInsertion(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindHeteroInsertion<span>::update()
{
    if(!this->_find->homo_only() )
    {
        //cout << this->_find->model().toString(this->_find->current_info().kmer) << endl;
        //if (this->_find->current_info().nb_out == 2 ) { cout << " NB_OUT = 2 :" <<  this->_find->model().toString(this->_find->current_info().kmer) << "POSITION : " << this->_find->position() << endl;}
        //if (this->_find->current_info().nb_in == 2 ) { cout <<  " NB_in = 2 :" <<  this->_find->model().toString(this->_find->current_info().kmer) << "POSITION : " << this->_find->position() << endl; }
        //cout << "Gap size" << this->_find->gap_stretch_size() <<"snp_near " << this->_find->recent_snp() <<  endl;
        // hetero site detection
        if(!this->_find->kmer_end_is_repeated() && this->_find->current_info().nb_in == 2 && !this->_find->recent_hetero() && this->_find->gap_stretch_size()==0)
        {
            //std::cout<<"\n Heterozygote suspect " << this->_find->position()-1  << endl;
            //std::cout << " gap stretch size " <<  this->_find->gap_stretch_size() << endl;
            //loop over putative repeat size (0=clean, >0 fuzzy), reports only the smallest repeat size found.

            for(int i = 0; i <= this->_find->max_repeat(); i++)
            {
        //cout<<" \n" << i << this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index() +i).kmer)<< endl;
                if(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 2 && !this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).is_repeated) //this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 2 &&
                {
                    //std::cout<<"\n Heterozygote suspect " << this->_find->position()-1  << endl;
                    //hetero breakpoint found
                    string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).kmer);
                    //string kmer_end_str = this->_find->model().toString(this->_find->current_info().kmer);
                    //modif 15/06/2018 to check !!! (before in case of fuzzy>0, the end and right kmers overlapped, => insertion of wrong size (- fuzzy), missing the repeat + loss of recall if insertion of size < repeat)
                    string kmer_end_str = string(&(this->_find->chrom_seq()[this->_find->position() + i]), this->_find->kmer_size());
                    string found_snp="none";

                    this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position()-1+i, kmer_begin_str, kmer_end_str,i, STR_HET_TYPE,found_snp,  this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).is_repeated,this->_find->kmer_end_is_repeated() );

                    this->_find->breakpoint_id_iterate();

                    if(i==0)
                    {
                        this->_find->hetero_clean_iterate();
                    }
                    else
                    {
                        this->_find->hetero_fuzzy_iterate();
                    }

                    this->_find->recent_hetero(this->_find->max_repeat()); // we found a breakpoint, the next hetero one mus be at least _max_repeat apart from this one.
                    return true; //reports only the smallest repeat size found.
                }
            }
        }
            // Find heterozygote insertion with SNP upstream the insertion (<k-1 distance)
           /*if(!this->_find->kmer_end_is_repeated() && this->_find->current_info().nb_in == 2 && !this->_find->recent_hetero() )
            {

               if(this->_find->recent_snp()>0 || this->_find->gap_stretch_size()>0 )
               {
                    for(size_t i =0; i <= (this->_find->kmer_size())*2; i++)
                    {


                       //cout << "\n I " << i << endl;
                       //cout << "\n " << " POSITION " <<  i << "\n SEQU : " << this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).kmer) <<" index " << this->_find->het_kmer_begin_index()+i << endl;
                       //cout << "\n " << " POSITION " <<  i << "\n SEQU : " << this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).kmer)<<" index " << this->_find->het_kmer_begin_index()-i << endl;
                       if(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).nb_out==1 && !this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).is_repeated)
                       {
                           for(int i = 0; i <= this->_find->max_repeat(); i++)
                           {}
                           string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i-1).kmer);
                           string kmer_end_str = this->_find->model().toString(this->_find->current_info().kmer);
                           //string kmer_end_str =string(&(this->_find->chrom_seq()[this->_find->position() - i]), this->_find->kmer_size());
                           string found_snp="up";

                           this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position()-1, kmer_begin_str, kmer_end_str,i, STR_HET_TYPE,found_snp,  this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).is_repeated,this->_find->kmer_end_is_repeated() );

                           this->_find->breakpoint_id_iterate();
                           // TODO FUZZY INSERTION
                           this->_find->hetero_clean_iterate();
                           this->_find->recent_hetero(this->_find->max_repeat());
                           this->_find->recent_snp(0);
                           return true;
                       }
                    }
                }
            }*/


        //Find heterozygote insertion with SNP downstream the insertion (<k-1 distance) TODO NOT DONE
        /*if(!this->_find->kmer_end_is_repeated() && this->_find->current_info().nb_in == 1 && !this->_find->recent_hetero() )
        {
           if(this->_find->recent_snp()>0 || this->_find->gap_stretch_size()>0 )
           {
           // cout << "\n POSITION INIT : " << this->_find->position() <<endl;
           // cout << "\n seq :" << this->_find->model().toString(this->_find->current_info().kmer) << endl;
           // cout << " \n GAP : " << this->_find->gap_stretch_size() << endl;
               // cout << "Gap size" << this->_find->gap_stretch_size() <<"snp_near " << this->_find->recent_snp() << "recent hetero" << this->_find->recent_hetero() <<  endl;
                for(size_t i =0; i <= (this->_find->kmer_size())*2; i++)
                {
                    //cout << "\n I " << i << endl;
                    //cout << "\n " << " POSITION " <<  this->_find->position()-1-i-this->_find->gap_stretch_size() << "\n SEQU : " <<this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i-+this->_find->gap_stretch_size()).kmer) << endl;
                   if(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).nb_out==2 && !this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).is_repeated)
                   {
                       //cout << "\n I " << i << endl;
                       //cout << "\n " << " POSITION " <<  this->_find->position()-1-i << "\n SEQU : " <<this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).kmer) << "\n nb out" << this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).nb_out << endl;
                       //cout << " \n GAP : " << this->_find->gap_stretch_size() << endl;
                       string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()-i).kmer);
                       string kmer_end_str = this->_find->model().toString(this->_find->current_info().kmer);
                       //string kmer_end_str =string(&(this->_find->chrom_seq()[this->_find->position() - i]), this->_find->kmer_size());
                       string found_snp="down";
                       //if ((kmer_end_str.find("N")) | (kmer_end_str.find("n"))) return false;

                       this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position()-1-i, kmer_begin_str, kmer_end_str,i, STR_HET_TYPE,found_snp,  this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).is_repeated,this->_find->kmer_end_is_repeated() );

                       this->_find->breakpoint_id_iterate();
                       // TODO FUZZY INSERTION
                       this->_find->hetero_clean_iterate();
                       this->_find->recent_hetero(this->_find->max_repeat());
                       this->_find->recent_snp(0);
                       return true;
                   }
                }
            }
        }*/


        //cout << "\n OK" <<  this->_find->m_recent_hetero() << endl;

        this->_find->recent_hetero(max(0, this->_find->recent_hetero() - 1));  // when recent_hetero=0 : we are sufficiently far from the previous hetero-site

    }
    return false;
}

#endif /* _TOOL_FindHetero_HPP_ */
