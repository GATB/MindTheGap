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
	typedef typename gatb::core::kmer::impl::Kmer<span> Kmer;
	typedef typename Kmer::ModelCanonical KmerModel;
	typedef typename KmerModel::Iterator KmerIterator;
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
	if(!this->_find->homo_only())
	{
        // branching filter parameters
        int branching_threshold = this->_find->branching_threshold(); //max number of branching kmers in the 100 bp window of previous kmers
        int max_branching_kmers = branching_threshold;
        bool filtering = true;
        if (branching_threshold<0){
            filtering = false;
            max_branching_kmers = 100;
        }
        int filter_window_size = 100 ; //should not be larger than the size of het_kmer_history = 256
  
        
		// hetero site detection
		if(!this->_find->kmer_end_is_repeated() && this->_find->current_info().nb_in == 2 && !this->_find->recent_hetero())
		{
			//loop over putative repeat size (0=clean, >0 fuzzy), reports only the smallest repeat size found.
			for(int i = 0; i <= this->_find->max_repeat(); i++)
			{
				bool found_base_one = false;
				if(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 2 && !this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).is_repeated)
				{
					//hetero breakpoint found
					string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).kmer);
					//string kmer_end_str = this->_find->model().toString(this->_find->current_info().kmer);
                    //modif 15/06/2018 to check !!! (before in case of fuzzy>0, the end and right kmers overlapped, => insertion of wrong size (- fuzzy), missing the repeat + loss of recall if insertion of size < repeat)
                    string kmer_end_str = string(&(this->_find->chrom_seq()[this->_find->position() + i]), this->_find->kmer_size());
					string ref = kmer_begin_str.substr(kmer_begin_str.size() - 1 - i, 1);
                    
                    //Tests if this can be a small (1-2 bp) insertion
					char nucleo[20][6] = {"A", "C", "G", "T", "AA", "AC", "AG", "AT", "CA", "CC", "CG", "CT", "GA", "GC", "GG", "GT", "TA", "TC", "TG", "TT"};
					KmerModel local_m(this->_find->kmer_size());
					KmerIterator local_it(local_m);
					std::string seq;
					string inser_base_one;
					if (!this->_find->model().codeSeed(&(this->_find->chrom_seq()[this->_find->position() +i]),Data::ASCII).isValid())
                    {
                               return false;
                    }

					for (int a = 0; a < 20; a++) // for all possible 1-2 bp insertions, perform a micro-assembly
					{
						seq = kmer_begin_str + nucleo[a] + kmer_end_str;
						Data local_d(const_cast<char *>(seq.c_str()));
						int sum_valid = 0;
						//        // Init this variable
						local_d.setRef(const_cast<char *>(seq.c_str()), (size_t)seq.length());
						local_it.setData(local_d);
						for (local_it.first(); !local_it.isDone(); local_it.next())
						{
							if (this->contains(local_it->forward()))
							{
								sum_valid++;
							}
							else
							{
								break;
							}
							if (sum_valid == this->_find->kmer_size())
							{
								inser_base_one = ref + nucleo[a];
								found_base_one = true;
							}
						}
						if (found_base_one == true)
							break;
					}
					if (found_base_one)
					{
						this->_find->writeIndel(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1, ref, inser_base_one, i, STR_HET_TYPE);
						this->_find->hetero_indel_iterate();
						this->_find->breakpoint_id_iterate();
						return true;
					}
					else
					{
                        
                        //this may be a large insertion
                        
                         int nb_branching = 0;
                        //Applying the branching-filter :
                        if (filtering){
                            //counts the number of branching-kmers among the 100 previous ones
                            int nb_prev = 0;
                            unsigned char begin_index = this->_find->het_kmer_begin_index()-1;
                            while ((nb_branching <= max_branching_kmers) && (nb_prev<filter_window_size)){
                                //cout << "in loop" << nb_prev << "  " << begin_index-nb_prev << endl;
                                if(this->_find->het_kmer_history(begin_index-nb_prev).nb_out >1 || this->_find->het_kmer_history(begin_index-nb_prev).nb_in >1 ){
                                    nb_branching ++;
                                }
                                nb_prev++;
                            }
                        }
                        
                        if(nb_branching <= max_branching_kmers){
                            this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1 + i, kmer_begin_str, kmer_end_str, i, STR_HET_TYPE, this->_find->het_kmer_history(this->_find->het_kmer_begin_index() + i).is_repeated, this->_find->kmer_end_is_repeated());

                            this->_find->breakpoint_id_iterate();

                            if (i == 0)
                            {
                                this->_find->hetero_clean_iterate();
                            }
                            else
                            {
                                this->_find->hetero_fuzzy_iterate();
                            }

                            this->_find->recent_hetero(this->_find->max_repeat()); // we found a breakpoint, the next hetero one mus be at least _max_repeat apart from this one.
                            return true;										   //reports only the smallest repeat size found.
                        }
                        else{ // stop the loop over fuzzy size, because the branching context will remain not good for other fuzzy sizes
                            this->_find->recent_hetero(max(0, this->_find->recent_hetero() - 1)); // when recent_hetero=0 : we are sufficiently far from the previous hetero-site
                            return false;
                        }
                    }
				}
			}
		}

		this->_find->recent_hetero(max(0, this->_find->recent_hetero() - 1)); // when recent_hetero=0 : we are sufficiently far from the previous hetero-site
	}

	return false;
}

#endif /* _TOOL_FindHetero_HPP_ */
