/*****************************************************************************
*   MindTheGap: Integrated detection and assembly of insertion variants
*   A tool from the GATB (Genome Assembly Tool Box)
*   Copyright (C) 2022  INRIA
*   Authors: C. Lemaitre, G. Rizk, P. Marijon, W. Delage
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

#ifndef FINDSMALLINSERTION_HPP_
#define FINDSMALLINSERTION_HPP_
//**********************************
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>


template<size_t span>
class FindSmallCleanInsertion : public IFindObserver<span>
{
public :

    typedef typename gatb::core::kmer::impl::Kmer<span> Kmer;

    typedef typename Kmer::ModelCanonical KmerModel;
    typedef typename KmerModel::Iterator KmerIterator;

public:

    /** \copydoc IFindObserver<span>
     */
    /** \copydoc IFindObserver::IFindObserver
     */
    FindSmallCleanInsertion(FindBreakpoints<span> * find);


    /** \copydoc IFindObserver::IFindObserver
     */
    bool update();
};

template<size_t span>
FindSmallCleanInsertion<span>::FindSmallCleanInsertion(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindSmallCleanInsertion<span>::update()
{

    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
        return false;
    }

    if(this->_find->gap_stretch_size() == (this->_find->kmer_size()-1)) //Check size of gap
    {
        // obtains the kmer sequence
        string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
        string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());
        string ref = kmer_begin_str.substr(kmer_begin_str.size()-1,1);

        //All possible insertions of size 1 and 2
        char nucleo[20][6] = {"A","C","G","T","AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};

        KmerModel local_m(this->_find->kmer_size());
        KmerIterator local_it(local_m);
        std::string seq;
        string inser_base_one;
        bool found_base_one=false;
        
        //Test all possible insertions, by performing a micro-guided-assembly, ie. checks if all kmers of the insertion are present in the graph
        for (int i=0; i<20; i++)
        {
            seq = kmer_begin_str+ nucleo[i] + kmer_end_str;
            Data local_d(const_cast<char*>(seq.c_str()));
            int sum_valid=0;
            // Init this variable
            local_d.setRef(const_cast<char*>(seq.c_str()), (size_t)seq.length());
            local_it.setData(local_d);
            for(local_it.first(); !local_it.isDone(); local_it.next())
            {
                if(this->contains(local_it->forward()))
                {
                    sum_valid++;
                }
                else
                {
                    break;
                }
                if (sum_valid==this->_find->kmer_size())
                {
                    inser_base_one=ref+nucleo[i];
                    found_base_one=true;
                }
            }
            if (found_base_one==true) break;
        }
        if (!found_base_one) return false;
        
        this->_find->writeIndel(this->_find->breakpoint_id(),this->_find->chrom_name(),this->_find->position()-2, ref, inser_base_one, 0, STR_HOM_TYPE);
        this->_find->homo_clean_indel_iterate();
        this->_find->breakpoint_id_iterate();
        
        return true;
    }
    return false;
}

///*
template<size_t span>
class FindSmallFuzzyInsertion : public IFindObserver<span>
{
public :

    typedef typename gatb::core::kmer::impl::Kmer<span> Kmer;

    typedef typename Kmer::ModelCanonical KmerModel;
    typedef typename KmerModel::Iterator KmerIterator;

public:

    /** \copydoc IFindObserver<span>
     */
    /** \copydoc IFindObserver::IFindObserver
     */
    FindSmallFuzzyInsertion(FindBreakpoints<span> * find);


    /** \copydoc IFindObserver::IFindObserver
     */
    bool update();
};

template<size_t span>
FindSmallFuzzyInsertion<span>::FindSmallFuzzyInsertion(FindBreakpoints<span> * find):IFindObserver<span>(find){}

template<size_t span>
bool FindSmallFuzzyInsertion<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
        return false;
    }

    if(this->_find->gap_stretch_size() < this->_find->kmer_size() - 1 && this->_find->gap_stretch_size() >= this->_find->kmer_size() - 1 - this->_find->max_repeat())
    {
        int repeat_size = this->_find->kmer_size() - 1 - this->_find->gap_stretch_size();
        // obtains the kmer sequence
        string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
        string kmer_end_str = string(&(this->_find->chrom_seq()[this->_find->position() - 1 + repeat_size]), this->_find->kmer_size());
        if ((this->nb_out_branch(this->_find->kmer_begin().forward())==0) || (this->nb_in_branch(this->_find->kmer_end().forward())==0) || (!this->_find->model().codeSeed(&(this->_find->chrom_seq()[this->_find->position() - 1 + repeat_size]),Data::ASCII).isValid()))
        {
                   return false;
        }
        else
        {
            string ref = kmer_begin_str.substr(kmer_begin_str.size()-1-repeat_size,1);
            
            //All possible insertions of size 1 and 2
            char nucleo[20][6] = {"A","C","G","T","AA","AC","AG","AT","CA","CC","CG","CT","GA","GC","GG","GT","TA","TC","TG","TT"};
            KmerModel local_m(this->_find->kmer_size());
            KmerIterator local_it(local_m);
            std::string seq;
            string inser_base_one;
            bool found_base_one=false;
            //std::list<char> fourth (nucleo, nucleo + sizeof(nucleo) / sizeof(char) );
            for (int i=0; i<20; i++)
            {
                seq = kmer_begin_str+ nucleo[i] + kmer_end_str;
                //std::cout << seq << endl;
                Data local_d(const_cast<char*>(seq.c_str()));
                int sum_valid=0;
        //        // Init this variable
                local_d.setRef(const_cast<char*>(seq.c_str()), (size_t)seq.length());
                local_it.setData(local_d);
            for(local_it.first(); !local_it.isDone(); local_it.next())
               {
                    if(this->contains(local_it->forward()))
                    {
                    sum_valid++;
                    }
                    else
                    {
                        break;
                    }
                    if (sum_valid==this->_find->kmer_size())
                    {
                        inser_base_one=ref+nucleo[i];
                        found_base_one=true;
                    }
               }
            if (found_base_one==true) break;
            }
            if (!found_base_one) return false;
            this->_find->writeIndel(this->_find->breakpoint_id(),this->_find->chrom_name(),this->_find->position()- 2, ref, inser_base_one, repeat_size, STR_HOM_TYPE);
            this->_find->homo_clean_indel_iterate();
            this->_find->breakpoint_id_iterate();

            return true;
        }
        }
            return false;
    }


#endif // FINDSMALLINSERTION_HPP_

