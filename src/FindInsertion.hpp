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

#ifndef _TOOL_FindInsert_HPP_
#define _TOOL_FindInsert_HPP_

/*******************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>

template<size_t span>
class FindCleanInsertion : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindCleanInsertion(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindCleanInsertion<span>::FindCleanInsertion(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindCleanInsertion<span>::update()
{

//cout<<" \n" << i << this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index() +i).kmer)<< endl;

    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
        return false;
    }

    if(this->_find->gap_stretch_size() == (this->_find->kmer_size()-1)) //Check size of gap
    {
        // obtains the kmer sequence
        string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
        string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());
        // Check that kmer_begin has neighbor, if not the breakpoint is not valid
        for(int i = -1; i <= this->_find->max_repeat()+1; i++)
        {
        string kmer_begin_str_1 = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).kmer); //this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 2 &&
        if (!kmer_begin_str_1.compare(kmer_begin_str))
        {
            //cout << kmer_begin_str << "kmer begin str" << kmer_begin_str_1 << "kmer searched in history" << endl;
            //cout << "branchig out" << this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out << endl;
            if (this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 0)
            {
               // cout << " bkpt  nb_out = 0" << endl;
                return false;
            }
        }
    }
    /*if (check==false)
    {
        cout << "error not found kmer begin in history" << endl;
    }
    else
    {
        cout << "Good" << endl;
    }*/


        //position : this->_find->position() is the beginning of the second found kmer after the gap : -2 ie position of the last 0, ie position just before (at the left of) the insertion site (0-based)
        this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 2, kmer_begin_str, kmer_end_str, 0,STR_HOM_TYPE,"none",  this->_find->kmer_begin_is_repeated() ,this->_find->kmer_end_is_repeated()  );

        // iterate counter
        this->_find->breakpoint_id_iterate();
        this->_find->homo_clean_iterate();
        return true;
    }
    /*if(this->_find->gap_stretch_size() > (this->_find->kmer_size()-1)) //Check size of gap
    {
        // obtains the kmer sequence
        string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
        string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());

        //position : this->_find->position() is the beginning of the second found kmer after the gap : -2 ie position of the last 0, ie position just before (at the left of) the insertion site (0-based)
        this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 2, kmer_begin_str, kmer_end_str, 0,STR_HOM_TYPE,"none",  this->_find->kmer_begin_is_repeated() ,this->_find->kmer_end_is_repeated()  );

        // iterate counter
        this->_find->breakpoint_id_iterate();
        this->_find->homo_clean_iterate();
        return true;
    }*/

    return false;
}

template<size_t span>
class FindFuzzyInsertion : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindobserver
     */
    FindFuzzyInsertion(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindFuzzyInsertion<span>::FindFuzzyInsertion(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindFuzzyInsertion<span>::update()
{

    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
        return false;
    }

    if(this->_find->gap_stretch_size() < this->_find->kmer_size() - 1 && this->_find->gap_stretch_size() >= this->_find->kmer_size() - 1 - this->_find->max_repeat())
    {
        // Fuzzy site, position and kmer_end are impacted by the repeat
        int repeat_size = this->_find->kmer_size() - 1 - this->_find->gap_stretch_size();

        // obtains the kmer sequence
        string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
        string kmer_end_str = string(&(this->_find->chrom_seq()[this->_find->position() - 1 + repeat_size]), this->_find->kmer_size());

        // Check that kmer_begin has neighbor, if not the breakpoint is not valid
        for(int i = -1; i <= this->_find->max_repeat()+1; i++)
        {
    //cout<<" \n" << i << this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index() +i).kmer)<< endl;
            string kmer_begin_str_1 = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).kmer); //this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 2 &&
            if (!kmer_begin_str_1.compare(kmer_begin_str))
            {
                //cout << kmer_begin_str << "kmer begin str" << kmer_begin_str_1 << "kmer searched in history" << endl;
                //cout << "branchig out" << this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out << endl;
                if (this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 0)
                {
                    //cout << " bkpt  nb_out = 0" << endl;
                    return false;
                }

            }
        }

        /*if (check==false)
        {
            cout << "error not found kmer begin in history" << endl;
        }
        else
        {
            cout << "Good" << endl;
        }
        */
        string found_snp="none";
        //position : this->_find->position() is the beginning of the second found kmer after the gap : -2 ie position of the last 0, ie position just before (at the left of) the insertion site (0-based)
        this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 2 + repeat_size, kmer_begin_str, kmer_end_str, repeat_size, STR_HOM_TYPE,found_snp,   this->_find->kmer_begin_is_repeated() , this->_find->kmer_end_is_repeated());

        //iterate counter
        this->_find->breakpoint_id_iterate();
        this->_find->homo_fuzzy_iterate();

        return true;
    }

    return false;
}

#endif /* _TOOL_FindInsert_HPP_ */
