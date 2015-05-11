/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G.Rizk
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

#ifndef _TOOL_FindObserver_HPP_
#define _TOOL_FindObserver_HPP_

/*******************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>

template<size_t span>
class FindCleanInsert : public IFindObserver<span>
{
public :
    
    /** \copydoc IFindObserver::IFindObserver
     */
    FindCleanInsert(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindCleanInsert<span>::FindCleanInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindCleanInsert<span>::update()
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

	this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1, kmer_begin_str, kmer_end_str, 0,STR_HOM_TYPE);

	// iterate counter
	this->_find->breakpoint_id_iterate();
	this->_find->homo_clean_iterate();
	return true;
    }
    
    return false;
}

template<size_t span>
class FindFuzzyInsert : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindobserver
     */
    FindFuzzyInsert(FindBreakpoints<span> * find);
    
    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindFuzzyInsert<span>::FindFuzzyInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindFuzzyInsert<span>::update()
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
	    
	this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1 + repeat_size, kmer_begin_str, kmer_end_str, repeat_size, STR_HOM_TYPE);

	//iterate counter
	this->_find->breakpoint_id_iterate();
	this->_find->homo_fuzzy_iterate();

	return true;
    }

    return false;
}

template<size_t span>
class FindSNP : public IFindObserver<span>
{
public :
    
    typedef typename gatb::core::kmer::impl::Kmer<span>::Type KmerType;

public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindSNP(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    virtual bool update() = 0;

protected :

    bool contains(KmerType kmer);
    KmerType correct(KmerType& kmer, KmerType& nuc, size_t pos);
    KmerType mutate_kmer(KmerType& kmer, KmerType& nuc, size_t pos);
    void remove_nuc(std::map<KmerType, unsigned int>& nuc, size_t pos);
};

template<size_t span>
FindSNP<span>::FindSNP(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindSNP<span>::contains(KmerType kmer)
{
    Node node = Node(Node::Value(kmer));
    return this->_find->graph_contains(node);
}

template<size_t span>
typename FindSNP<span>::KmerType FindSNP<span>::correct(KmerType& kmer, KmerType& nuc, size_t pos)
{
    KmerType mutate = this->mutate_kmer(kmer, nuc, pos);
    return min(mutate, revcomp(mutate, this->_find->kmer_size()));
}

//// mutate_kmer
// fonction to mutate kmer : takes kmer, pos (de la fin ,0-based), and nt (nt = 0,1,2ou 3)
// par exemple  kmer ,2 , C    avec kmer =  AAAAAAAAAA
//
// return :
// AAAAAAACAA
template<size_t span>
typename FindSNP<span>::KmerType FindSNP<span>::mutate_kmer(KmerType& kmer, KmerType& nuc, size_t pos)
{
    KmerType reset_mask = ~((KmerType)3 << (pos*2));
    KmerType set_mask = nuc << (pos*2);
    return (kmer & reset_mask ) | set_mask;
}

template<size_t span>
void FindSNP<span>::remove_nuc(std::map<KmerType, unsigned int>& nuc, size_t pos)
{
    // pos in history = integer part of pos / kmer_size + 1
    unsigned int index_swip = (pos / this->_find->kmer_size()) + 1;
    unsigned int bin_cache = 3 << (pos % this->_find->kmer_size());

    //Get the mutate nuc and remove this in map
    KmerType snp_nuc = this->_find->het_kmer_history(this->_find->het_kmer_begin_index() - index_swip).kmer & bin_cache >> (pos % this->_find->kmer_size()); // obtain the nuc

    nuc.erase(snp_nuc);
}

template<size_t span>
class FindSoloSNP : public FindSNP<span>
{
public :

    typedef typename FindSNP<span>::KmerType KmerType;

public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindSoloSNP(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindSoloSNP<span>::FindSoloSNP(FindBreakpoints<span> * find) : FindSNP<span>(find){}

template<size_t span>
bool FindSoloSNP<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }

    if(this->_find->gap_stretch_size() == this->_find->kmer_size())
    {
	// Create map with all nuc A = 0, C = 1, T = 2, G = 3
	std::map<KmerType, unsigned int> nuc;
	nuc[0] = 0;
	nuc[1] = 0;
	nuc[2] = 0;
	nuc[3] = 0;
	
	this->remove_nuc(nuc, this->_find->kmer_size());

	// iterate one possible new value
	for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end(); nuc_it++)
	{
	    for(size_t i = this->_find->het_kmer_begin_index() - 1, j = 0; i != this->_find->het_kmer_begin_index() + this->_find->kmer_size() - 1; i++, j++)
	    {
	        KmerType tmp = nuc_it->first;
		if(this->contains(this->correct(this->_find->het_kmer_history(i).kmer, tmp, j)))
		{
		    nuc_it->second++;
		}
	    }
	}

	for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end(); nuc_it++)
	{
	    if(nuc_it->second == this->_find->kmer_size())
	    {
		string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
		string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());

		this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1, kmer_begin_str, kmer_end_str, 0,STR_SNP_TYPE);
		this->_find->breakpoint_id_iterate();

		return true;
	    }
	}
    }

    return false;
}

template<size_t span>
class FindDuoSNP : public FindSNP<span>
{
public :

    typedef typename FindSNP<span>::KmerType KmerType;

public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindDuoSNP(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindDuoSNP<span>::FindDuoSNP(FindBreakpoints<span> * find) : FindSNP<span>(find){}

template<size_t span>
bool FindDuoSNP<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }

    if(this->_find->gap_stretch_size() > this->_find->kmer_size())
    {
	unsigned int snp_dist = this->_find->gap_stretch_size() - this->_find->kmer_size();

	// Create map with all nuc A = 0, C = 1, T = 2, G = 3
	std::map<KmerType, unsigned int> nuc_first;
	nuc_first[0] = 0;
	nuc_first[1] = 0;
	nuc_first[2] = 0;
	nuc_first[3] = 0;

	std::map<KmerType, unsigned int> nuc_second;
	nuc_second[0] = 0;
	nuc_second[1] = 0;
	nuc_second[2] = 0;
	nuc_second[3] = 0;

	this->remove_nuc(nuc_first, this->_find->kmer_size());
	this->remove_nuc(nuc_second, this->_find->kmer_size() + snp_dist);

	// iterate one possible new value
	for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc_first.begin(); nuc_it != nuc_first.end(); nuc_it++)
	{
	    for(size_t i = this->_find->het_kmer_begin_index() - 1, j = 0; i != this->_find->het_kmer_begin_index() + this->_find->kmer_size() - 1; i++, j++)
	    {
		KmerType tmp = nuc_it->first;
		if(this->contains(this->correct(this->correct(this->_find->het_kmer_history(i).kmer, tmp, j), tmp, j + snp_dist)))
		{
		    nuc_it->second++;
		}
	    }
	}

	for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc_first.begin(); nuc_it != nuc_first.end(); nuc_it++)
	{
	    if(nuc_it->second == this->_find->kmer_size())
	    {
		string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
		string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());

		this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1, kmer_begin_str, kmer_end_str, 0,STR_MSNP_TYPE);
		this->_find->breakpoint_id_iterate();

		return true;
	    }
	}
    }

    return false;
}

template<size_t span>
class FindBackup : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindBackup(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindBackup<span>::FindBackup(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindBackup<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }

    if(this->_find->gap_stretch_size() > (this->_find->kmer_size() / 2)) {
	string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
	string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());
	string chrom_name_bak = this->_find->chrom_name()+"_backup";

	this->_find->writeBreakpoint(this->_find->breakpoint_id(), chrom_name_bak, this->_find->position() - 1, kmer_begin_str, kmer_end_str, 0, STR_BKP_TYPE);

	return true;
    }

    return false;
}

template<size_t span>
class FindHeteroInsert : public IFindObserver<span>
{
public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindHeteroInsert(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindHeteroInsert<span>::FindHeteroInsert(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindHeteroInsert<span>::update()
{
    if(!this->_find->homo_only())
    {
        // hetero site detection
	if(!this->_find->kmer_end_is_repeated() && this->_find->current_info().nb_in == 2 && !this->_find->recent_hetero())
	{
	    //loop over putative repeat size (0=clean, >0 fuzzy), reports only the smallest repeat size found.
	    for(int i = 0; i <= this->_find->max_repeat(); i++)
	    {
		if(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).nb_out == 2 && !this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).is_repeated)
		{
		    //hetero breakpoint found
		    string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index()+i).kmer);
		    string kmer_end_str = this->_find->model().toString(this->_find->current_info().kmer);
		    this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position()-1+i, kmer_begin_str, kmer_end_str,i, STR_HET_TYPE);

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
    }

    return false;
}


#endif /* _TOOL_FindObserver_HPP_ */
