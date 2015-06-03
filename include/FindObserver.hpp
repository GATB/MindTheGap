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
    KmerType mutate_kmer(KmerType& kmer, KmerType& nuc, size_t pos);
    void remove_nuc(std::map<KmerType, unsigned int>& nuc, size_t pos);
    bool snp_at_end(unsigned char* beginpos, size_t limit, KmerType* ret_nuc);
};

template<size_t span>
FindSNP<span>::FindSNP(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

template<size_t span>
bool FindSNP<span>::contains(KmerType kmer)
{
    kmer = std::min(kmer, revcomp(kmer, this->_find->kmer_size()));
    Node node = Node(Node::Value(kmer));
    return this->_find->graph_contains(node);
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
    size_t p = this->_find->kmer_size() - pos;
    KmerType reset_mask = ~((KmerType)3 << (p*2));
    KmerType set_mask = nuc << (p*2);

    return (kmer & reset_mask ) | set_mask;
}

template<size_t span>
void FindSNP<span>::remove_nuc(std::map<KmerType, unsigned int>& nuc, size_t pos)
{
    // pos in history = integer part of pos / kmer_size + 1
    unsigned int index_swip = (pos / this->_find->kmer_size());
    unsigned int bin_cache = 3 << (pos % this->_find->kmer_size());

    //Get the mutate nuc and remove this in map
    KmerType snp_nuc = this->_find->het_kmer_history(this->_find->het_kmer_begin_index() - index_swip).kmer & bin_cache >> (pos % this->_find->kmer_size()); // obtain the nuc

    nuc.erase(snp_nuc);
}

template<size_t span>
bool FindSNP<span>::snp_at_end(unsigned char* beginpos, size_t limit, KmerType* ret_nuc)
{
    // Create map with all nuc A = 0, C = 1, T = 2, G = 3
    std::map<KmerType, unsigned int> nuc;
    nuc[0] = 0;
    nuc[1] = 0;
    nuc[2] = 0;
    nuc[3] = 0;

    // Compute index_begin and index_end
    unsigned char endpos = (*beginpos + limit) % 256;

    this->remove_nuc(nuc, this->_find->kmer_size());

    // iterate one possible new value
    bool end = false;

    for(unsigned char j = 0; !end && j != this->_find->kmer_size(); (*beginpos)++, j++)
    {
	for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end();)
	{
	    KmerType const_fix = nuc_it->first;
	    KmerType correct_kmer = this->mutate_kmer(this->_find->het_kmer_history(*beginpos).kmer, const_fix, this->_find->kmer_size() - j);
	    if(this->contains(correct_kmer))
	    {
		nuc[nuc_it->first]++;
		++nuc_it;
	    }
	    else
	    {
        	if(nuc.size() == 1)
		{
		    end = true;
		    (*beginpos) -= 1; // Last iteration didn't create valid kmer we need decrement value
		    break;
		}
		nuc.erase(nuc_it++);
	    }
	}
    }

    KmerType max = nuc.begin()->first;
    for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end(); nuc_it++)
    {
	if(nuc_it->second > nuc[max])
	{
	    max = nuc_it->first;
	}
    }

    if((unsigned int)nuc[max] >= limit)
    {
	*ret_nuc = max;
	return true;
    }

    return false;
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
	KmerType nuc;
	unsigned char pos = this->_find->het_kmer_begin_index() - 1;
	if(this->snp_at_end(&pos, this->_find->kmer_size(), &nuc))
	{
	    string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
	    string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());

	    this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 1, kmer_begin_str, kmer_end_str, 0, STR_SNP_TYPE);
	    this->_find->breakpoint_id_iterate();

	    return true;
	}
    }

    return false;
}

template<size_t span>
class FindMultiSNP : public FindSNP<span>
{
public :

    typedef typename FindSNP<span>::KmerType KmerType;

public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindMultiSNP(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();

    void correct_history(unsigned char pos, KmerType nuc);
};

template<size_t span>
FindMultiSNP<span>::FindMultiSNP(FindBreakpoints<span> * find) : FindSNP<span>(find){}

template<size_t span>
bool FindMultiSNP<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }

    int kmer_threshold = 2;
    // Not content 2 snp with minimal distance
    if(this->_find->gap_stretch_size() > this->_find->kmer_size() + kmer_threshold)
    {
	int nb_snp = 0;
	// - 1 because pos is upper 1 at the end of gap
	size_t begin_pos = this->_find->position() - this->_find->gap_stretch_size() + this->_find->kmer_size() - 1;
	size_t end_pos = this->_find->position() - 1;

	// % 256 because buffer history size is equale to 256 
	unsigned char index_end = this->_find->het_kmer_begin_index() + this->_find->kmer_size() - 1;
	unsigned char index_pos = index_end - this->_find->gap_stretch_size();

	// We read all kmer in gap
	while(index_pos != index_end)
	{
	    unsigned char save_index = index_pos;
	    KmerType nuc;

	    // if detect snp at end
	    if(this->snp_at_end(&index_pos, kmer_threshold, &nuc))
	    {
		this->correct_history(save_index, nuc);
		nb_snp++;
		string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(save_index-1).kmer);
		string kmer_end_str = this->_find->model().toString(this->_find->het_kmer_history(save_index+this->_find->kmer_size()).kmer);
       
		this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), begin_pos - 1, kmer_begin_str, kmer_end_str, 0, STR_MSNP_TYPE);
		this->_find->breakpoint_id_iterate();

	        begin_pos += (index_pos - save_index);
	    }
	    // else return false
	    else
	    {
		break;
	    }
        }

	//Set value for future detection
	int dist_to_end = end_pos - begin_pos;
	if(dist_to_end > 0)
	{
	    this->_find->m_gap_stretch_size -= dist_to_end;
	    this->_find->m_solid_stretch_size += dist_to_end;
	    this->_find->m_kmer_begin.set(this->_find->het_kmer_history(index_pos).kmer, revcomp(this->_find->het_kmer_history(index_pos).kmer, this->_find->kmer_size()));

	    return false;
	}

	return true;
    }

    return false;
}

template<size_t span>
void FindMultiSNP<span>::correct_history(unsigned char pos, KmerType nuc)
{
    for(unsigned int i = 0; i != this->_find->kmer_size(); i++)
    {
	unsigned char index = (i + pos) % 256;
        this->_find->het_kmer_history(index).kmer = this->mutate_kmer(this->_find->het_kmer_history(index).kmer, nuc, this->_find->kmer_size() - i);
    }
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

	this->_find->recent_hetero(max(0, this->_find->recent_hetero() - 1));  // when recent_hetero=0 : we are sufficiently far from the previous hetero-site
    }

    return false;
}


#endif /* _TOOL_FindObserver_HPP_ */
