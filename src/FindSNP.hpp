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

#ifndef _TOOL_FindSNP_HPP_
#define _TOOL_FindSNP_HPP_

/*******************************************************************************/
#include <IFindObserver.hpp>
#include <FindBreakpoints.hpp>

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

    /** Change nucleotide at position in kmer with new nucleotide
     * \param[in] kmer The kmer you want modificate
     * \param[in] nuc The new nucleotide
     * \param[in] pos The position where you want place the new nucleotide
     * \return The new kmer
     */
    KmerType mutate_kmer(KmerType& kmer, KmerType& nuc, size_t pos);

    /** Find the last nucleotide of a kmer in history and remove this nucleotide in the hash map
     * \param[in] nuc nucleotide hash map
     * \param[in] pos position of kmer in history
     */
    void remove_nuc(std::map<KmerType, unsigned int>& nuc, size_t pos);

    /** Find if we can find a snp at the end of kmer at a position in history
     * \param[in] The position of the first kmer in history
     * \param[in] The number of kmer need to be validate
     * \param[out] The nucleotide in reads
     * \param[out] The nucleotide in reference
     */
    bool snp_at_end(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc);

    bool snp_at_end(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc, unsigned int* nb_kmer_val);
    
    bool snp_at_begin(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc);

    bool snp_at_begin(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc, unsigned int* nb_kmer_val);

    char nuc_to_char(KmerType nuc);
};

template<size_t span>
FindSNP<span>::FindSNP(FindBreakpoints<span> * find) : IFindObserver<span>(find){}

//// mutate_kmer
// fonction to mutate kmer : takes kmer, pos (du debut ,1-based), and nt (nt = 0,1,2ou 3)
// par exemple  kmer ,2 , C    avec kmer =  AAAAAAAAAA
//
// return :
// ACAAAAAAAA
template<size_t span>
typename FindSNP<span>::KmerType FindSNP<span>::mutate_kmer(KmerType& kmer, KmerType& nuc, size_t pos)
{
    size_t p = this->_find->kmer_size() - pos;
    KmerType trois; trois.setVal(3);
    KmerType reset_mask = ~(trois << (p*2));
    KmerType set_mask = nuc << (p*2);

    return (kmer & reset_mask ) | set_mask;
}

template<size_t span>
char FindSNP<span>::nuc_to_char(KmerType nuc)
{
    if(nuc==0){
    	return 'A';
    }
    else{
    	if(nuc==1){
    		return 'C';
    	}
    	else{
    		if(nuc==2){
    			return 'T';
    		}
    		else{
    			return 'G';
    		}
    	}
    }
}

//NO LONGER used
template<size_t span>
void FindSNP<span>::remove_nuc(std::map<KmerType, unsigned int>& nuc, size_t pos)
{
    //Get the mutate nuc and remove this in map
    KmerType snp_nuc = this->_find->het_kmer_history(pos).kmer & 3; // obtain the nuc

    nuc.erase(snp_nuc);
}

template<size_t span>
bool FindSNP<span>::snp_at_end(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc)
{
    unsigned int tmp;
    return snp_at_end(beginpos, limit, ret_nuc, ref_nuc, &tmp);
}


// when exiting, beginpos always point to first non solid kmer
template<size_t span>
bool FindSNP<span>::snp_at_end(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc, unsigned int* nb_kmer_val) //mute le dernier nt du kmer a position beginpos de listorique puis avance
{
    // Create map with all nuc A = 0, C = 1, T = 2, G = 3
    std::map<KmerType, unsigned int> nuc;
    KmerType zero, un, deux, trois;
    zero.setVal(0); un.setVal(1); deux.setVal(2); trois.setVal(3);
    nuc[zero] = 0;
    nuc[un] = 0;
    nuc[deux] = 0;
    nuc[trois] = 0;
	
    unsigned char endpos = (*beginpos + limit) % 256;
	
    unsigned char  beginpos_init = (*beginpos);
    //this->remove_nuc(nuc, *beginpos);
    *ref_nuc = this->_find->het_kmer_history(*beginpos).kmer & 3; // obtain the reference nuc
    nuc.erase(*ref_nuc);
	
    // if end is false or if didn't read all kmer loop
    bool end = false;
    for(unsigned char j = 0; !end && j != this->_find->kmer_size(); (*beginpos)++, j++)
	{
		// for each nucleotide of nuc map
		for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end();)
		{
			KmerType const_fix = nuc_it->first; // fix conversion error
			KmerType correct_kmer = this->mutate_kmer(this->_find->het_kmer_history(*beginpos).kmer, const_fix, this->_find->kmer_size() - j);
			if(this->contains(correct_kmer))
			{
				nuc[nuc_it->first]++;
				++nuc_it;
			}
			else
			{
				
				if(nuc.size() == 1) //Is the last nucleotide and the last iteration
				{
					end = true;
					(*beginpos) -= 1; // Last iteration didn't create valid kmer we need decrement value //
					//will still be incr by end of upper for loop
					
					break;
				}
				nuc.erase(nuc_it++); // This nucleotide didn't valid kmer we remove it
			}
		}
	}
	
    //Find the max nucleotide correct most kmer
    KmerType max = nuc.begin()->first;
    for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end(); nuc_it++)
	{
		if(nuc_it->second > nuc[max])
		{
			max = nuc_it->first;
		}
	}
	
    // If nuc max is upper or equale limit we find a snp
    if((unsigned int)nuc[max] >= limit)
	{
		*ret_nuc = max;
		*nb_kmer_val = nuc[max];
		return true;
	}
    else
	{
		*beginpos = beginpos_init;
		return false;
	}
    return false;
}

template<size_t span>
bool FindSNP<span>::snp_at_begin(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc)
{
    unsigned int tmp;
    return snp_at_begin(beginpos, limit, ret_nuc, ref_nuc, &tmp);
}

// when exiting, beginpos always point to first non solid kmer
template<size_t span>
bool FindSNP<span>::snp_at_begin(unsigned char* beginpos, size_t limit, KmerType* ret_nuc, KmerType* ref_nuc, unsigned int* nb_kmer_val) //mute le premier nt du kmer a position beginpos de listorique puis recule
{
	// Create map with all nuc A = 0, C = 1, T = 2, G = 3
	std::map<KmerType, unsigned int> nuc;
	KmerType zero, un, deux, trois;
	zero.setVal(0); un.setVal(1); deux.setVal(2); trois.setVal(3);
	nuc[zero] = 0;
	nuc[un] = 0;
	nuc[deux] = 0;
	nuc[trois] = 0;
	
	
	unsigned char  beginpos_init = (*beginpos);
	//this->remove_nuc(nuc, *beginpos - (this->_find->kmer_size()-1));
	*ref_nuc = (this->_find->het_kmer_history(*beginpos).kmer) >>  (2*(this->_find->kmer_size()-1)) & 3; // obtain the reference nuc
	//bug : should take first nt here, not last one
	nuc.erase(*ref_nuc);
	
	//printf("snp at begin : bpos %i  : %c \n",*beginpos, this->nuc_to_char(*ref_nuc));
	
	// if end is false or if didn't read all kmer loop
	bool end = false;
	for(unsigned char j = 0; !end && j != this->_find->kmer_size(); (*beginpos)--, j++)
	{
		// for each nucleotide of nuc map
		for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end();)
		{
			KmerType const_fix = nuc_it->first; // fix conversion error
			KmerType correct_kmer = this->mutate_kmer(this->_find->het_kmer_history(*beginpos).kmer, const_fix,  j+1);
			
			
			if(this->contains(correct_kmer))
			{
				nuc[nuc_it->first]++;
				++nuc_it;
			}
			else
			{
				
				if(nuc.size() == 1) //Is the last nucleotide and the last iteration
				{
					end = true;
					(*beginpos) += 1; // Last iteration didn't create valid kmer we need decrement value //
					//will still be incr by end of upper for loop
					break;
				}
				nuc.erase(nuc_it++); // This nucleotide didn't valid kmer we remove it
			}
		}
	}
	
	//Find the max nucleotide correct most kmer
	KmerType max = nuc.begin()->first;
	for(typename std::map<KmerType, unsigned int>::iterator nuc_it = nuc.begin(); nuc_it != nuc.end(); nuc_it++)
	{
		if(nuc_it->second > nuc[max])
		{
			max = nuc_it->first;
		}
	}
	
	// If nuc max is upper or equale limit we find a snp
	if((unsigned int)nuc[max] >= limit)
	{
		*ret_nuc = max;
		*nb_kmer_val = nuc[max];
		return true;
	}
	else
	{
		*beginpos = beginpos_init;
		return false;
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
		KmerType ref_nuc; // reference nucleotide
		KmerType nuc; // alternative nucleotide
		unsigned char pos = this->_find->het_kmer_begin_index() - 1;
		if(this->snp_at_end(&pos, this->_find->kmer_size(), &nuc, &ref_nuc))
		{
			//string kmer_begin_str = this->_find->model().toString(this->_find->kmer_begin().forward());
			//string kmer_end_str = this->_find->model().toString(this->_find->kmer_end().forward());
			//this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 2, kmer_begin_str, kmer_end_str, 0, STR_SNP_TYPE);
			
			char ref_char [2];
			ref_char[0] = this->nuc_to_char(ref_nuc);
			ref_char[1] = '\0';
			char alt_char [2];
			alt_char[0] = this->nuc_to_char(nuc);
			alt_char[1] = '\0';
			
			this->_find->writeVcfVariant(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - 2, ref_char, alt_char, 0, STR_SNP_TYPE);
			this->_find->breakpoint_id_iterate();
			this->_find->solo_snp_iterate();
			return true;
		}
	}
	
	return false;
}

template<size_t span>
class FindFuzzySNP : public FindSNP<span>
{
public :

    typedef typename FindSNP<span>::KmerType KmerType;

public :

    /** \copydoc IFindObserver::IFindObserver
     */
    FindFuzzySNP(FindBreakpoints<span> * find);

    /** \copydoc IFindObserver::update
     */
    bool update();
};

template<size_t span>
FindFuzzySNP<span>::FindFuzzySNP(FindBreakpoints<span> * find) : FindSNP<span>(find){}

template<size_t span>
bool FindFuzzySNP<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }

    if(this->_find->gap_stretch_size() >= this->_find->kmer_size() - this->_find->max_repeat())
    {
	int delta = this->_find->kmer_size() - this->_find->gap_stretch_size();
	KmerType ref_nuc;
	KmerType nuc;
	unsigned char pos = this->_find->het_kmer_begin_index() - 1;
	if(this->snp_at_end(&pos, this->_find->kmer_size(), &nuc, &ref_nuc))
	{
	    string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(this->_find->het_kmer_begin_index() - delta - 1).kmer);
	    string kmer_end_str = this->_find->model().toString(this->_find->het_kmer_history(pos).kmer);

		
	    this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), this->_find->position() - delta, kmer_begin_str, kmer_end_str, 0, STR_SNP_TYPE);
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
	int kmer_threshold = this->_find->snp_min_val();
	// Not content 2 snp with minimal distance
	if(this->_find->gap_stretch_size() > this->_find->kmer_size() + kmer_threshold)
	{
		int nb_snp = 0;
		// - 1 because pos is upper 1 at the end of gap
		size_t begin_pos = this->_find->position() -1  - this->_find->gap_stretch_size() + this->_find->kmer_size() - 1;//position dans le genome du snp
		size_t begin_pos_init = begin_pos;

		// % 256 because buffer history size is equale to 256
		unsigned char index_end = this->_find->het_kmer_begin_index() + this->_find->kmer_size() - 1; // premier kmer solide
		unsigned char index_pos = index_end - this->_find->gap_stretch_size(); //premier kmer non solide


		// We read all kmer in gap
		while(index_pos != index_end)
		{

			unsigned char save_index = index_pos;
			unsigned int nb_kmer_val = 0;
			KmerType ref_nuc;
			KmerType nuc;

			// if detect snp at end
			if(this->snp_at_end(&index_pos, kmer_threshold, &nuc, &ref_nuc, &nb_kmer_val)) // avance tant que au moins un kmer solide
			{
				if(begin_pos + nb_kmer_val - begin_pos_init > this->_find->m_gap_stretch_size){
					// verifying that we did not go beyond the gap
					break;
				}
				
				this->correct_history(save_index, nuc);
				nb_snp++;
				//string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(save_index-1).kmer);
				//string kmer_end_str = this->_find->model().toString(this->_find->het_kmer_history(save_index+this->_find->kmer_size()).kmer);
				//this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), begin_pos , kmer_begin_str, kmer_end_str, 0, STR_MSNP_TYPE);

				char ref_char [2];
				ref_char[0] = this->nuc_to_char(ref_nuc);
				ref_char[1] = '\0';
				char alt_char [2];
				alt_char[0] = this->nuc_to_char(nuc);
				alt_char[1] = '\0';


				this->_find->writeVcfVariant(this->_find->breakpoint_id(), this->_find->chrom_name(), begin_pos, ref_char, alt_char, 0, STR_SNP_TYPE);
				//in vcf : type = SNP and not MSNP, not to conduse the user
				this->_find->breakpoint_id_iterate();
				this->_find->multi_snp_iterate();

				begin_pos += nb_kmer_val;
			}
			// else return false
			else
			{
				break;
			}
		}

		//Set value for future detection
		unsigned int nb_kmer_correct = begin_pos - begin_pos_init;
		if(nb_kmer_correct == 0)
		{
			return false;
		}

		if(nb_kmer_correct != this->_find->gap_stretch_size())
		{
			this->_find->m_gap_stretch_size -= nb_kmer_correct;
			this->_find->m_solid_stretch_size += nb_kmer_correct;
			this->_find->m_kmer_begin.set(this->_find->het_kmer_history(index_pos-1).kmer, revcomp(this->_find->het_kmer_history(index_pos-1).kmer, this->_find->kmer_size()));


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



//find one or more SNP starting from the right side of the gap, then going left
template<size_t span>
class FindMultiSNPrev : public FindSNP<span>
{
public :
	
    typedef typename FindSNP<span>::KmerType KmerType;
	
public :
	
    /** \copydoc IFindObserver::IFindObserver
     */
    FindMultiSNPrev(FindBreakpoints<span> * find);
	
    /** \copydoc IFindObserver::update
     */
    bool update();
	
    void correct_history(unsigned char pos, KmerType nuc);
};

template<size_t span>
FindMultiSNPrev<span>::FindMultiSNPrev(FindBreakpoints<span> * find) : FindSNP<span>(find){}

template<size_t span>
bool FindMultiSNPrev<span>::update()
{
    if((this->_find->kmer_begin().isValid() && this->_find->kmer_end().isValid()) == false)
    {
	return false;
    }
	
    int kmer_threshold = this->_find->snp_min_val();

    if(this->_find->gap_stretch_size() > this->_find->kmer_size() + kmer_threshold)
    {
	int nb_snp = 0;

	size_t begin_pos = this->_find->position() - 2;//position dans le genome du  dernier snp du trou (pos du dernier 0)
	size_t begin_pos_init = begin_pos;
		
	// % 256 because buffer history size is equal to 256
	unsigned char index_limit = this->_find->het_kmer_end_index() - 2 - this->_find->gap_stretch_size(); // dernier kmer solide avant trou
	unsigned char index_pos = this->_find->het_kmer_end_index() - 2; //dernier kmer non solide
		
//		//debug
//		printf("histo pos %i \n",index_pos);
//	KmerType tt = 	 this->_find->het_kmer_history(index_pos).kmer;
//		cout << tt.toString(this->_find->kmer_size()) << endl;
///
		
	// We read all kmer in gap
	while(index_pos != index_limit)
	{
		unsigned char save_index = index_pos;
		unsigned int nb_kmer_val = 0;
		KmerType ref_nuc;
		KmerType nuc;

		// if detect snp at end
		if(this->snp_at_begin(&index_pos, kmer_threshold, &nuc, &ref_nuc, &nb_kmer_val)) // recule tant que au moins un kmer solide
		{
			if(begin_pos_init - (begin_pos - nb_kmer_val) > this->_find->m_gap_stretch_size){
				// verifying that we did not go beyond the gap
				break;
			}

			this->correct_history(save_index-( this->_find->kmer_size()-1), nuc);
			nb_snp++;
			//string kmer_begin_str = this->_find->model().toString(this->_find->het_kmer_history(save_index- this->_find->kmer_size() ).kmer);
			//string kmer_end_str = this->_find->model().toString(this->_find->het_kmer_history(save_index+1).kmer);
			//this->_find->writeBreakpoint(this->_find->breakpoint_id(), this->_find->chrom_name(), begin_pos , kmer_begin_str, kmer_end_str, 0, STR_MSNP_TYPE);

			char ref_char [2];
			ref_char[0] = this->nuc_to_char(ref_nuc);
			ref_char[1] = '\0';
			char alt_char [2];
			alt_char[0] = this->nuc_to_char(nuc);
			alt_char[1] = '\0';

			this->_find->writeVcfVariant(this->_find->breakpoint_id(), this->_find->chrom_name(), begin_pos, ref_char, alt_char, 0, STR_SNP_TYPE);
			//in vcf : type = SNP and not MSNP, not to conduse the user

			this->_find->breakpoint_id_iterate();
			this->_find->multi_snp_iterate();

			begin_pos -= nb_kmer_val;
		}
		// else return false
		else
		{
			break;
		}
	}

	//Set value for future detection
	unsigned int nb_kmer_correct = begin_pos_init - begin_pos ;
	if(nb_kmer_correct == 0)
	{
		return false;
	}

	if(nb_kmer_correct != this->_find->gap_stretch_size())
	{
		this->_find->m_position -=  nb_kmer_correct;
        this->_find->m_het_kmer_end_index -= nb_kmer_correct;
        this->_find->m_het_kmer_begin_index -= nb_kmer_correct;
        
		//cout << "before gap_stretch_size = " << this->_find->m_gap_stretch_size;
		this->_find->m_gap_stretch_size -= nb_kmer_correct;
		//cout << " after gap_stretch_size = " << this->_find->m_gap_stretch_size << endl;
		//this->_find->m_solid_stretch_size += nb_kmer_correct;
		this->_find->m_kmer_end.set(this->_find->het_kmer_history(index_pos+1).kmer, revcomp(this->_find->het_kmer_history(index_pos+1).kmer, this->_find->kmer_size()));

		return false;
	}

	return true;
    }

    return false;
}

//same as findmultisnp
template<size_t span>
void FindMultiSNPrev<span>::correct_history(unsigned char pos, KmerType nuc)
{
    for(unsigned int i = 0; i != this->_find->kmer_size(); i++)
    {
	unsigned char index = (i + pos) % 256;
	this->_find->het_kmer_history(index).kmer = this->mutate_kmer(this->_find->het_kmer_history(index).kmer, nuc, this->_find->kmer_size() - i);
    }
}

#endif /* _TOOL_FindSNP_HPP_ */
