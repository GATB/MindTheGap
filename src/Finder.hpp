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

#ifndef _TOOL_Finder_HPP_
#define _TOOL_Finder_HPP_

/********************************************************************************/
#include <gatb/gatb_core.hpp>
using namespace std;

/********************************************************************************/

static const char* STR_URI_REF = "-ref";
static const char* STR_MAX_REPEAT = "-max-rep";;
static const char* STR_HET_MAX_OCC = "-het-max-occ";
static const char* STR_SNP_MIN_VAL = "-snp-min-val";
static const char* STR_BRANCHING_FILTER = "-branching-filter";

static const char* STR_HOMO_ONLY = "-homo-only";
static const char* STR_INSERT_ONLY = "-insert-only";
static const char* STR_SNP_ONLY = "-snp-only";
static const char* STR_DELETION_ONLY = "-deletion-only";
static const char* STR_HETERO_ONLY = "-hete-only";
static const char* STR_NO_BACKUP = "-no-backup";
static const char* STR_WITH_BACKUP = "-backup";
static const char* STR_NO_SNP = "-no-snp";
static const char* STR_NO_INSERT = "-no-insert";
static const char* STR_NO_DELETION = "-no-deletion";
static const char* STR_NO_HETERO = "-no-hetero";
static const char* STR_BED = "-bed";

static const char* STR_HOM_TYPE = "HOM";
static const char* STR_HET_TYPE = "HET";
static const char* STR_SNP_TYPE = "SNP";
static const char* STR_MSNP_TYPE = "MSNP";
static const char* STR_DEL_TYPE = "DEL";
static const char* STR_BKP_TYPE = "BACKUP";

class Finder : public Tool
{
public:

    // Constructor
    Finder ();
    ~Finder ();
	
	void FinderHelp();

    const char* _mtg_version;
    size_t _kmerSize;
    Graph _graph;

    //parameters
    int _max_repeat;
    int _het_max_occ;
    int _snp_min_val;
    int _branching_threshold;
    int _nbCores;
    bool _homo_only;
    bool _homo_insert;
    bool _hete_insert;
    bool _snp;
    bool _backup;
    bool _deletion;
    bool _small_homo;
    bool _small_hetero;
    
    //input/output files
    IBank* _refBank;
    string _breakpoint_file_name;
    FILE * _breakpoint_file;
    string _vcf_file_name;
    FILE * _vcf_file;

    string _bed_file_name;

    //results statistics
    int _nb_homo_clean;
    int _nb_homo_fuzzy;
    int _nb_hetero_clean;
    int _nb_hetero_fuzzy;
    int _nb_fuzzy_deletion;
    int _nb_clean_deletion;
    int _nb_solo_snp;
    int _nb_multi_snp;
    int _nb_backup;
    int _nb_homo_clean_indel;
    int _nb_homo_fuzzy_indel;
    int _nb_hetero_indel;
    // Actual job done by the tool is here
    void execute ();

    
private:
    
    /** fills getInfo() with parameters informations
     */
    void resumeParameters();

    /** fills getInfo() with results informations
     * arg seconds running time
     */
    void resumeResults(double seconds);

    /** writes the header of the vcf file
         */
    void writeVcfHeader();

    /** Create and use FindBreakpoints class to find gaps in the reference genome
     */
    template<size_t span>
    struct runFindBreakpoints {  void operator ()  (Finder* object); };
};

/********************************************************************************/

#endif /* _TOOL_Finder_HPP_ */

