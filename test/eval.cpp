#include <stdio.h>
#include <string.h>
#include <math.h>
#include <sys/time.h>
#include <stdlib.h>
#include <unistd.h>
#include <stdint.h>
#include <vector>
#include <array>
#include <map>
#include <string>
#include <iostream>

#define MISMATCH_COST 1
#define GAP_COST 1
#define MATCH_COST 0


#define gmatrix(i,j)  matrix[(j)+ (i)*(size2+1)]


//returns nb of errors (mismatch or gap = 1 err)
//seq 1 is ref
//seq 2 is query
//peut etre donner qq char en plus pour ref au bout , permet meilelur score
//
//
//
//
//
//
//
//

//align global-global
int compare_WN(const char * s1, const char* s2, int size1, int size2, int max_gaps)
{
	max_gaps = std::max(size1,size2); // compute full matrix
	//max_gaps = size;
	int * matrix = (int *) calloc((size1+1)*(size2+1),sizeof(int));
	//int matrix[size+1][size+1];
	//matrix[0][0] = 0;
	
	//LOG("-------- WN --------\n");
	//		LOG("%.*s\n",size,s1);
	//		LOG("%.*s\n",size,s2);
	
	//		for(int ii=0; ii< size+1; ii++)
	//		{
	//			for(int jj=0; jj< size+1; jj++)
	//			{
	//				matrix[ii][jj] = 0;
	//			}
	//		}
	
	//free gaps here
	for(int i=1; i <= std::min( max_gaps,size1); ++i)
	{
		gmatrix(i,0) = i * GAP_COST; //ou 0 pour free gaps au debut
	}
	
	//first line
	for(int j=1; j <= std::min(max_gaps,size2); ++j)
	{
		gmatrix(0,j) =  j * GAP_COST ; //j * GAP_COST;
	}
	int binf, bsup;
	
	int dim1 = size1+1;
	int dim2 = size2+1;
	
	for(int i = 1; i < dim1; ++i)
	{
		binf = std::max(i-  max_gaps,1);//  was 2 * max_gaps ?
		bsup = std::min(i + max_gaps,dim2-1); //was i ?
		
		for(int j = binf; j <= bsup; ++j)
		{
			gmatrix(i,j) = 10000 ;
			
			gmatrix(i,j) = gmatrix(i-1,j-1)  +  ((s1[i-1] != s2[j-1]) ? MISMATCH_COST : MATCH_COST);
			
			if( (i>(j-max_gaps)) && ( gmatrix(i-1,j) + GAP_COST < gmatrix(i,j))   )
			{
				gmatrix(i,j) =  gmatrix(i-1,j) + GAP_COST;
			}
			if((j>(i-max_gaps)) && ( gmatrix(i,j-1) + GAP_COST < gmatrix(i,j)) )
			{
				gmatrix(i,j) = gmatrix(i,j-1)  + GAP_COST;
			}
		}
	}
	
//	
//			for(int ii=0; ii< size1+1; ii++)
//			{
//				for(int jj=0; jj< size2+1; jj++)
//				{
//					printf("%2i ",gmatrix(ii,jj));
//				}
//				printf("\n");
//			}
	
	int resu = gmatrix(size1,size2); // pas de free gaps
//	int resu = 10000;
//	//get max score last column with max gaps
//	for(int ii=size1; ii>= size1 - max_gaps; ii--)
//	{
//		//printf("cell i j  %i \n",matrix[ii,size]);
//		resu = std::min(resu,gmatrix(ii,size2));
//	}
//	

	free(matrix);
	return resu;
}


typedef struct
{
	int pos;
	int cid;
	std::string seq;
	bool truei;
	
} insert_info_t;



typedef struct
{
	int pos;
	int cid;
	bool truei;
	
} bkpt_info_t;




int main(int argc,char * argv[]){
	
	//debug
	
//	std::string s1 = "TTCAGATTCAGCGTCATTTTTTACCCTGAAACTGAATTTTAAA";
//	std::string s2 = "TTCAGATTCAGCGTCATTTTTTACCCTGAAACTGAATTTTAAAGCTTTTTGTGCAAAATTCCGTAAGATTCAGTCTTTTGGGTTACTGTAGTCGGTGTTTGCACAGAGATAGAGCTACAGTACCCAAAATCGTCATTTTTCATGAAATTCTCCATTTTGAAGCAAA";
////	
//	int  nberrs = compare_WN( s1.c_str(),  s2.c_str(), s1.size(), s2.size(), 10);
//	std::cout << s1<< std::endl;
//	std::cout << s2<< std::endl;
//
//	float pid = 1.0 - ( nberrs / (float)std::max(s1.size(),s2.size()));
//
//	
//	printf("pid %f  ( %lu %lu ) nbdiff %i\n",pid,s1.size(),s2.size(),nberrs);
//	
//	
//	printf("nb errs %i \n",nberrs);
//	
//	
//	exit(0);
////
//
	
	//
	
	
	if(argc<4)
	{
		printf("eval ref_fasta  breakpoint_file insert_fasta\n");
		exit(1);
	}
	
	int nw = 90;
	if(argc==5)
	{
		 nw = atoi(argv[4]);
	}
	
	float nw_pass= nw/(float) 100;

	
	FILE * log_err = fopen("log_err","w");
	FILE * log_true = fopen("log_true","w");
	
	

	
	FILE *ref_fasta = fopen(argv[1], "r");
	fseek(ref_fasta, 0, SEEK_END);
	long fsize = ftell(ref_fasta);
	fseek(ref_fasta, 0, SEEK_SET);
	
	char *ref_raw = (char *)malloc(fsize + 1);
	fread(ref_raw, fsize, 1, ref_fasta);
	fclose(ref_fasta);
	
	ref_raw[fsize] = 0;

	
	//map pos --> insert info
	std::map< int , insert_info_t > rmap;

	
	std::vector< std::string  > refdata;
	std::vector< std::string > refheaders;
	int i= 0;
	int j;
	int ncid = -1 ;
	char tempheader [2000];
	char * tempseq  = (char *) malloc(sizeof(char) * 1000000) ;

	char hid [2000];
	int crid, gstart,gend,glen,estart,eend,elen;
	char str;
	
	int did =0;
	int cid=0;
	int pos=0;
	
	///////////////////////////////////////////////////////////////////////////////////
	//////////////////////////// map  reference true insert in rmap ///////////////////
	///////////////////////////////////////////////////////////////////////////////////
	while (i < fsize)
	{
		//read header
		if (ref_raw[i]=='>')
		{
			ncid++;
			refheaders.push_back(std::string());
			int pp=0;
			while (ref_raw[i]!='\n')
			{
				tempheader[pp++] = ref_raw[i];
				++i;
			}
			++i;
			
			tempheader[pp] = '\0';
			
			
			
			//>deletion_1 : chr1_69719
			sscanf(tempheader,">deletion_%i : chr%i_%i ",&did,&cid,&pos);

			refheaders.back() = std::string(tempheader);

			refdata.push_back(std::string());
			
		}
		else
		{
			/// read nt sequence , skip  end line
			j=0;
			int llen = 0;
			int pp=0;
			/// On boucle jusqu'a la fin de la partition ou de la sequence
			while ((i<fsize)&&(ref_raw[i]!='>'))
			{
				if (((ref_raw[i]>='a')&&(ref_raw[i]<'z'))||((ref_raw[i]>='A')&&(ref_raw[i]<'Z')))
				{
					++j;
					tempseq[pp++] = ref_raw[i];
					llen++;
				}
				++i;
			}
			
			tempseq[pp] = '\0';

			insert_info_t insert_info;
			insert_info.pos = pos;
			insert_info.cid = cid;
			insert_info.seq =std::string(tempseq);
			
			if(rmap.find(pos)!=rmap.end())
			{
				printf("-----two insert at same pos (maybe diff chrom), contact dev about this, not yet supoorted in this eval script-----\n");
			}
			
			rmap[pos] = insert_info;
			//printf("%i\n%s \n",insert_info.pos,insert_info.seq.c_str());
		}
	}
	
	//rmap pos -> insert
	

	FILE *insert_fasta = fopen(argv[3], "r");
	fseek(insert_fasta, 0, SEEK_END);
	 fsize = ftell(insert_fasta);
	fseek(insert_fasta, 0, SEEK_SET);
	
	char *insert_raw = (char *)malloc(fsize + 1);
	fread(insert_raw, fsize, 1, insert_fasta);
	fclose(insert_fasta);
	
	insert_raw[fsize] = 0;
	std::map< int , std::vector<insert_info_t> > imap;
	
	
	 i= 0;
	 ncid = -1 ;
	
	 did =0;
	 cid=0;
	 pos=0;
	//insert
	
	////

	
	
	///////////////////////////////////////////////////////////////////////////////////
	//////////////////////////// map  insert assmbled by mtg in imap //////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	while (i < fsize)
	{
		//read header
		if (insert_raw[i]=='>')
		{
			ncid++;
			int pp=0;
			while (insert_raw[i]!='\n')
			{
				tempheader[pp++] = insert_raw[i];
				++i;
			}
			++i;
			
			tempheader[pp] = '\0';
			
			
			
			//>deletion_1 : chr1_69719
			//idem avec >bkpt2 insertion_len_18_chr1_pos_290291_repeat_2_HOM
		//	>bkpt7_chr1_pos_931449_fuzzy_0_HOM_len_65

			sscanf(tempheader,">bkpt%i_chr%i_pos_%i",&did,&cid,&pos);

			//sscanf(tempheader,">bkpt%i insertion_len_%i_chr%i_pos_%i",&did,&ilen,&cid,&pos);
			//printf("%i %i %i  \n",did,cid,pos);
			
			
		//	refdata.push_back(std::string());
			
		}
		else
		{
			/// read nt sequence , skip  end line
			j=0;
			int llen = 0;
			int pp=0;
			/// On boucle jusqu'a la fin de la partition ou de la sequence
			while ((i<fsize)&&(insert_raw[i]!='>'))
			{
				if (((insert_raw[i]>='a')&&(insert_raw[i]<'z'))||((insert_raw[i]>='A')&&(insert_raw[i]<'Z')))
				{
					++j;
					tempseq[pp++] = insert_raw[i];
					llen++;
				}
				++i;
			}
			
			tempseq[pp] = '\0';
			
			insert_info_t insert_info;
			insert_info.pos = pos;
			insert_info.cid = cid;
			insert_info.seq =std::string(tempseq);
			insert_info.truei = false;
			imap[pos].push_back(insert_info);
			//printf("%i\n%s \n",insert_info.pos,insert_info.seq.c_str());
		}
	}
	
	
	
	
	///////////////////////////////////////////////////////////////////////////////////
	//////////////////////////// map  bkpt in bmap                   //////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	FILE *bkpt = fopen(argv[2], "r");
	fseek(bkpt, 0, SEEK_END);
	fsize = ftell(bkpt);
	fseek(bkpt, 0, SEEK_SET);
	
	char *bkpt_raw = (char *)malloc(fsize + 1);
	fread(bkpt_raw, fsize, 1, bkpt);
	fclose(bkpt);
	
	bkpt_raw[fsize] = 0;
	std::map< int , bkpt_info_t > bmap;
	
	
	i= 0;
	ncid = -1 ;
	
	did =0;
	cid=0;
	pos=0;
	while (i < fsize)
	{
		//read header
		if (bkpt_raw[i]=='>')
		{
			ncid++;
			int pp=0;
			while (bkpt_raw[i]!='\n')
			{
				tempheader[pp++] = bkpt_raw[i];
				++i;
			}
			++i;
			
			tempheader[pp] = '\0';
			
			
			//>bkpt1_chr1_pos_15777_fuzzy_1_HET left_kmer
			
			
			//>deletion_1 : chr1_69719
			//idem avec >bkpt2 insertion_len_18_chr1_pos_290291_repeat_2_HOM
			//	>bkpt7_chr1_pos_931449_fuzzy_0_HOM_len_65
			
			sscanf(tempheader,">bkpt%i_chr%i_pos_%i",&did,&cid,&pos);
			
			bkpt_info_t bkpt_info;

			bkpt_info.pos = pos;
			bkpt_info.cid = cid;
			
		//	printf("did %i cid %i pos %i    %s\n",did,cid,pos,tempheader);
			if(bmap.find(pos)!=bmap.end())
			{
				bkpt_info_t bold = bmap[pos];
				if(bold.cid != cid)
					printf("-----two bkpt at same pos (with diff chrom), contact dev about this, not yet supoorted in this eval script-----\n");

			}
			bmap[pos] = bkpt_info;
			
		}
		else
		{
			j=0;
			int llen = 0;
			int pp=0;
			/// On boucle jusqu'a la fin de la partition ou de la sequence
			while ((i<fsize)&&(bkpt_raw[i]!='>'))
			{
				if (((bkpt_raw[i]>='a')&&(bkpt_raw[i]<'z'))||((bkpt_raw[i]>='A')&&(bkpt_raw[i]<'Z')))
				{
					++j;
					tempseq[pp++] = bkpt_raw[i];
					llen++;
				}
				++i;
			}
			tempseq[pp] = '\0';
		//	printf("%s\n",tempseq);
		}
	}
	
	
	int ll = 5;

	
	///////////////////////////////////////////////////////////////////////////////////
	//////////////////////////// compute recall    du find           //////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	int true_bkpt =0;
	for(auto iterator = rmap.begin(); iterator != rmap.end(); iterator++) {
		int tpos = iterator->first;
		int r_cid = iterator->second.cid;
		std::string refseq = iterator->second.seq;
		for(int ii= -ll; ii<=ll; ii++)
		{
			if(bmap.find(tpos+ii)!=bmap.end())
			{
				bkpt_info_t bfound = bmap[tpos+ii];
				if(bfound.cid == r_cid)
				{
					true_bkpt++;
				}
			}
		}
	}
	
	
	///////////////////////////////////////////////////////////////////////////////////
	//////////////////////////// compute recall                      //////////////////
	///////////////////////////////////////////////////////////////////////////////////
	
	
		//printf("-----------------\n");
	

	
	int nb_insert =0;
	int tp = 0;
	int fn = 0;
	int fp = 0;
	int good_pos = 0;
	bool found_good_pos = false;
	std::vector<insert_info_t>  * vec_insert;
	for(auto iterator = rmap.begin(); iterator != rmap.end(); iterator++) {
		//printf("%i ---> %s \n",iterator->first, iterator->second.seq.c_str() );
		int tpos = iterator->first;
		int r_cid = iterator->second.cid;
		std::string refseq = iterator->second.seq;
		
		bool found = false;
		found_good_pos = false;
		//printf("ref %i tpos %i r_cid %i \n",nb_insert,tpos,r_cid);
		for(int ii= -ll; ii<=ll; ii++)
		{
			if(imap.find(tpos+ii)!=imap.end())
			{
				vec_insert =  & (imap[tpos+ii]);
				//std::map< int , std::vector<insert_info_t> > imap;
				for(int jj=0; jj < vec_insert->size(); jj++) // loop on different insert at this pos
				{
					insert_info_t  & insert_info = (*vec_insert)[jj];
					int i_cid = insert_info.cid;
					std::string iseq = insert_info.seq;
					int nberrs =   compare_WN(refseq.c_str(), iseq.c_str(), refseq.size(),iseq.size(), 10);
					float pid = 1.0 - ( nberrs / (float)std::max(refseq.size(),iseq.size()));
					
					//std::cout << refseq << std::endl;
					//std::cout << iseq << std::endl;
					//printf("pid %f  ( %lu %lu ) nbdiff %i\n",pid,refseq.size(),iseq.size(),nberrs);
					//printf("--------------------------------------------\n");
					
					
					if(i_cid == r_cid ) //
					{
						if(!found_good_pos)
						{
							good_pos++;
							found_good_pos = true;
						}
					}
					if(i_cid == r_cid  &&  pid > nw_pass) //
					{
						tp++;
						found = true;
						insert_info.truei = true;
						
						fprintf(log_true,"%s\n",refseq.c_str());
						fprintf(log_true,"%s\n",iseq.c_str());
						fprintf(log_true,"pid %f  ( %lu %lu ) nbdiff %i  pos %i  %i/%lu \n",pid,refseq.size(),iseq.size(),nberrs,i_cid,jj+1,vec_insert->size());
						break;
					}

				}


			}
			if(found)
				break;
		}
		
		
		//debug
		//debug
		if( found_good_pos &&  !found )
		{
			fprintf(log_err,"----------- Good pos seq diff-------------\n");
			fprintf(log_err,"%s\n",refseq.c_str());
			fprintf(log_err,"------------------------------------------\n");

			for(int jj=0; jj < vec_insert->size(); jj++) // loop on different insert at this pos
			{
				insert_info_t insert_info = (*vec_insert)[jj];
				int i_cid = insert_info.cid;
				std::string iseq = insert_info.seq;
				int nberrs =   compare_WN(refseq.c_str(), iseq.c_str(), refseq.size(),iseq.size(), 10);
				float pid = 1.0 - ( nberrs / (float)std::max(refseq.size(),iseq.size()));
				
				fprintf(log_err,"%s\n",iseq.c_str());
				fprintf(log_err,"pid %f  ( %lu %lu ) nbdiff %i  pos %i  %i/%lu \n",pid,refseq.size(),iseq.size(),nberrs,i_cid,jj+1,vec_insert->size());
				fprintf(log_err,"--------------------------------------------\n");
				
			}
		}
		
		
		if(!found && !found_good_pos)
		{
			fprintf(log_err,"----------- Not found-------------\n");
			fprintf(log_err,"%s\n",refseq.c_str());
			fprintf(log_err,"----------------------------------\n");

		}
			
		
		
//		if(!found)
//		{
//			//printf("not found : %i cr %i  ---> %s \n",iterator->first,iterator->second.cid, iterator->second.seq.c_str() );
//			//fn++;
//		}
		
		nb_insert++;
	}
	
	
	/////////////////compute precision fill
	int nb_insert_filled =0;
	int nb_true_insert =0;

	
	for(auto iterator = imap.begin(); iterator != imap.end(); iterator++) {

		std::vector<insert_info_t>  & vec_insert = iterator->second;
		bool found = false;
		for(int jj=0; jj < vec_insert.size(); jj++)
		{
			insert_info_t insert_info = vec_insert[jj];

			if(insert_info.truei)
			{
				nb_true_insert ++;
				break;
			}
		}

		nb_insert_filled++;
	}
	
	
	
	
	
	printf("Find recall         %i / %lu  : %.3f\n", true_bkpt,rmap.size(), true_bkpt/(float)rmap.size()  );
	printf("Fill good loc       %i / %i  : %.3f \n",good_pos,nb_insert,good_pos/(float)nb_insert);
	printf("Recall (> %.2f)     %i / %i  : %.3f \n",nw_pass,tp,nb_insert,tp/(float)nb_insert);
	printf("Fill prec           %i / %i  : %.3f \n",nb_true_insert,nb_insert_filled,  nb_true_insert/(float)nb_insert_filled );

	fclose(log_err);
	fclose(log_true);
	

}