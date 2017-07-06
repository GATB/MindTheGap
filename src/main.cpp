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

// We include the header file for the tool
#include <Finder.hpp>
#include <Filler.hpp>

/********************************************************************************/

using namespace std;

static const char* MTG_VERSION = "2.0.2";

static const char* STR_FIND        = "find";
static const char* STR_FILL = "fill";

void displayVersion(std::ostream& os){

	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
	os << "* MindTheGap version "<< MTG_VERSION << "                *" << endl; //<< " AGPL licence" <<endl;
	os << "* Using gatb-core version " << System::info().getVersion() <<  "           *" << endl;
	os << "* Supported kmer sizes <" << KSIZE_STRING << "   *" << endl;
	os << "* * * * * * * * * * * * * * * * * * * * * *" << endl;
}

void displayHelp(std::ostream& os){

	os << endl <<"MindTheGap version "<< MTG_VERSION << endl << endl;
	os << "Usage: MindTheGap <module> [module options]" <<endl << endl;
	os << "[MindTheGap modules]" << endl;
	os << "    find     :    insertion breakpoint detection" << endl;
	os << "                  usage: MindTheGap find (-in <reads.fq> | -graph <graph.h5>) -ref <reference.fa> [options]" << endl;
	os << "                  help: MindTheGap find -help"<< endl;
	os << "    fill     :    gap-filler or insertion assembly"<< endl;
	os << "                  usage: MindTheGap fill (-in <reads.fq> | -graph <graph.h5>) -bkpt <breakpoints.fa> [options]" << endl;
	os << "                  help: MindTheGap fill -help"<< endl;
	os << "[Common options]" << endl;
	os << "    -help    :    display this help menu" << endl;
	os << "    -version :    display current version" << endl;
	os << endl;

}


int main (int argc, char* argv[])
{

    
    if(argc<2){
    	displayHelp(cout);
    	return EXIT_FAILURE;
    }

    if(strcmp(argv[1],STR_VERSION)==0  ||  strcmp(argv[1],"-v")==0 ){
    	displayVersion(cout);
        return EXIT_FAILURE;
    }

    if(strcmp(argv[1],STR_HELP)==0){
    	displayHelp(cout);
    	return EXIT_FAILURE;
    }

    if ((strcmp(argv[1],STR_FIND) != 0 && strcmp(argv[1],STR_FILL) != 0 ) || (strcmp(argv[1],STR_FIND) == 0 && strcmp(argv[1],STR_FILL) == 0 ))
    {
        cerr << "options find and fill are incompatible, but at least one of these is mandatory" << endl;
        return EXIT_FAILURE;

    }
    
    if (strcmp(argv[1],STR_FIND) == 0)
    {
        try
        {
        	Finder finder = Finder();
        	finder._mtg_version = MTG_VERSION;
            finder.run (argc-1, argv+1);
        }
        catch (Exception& e)
        {
        	if(strcmp(e.getMessage(),"")!=0){
        		std::cout << std::endl << "EXCEPTION: " << e.getMessage() << std::endl;
        	}
            return EXIT_FAILURE;
        }
    }

    if (strcmp(argv[1],STR_FILL) == 0)
        {
            try
            {
				Filler filler = Filler();
				filler._mtg_version = MTG_VERSION;
                filler.run (argc-1, argv+1);
            }
            catch (Exception& e)
            {
            	if(strcmp(e.getMessage(),"")!=0){
            		std::cout << std::endl << "EXCEPTION: " << e.getMessage() << std::endl;
            	}
                return EXIT_FAILURE;
            }
        }

    return EXIT_SUCCESS;
    
}

