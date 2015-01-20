/*****************************************************************************
 *   MindTheGap: Integrated detection and assembly of insertion variants
 *   A tool from the GATB (Genome Assembly Tool Box)
 *   Copyright (C) 2014  INRIA
 *   Authors: C.Lemaitre, G. Rizk
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

/********************************************************************************/

using namespace std;


static const char* STR_FIND        = "find";
static const char* STR_FILL = "fill";

int main (int argc, char* argv[])
{

	/** We create a command line parser. */
//	 OptionsParser parser ("MindTheGap");
//	 parser.push_back (new OptionOneParam (STR_NB_CORES,       "number of cores",                      false, "0"  ));
//	 parser.push_back (new OptionOneParam (STR_VERBOSE,        "verbosity level",                      false,  "1"));
//	 parser.push_back (new OptionNoParam  (STR_HELP,           "display help about possible options",  false       ));
//
//	 OptionsParser tmp = Graph::getOptionsParser(false);
//	 parser.add(tmp);
//	 parser.push_front (new OptionNoParam (STR_FIND, "find module", false));
//	 parser.push_front (new OptionNoParam (STR_FILL, "find module", false));

	 // parser.displayHelp(stdout);

//    try
//    {
//        /** We parse the user options. */
//        IProperties* options = parser.parse (argc, argv);
//        
//        if ((options->get(STR_FIND) != 0 && options->get(STR_FILL) != 0) || (options->get(STR_FIND) == 0 && options->get(STR_FILL) == 0))
//        {
//            throw Exception("options find and fill are incompatible, but at least one of these is mandatory");
//        }
//        
//        if (options->get(STR_FIND) != 0){
//            Finder().run (argc-1, argv+1);
//        }
//        
//    }
//    catch (OptionFailure& e)
//    {
//        e.getParser().displayErrors   (stdout);
//        e.getParser().displayWarnings (stdout);
//        e.getParser().displayHelp     (stdout);
//        e.getParser().displayVersion  (stdout);
//        return EXIT_FAILURE;
//    }
//    catch (Exception& e)
//    {
//        std::cout << std::endl << "EXCEPTION: " << e.getMessage() << std::endl;
//        return EXIT_FAILURE;
//    }
    
    if(argc<2){
        //TODO display help
        cerr << "help" << endl;
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
            Finder().run (argc-1, argv+1);
        }
        catch (OptionFailure& e) //ces 2 catch ne marchent pas si option inconnue (pb de Tool)
        {
            std::cout << "essai" << std::endl;
//            e.getParser().displayErrors   (stdout);
//            e.getParser().displayWarnings (stdout);
//            e.getParser().displayHelp     (stdout);
//            e.getParser().displayVersion  (stdout);
            return EXIT_FAILURE;
        }
        catch (Exception& e)
        {            
            std::cout << std::endl << "EXCEPTION: " << e.getMessage() << std::endl;
            return EXIT_FAILURE;
        }
    }


    return EXIT_SUCCESS;
    
}

