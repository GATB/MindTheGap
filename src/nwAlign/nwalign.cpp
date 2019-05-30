/*********************************************************************
Minimalist utility to use perform fast Needleman_Wunsch alignment outside of MindTheGap

Usage : nwalign < infile
Where infile is a two lines file with the two sequences to compare
Outputs identity score in stdout
*********************************************************************/

#include <gatb/gatb_core.hpp>
#include <Utils.cpp>
#include <iostream>

using namespace std;

int main (int argc, char* argv[])
{
    // We use a try/catch block since GATB functions may throw exceptions
    try
    {
        int nbLine = 0;
        string seq1;
        string seq2;

        float score;

        for (std::string line; std::getline(std::cin, line);) {
            nbLine += 1;
            if (nbLine == 1){
                seq1 = line;
            } else if (nbLine == 2){
                seq2 = line;
            } else{
                cout << "Only two lines expected" << endl;
                break;
            }
        }
        score = needleman_wunsch(seq1,seq2, NULL, NULL, NULL);
        cout << score << endl;
        return 0;

    }
    catch (Exception& e)
    {
        std::cout << "EXCEPTION: " << e.getMessage() << std::endl;
        return EXIT_FAILURE;
    }
}