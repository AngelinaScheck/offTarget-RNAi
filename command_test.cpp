#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"
#include "expLev.h"
#include "matchMRNA.h"
#include "inputTransc.h"




int main(int argc, char const ** argv)
{
   //--------------for when option to build input file from scratch is set--------------------------
    // Parse the command line.
    ModifyStringOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);
// 
       //--------------for when option to build input file from scratch is set--------------------------

//     // If parsing was not successful then exit with code 1 if there were errors.
//     // Otherwise, exit with code 0 (e.g. help was printed).
//     if (res != seqan::ArgumentParser::PARSE_OK)
//         return res == seqan::ArgumentParser::PARSE_ERROR;
//     
//     
//     //read expression file
//     //get downregulated mRNAs
//     std::vector<long> IDsDown;
//     std::vector<float> ExpValuesDown;
//     
//     getExpLev (options, IDsDown, ExpValuesDown);
// 
//     
//     //==========================================================
//     
//     //read mRNA file into String Set
//     seqan::StringSet<seqan::DnaString> mRNAset;
//     //Entries from Expression List removed, if no matching sequence
//     std::vector<long> IDsDownInList;
//     std::vector<float> ExpValuesDownInList;
//     matchMRNA (options, IDsDown, ExpValuesDown, mRNAset,IDsDownInList,ExpValuesDownInList);
    
    //-----------------------------------------------------------------------------------------------------
    
    
    //declarations for storing transcriptome
    std::vector<std::string>  Ids;
    std::vector<std::string> Genes;
    std::vector<float> Values;
    seqan::StringSet<seqan::DnaString> mRNAset;
    //read transcriptome data from table
    getTransc (options, Ids, Genes, Values, mRNAset);
    


    return 0;
    

}