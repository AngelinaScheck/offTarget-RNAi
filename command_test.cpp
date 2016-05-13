#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"
#include "expLev.h"
#include "matchMRNA.h"
#include "inputTransc.h"
#include "kmers.h"
#include "countFinds.h"




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
    
    

    //read transcriptome data from table
    Transcriptome transcripts;
    getTransc (options, transcripts);
    
    //mRNA affected by RNAi
    Transcriptome regulated;
    //mRNA not affected by RNAi
    Transcriptome notRegulated;
    //decide, if mRNA is up-/downregulated or not affected bei RNAi according to cutoff;
    sortMRNA (options, transcripts, regulated, notRegulated);
    std::cout<< regulated.ids.size()<< '\n';
    
    //generate kmers
    seqan::StringSet<seqan::DnaString> kmers;
    kmers=makeKmer(options.k);
    std::cout<< length(kmers)<< '\n';
    
    //find kmers
    //initialization of Contingency Table
    std::vector<Contingency> contingencies;
    initializeCont(kmers, contingencies);
    //fill counters
    countFinds (contingencies, kmers, regulated, notRegulated);
    std::cout<< contingencies.size()<< '\n';


    return 0;
}