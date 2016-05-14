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
#include "fisher.h"



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
    //number of mRNAs in total
    unsigned nMRNAs;
    nMRNAs=length(transcripts.mRNAset);
    
    //mRNA affected by RNAi
    Transcriptome regulated;
    //mRNA not affected by RNAi
    Transcriptome notRegulated;
    //decide, if mRNA is up-/downregulated or not affected bei RNAi according to cutoff;
    sortMRNA (transcripts, regulated, notRegulated);
    //number of downregulated mRNAs
    unsigned nReg;
    nReg= length(regulated.mRNAset);
    
    //generate kmers
    seqan::StringSet<seqan::DnaString> kmers;
    kmers=makeKmer(options.k);
    
    //find kmers
    //initialization of Contingency Table
    std::vector<Contingency> contingencies;
    contingencies.reserve(length(kmers));
    initializeCont(kmers, contingencies);
    //fill counters
    countFinds (contingencies, kmers, transcripts);
    //fill empty fields of contingency tables
    fillFields (contingencies, regulated, notRegulated);
    
    //decide if enrichment of kmers in downregulated mRNA is significant, create results
    Results results;
    significant(contingencies, 0.05, nReg, nMRNAs, results, transcripts);


    return 0;
}