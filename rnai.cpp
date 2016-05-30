#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>
#include <algorithm>
#include <numeric>
#include <sys/time.h>


#include <seqan/sequence.h>
#include <seqan/stream.h>
#include <seqan/align.h>
#include <seqan/index.h>

#include "parse.h"
#include "expLev.h"
#include "matchMRNA.h"
#include "inputTransc.h"
#include "kmers.h"
#include "countFinds.h"
#include "fisher.h"


template <typename TStringSet, typename TIndexSpec>
void qgramCounting(Transcriptome & transcripts, TStringSet & kmers, Contingency & contingencies, TIndexSpec)
{
    typedef seqan::Index<TStringSet, TIndexSpec> TIndex;

    // build Index over mRNAs
    timeval indexTime;
    TStringSet set;
    set=transcripts.mRNAset;
    
    gettimeofday(&indexTime, NULL);
    double startIndex= (indexTime.tv_sec *1000000) + indexTime.tv_usec;
    TIndex index(set);
    indexRequire(index, seqan::QGramCounts());
    gettimeofday(&indexTime, NULL);
    double endIndex= (indexTime.tv_sec *1000000) + indexTime.tv_usec;
    
    std::cout << "Index for "<< length(set) << " mRNAs build in " << (endIndex-startIndex)/1000000  << " seconds" << '\n';

    ////////////////////////////////////////////////////////////////////
    //Indexed Pattern Matching
    seqan::Finder<TIndex> finder(set);
    for(unsigned int i=0; i < length(kmers); ++i){
        clear(finder);
        while (find(finder, kmers[i] ))
        {
            // getValueI1 returns index of a string in mRNAset, if mRNA is downregulated, the index is appended to idDN in the contingencies table at the position of the kmer
            //if the mRNA is not downregulated, the id is appended to idNoDN at the position of the kmer in the contingencies
            if (transcripts.isReg[getValueI1(beginPosition(finder))]){
                contingencies.idDN[i].push_back(getValueI1(beginPosition(finder)));
            }
            else {
                contingencies.idNoDN[i].push_back(getValueI1(beginPosition(finder)));
            }
        }
    }
    gettimeofday(&indexTime, NULL);
    double endFind= (indexTime.tv_sec *1000000) + indexTime.tv_usec;
    std::cout << length(kmers) << " kmers checked in " << (endFind-endIndex)/1000000  << " seconds with indexed pattern matching" << '\n';
    
}






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
    if(transcripts.ids.size()==0){
        return 1;
    }
    //number of mRNAs in total
    unsigned nMRNAs;
    nMRNAs=length(transcripts.mRNAset);
    
    //mRNA affected by RNAi
    Transcriptome regulated;
    //mRNA not affected by RNAi
    Transcriptome notRegulated;
    //decide, if mRNA is up-/downregulated or not affected bei RNAi according to cutoff;
    
    sortMRNA (options, transcripts, regulated, notRegulated);
    //number of downregulated mRNAs
    unsigned nReg;
    nReg= length(regulated.mRNAset);
    
    //generate kmers
    seqan::StringSet<seqan::DnaString> kmers;
    kmers=makeKmer(options.k);
    
    //find kmers
    //initialization of Contingency Table
    Contingency contingencies;
    //contingencies.reserve(length(kmers));
    initializeCont(kmers, contingencies);
    
    //indexed pattern matching
    qgramCounting(transcripts, kmers, contingencies, seqan::IndexQGram<seqan::UngappedShape<8> >() );
    
    //fill counters
   //countFinds (contingencies, kmers, transcripts);
    //fill empty fields of contingency tables
    fillFields (contingencies, regulated, notRegulated);
    //decide if enrichment of kmers in downregulated mRNA is significant, create results
    Results results;
    significant(contingencies, options, nReg, nMRNAs, results, transcripts);
    //multiple testing correction (Benjamin Hochberg)
    std::cout <<length(results.signfKmers)<<'\t' << "suspicous kmers" << '\n';
    //benjHoch (results, options);
    
    //print results to output file
    std::ofstream output;
    output.open (options.output);
    
    output << "rank" << '\t' << "kmer" <<'\t' << "score" <<'\n';
    for (unsigned i=0; i<results.kmerDN.size(); i++){
        output << "rank " << i+1 << '\t' << results.signfKmers[i] << '\t' << results.enrichment[i] << '\n';
            
    }
    output.close();


    return 0;
}