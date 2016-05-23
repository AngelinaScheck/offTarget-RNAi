#include <iostream>
#include <vector>
#include <algorithm>
#include <sys/time.h>

#include <seqan/sequence.h>
#include <seqan/find.h>
#include <seqan/index.h> 

#include "countFinds.h"
#include "inputTransc.h"
#include <math.h>



//translate kmer to a number to later hash the kmer in pattern matching
unsigned kmerToID(seqan::DnaString kmer){
//     total number of kmers is divided into quarters:
//     first quarter: kmer starts with A
//     second: kmer starts with C
//     third: kmer starts with G
//     fourth: kmer starts with T
//     for next position A: 1st quarter of 1st quarter, etc...
    
    //length of a single kmer
    unsigned l;
    l=length(kmer);
    //temporary variables for the loop
    unsigned add;
    unsigned position;
    position=0;

    for (unsigned i=0; i<l; i++){
        if(getValue(kmer, i)=='A') {add=0;}
        else if (getValue(kmer, i)=='C') {add=pow(4, (l-1-i));}
        else if (getValue(kmer, i)=='G') {add=(pow(4, l-1-i))*2;}
        else {add=(pow(4, l-1-i))*3;}
        position=position + add;
    }
    
    return position;
}




//_______________________________________Counting with indexed pattern matching________________________________________________
void countFindsIndex (Contingency & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome transcripts){
    
    //initialize finders
    seqan::Index<seqan::StringSet<seqan::DnaString>, seqan::IndexQGram<seqan::UngappedShape<6>, seqan::OpenAddressing>> index(transcripts.mRNAset);
    seqan::Finder<seqan::Index<seqan::StringSet<seqan::DnaString>> > finder(transcripts.mRNAset);
    //kmer temporary saved in pattern
    seqan::DnaString pattern;
    
    for (unsigned kN=0; kN < length(kmers); kN++){
        clear(finder);
        pattern=getValue(kmers, kN);
        std::cout << "pattern matching kmer" <<'\t'<<kN << '\n';
        while(find(finder, pattern)){
            //assign the hit to the correct kmer sequence by converting the pattern sequence into an ID
            
            //IMPROVE HERE!!!!!!!!!!!!!!!!!!!!!!!
            // unsigned id= kmerToID(pattern);
//             if(transcripts.isReg[i]){
//                 allContigs.idDN[id].push_back(i);
//             }
//             else {
//                 allContigs.idNoDN[id].push_back(i);
//             }
        }
    }
}


//initialize List of contingency tables (one for every kmer)
void initializeCont (seqan::StringSet<seqan::DnaString> kmers, Contingency & allContigs){
//     //Contingeny table for finding of one single kmer to store intermediate results
//     Contingency oneContig;
    unsigned l;
    l=length(kmers);
    //mock vector for mRNA IDs
    std::vector<int> mock;
    mock.push_back(-1);
    //reserve space
    allContigs.kmerDN.reserve(l);
    allContigs.kmerNoDN.reserve(l);
    allContigs.noKmerDN.reserve(l);
    allContigs.noKmerNoDN.reserve(l);
    allContigs.idNoDN.reserve(l);
    allContigs.idDN.reserve(l);
    
    //Initialize list of Contingency table for all the kmers
    for (unsigned i=0; i<l; i++){
        appendValue(allContigs.kmerSeq, "");
        allContigs.kmerDN.push_back(0);
        allContigs.kmerNoDN.push_back(0);
        allContigs.noKmerDN.push_back(0);
        allContigs.noKmerNoDN.push_back(0);
        //initialize ID Lists for the mRNAs, in which the kmer was found with -1, since there is no mRNA with rank -1
        allContigs.idNoDN.push_back(mock);
        allContigs.idDN.push_back(mock);
    }
    //assign the kmer to the position in the Contingencies Table
    for (unsigned j=0; j<l; j++){
        if (length(getValue(allContigs.kmerSeq, kmerToID(getValue(kmers, j))))>0){
            std::cerr << "kmer already set" << '\n';
            std::cerr << getValue(allContigs.kmerSeq, kmerToID(getValue(kmers, j)))<< '\n';
            return;
        }
        else {assignValue(allContigs.kmerSeq, kmerToID(getValue(kmers, j)), getValue(kmers, j));}
    }
}

//______________________________________________________________________________________________________________________________________________

//fills contingency tables with online pattern matching and kmerToID- method (see above)
void  countFinds (Contingency & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome transcripts){
    timeval online;
    gettimeofday(&online, NULL);
    double startOnline= (online.tv_sec *1000000) + online.tv_usec;
    // Database are the mRNAs from Transcriptome;
    
    // Define the Aho-Corasick pattern over the kmers
    seqan::Pattern<seqan::StringSet<seqan::DnaString>, seqan::AhoCorasick> pattern(kmers);
    // Search for the kmers database sequence by database sequence.
    //search kmers in down-/upregulated genes
    unsigned id;
    std::cout <<"processing mRNAs"<<'\n';
    for (unsigned i = 0; i < length(transcripts.mRNAset); i++){
        //std::cout <<"processing mRNA number"<<'\t'<< i<<'\n';
        seqan::Finder<seqan::DnaString> finder(transcripts.mRNAset[i]);  // new finder for each seq
        while(find(finder, pattern)){
            
            //assign the hit to the correct kmer sequence by converting the pattern sequence into an ID
            id= kmerToID(infix(finder));
            
            if(transcripts.isReg[i]){
                allContigs.idDN[id].push_back(i);
            }
            else {
                allContigs.idNoDN[id].push_back(i);
            }
        }
    }
    gettimeofday(&online, NULL);
    double endOnline= (online.tv_sec *1000000) + online.tv_usec;
    std::cout << length(kmers) << " kmers checked in " << (endOnline-startOnline)/1000000  << " seconds with Aho-Corasick pattern matching" << '\n';
}

//_________________________________________________________________________________________________________________________________________
//calculate #no hits in affected and unaffected genes from results of countFinds, correct: every mRNA is only counted as one hit, even if kmer appears multiple times
void fillFields (Contingency & allContigs, Transcriptome regulated, Transcriptome notRegulated){
    //total number of down-/upregulated mRNAs
     int nReg=regulated.ids.size()<< '\n';
    //total number of unaffected mRNAs
    int nNoReg=notRegulated.ids.size();
    //iterator to delete duplicates
    //std::vector<int>::iterator it;
    for(unsigned i=0; i<allContigs.kmerDN.size(); i++){
        //correct counters, for multiple hits in one mRNA (list is sorted--> we can use std::unique)
        //for affected mRNAs that contain the kmer (sort, iterator for duplicates, delete duplicates)
        std::sort(allContigs.idDN[i].begin(), allContigs.idDN[i].end());
        auto it = std::unique (allContigs.idDN[i].begin(), allContigs.idDN[i].end());
        allContigs.idDN[i].resize( std::distance(allContigs.idDN[i].begin(), it));
        //-1 because of the mock place holder at initialization of idDN
        allContigs.kmerDN[i]=allContigs.idDN[i].size()-1;
        //for unaffected mRNAs that contain the kmer
        std::sort(allContigs.idNoDN[i].begin(), allContigs.idNoDN[i].end());
        auto it2 = std::unique (allContigs.idNoDN[i].begin(), allContigs.idNoDN[i].end());
        allContigs.idNoDN[i].resize( std::distance(allContigs.idNoDN[i].begin(),it2));
        allContigs.kmerNoDN[i]=allContigs.idNoDN[i].size()-1;
        //#kmer not contained=#mRNA-#kmer contained
        allContigs.noKmerDN[i]=nReg-allContigs.kmerDN[i];
        allContigs.noKmerNoDN[i]=nNoReg-allContigs.kmerNoDN[i];
    }
}