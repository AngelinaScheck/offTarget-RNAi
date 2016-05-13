#include <iostream>
#include <vector>
#include <algorithm>

#include <seqan/sequence.h>
#include <seqan/find.h>

#include "countFinds.h"
#include "inputTransc.h"

//initialize List of contingency tables (one for every kmer)
void initializeCont (seqan::StringSet<seqan::DnaString> kmers, std::vector<Contingency> & allContigs){
    //Contingeny table for finding of one single kmer to store intermediate results
    Contingency oneContig;
    //Initialize list of Contingency table for all the kmers
    for (unsigned i=0; i<length(kmers); i++){
        oneContig.kmerSeq=getValue(kmers, i);
        oneContig.kmerDN=0;
        oneContig.kmerNoDN=0;
        oneContig.noKmerDN=0;
        oneContig.noKmerNoDN=0;
    
        allContigs.push_back(oneContig);
    }
}

//______________________________________________________________________________________________________________________________________________

//A new contingency table is created for each occucance of each kmer, later the contingency tables for same kmer are combined
void  countFinds (std::vector<Contingency> & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome transcripts){
    // Database are the mRNAs from Transcriptome;
    
    // Define the Aho-Corasick pattern over the kmers
    seqan::Pattern<seqan::StringSet<seqan::DnaString>, seqan::AhoCorasick> pattern(kmers);

    // Search for the kmers database sequence by database sequence.
    //search kmers in down-/upregulated genes
    for (unsigned i = 0; i < length(transcripts.mRNAset); i++){
        std::cout <<"processing mRNA number"<<'\t'<< i<<'\n';
        seqan::Finder<seqan::DnaString> finder(transcripts.mRNAset[i]);  // new finder for each seq
        while(find(finder, pattern)){
            //assign the hit to the correct kmer sequence by looping over the kmer sequences and comparing
            for(unsigned j=0; j<allContigs.size();j++){
                if(allContigs[j].kmerSeq==infix(finder)){
                    //if mRNA is up-/downregulated enhance counter for kmerDN
                    if(transcripts.isReg[i]){
                        allContigs[j].idDN.push_back(i);
                    }
                    else {
                        allContigs[j].idNoDN.push_back(i);
                    }
                    break;
                }
            }
        }
        //to check just the first mRNA (remove later)
        break;
    }
}

//_________________________________________________________________________________________________________________________________________
//calculate #no hits in affected and unaffected genes from results of countFinds, correct: every mRNA is only counted as one hit, even if kmer appears multiple times
void fillFields (std::vector<Contingency> & allContigs, Transcriptome regulated, Transcriptome notRegulated){
    //total number of down-/upregulated mRNAs
     int nReg=regulated.ids.size()<< '\n';
    //total number of unaffected mRNAs
    int nNoReg=notRegulated.ids.size();
    //iterator to delete duplicates
    //std::vector<int>::iterator it;
    for(unsigned i=0; i<allContigs.size(); i++){
        //correct counters, for multiple hits in one mRNA (list is sorted--> we can use std::unique)
        //for affected mRNAs that contain the kmer
        auto it = std::unique (allContigs[i].idDN.begin(), allContigs[i].idDN.end());
        allContigs[i].idDN.resize( std::distance(allContigs[i].idDN.begin(),it));
        allContigs[i].kmerDN=allContigs[i].idDN.size();
        //for unaffected mRNAs that contain the kmer
        auto it2 = std::unique (allContigs[i].idNoDN.begin(), allContigs[i].idNoDN.end());
        allContigs[i].idNoDN.resize( std::distance(allContigs[i].idNoDN.begin(),it2));
        allContigs[i].kmerNoDN=allContigs[i].idNoDN.size();
        //#kmer not contained=#mRNA-#kmer contained
        allContigs[i].noKmerDN=nReg-allContigs[i].kmerDN;
        allContigs[i].noKmerNoDN=nNoReg-allContigs[i].kmerNoDN;
    }
}