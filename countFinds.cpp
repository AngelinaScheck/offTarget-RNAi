#include <iostream>
#include <vector>

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



//A new contingency table is created for each occucance of each kmer, later the contingency tables for same kmer are combined
void  countFinds (std::vector<Contingency> & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome regulated, Transcriptome notRegulated){
    // Database for down-/upregulated genes in regulated.mRNAset;
    // Database for down-/upregulated genes in notRegulated.mRNAset;
    
    // Define the Aho-Corasick pattern over the kmers
    seqan::Pattern<seqan::StringSet<seqan::DnaString>, seqan::AhoCorasick> pattern(kmers);

    // Search for the kmers in the databases.  We have to search database sequence by database sequence.
    //search kmers in down-/upregulated genes
    for (unsigned i = 0; i < length(regulated.mRNAset); i++){
        seqan::Finder<seqan::DnaString> finder(regulated.mRNAset[i]);  // new finder for each seq
        while (find(finder, pattern)){
            //assign the hit to the correct kmer sequence
            for(unsigned j=0; j<allContigs.size();j++){
                if(allContigs[j].kmerSeq==infix(finder)){
                    allContigs[j].kmerDN++;
                    allContigs[j].idDN.push_back(i);
                }
            }
        }
    }
     //search kmers in not affected genes, like above
     for (unsigned k = 0; k < length(notRegulated.mRNAset); k++){
       seqan::Finder<seqan::DnaString> finder(notRegulated.mRNAset[k]);  // new finder for each seq
       while (find(finder, pattern)){
            //if the kmer was found before, combine the contingency tables
            for(unsigned l=0; l<allContigs.size();l++){
            //assign the hit to the correct kmer sequence
                if(allContigs[l].kmerSeq==infix(finder)){
                    allContigs[l].kmerNoDN++;
                    //allContigs[j].idDN.push_back(i);
                }
            }
         }
   }


    
 //   return allContigs;
}