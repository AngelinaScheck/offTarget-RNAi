#ifndef COUNTFIND_H
#define COUNTFIND_H

#include <iostream>
#include <vector>

#include <seqan/sequence.h>
#include <seqan/find.h>

#include "inputTransc.h"

struct Contingency
{
    //Sequence of the kmer
    seqan::StringSet<seqan::DnaString> kmerSeq;
    //number of regulated mRNAs, in which the kmer is contained
    std::vector<unsigned> kmerDN;
    //numbers in transcriptome of the regulated mRNAs that contain the kmer
    std::vector<std::vector<int>> idDN;
    //number of not regulated mRNAs, in which the kmer is contained
    std::vector<unsigned> kmerNoDN;
    //number in the transcriptome structure of the not regulated mRNAs that contain the kmer
    std::vector<std::vector<int>> idNoDN;
    //number of regulated mRNAs, that do not contain the kmer
    std::vector<unsigned> noKmerDN;
    //number of not regulated mRNAs, that do not contain the kmer
    std::vector<unsigned> noKmerNoDN;
};

// //qgram
// void qgramCounting(seqan::StringSet<seqan::DnaString> & kmers, Transcriptome & transcripts, Contingency & allContigs );

//indexed pattern matching
void countFindsIndex (Contingency & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome transcripts);

//__________________________________________________________________________________________
//pattern matching of kmers and mRNAs--> count #hits in affected and unaffected genes
void countFinds (Contingency & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome transcripts);

//initialization of Contingency List
void initializeCont (seqan::StringSet<seqan::DnaString> kmers, Contingency & allContigs);

//calculate #no hits in affected and unaffected genes from results of countFinds
void fillFields (Contingency & allContigs, Transcriptome regulated, Transcriptome notRegulated);

//translate kmer sequence to an id
unsigned kmerToID(seqan::DnaString kmer);


#endif