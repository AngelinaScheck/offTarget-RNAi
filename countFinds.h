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
    seqan::DnaString kmerSeq;
    //number of regulated mRNAs, in which the kmer is contained
    int kmerDN;
    //numbers in transcriptome of the regulated mRNAs that contain the kmer
    std::vector<unsigned> idDN;
    //number of not regulated mRNAs, in which the kmer is contained
    int kmerNoDN;
    //number in the transcriptome structure of the not regulated mRNAs that contain the kmer
    std::vector<unsigned> idNoDN;
    //number of regulated mRNAs, that do not contain the kmer
    int noKmerDN;
    //number of not regulated mRNAs, that do not contain the kmer
    int noKmerNoDN;
};
//pattern matching of kmers and mRNAs--> count #hits in affected and unaffected genes
void countFinds (std::vector<Contingency> & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome transcripts);

//initialization of Contingency List
void initializeCont (seqan::StringSet<seqan::DnaString> kmers, std::vector<Contingency> & allContigs);

//calculate #no hits in affected and unaffected genes from results of countFinds
void fillFields (std::vector<Contingency> & allContigs, Transcriptome regulated, Transcriptome notRegulated);


#endif