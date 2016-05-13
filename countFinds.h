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
    //numbers in the regulated transcriptome structure of the regulated mRNAs that contain the kmer
    std::vector<unsigned> idDN;
    //number of not regulated mRNAs, in which the kmer is contained
    int kmerNoDN;
    //number in the not regulated transcriptome structure of the not regulated mRNAs that contain the kmer
    std::vector<std::string> idNoDN;
    //number of regulated mRNAs, that do not contain the kmer
    int noKmerDN;
    //number of not regulated mRNAs, that do not contain the kmer
    int noKmerNoDN;
};

void countFinds (std::vector<Contingency> & allContigs, seqan::StringSet<seqan::DnaString> kmers, Transcriptome regulated, Transcriptome notRegulated);

void initializeCont (seqan::StringSet<seqan::DnaString> kmers, std::vector<Contingency> & allContigs);


#endif