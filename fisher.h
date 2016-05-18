#ifndef FISHER_H
#define FISHER_H

#include <iostream>
#include <vector>
#include <boost/math/distributions/hypergeometric.hpp>

#include "inputTransc.h"
#include "countFinds.h"

struct Results{
    //kmer sequence
    seqan::StringSet<seqan::DnaString> signfKmers;
    //counter:
    std::vector<int> kmerDN;
    std::vector<int> kmerNoDN;
    std::vector<int> noKmerDN;
    std::vector<int> noKmerNoDN;
    //p-value, q-value=p-Value after multiple testing correction
    std::vector<double>pValue;
    std::vector<double>qValue;
    //IDs of the downregulated mRNAs, where the kmer was found
    std::vector<std::string> mRNAIDs;
};

//calculates probability of the contingency table, its p-value and sorts entry according to level of significance
void significant(Contingency & allContigs, double alpha, unsigned nReg, unsigned nMRNAs, Results & results, Transcriptome & transcripts);

//multiple testing correction
void benjHoch (Results & results, double alpha);

//quicksort for Results Datastructure
void quickSort(Results & results, unsigned left, unsigned right);

#endif