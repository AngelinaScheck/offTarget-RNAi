#ifndef FISHER_H
#define FISHER_H

#include <iostream>
#include <vector>
#include <boost/math/distributions/hypergeometric.hpp>

#include "inputTransc.h"
#include "countFinds.h"

struct Results{
    seqan::StringSet<seqan::DnaString> signfKmers;
    std::vector<int> kmerDN;
    std::vector<int> kmerNoDN;
    std::vector<int> noKmerDN;
    std::vector<int> noKmerNoDN;
    std::vector<double>pValue;
    std::vector<std::string> mRNAIDs;
};

//calculates probability of the contingency table, its p-value and sorts entry according to level of significance
void significant(Contingency & allContigs, double alpha, unsigned nReg, unsigned nMRNAs, Results & results, Transcriptome & transcripts);


#endif