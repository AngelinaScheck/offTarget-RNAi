#include <iostream>
#include <vector>
#include <boost/math/distributions/hypergeometric.hpp>

#include "inputTransc.h"
#include "countFinds.h"
#include "fisher.h"

//check also documentation on boost.org on Fishers test
void significant(Contingency & allContigs, double alpha, unsigned nReg, unsigned nMRNAs, Results & results, Transcriptome & transcripts){
    std::string mRNAIDs;
    double pvalue;
    //how often was the kmer found?
    unsigned kmerIn;
    //Contingency Table Size
    unsigned s;
    s=allContigs.kmerDN.size();
    //reserve space
    results.kmerDN.reserve(s);
    results.kmerNoDN.reserve(s);
    results.noKmerDN.reserve(s);
    results.noKmerNoDN.reserve(s);
    results.pValue.reserve(s);
    results.mRNAIDs.reserve(s);
    //perform fisher's test on every Contingency Table
    for (unsigned i=0; i<s; i++){
    //Nullhypothesis: kmer is equally likely to occur in up-/downregulated mRNAs as in non-affected mRNAs
    //hypergeometric distribution from boost: hypergeometric_distribution(unsigned r, unsigned n, unsigned N) with r=number of mRNAs that contain the kmer, n=number of down-/upregulated mRNAs, N=total number of mRNAs
    
    //for the one tailed test, we are interested in cumulative probability of the distribution (cdf in boost), that kmerDN reaches a value higher than the observed. (complement, because otherwise it would sum up propabilties for values equal or lower, kmerDN-1, because otherwise the value equal to the observed would not be included)
        
        //how often was the kmer found?
        kmerIn=allContigs.kmerDN[i] + allContigs.kmerNoDN[i];
        
        //skip if the kmer was not found at all
        if(kmerIn>0){
            boost::math::hypergeometric_distribution<double> hg_dist (kmerIn, nReg, nMRNAs);
            pvalue=boost::math::cdf(boost::math::complement(hg_dist, allContigs.kmerDN[i]));
            //add the kmer, its pvalue,IDs of mRNAs, where it was found to list of results if p-value is equal or smaller alpha/2
            if(pvalue<alpha/2){
                //convert vector<unsigned> of positions of the mRNA in the transcriptome into string of the mRNA IDs 
                //starting from 1 because of the -1 place holder
                for(unsigned j=1; j<allContigs.idDN[i].size(); j++){
                    unsigned id= allContigs.idDN[i][j];
                    mRNAIDs.append(transcripts.ids[id] + " ");                    
                }
                appendValue(results.signfKmers, allContigs.kmerSeq[i]);
                results.kmerDN.push_back(allContigs.kmerDN[i]);
                results.kmerNoDN.push_back(allContigs.kmerNoDN[i]);
                results.noKmerDN.push_back(allContigs.noKmerDN[i]);
                results.noKmerNoDN.push_back(allContigs.noKmerNoDN[i]);
                results.pValue.push_back(pvalue);
                results.mRNAIDs.push_back(mRNAIDs);
                mRNAIDs.clear();
            }
        }
    }
}