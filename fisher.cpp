#include <iostream>
#include <vector>
#include <boost/math/distributions/hypergeometric.hpp>

#include "inputTransc.h"
#include "countFinds.h"
#include "fisher.h"

//check also documentation on boost.org on Fishers test
void significant(std::vector<Contingency> & allContigs, double alpha, unsigned nReg, unsigned nMRNAs, Results & results, Transcriptome & transcripts){
    std::string mRNAIDs;
    double pvalue;
    //how often was the kmer found?
    unsigned kmerIn;
    //perform fisher's test on every Contingency Table
    for (unsigned i=0; i<allContigs.size(); i++){
    //Nullhypothesis: kmer is equally likely to occur in up-/downregulated mRNAs as in non-affected mRNAs
    //hypergeometric distribution from boost: hypergeometric_distribution(unsigned r, unsigned n, unsigned N) with r=number of mRNAs that contain the kmer, n=number of down-/upregulated mRNAs, N=total number of mRNAs
    
    //for the one tailed test, we are interested in cumulative probability of the distribution (cdf in boost), that kmerDN reaches a value higher than the observed. (complement, because otherwise it would sum up propabilties for values equal or lower, kmerDN-1, because otherwise the value equal to the observed would not be included)
        
        //how often was the kmer found?
        kmerIn=allContigs[i].kmerDN + allContigs[i].kmerNoDN;
        //skip if the kmer was not found at all
        if(kmerIn<0){
            boost::math::hypergeometric_distribution<double> hg_dist (kmerIn, nReg, nMRNAs);
            pvalue=boost::math::cdf(boost::math::complement(hg_dist, allContigs[i].kmerDN-1));
            //add the kmer, its pvalue,IDs of mRNAs, where it was found to list of results if p-value is equal or smaller alpha/2
            if(pvalue<alpha/2){
                //convert vector<unsigned> of positions of the mRNA in the transcriptome into string of the mRNA IDs 
                for(unsigned j=0; j<allContigs[i].idDN.size(); j++){
                mRNAIDs= mRNAIDs + transcripts.ids[allContigs[i].idDN[j]] + "";
                }
                appendValue(results.signfKmers, allContigs[i].kmerSeq);
                results.kmerDN.push_back(allContigs[i].kmerDN);
                results.kmerNoDN.push_back(allContigs[i].kmerNoDN);
                results.noKmerDN.push_back(allContigs[i].noKmerDN);
                results.noKmerNoDN.push_back(allContigs[i].noKmerNoDN);
                results.pValue.push_back(pvalue);
                results.mRNAIDs.push_back(mRNAIDs);
            }
        }
    }
}