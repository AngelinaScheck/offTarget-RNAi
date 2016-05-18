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
    results.qValue.reserve(s);
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
                results.qValue.push_back(pvalue);
                results.mRNAIDs.push_back(mRNAIDs);
                mRNAIDs.clear();
            }
        }
    }
}

//______________________________________________________________________________________________________________________________________________________________

void benjHoch (Results & results, double alpha){
 
    //sort results by p-value (quicksort)
    quickSort(results, 0, results.kmerDN.size()-1);
    //according to Benjamin Hochber all kmers with qValue=alpha*rank(pValue)/#entries<= pVale are significant
    
    //calculate qValue for kmers, delete if qValue>pValue
    unsigned nEntries = results.kmerDN.size();
    for(unsigned i=0; i<results.kmerDN.size(); i++){
        results.qValue[i] = alpha * (i+1) / nEntries;
        if (results.qValue[i] > results.pValue[i]) {
            erase(results.signfKmers, i);
            results.kmerDN.erase(results.kmerDN.begin() + i);
            results.kmerNoDN.erase(results.kmerNoDN.begin() + i);
            results.noKmerDN.erase (results.noKmerDN.begin() + i);
            results.noKmerNoDN.erase (results.noKmerNoDN.begin() + i);
            results.pValue.erase (results.pValue.begin() + i);
            results.qValue.erase (results.qValue.begin() + i);
            results.mRNAIDs.erase (results.mRNAIDs.begin() + i);
        }
    }
    
}

//______________________________________________________________________________________________________________________________________________________________

void quickSort(Results & results, unsigned left, unsigned right){
    //intialize counter for left and right sublists
    unsigned i = left; 
    unsigned j = right;
    //initialize pivot element (middle)
    unsigned mid= left + (right - left) / 2;
    double pivot = results.pValue[mid];
    //temporary values
    seqan::DnaString tempSignfKmer;
    int tempKmerDN;
    int tempKmerNoDN;
    int tempNoKmerDN;
    int tempNoKmerNoDN;
    double tempPValue;
    double tempQValue;
    std::string tempMRNAID;
    
    //divide and conq.
    while (i<=j) {
        while (results.pValue[i] <= pivot){
            //all elements in left sublist smaller/equal than the pivot element remain untouched
            if (i>=j){
                break;
            }
            i++;
        }
        while (results.pValue[j] > pivot){
            //all elements in right sublist bigger than the pivot remain untouched  
            if (j<=i){
                break;
            }
            j--;
        }
        //if an element in the left sublist bigger pivot or an element in the right sublist smaller/equal the pivot element was find, swap the entries
        if (i <= j) {
            //assign the temporary variables
            tempSignfKmer=getValue(results.signfKmers, i);
            tempKmerDN= results.kmerDN[i];
            tempKmerNoDN = results.kmerNoDN[i];
            tempNoKmerDN = results.noKmerDN[i];
            tempNoKmerNoDN = results.noKmerNoDN[i];
            tempPValue = results.pValue[i];
            tempQValue = results.qValue[i];
            tempMRNAID = results.mRNAIDs[i];
            
            //swap
            assignValue(results.signfKmers, i, getValue(results.signfKmers, j));
            results.kmerDN[i] = results.kmerDN[j];
            results.kmerNoDN[i] = results.kmerNoDN[j];
            results.noKmerDN[i] = results.noKmerDN[j];
            results.noKmerNoDN[i] = results.noKmerNoDN[j];
            results.pValue[i] = results.pValue[j];
            results.qValue[i] = results.qValue[j];
            results.mRNAIDs[i] = results.mRNAIDs[j];
            
            assignValue(results.signfKmers, j, tempSignfKmer);
            results.kmerDN[j] = tempKmerDN;
            results.kmerNoDN[j] = tempKmerNoDN;
            results.noKmerDN[j] = tempNoKmerDN;
            results.noKmerNoDN[j] = tempNoKmerNoDN;
            results.pValue[j] = tempPValue;
            results.qValue[j] = tempQValue;
            results.mRNAIDs[j] = tempMRNAID;
            
            //continue
            i++;
            j--;
            }
    
        //if sublist already sorted got one level deeper    
        if(left<j){
            quickSort(results, left, j);
        }
        if(i<right){
            quickSort(results,i,right);
        }
    }
}