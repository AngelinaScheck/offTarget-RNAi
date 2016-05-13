#include <iostream>

#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "kmers.h"


//kmers through permutation (takes long)
seqan::StringSet<seqan::DnaString> makeKmer(unsigned int k){
    
    //Idea: copy every precursor 4 times, append A, C, G or T, appended string is new precursor
    //anchor

    seqan::StringSet<seqan::DnaString> precursors;
    seqan::StringSet<seqan::DnaString> kmers;
    
    appendValue(precursors, "");
    int oldLength=1;
    
    for (unsigned j=0; j<k;j++) {
        
        //remove old precursors
        clear(kmers);
        
        //typedef for Iterator in kmers.h
        TStringSetIterator it;
        for (TStringSetIterator it = begin(precursors); it != seqan::end(precursors); ++it){
            
            //copy and append
            seqan::DnaString addA = *it;
            appendValue(addA, 'A');
            seqan::DnaString addC = *it;
            appendValue(addC, 'C');
            seqan::DnaString addG = *it;
            appendValue(addG, 'G');
            seqan::DnaString addT = *it;
            appendValue(addT, 'T');
            
            appendValue(kmers, addA);
            appendValue(kmers, addC);
            appendValue(kmers, addG);
            appendValue(kmers, addT);
        }
        precursors=kmers;

    }
    std::cout<<length(kmers)<<'\t'<< "kmers generated"<<'\n';
    return kmers;
}
