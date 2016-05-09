#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"
#include "expLev.h"
#include "matchMRNA.h"

//Match mRNA sequences to expression data
void matchMRNA (ModifyStringOptions & options, std::vector<long> & IDsDown, std::vector<float> & ExpValuesDown, seqan::StringSet<seqan::DnaString> & mRNAset, std::vector<long> & IDsDownInList, std::vector<float> & ExpValuesDownInList){
    
    std::ifstream inMRNA;
    
    //if mRNA file is incorrect
    inMRNA.open(options.sequences);
    if (! inMRNA) {
        std::cerr << "unable to open mRNA file" << '\n';
        return;
    }
    
    // ID comparison
    std::string mRNAidString;
    unsigned long mRNAid;
    
    //save mRNAs in String Set
    std::string mRNAstring;
   
    //read one line at the time
    std::stringstream mlinestream;
    
    for (std::string mline; std::getline(inMRNA, mline); ) {
        mlinestream.clear();
        mlinestream.str(mline);
        
        //read ID as String and translate
        mlinestream >> mRNAidString;
        mRNAid= std::stol (mRNAidString);
        
        //sequence in column 6, scip first 5 columns;
        std::string col2, col3, col4, col5;
        mlinestream >> col2 >> col3 >> col4 >> col5;
        if (mlinestream >> mRNAstring) {
            // if ID of mRNA is also in downregulated Expression List, append the mRNA
            //Maybe there are mRNAs in the expression list without correspondance in the mRNA list--> copy all IDs with known mRNA in new vector
            for (unsigned i=0; i < IDsDown.size(); i++){
                if(IDsDown[i]==mRNAid){
                    seqan::DnaString mRNA = mRNAstring;
                    appendValue(mRNAset, mRNA);
                    IDsDownInList.push_back(IDsDown[i]);
                    ExpValuesDownInList.push_back(ExpValuesDown[i]);
                }
            }
        }
        
    }
    
    if (length(mRNAset)==0) {
        std::cerr << "no sequences that match IDs" << '\n';
        return;
    }
    
    inMRNA.close();
}