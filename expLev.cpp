#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"
#include "expLev.h"


void getExpLev (ModifyStringOptions & options, std::vector<long> & IDsDown, std::vector<float> & ExpValuesDown){
    std::ifstream inExpression;
    
    //if expression data file is incorrect
    inExpression.open(options.expression);
    if (! inExpression) {
        std::cerr << "unable to open expression file" << '\n';
        return;
    }
    
    //save mRNA IDs in vector of Integers
    unsigned long id;
    std::vector<long> IDs;
    //save expression level in vector of Floats
    float expLev;
    std::vector<float> ExpValues;
    
    //read one line at the time, jump to nex entry
    std::stringstream linestream;
    for (std::string line; std::getline(inExpression, line); ) {
        linestream.clear();
        linestream.str(line);
        
        //ignore header
        if (linestream >> id >> expLev) {
            IDs.push_back(id);
            ExpValues.push_back(expLev);
        }
    }
    
    inExpression.close();
    
    
    //copy entries, which are downregulated
    for (int i=0; i < IDs.size(); i++){
        if (ExpValues[i]>=options.cutoff){
            IDsDown.push_back(IDs[i]);
            ExpValuesDown.push_back(ExpValues[i]);
        }
    }
    
    if (IDsDown.size()==0){
        std::cerr << "no expression level bigger cutoff" << '\n';
        return;
    }
    
}