#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"
#include "inputTransc.h"

void getTransc (ModifyStringOptions & options, std::vector<std::string> & Ids, std::vector<std::string> & Genes, std::vector<float> & Values, seqan::StringSet<seqan::DnaString> & mRNAset){
    Ids.clear();
    Values.clear();
    clear(mRNAset);
    
    //build and check stream
    std::ifstream inTransc;
    inTransc.open(options.transcriptome);
    if (! inTransc) {
        std::cerr << "unable to open transcriptome table" << '\n';
        return;
    }
    
    //read one line at the time, jump to nex entry
    //declare variables for intermediate storage
    std::string id;
    std::string gene;
    float value;
    std::string mRNA;
    
    std::stringstream linestream;
    for (std::string line; std::getline(inTransc, line); ) {
        linestream.clear();
        linestream.str(line);
        
        //ignore header
        if (linestream >> id >> gene >> value >> mRNA) {
            //if cutoff matched
            if (options.cutoff<=value){
                Ids.push_back(id);
                Genes.push_back(gene);
                Values.push_back(value);
                appendValue(mRNAset, mRNA);
            }
        }
    }
    
    inTransc.close();
    
    //error when no up/downregulated genes
    if(Ids.size()==0){
            std::cerr << "no expression level meets cutoff" << '\n';
            return;
        }
    
}