#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"
#include "inputTransc.h"

void getTransc (ModifyStringOptions & options, Transcriptome & transcripts){
    transcripts.ids.clear();
    transcripts.values.clear();
    transcripts.genes.clear();
    clear(transcripts.mRNAset);
    
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
    //Stream
    std::stringstream linestream;
    for (std::string line; std::getline(inTransc, line); ) {
        linestream.clear();
        linestream.str(line);
        
        //ignore header, read into Transcriptome datastructure
        if (linestream >> id >> gene >> value >> mRNA) {
            transcripts.ids.push_back(id);
            transcripts.genes.push_back(gene);
            transcripts.values.push_back(value);
            appendValue(transcripts.mRNAset, mRNA);
        }
    }
    inTransc.close();
}


// split whole transcriptome: if mRNA expression meets cutoff -->regulated, else notRegulated
void sortMRNA (ModifyStringOptions & options, Transcriptome & transcripts, Transcriptome & regulated, Transcriptome & notRegulated){
    //clear
    regulated.ids.clear();
    regulated.values.clear();
    regulated.genes.clear();
    clear(regulated.mRNAset);
    
    notRegulated.ids.clear();
    notRegulated.values.clear();
    notRegulated.genes.clear();
    clear(notRegulated.mRNAset);
    
    //sort according to cutoff
    for (int i=0; i<transcripts.ids.size(); i++){
        if(transcripts.values[i] < options.cutoff){
            notRegulated.ids.push_back(transcripts.ids[i]);
            notRegulated.genes.push_back(transcripts.genes[i]);
            notRegulated.values.push_back(transcripts.values[i]);
            appendValue(notRegulated.mRNAset, (getValue(transcripts.mRNAset, i)));
        }
        
        else{
            regulated.ids.push_back(transcripts.ids[i]);
            regulated.genes.push_back(transcripts.genes[i]);
            regulated.values.push_back(transcripts.values[i]);
            appendValue(regulated.mRNAset, (getValue(transcripts.mRNAset, i)));
        }
    }
    //error when no up/downregulated genes that fullfill cutoff criteria
    if(regulated.ids.size()==0){
        std::cerr << "no expression level meets cutoff" << '\n';
        return;
        }
}
