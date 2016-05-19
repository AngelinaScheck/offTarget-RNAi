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
    inTransc.open(options.preprocessed);
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
            transcripts.isReg.push_back(false);
        }
    }
    inTransc.close();
}
//___________________________________________________________________________________________________
//sorts transcriptome data according to expression level using the quicksort algorithm (see also quickSort in fisher.cpp)
void sortTransc (Transcriptome & transcripts, unsigned left, unsigned right){
   
// Reminder Transcriptome Data Structure
//     std::vector<std::string> ids;
//     std::vector<std::string> genes;
//     std::vector<float> values;
//     seqan::StringSet<seqan::DnaString> mRNAset;
//     std::vector<bool> isReg;
    
    
 //intialize counter for left and right sublists
    unsigned i = left; 
    unsigned j = right;
    //initialize pivot element (middle)
    unsigned mid= left + (right - left) / 2;
    double pivot = transcripts.values[mid];
    //temporary values
    std::string tempId;
    std::string tempGene;
    float tempValue;
    seqan::DnaString tempMRNA;
    bool tempIsReg;
    //divide and conq.
    while (i<=j) {
        while (transcripts.values[i] <= pivot){
            //all elements in left sublist smaller/equal than the pivot element remain untouched
            if (i>=j){
                break;
            }
            i++;
        }
        while (transcripts.values[j] > pivot){
            //all elements in right sublist bigger than the pivot remain untouched  
            if (j<=i){
                break;
            }
            j--;
        }
        //if an element in the left sublist bigger pivot or an element in the right sublist smaller/equal the pivot element was find, swap the entries
        if (i <= j) {
            //assign the temporary variables
            tempId = transcripts.ids[i];
            tempGene = transcripts.genes[i];
            tempValue = transcripts.values[i];
            tempMRNA = getValue(transcripts.mRNAset, i);
            tempIsReg = transcripts.isReg[i];
            
            //swap
            transcripts.ids[i] = transcripts.ids[j];
            transcripts.genes[i] = transcripts.genes[j];
            transcripts.values[i] = transcripts.values[j];
            assignValue(transcripts.mRNAset, i, getValue(transcripts.mRNAset, j));
            transcripts.isReg[i] = transcripts.isReg[j];
            
            transcripts.ids[j] = tempId;
            transcripts.genes[j] = tempGene;
            transcripts.values[j] = tempValue;
            assignValue(transcripts.mRNAset, j, tempMRNA);
            transcripts.isReg[j] = tempIsReg;
            
            //continue
            i++;
            j--;
            }
    
        //if sublist already sorted got one level deeper    
        if(left<j){
            sortTransc(transcripts, left, j);
        }
        if(i<right){
            sortTransc(transcripts,i,right);
        }
    }   
    
}

//____________________________________________________________________________________________________
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
    
    //sort transcripts by expression level, rank= place in vector +1
    sortTransc (transcripts, 0 , transcripts.ids.size()-1);
    
    unsigned absCutoff = options.cutoff * transcripts.ids.size();
    //sort according to bool value in transcripts which was set according to cutoff
    for (int i=0; i<transcripts.ids.size(); i++){
        if((i+1)<= absCutoff){
            transcripts.isReg[i]=true;
            regulated.ids.push_back(transcripts.ids[i]);
            regulated.genes.push_back(transcripts.genes[i]);
            regulated.values.push_back(transcripts.values[i]);
            appendValue(regulated.mRNAset, (getValue(transcripts.mRNAset, i)));
            notRegulated.ids.push_back(transcripts.ids[i]);
            notRegulated.genes.push_back(transcripts.genes[i]);
            notRegulated.values.push_back(transcripts.values[i]);
            appendValue(notRegulated.mRNAset, (getValue(transcripts.mRNAset, i)));
        }
        
        else{
            notRegulated.ids.push_back(transcripts.ids[i]);
            notRegulated.genes.push_back(transcripts.genes[i]);
            notRegulated.values.push_back(transcripts.values[i]);
            appendValue(notRegulated.mRNAset, (getValue(transcripts.mRNAset, i)));
        }
    }
    //error when no up/downregulated genes that fullfill cutoff criteria
    if(regulated.ids.size()==0){
        std::cerr << "no expression level meets cutoff" << '\n';
        return;
        }
    std::cout << regulated.ids.size() << '\t' << "mRNAs affected by RNAi"<< '\n';
}
