#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"


//==============declarations===================================
//Options
struct ModifyStringOptions
{
    //Expression Data File
    std::string expression;
    //mRNA Sequences File
    std::string sequences;
    double cutoff;
};


//==========================methods=============================
//parser-Function
seqan::ArgumentParser::ParseResult
parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
    //ArgumentParser for reading the expression level table
    seqan::ArgumentParser parser("readExp");

    //filename for expression.txt
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "EXPRESSION"));
    //filename for mRNA sequences.txt
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "SEQUENCES"));
    //cutoff
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::DOUBLE, "CUTOFF"));

    //Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    //If parsing was not successful then exit with code 1 if there were errors.
    //Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
   
    //assign parsed arguments
    getArgumentValue(options.expression, parser, 0);
    getArgumentValue(options.sequences, parser, 1);
    getArgumentValue(options.cutoff, parser, 2);
    
    return seqan::ArgumentParser::PARSE_OK;
}
//------------------------------------------------------------------------------------

//downregulated mRNA IDS and expression levels from file
    //read expression file
   
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

//-------------------------------------------------------------------

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
//---------------------------------------------------

int main(int argc, char const ** argv)
{
    // Parse the command line.
    ModifyStringOptions options;
    seqan::ArgumentParser::ParseResult res = parseCommandLine(options, argc, argv);

    // If parsing was not successful then exit with code 1 if there were errors.
    // Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res == seqan::ArgumentParser::PARSE_ERROR;
    
    
    //read expression file
    //get downregulated mRNAs
    std::vector<long> IDsDown;
    std::vector<float> ExpValuesDown;
    
    getExpLev (options, IDsDown, ExpValuesDown);

    
    //==========================================================
    
    //read mRNA file into String Set
    seqan::StringSet<seqan::DnaString> mRNAset;
    //Entries from Expression List removed, if no matching sequence
    std::vector<long> IDsDownInList;
    std::vector<float> ExpValuesDownInList;
    matchMRNA (options, IDsDown, ExpValuesDown, mRNAset,IDsDownInList,ExpValuesDownInList);


    return 0;
    

}