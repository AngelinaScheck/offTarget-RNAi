#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <seqan/arg_parse.h>

struct ModifyStringOptions
{
    //Input
    //name of preprocessed mRNA-expression level table
    std::string preprocessed;
    //names of unprocessed expression level list and unprocessed mRNA list
    std::string expression;
    std::string sequences;
    //cutoff for expression table
    double cutoff;
    //length of the kmers
    int k;
    //significance level alpha
    double signf;
    //Output
    std::string output;
};

seqan::ArgumentParser::ParseResult
parseCommandLine(struct ModifyStringOptions & options, int argc, char const ** argv);



#endif