#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <seqan/arg_parse.h>

struct ModifyStringOptions
{
    //Expression Data File
    std::string expression;
    //mRNA Sequences File
    std::string sequences;
    double cutoff;
    std::string transcriptome;
};

seqan::ArgumentParser::ParseResult
parseCommandLine(struct ModifyStringOptions & options, int argc, char const ** argv);



#endif