#include <iostream>

#include <seqan/arg_parse.h>

#include "parse.h"

struct ModifyStringOptions
{
    //Expression Data File
    std::string expression;
    //mRNA Sequences File
    std::string sequences;
    double cutoff;
};

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