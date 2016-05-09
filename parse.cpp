#include <iostream>

#include <seqan/arg_parse.h>

#include "parse.h"


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
    
    //filename expression List with attached mRNA sequences
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "TRANSCRIPTOME"));

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
    getArgumentValue(options.transcriptome, parser, 3);
    
    return seqan::ArgumentParser::PARSE_OK;
}