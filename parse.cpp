#include <iostream>

#include <seqan/arg_parse.h>

#include "parse.h"


//parser-Function
seqan::ArgumentParser::ParseResult
parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
        //ArgumentParser
    seqan::ArgumentParser parser("input");
    
        //filename expression List with attached mRNA sequences
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::STRING, "TRANSCRIPTOME"));
    
    //cutoff
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::DOUBLE, "CUTOFF"));
    
    //k for kmers
    addArgument(parser, seqan::ArgParseArgument(
        seqan::ArgParseArgument::INTEGER, "k"));
    
    
//         //filename for expression.txt
//     addArgument(parser, seqan::ArgParseArgument(
//         seqan::ArgParseArgument::STRING, "EXPRESSION"));
//     
//     //filename for mRNA sequences.txt
//     addArgument(parser, seqan::ArgParseArgument(
//         seqan::ArgParseArgument::STRING, "SEQUENCES"));
    

    //Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    //If parsing was not successful then exit with code 1 if there were errors.
    //Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK)
        return res;
   
    //assign parsed arguments
    getArgumentValue(options.transcriptome, parser, 0);
    getArgumentValue(options.cutoff, parser, 1);
    getArgumentValue(options.k, parser, 2);
//  getArgumentValue(options.expression, parser, 3);
//  getArgumentValue(options.sequences, parser, 4);
    
    return seqan::ArgumentParser::PARSE_OK;
}