#include <iostream>

#include <seqan/arg_parse.h>

#include "parse.h"


//parser-Function
seqan::ArgumentParser::ParseResult
parseCommandLine(ModifyStringOptions & options, int argc, char const ** argv)
{
        //ArgumentParser
    seqan::ArgumentParser parser("input");
    
    //output filename
    addOption(parser, seqan::ArgParseOption(
        "o", "output", "Name of output file",
        seqan::ArgParseArgument::STRING, "OUTPUT"));
    setDefaultValue(parser, "output", "output");
    
    //filename preprocessed expression list with mRNA sequences
    addOption(parser, seqan::ArgParseOption(
        "pre", "preprocessed", "Name of the preprocessed transcriptome data",
        seqan::ArgParseArgument::STRING, "PREPROCESSED"));
    
    //ONLY in case there is no preprocessed transcriptome table, one can also give the expression level list and the mRNA sequences seperately
    //not processed expression level list
    addOption(parser, seqan::ArgParseOption(
        "exp", "expression", "Name of the unprocessed expression level list",
        seqan::ArgParseArgument::STRING, "EXPRESSION"));
    addOption(parser, seqan::ArgParseOption(
        "seq", "sequences", "Name of the unprocessed list of mRNAs",
        seqan::ArgParseArgument::STRING, "SEQUENCES"));
    
    //cutoff
    addOption(parser, seqan::ArgParseOption(
        "c", "cutoff", "Cutoff for the expression level",
        seqan::ArgParseArgument::DOUBLE, "CUTOFF"));
    //default cutoff is 10 %
    setDefaultValue(parser, "cutoff", "0.1");
    
    //k for kmers
    addOption(parser, seqan::ArgParseOption(
        "k", "kmer", "Length of a kmer",
        seqan::ArgParseArgument::INTEGER, "KMER"));
    setDefaultValue(parser, "kmer", "6");
    
    //significance level
    addOption(parser, seqan::ArgParseOption(
        "a", "alpha", "Significance threshold",
        seqan::ArgParseArgument::DOUBLE, "SIGNF"));
    setDefaultValue(parser, "alpha", "0.05");

    //Parse command line.
    seqan::ArgumentParser::ParseResult res = seqan::parse(parser, argc, argv);

    //If parsing was not successful then exit with code 1 if there were errors.
    //Otherwise, exit with code 0 (e.g. help was printed).
    if (res != seqan::ArgumentParser::PARSE_OK){
        return res;}
    
    //The program need ether a processed table or two files: one for the expression levels, one for the mRNA sequences
    //If both (processed + expression/sequences) or expression but not sequences (or vice versa) are set, it's an error
    if ((isSet(parser, "expression") && !isSet(parser, "sequences")) || (isSet(parser, "sequences") && !isSet(parser, "expression"))){
        std::cerr<< "Program requires both: list of expression levels and list of mRNA sequences" << '\n';
        res= seqan::ArgumentParser::PARSE_ERROR;
        return res;
    }
    if (isSet(parser, "preprocessed") && isSet(parser, "sequences")){
        std::cerr<< "If a preprocessed transcriptome table is set, extra lists for expression levels or mRNA sequences are not allowed" << '\n';
        res= seqan::ArgumentParser::PARSE_ERROR;
        return res;
    }
   
    //assign parsed arguments
    getOptionValue(options.output, parser, "output");
    getOptionValue(options.preprocessed, parser, "preprocessed");
    getOptionValue(options.expression, parser, "expression");
    getOptionValue(options.sequences, parser, "sequences");
    getOptionValue(options.cutoff, parser, "cutoff");
    getOptionValue(options.k, parser, "kmer");
    getOptionValue(options.signf, parser, "alpha");
    
    return seqan::ArgumentParser::PARSE_OK;
}