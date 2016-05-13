#ifndef TRANSC_H
#define TRANSC_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"

struct Transcriptome
{
    std::vector<std::string> ids;
    std::vector<std::string> genes;
    std::vector<float> values;
    seqan::StringSet<seqan::DnaString> mRNAset;
    std::vector<bool> isReg;
    
};

void getTransc (ModifyStringOptions & options, Transcriptome & transcripts);

void sortMRNA (Transcriptome & transcripts, Transcriptome & regulated, Transcriptome & notRegulated);

#endif