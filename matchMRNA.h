#ifndef MATCHMRNA_H
#define MATCHMRNA_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"
#include "expLev.h"


void matchMRNA (ModifyStringOptions & options, std::vector<long> & IDsDown, std::vector<float> & ExpValuesDown, seqan::StringSet<seqan::DnaString> & mRNAset, std::vector<long> & IDsDownInList, std::vector<float> & ExpValuesDownInList);

#endif