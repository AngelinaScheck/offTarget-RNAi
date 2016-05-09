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

void getTransc (ModifyStringOptions & options, std::vector<std::string> & Ids, std::vector<std::string> & Genes, std::vector<float> & Values, seqan::StringSet<seqan::DnaString> & mRNAset);

#endif