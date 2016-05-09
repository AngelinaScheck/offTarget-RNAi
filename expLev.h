#ifndef EXPLEV_H
#define EXPLEV_H

#include <vector>
#include <string>
#include <fstream>
#include <sstream>
#include <iostream>

#include <seqan/arg_parse.h>
#include <seqan/sequence.h>
#include <seqan/stream.h>

#include "parse.h"

void getExpLev (ModifyStringOptions & options, std::vector<long> & IDsDown, std::vector<float> & ExpValuesDown);

#endif