#ifndef PARSE_H
#define PARSE_H

#include <iostream>
#include <seqan/arg_parse.h>

seqan::ArgumentParser::ParseResult
parseCommandLine(struct ModifyStringOptions & options, int argc, char const ** argv);



#endif