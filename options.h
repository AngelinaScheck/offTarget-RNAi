#ifndef OPTIONS_H
#define OPTIONS_H

#include <iostream>

struct ModifyStringOptions
{
    //Expression Data File
    std::string expression;
    //mRNA Sequences File
    std::string sequences;
    double cutoff;
};

#endif