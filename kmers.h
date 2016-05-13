#ifndef KMERS_H
#define KMERS_H

#include <iostream>
#include <seqan/arg_parse.h>


typedef seqan::Iterator<seqan::StringSet<seqan::DnaString> >::Type TStringSetIterator;

seqan::StringSet<seqan::DnaString> makeKmer(unsigned int k);


#endif