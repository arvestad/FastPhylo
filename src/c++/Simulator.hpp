//
//
// Functions for generating trees and sequences with the help of Beep, SeqGen, and ROSE.
//
//

#ifndef SIMULATOR_HPP
#define SIMULATOR_HPP


#include "stl_utils.hpp"
#include <stdio.h>
#include <iostream>
#include <fstream>
#include <stdlib.h>
#include <string>
#include "SequenceTree.hpp"



//ultrametricDeviation=the ultrametric tree outputed from beep is made
//non-ultrametric by multiplying each branch by random number in
//[1/ultraD,ultraD]
void
createRandomTreeUsingBeep(SequenceTree &tree, int numLeafs, int ultrametricDeviation);

//
//diameterFactor = each branch is multiplied by this number in the seq gen program
void
evolveSequencesUsingSeqGen(SequenceTree &tree, int seqlen, float diameterFactor, bool writeAncestral=false);

//exp_sub = the expected number of substitutions IN ONE POSITION from root down to leaf
//exp_indel = expected number of indels IN ONE POSITION from the root down to leaf
void
evolveSequencesUsingROSE(SequenceTree &tree,int seqlen, double exp_sub, double exp_indel);

void
createFileName(std::string &name, int numLeafs, int ultrametricDeviation, int seqlen, float diameterFactor);

void
createGnuplotFromDatFile(const char *datfilename, 
			 const std::string xlabel, const std::string ylabel,
			 const std::vector<int> columns, 
			 const std::vector<std::string> names, 
			 const std::string title);

#endif
