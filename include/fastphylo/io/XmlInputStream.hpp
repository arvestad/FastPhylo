#pragma once

#include <cstdio>
#include <iostream>
#include <fstream>
#include <libxml/xmlreader.h>
#include "fastphylo/core/Sequence.hpp"
#include "fastphylo/io/Extrainfos.hpp"

// Layout Phase C (input side): shared XML sequence-file parsing for
// fastdist and fastprot. readSequences()/the constructor were ~85-90%
// identical (confirmed by diff before merging) - genuine differences
// were the RelaxNG validation schema (DNA vs protein, now a constructor
// parameter) and a real bug: fastprot's constructor called
// strncmp(filename, "-", 1) *before* checking filename for null, unlike
// fastdist's correctly-guarded `filename && !strncmp(...)` - a crash on
// `fastprot -I xml` reading from stdin (no filename given), not covered
// by any RunExamples.sh fixture. Fixed here using fastdist's correct
// guard, not silently carried forward.
enum streamstatus { RUN_NOT_FINISHED = 0, RUN_FINISHED = 1 };

struct locator_t {
  int in_root;
  int in_runs;
  int in_run;
  int in_seq;
};

class XmlSequenceReader {
public:
  XmlSequenceReader(char *filename, const char *relaxngSchemaStr);
  ~XmlSequenceReader();

  bool readSequences(std::vector<Sequence> &seqs, std::string &runId, Extrainfos &extrainfos);

protected:
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
};
