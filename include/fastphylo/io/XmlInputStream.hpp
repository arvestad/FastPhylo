#pragma once

#include <cstdio>
#include <iostream>
#include <fstream>
#include <optional>
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
  bool in_root;
  bool in_runs;
  bool in_run;
  bool in_seq;
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

private:
  // readSequences()'s state machine, one handler per XML element it
  // recognizes at a given nesting depth/parent state - same technique
  // as fnj/XmlInputStream.cpp's readDM() split (see that file's
  // handleXNode() methods for the full rationale). Each handler
  // re-checks its own guard condition and returns std::nullopt when it
  // doesn't apply; a non-nullopt return is readSequences()'s cue to
  // return that value immediately (only handleRunNode() ever does,
  // for `return true;` at the end of a run).
  //
  // Lint (bugprone-easily-swappable-parameters): `type` is
  // libxml2's own `xmlReaderTypes` rather than the bare `int`
  // xmlTextReaderNodeType() returns, so it can't be silently
  // transposed with the adjacent `depth` at a call site - the
  // compiler rejects it instead of a lint check having to notice.
  std::optional<bool> handleSeqNode(int depth, xmlReaderTypes type, const xmlChar *name, std::vector<Sequence> &seqs, Extrainfos &extrainfos, int &numSequences);
  std::optional<bool> handleExtrainfoNode(int depth, xmlReaderTypes type, const xmlChar *name, Extrainfos &extrainfos);
  std::optional<bool> handleRootNode(int depth, xmlReaderTypes type, const xmlChar *name);
  std::optional<bool> handleRunsNode(int depth, xmlReaderTypes type, const xmlChar *name);
  std::optional<bool> handleRunNode(int depth, xmlReaderTypes type, const xmlChar *name, std::string &runId, Extrainfos &extrainfos);
};
