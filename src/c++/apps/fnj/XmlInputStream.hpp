#pragma once

#include <cstdio>
#include <iostream>
#include <fstream>
#include <optional>
#include <libxml/xmlreader.h>
#include "DataInputStream.hpp"
#include "fastphylo/io/fileFormatSchema.hpp"

struct locator_t {
  bool in_root;
  bool in_runs;
  bool in_run;
  bool in_identities;
  bool in_identity;
  bool in_dms;
  bool in_dm;
  bool in_row;
  int row_nr;
  int entry_nr;
};

class XmlInputStream : public DataInputStream {
public:
  XmlInputStream(char *filename);
  ~XmlInputStream() override;
  readstatus readDM( StrDblMatrix & dm, std::vector<std::string> & names, std::string & runId, Extrainfos & extrainfos) override;
protected:
  xmlTextReaderPtr reader;
  locator_t l;
  int fd;
  int dmSize;

private:
  // readDM(StrDblMatrix, ...)'s state machine, one handler per XML
  // element it recognizes at a given nesting depth/parent state
  // (readDM()'s own cognitive complexity was 82 before this split - a
  // chain of ten near-identical "if (right depth/parents/name) {
  // switch(type) {...} }" blocks). Each handler re-checks its own
  // guard condition (the same one the corresponding block in readDM()
  // used to check) and returns std::nullopt when it doesn't apply -
  // exactly equivalent to the original's "fall through to the next
  // if" behavior, since only one handler's name check can ever match
  // a given node. A non-nullopt return is readDM()'s cue to return
  // immediately with that status, matching the original's inline
  // `return DM_READ/END_OF_RUN/END_OF_RUNS;` statements.
  //
  // Lint (bugprone-easily-swappable-parameters): `type` is libxml2's
  // own `xmlReaderTypes` rather than the bare `int`
  // xmlTextReaderNodeType() returns, so it can't be silently
  // transposed with the adjacent `depth` at a call site - the
  // compiler rejects it instead of a lint check having to notice.
  std::optional<readstatus> handleEntryNode(int depth, xmlReaderTypes type, const xmlChar *name, StrDblMatrix &dm);
  std::optional<readstatus> handleRowNode(int depth, xmlReaderTypes type, const xmlChar *name);
  std::optional<readstatus> handleDmNode(int depth, xmlReaderTypes type, const xmlChar *name, std::vector<std::string> &names, StrDblMatrix &dm);
  std::optional<readstatus> handleIdentityNode(int depth, xmlReaderTypes type, const xmlChar *name, std::vector<std::string> &names, Extrainfos &extrainfos, int &nr_of_ids);
  std::optional<readstatus> handleExtrainfoNode(int depth, xmlReaderTypes type, const xmlChar *name, Extrainfos &extrainfos);
  std::optional<readstatus> handleDmsNode(int depth, xmlReaderTypes type, const xmlChar *name);
  std::optional<readstatus> handleRunNode(int depth, xmlReaderTypes type, const xmlChar *name, std::string &runId);
  std::optional<readstatus> handleIdentitiesNode(int depth, xmlReaderTypes type, const xmlChar *name, std::vector<std::string> &names, Extrainfos &extrainfos);
  std::optional<readstatus> handleRunsNode(int depth, xmlReaderTypes type, const xmlChar *name);
  std::optional<readstatus> handleRootNode(int depth, xmlReaderTypes type, const xmlChar *name);
};

