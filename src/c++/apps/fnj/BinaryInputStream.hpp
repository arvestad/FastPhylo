/*
 * BinaryInputStream.hpp
 *
 * 	Created on: Dec 14, 2011
 *      Auther: Mehmood Alam Khan
 *       Email: malagori@kth.se
 */
#pragma once

#include "DataInputStream.hpp"
#include <cstdio>
#include <iostream>
#include <fstream>
#include <vector>
#include <string>

using namespace std;

class BinaryInputStream : public DataInputStream
{
  public:
    BinaryInputStream(char *filename);
    ~BinaryInputStream() override;
    readstatus readDM(StrDblMatrix &dm, std::vector<std::string> &names, std::string &runId,
                      Extrainfos &extrainfos) override;

  protected:
    istream *fp;
    ifstream fin;
    bool file_was_opened;
    int newSize;

  private:
    // binary_dm_format_plan.md: replaces the old single-shot
    // `bool input_was_read` (which could only ever fire once for the
    // whole stream, so a second run's header was misread as more matrix
    // floats) with an explicit per-run state machine, now that the wire
    // format records how many matrix bodies belong to each header.
    //
    // NeedHeader: the next readDM() call must parse a fresh header.
    // InRun: reuse the cached names; read one more body.
    // JustFinishedRun: the call right after a run's last body was
    // delivered - returns END_OF_RUN without touching the stream, then
    // goes back to NeedHeader for the call after that. Required by
    // processRuns()'s `for (; readDM(...)==DM_READ; )` loop shape, which
    // needs a dedicated call to learn "stop" before the next run's real
    // data starts - mirrors how XmlInputStream already signals run
    // boundaries via END_OF_RUN/END_OF_RUNS.
    enum class RunState
    {
        NeedHeader,
        InRun,
        JustFinishedRun
    };
    RunState state = RunState::NeedHeader;
    size_t bodiesRemainingInRun = 0;

    // Parses the next run's "FASTPHYLO 2" tag + node count +
    // matricesPerRun into newSize/bodiesRemainingInRun (validates the
    // tag, throws on a mismatch rather than silently misparsing).
    // Returns false if the stream has nothing left at all (true end of
    // file - the same "no more real data" signal PhylipDmInputStream
    // already uses, so fnj's existing -r-bounded outer loop needs no
    // changes).
    bool readNextHeader();
};
