// Regression test for BinaryInputStream::readDM(StrDblMatrix&, ...)
// (src/c++/apps/fnj/BinaryInputStream.cpp) - previously an unimplemented
// stub that printed "Not implemented!" and exit(-1)'d, meaning fnj could
// never read back the binary distance-matrix format fastprot/fastdist
// have long been able to write (-O binary). This is the overload fnj's
// main.cpp actually calls (see DataInputStream.hpp's 2016-06-14 comment:
// only the StrDblMatrix form is part of the polymorphic interface).
//
// Round-trips real StrDblMatrix data through BinaryDmOutputStream (the
// writer, shared library code) and BinaryInputStream (the reader,
// fnj-local) via a real file on disk, mirroring exactly how
// BinaryDmOutputStream::printHeader()/print() are actually called from
// fastprot/fastdist's main.cpp (printStartRun() once, printHeader()
// once, then one print() per matrix - original data plus any bootstrap
// replicates - all under the same header).
//
// binary_dm_format_plan.md: also covers the run-boundary fix -
// BinaryDmOutputStream's constructor now takes a matricesPerRun count
// (how many print() bodies follow each header), written into the
// "FASTPHYLO 2" header so BinaryInputStream can tell "one more body"
// from "a new run's header" apart. test_multiple_runs_are_correctly_
// separated() below is the exact scenario that used to silently
// corrupt (see project memory / this plan doc): two full runs, each
// with its own header, back to back in one stream.

#undef NDEBUG
#include <cassert>
#include <cmath>
#include <cstdio>
#include <iostream>
#include <limits>
#include <string>
#include <vector>

#include "fastphylo/core/DistanceMatrix.hpp"
#include "fastphylo/io/BinaryDmOutputStream.hpp"
#include "fastphylo/io/Extrainfos.hpp"
#include "BinaryInputStream.hpp"

using namespace std;

namespace
{

const char *kTestFile = "BinaryDmIO_test.bin";

StrDblMatrix makeMatrix(vector<string> &names, double lastEntry)
{
    StrDblMatrix dm(names.size());
    dm.setIdentifiers(names);
    double v = 1.0;
    for (size_t i = 0; i < names.size(); i++)
    {
        for (size_t j = i; j < names.size(); j++)
        {
            dm.setDistance(i, j, (i == names.size() - 1 && j == names.size() - 1) ? lastEntry : v);
            v += 1.0;
        }
    }
    return dm;
}

void assertMatrixEquals(StrDblMatrix &dm, vector<string> &names, double lastEntry)
{
    assert(dm.getSize() == names.size());
    double v = 1.0;
    for (size_t i = 0; i < names.size(); i++)
    {
        assert(dm.getIdentifier(i) == names[i]);
        for (size_t j = i; j < names.size(); j++)
        {
            double expected = (i == names.size() - 1 && j == names.size() - 1) ? lastEntry : v;
            assert(dm.getDistance(i, j) == expected);
            v += 1.0;
        }
    }
}

} // namespace

// A single matrix, no bootstrap replicates - the plain
// "fastprot -O binary" case (examples/RunExamples.sh's ex16). Also
// covers the -1.0-for-non-finite convention BinaryDmOutputStream::
// print() applies to infinite/NaN distances.
static void test_single_matrix_round_trip()
{
    vector<string> names{"Alpha", "Beta", "Gamma"};
    StrDblMatrix written = makeMatrix(names, numeric_limits<double>::infinity());

    {
        BinaryDmOutputStream out(const_cast<char *>(kTestFile), 1);
        string runId = "run1";
        Extrainfos extrainfos;
        out.printStartRun(names, runId, extrainfos);
        out.printHeader(names.size());
        out.print(written);
    }

    BinaryInputStream in(const_cast<char *>(kTestFile));
    StrDblMatrix readBack;
    vector<string> readNames;
    string readRunId;
    Extrainfos readExtrainfos;

    readstatus status = in.readDM(readBack, readNames, readRunId, readExtrainfos);
    assert(status == DM_READ);
    assertMatrixEquals(readBack, names, -1.0);

    // Nothing more in the file - the second call must signal end-of-run,
    // not misread trailing garbage or hang.
    status = in.readDM(readBack, readNames, readRunId, readExtrainfos);
    assert(status == END_OF_RUN);

    remove(kTestFile);
}

// One header, three matrix bodies - fastprot/fastdist's "-b 2" shape
// (original data + 2 bootstrap replicates all under one printHeader()
// call, see fastprot/main.cpp's processRuns()). Verifies the
// state-machine's header/identifier caching: only the first readDM()
// call for this run should consume the "FASTPHYLO 2"+size+
// matricesPerRun+names header, and later calls must reuse the cached
// names rather than re-parsing a header from mid-body bytes.
static void test_multiple_matrices_share_one_header()
{
    vector<string> names{"One", "Two"};
    StrDblMatrix original = makeMatrix(names, 42.0);
    StrDblMatrix boot1 = makeMatrix(names, 43.0);
    StrDblMatrix boot2 = makeMatrix(names, 44.0);

    {
        BinaryDmOutputStream out(const_cast<char *>(kTestFile), 3);
        string runId = "run1";
        Extrainfos extrainfos;
        out.printStartRun(names, runId, extrainfos);
        out.printHeader(names.size());
        out.print(original);
        out.print(boot1);
        out.print(boot2);
    }

    BinaryInputStream in(const_cast<char *>(kTestFile));
    StrDblMatrix readBack;
    vector<string> readNames;
    string readRunId;
    Extrainfos readExtrainfos;

    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == DM_READ);
    assertMatrixEquals(readBack, names, 42.0);

    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == DM_READ);
    assertMatrixEquals(readBack, names, 43.0);

    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == DM_READ);
    assertMatrixEquals(readBack, names, 44.0);

    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == END_OF_RUN);

    remove(kTestFile);
}

// Two full runs (each its own header + body) back to back in one
// stream - fastdist's "-r 2" shape. This is the exact scenario
// binary_dm_format_plan.md was written to fix: before the run-boundary
// count existed, the reader had no way to tell run 2's header bytes
// apart from more of run 1's matrix floats, and silently misread them
// as data (confirmed live: reported count=3 for run 1, nothing for
// run 2, exit code 0). Both runs use the same matricesPerRun (1) since
// that's the real constraint - it's fixed once per BinaryDmOutputStream
// instance/program invocation (matching fastprot/fastdist's actual
// (no_incl_orig?0:1)+bootstraps, which never varies per dataset - see
// the plan doc), not a per-header parameter. Run 2 deliberately uses
// different names than run 1, so this also proves the reader re-parses
// fresh identifiers rather than reusing anything cached from run 1.
static void test_multiple_runs_are_correctly_separated()
{
    vector<string> names1{"Alpha", "Beta", "Gamma"};
    StrDblMatrix run1matrix = makeMatrix(names1, 100.0);

    vector<string> names2{"X", "Y"};
    StrDblMatrix run2matrix = makeMatrix(names2, 200.0);

    {
        // One BinaryDmOutputStream across both runs, calling
        // printStartRun()/printHeader()/print() twice in sequence -
        // exactly how fastdist's processRuns() writes multiple datasets
        // through a single DataOutputStream.
        BinaryDmOutputStream out(const_cast<char *>(kTestFile), 1);
        string runId1 = "run1";
        Extrainfos extrainfos1;
        out.printStartRun(names1, runId1, extrainfos1);
        out.printHeader(names1.size());
        out.print(run1matrix);

        string runId2 = "run2";
        Extrainfos extrainfos2;
        out.printStartRun(names2, runId2, extrainfos2);
        out.printHeader(names2.size());
        out.print(run2matrix);
    }

    BinaryInputStream in(const_cast<char *>(kTestFile));
    StrDblMatrix readBack;
    vector<string> readNames;
    string readRunId;
    Extrainfos readExtrainfos;

    // Run 1: one matrix, then END_OF_RUN.
    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == DM_READ);
    assertMatrixEquals(readBack, names1, 100.0);
    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == END_OF_RUN);

    // Run 2: a fresh header (different names, size 2 not 3) and its one
    // matrix, then END_OF_RUN.
    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == DM_READ);
    assertMatrixEquals(readBack, names2, 200.0);
    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == END_OF_RUN);

    // True end of file - calling again still returns END_OF_RUN, not a
    // crash or a misread of nonexistent trailing data.
    assert(in.readDM(readBack, readNames, readRunId, readExtrainfos) == END_OF_RUN);

    remove(kTestFile);
}

int main()
{
    test_single_matrix_round_trip();
    test_multiple_matrices_share_one_header();
    test_multiple_runs_are_correctly_separated();
    cout << "BinaryDmIO_test: all tests passed" << endl;
    return 0;
}
