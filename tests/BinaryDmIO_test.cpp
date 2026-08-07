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

namespace {

const char *kTestFile = "BinaryDmIO_test.bin";

StrDblMatrix makeMatrix(vector<string> &names, double lastEntry) {
	StrDblMatrix dm(names.size());
	dm.setIdentifiers(names);
	double v = 1.0;
	for (size_t i = 0; i < names.size(); i++) {
		for (size_t j = i; j < names.size(); j++) {
			dm.setDistance(i, j, (i == names.size() - 1 && j == names.size() - 1) ? lastEntry : v);
			v += 1.0;
		}
	}
	return dm;
}

void assertMatrixEquals(StrDblMatrix &dm, vector<string> &names, double lastEntry) {
	assert(dm.getSize() == names.size());
	double v = 1.0;
	for (size_t i = 0; i < names.size(); i++) {
		assert(dm.getIdentifier(i) == names[i]);
		for (size_t j = i; j < names.size(); j++) {
			double expected = (i == names.size() - 1 && j == names.size() - 1) ? lastEntry : v;
			assert(dm.getDistance(i, j) == expected);
			v += 1.0;
		}
	}
}

}  // namespace

// A single matrix, no bootstrap replicates - the plain
// "fastprot -O binary" case (examples/RunExamples.sh's ex16). Also
// covers the -1.0-for-non-finite convention BinaryDmOutputStream::
// print() applies to infinite/NaN distances.
static void test_single_matrix_round_trip() {
	vector<string> names{"Alpha", "Beta", "Gamma"};
	StrDblMatrix written = makeMatrix(names, numeric_limits<double>::infinity());

	{
		BinaryDmOutputStream out(const_cast<char *>(kTestFile));
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
// input_was_read-guarded header/identifier caching: only the first
// readDM() call should consume the "FASTPHYLO 1"+size+names header, and
// later calls must reuse the cached names rather than misreading
// subsequent matrix bytes as a second header.
static void test_multiple_matrices_share_one_header() {
	vector<string> names{"One", "Two"};
	StrDblMatrix original = makeMatrix(names, 42.0);
	StrDblMatrix boot1 = makeMatrix(names, 43.0);
	StrDblMatrix boot2 = makeMatrix(names, 44.0);

	{
		BinaryDmOutputStream out(const_cast<char *>(kTestFile));
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

int main() {
	test_single_matrix_round_trip();
	test_multiple_matrices_share_one_header();
	cout << "BinaryDmIO_test: all tests passed" << endl;
	return 0;
}
