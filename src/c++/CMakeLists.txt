# Modernization Phase 5 (modernization_plan.md): was CMakeForceCompiler's
# CMAKE_FORCE_CXX_COMPILER(mpic++/mpiCC ...), a deprecated (CMake docs:
# "once intended for... cross-compiling toolchain files", superseded
# long ago) way to fake CMake's compiler identification. find_package's
# own MPI::MPI_CXX imported target is the standard modern replacement -
# it carries the right compile flags/include dirs/link flags without
# needing to lie to CMake about which compiler is in use. Unverified
# end-to-end (no MPI installation in this environment, and
# BUILD_WITH_MPI defaults OFF), but this is the documented CMake-
# upstream-recommended pattern.
if(BUILD_WITH_MPI)
  find_package(MPI REQUIRED)
endif()


set(ARCH_FLAGS "-msse2")

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_CURRENT_BINARY_DIR}/script.nsi" DEPENDS fastdist fnj fastprot fastprot_mpi)

  add_custom_target(win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe")

# gengetopt_migration_plan.md, Phase E: fastprot/fastprot_mpi (Phase D)
# were the last two apps still using gengetopt - every self-built
# gengetopt/FTP-fallback/generated-code machinery that used to live
# here is now fully removed, not just unreferenced.

# XML is an optional feature, auto-detected: if libxml2 isn't found,
# build without XML support (a warning, not a configure/build failure)
# rather than the previous behavior, where WITH_LIBXML=ON (the default)
# unconditionally compiled the Xml*.cpp files regardless of whether
# find_package(LibXml2) actually found anything - a missing libxml2
# used to surface as a confusing build-time "header not found" error
# deep in XmlInputStream.cpp, not a clear configure-time message.
# Overwriting the WITH_LIBXML variable itself (not a separate derived
# one) when it's not found means every later `if(WITH_LIBXML)` source-
# list guard in this file, and the `#cmakedefine WITH_LIBXML` in
# config.h.cmake (so every `#ifdef WITH_LIBXML` in the apps' own
# main.cpp files too), automatically reflect real availability with no
# further changes needed anywhere else.
if(WITH_LIBXML)
  find_package(LibXml2)
  if(TARGET LibXml2::LibXml2)
    set(FASTPHYLO_XML_LIBS LibXml2::LibXml2)
  else()
    message(WARNING "libxml2 not found - building without XML input/output "
      "support (fastprot/fastdist/fnj -I/-O xml will not be available). "
      "Install libxml2 development headers (e.g. 'apt install libxml2-dev') "
      "to enable it, or pass -DWITH_LIBXML=OFF to build without it "
      "intentionally and silence this warning.")
    set(WITH_LIBXML OFF)
  endif()
endif()


configure_file(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

configure_file(${CMAKE_CURRENT_SOURCE_DIR}/script.nsi.cmake ${CMAKE_CURRENT_BINARY_DIR}/script.nsi)

# DNA_b128/sse2_wrapper.h uses x86 SSE2 intrinsics directly. On non-x86
# platforms (e.g. Apple Silicon / other ARM), it falls back to simde
# (https://github.com/simd-everywhere/simde), a header-only library that
# translates the same _mm_* calls onto NEON. Not needed, and not searched
# for, on real x86 hardware.
set(FASTPHYLO_EXTRA_INCLUDE_DIRS "")
if(NOT (CMAKE_SYSTEM_PROCESSOR MATCHES "^(x86_64|amd64|i.86|AMD64)$"))
  find_path(SIMDE_INCLUDE_DIR simde/x86/sse2.h
    HINTS /opt/homebrew/include /usr/local/include)
  if(SIMDE_INCLUDE_DIR)
    list(APPEND FASTPHYLO_EXTRA_INCLUDE_DIRS ${SIMDE_INCLUDE_DIR})
  else()
    message(WARNING "simde not found (needed for SSE2-via-NEON on ${CMAKE_SYSTEM_PROCESSOR}); "
      "install it (e.g. 'brew install simde') or DNA_b128/sse2_wrapper.h will fail to compile.")
  endif()
endif()

# Modernization Phase B (project_layout_plan.md): sources now live under
# src/fastphylo/{core,dna}/, mirroring the public headers under
# include/fastphylo/{core,dna}/ - see that plan doc for the rationale
# (a discoverable public API, and a step toward de-duplicating the I/O
# stream classes still living per-program under apps/).
set(FASTPHYLO_CORE_DIR ${CMAKE_SOURCE_DIR}/src/fastphylo/core)
set(FASTPHYLO_DNA_DIR ${CMAKE_SOURCE_DIR}/src/fastphylo/dna)

set(FASTPHYLO_SRCS
${FASTPHYLO_CORE_DIR}/BitVector.cpp
${FASTPHYLO_CORE_DIR}/InitAndPrintOn_utils.cpp
${FASTPHYLO_CORE_DIR}/Sequence.cpp
${FASTPHYLO_CORE_DIR}/SequenceTree.cpp
${FASTPHYLO_CORE_DIR}/file_utils.cpp
${FASTPHYLO_CORE_DIR}/random_utils.cpp
${FASTPHYLO_DNA_DIR}/DNA_b128_String.cpp
${FASTPHYLO_DNA_DIR}/Sequences2DistanceMatrix.cpp
${FASTPHYLO_DNA_DIR}/NeighborJoining.cpp
${FASTPHYLO_DNA_DIR}/Kimura2parameter.cpp
${FASTPHYLO_DNA_DIR}/TamuraNei.cpp
${FASTPHYLO_DNA_DIR}/ambiguity_nucleotide.cpp
${FASTPHYLO_DNA_DIR}/dna_pairwise_sequence_likelihood.cpp
${FASTPHYLO_DNA_DIR}/string_compare.cpp
${FASTPHYLO_CORE_DIR}/DistanceMatrix.cpp
${FASTPHYLO_CORE_DIR}/DistanceRow.cpp
${FASTPHYLO_CORE_DIR}/xml_output_global.cpp
)

# C++ standard is now set project-wide via CMAKE_CXX_STANDARD in the
# root CMakeLists.txt (modernization_plan.md, Phase 0) - portable across
# compilers, unlike the hand-appended "--std=c++11" GCC/Clang-only flag
# string this used to be (which also silently never applied on Clang at
# all until a fix earlier on this branch's parent, speed2026a - see that
# branch's "Fix silent C++98 build on Clang" commit).

set(FASTPHYLO_SPECIAL_SRCS
${FASTPHYLO_DNA_DIR}/sse2_wrapper.c
${FASTPHYLO_DNA_DIR}/computeTAMURANEIDistance_DNA_b128_String.cpp
${FASTPHYLO_DNA_DIR}/computeDistance_DNA_b128_String.cpp)

# Layout Phase C (project_layout_plan.md): shared output-stream classes,
# consolidated out of fastdist's and fastprot's independently-diverged
# per-program copies (see that plan doc for the audit - some, like
# Extrainfos, were near-identical; PhylipDmOutputStream/XmlOutputStream/
# DataOutputStream were genuinely different, with fastdist still on the
# pre-speed2026a unbatched write path, now ported to fastprot's faster
# one). fnj's own DataOutputStream/XmlOutputStream stay local to
# apps/fnj/ - fnj writes trees built from a distance matrix, not a
# distance matrix itself, a genuinely different shape, not just naming
# drift.
set(FASTPHYLO_IO_DIR ${CMAKE_SOURCE_DIR}/src/fastphylo/io)
set(FASTPHYLO_IO_SRCS
${FASTPHYLO_IO_DIR}/DataOutputStream.cpp
${FASTPHYLO_IO_DIR}/PhylipDmOutputStream.cpp
${FASTPHYLO_IO_DIR}/XmlOutputStream.cpp
${FASTPHYLO_IO_DIR}/BinaryDmOutputStream.cpp
${FASTPHYLO_IO_DIR}/FastaInputStream.cpp
${FASTPHYLO_IO_DIR}/PhylipMaInputStream.cpp
)

# XmlInputStream.cpp (input side) uses real libxml2 xmlTextReader APIs,
# unlike XmlOutputStream.cpp above (hand-rolled fprintf, no libxml2
# dependency) - needs its own WITH_LIBXML guard so -DWITH_LIBXML=OFF
# builds (verified working in Phase 5) don't try to compile it.
if(WITH_LIBXML)
  set(FASTPHYLO_IO_XML_SRCS ${FASTPHYLO_IO_DIR}/XmlInputStream.cpp)
endif()

# fastdist is the second app migrated off gengetopt onto CLI11
# (gengetopt_migration_plan.md, Phase C) - no generated cmdline.c to
# list here anymore.
set(FASTDIST_SRCS apps/fastdist/main.cpp
apps/fastdist/PhylipMaInputStream.cpp
apps/fastdist/FastaInputStream.cpp
)

if(WITH_LIBXML)
 set(FASTDIST_XML_SRCS apps/fastdist/XmlInputStream.cpp)
endif(WITH_LIBXML)


# fnj was the first app migrated off gengetopt onto CLI11
# (gengetopt_migration_plan.md, Phase B) - no generated cmdline.c to
# list here anymore.
set(FNJ_SRCS apps/fnj/main.cpp
apps/fnj/DataInputStream.cpp
apps/fnj/DataOutputStream.cpp
apps/fnj/TreeTextOutputStream.cpp
apps/fnj/XmlOutputStream.cpp
apps/fnj/PhylipDmInputStream.cpp
apps/fnj/BinaryInputStream.cpp
)

if(WITH_LIBXML)
 set(FNJ_XML_SRCS apps/fnj/XmlInputStream.cpp)
endif(WITH_LIBXML)

# Layout Phase C2 (project_layout_plan.md): fastprot's protein-specific
# "library" pieces, moved into the shared fastphylo library (same
# treatment as dna/ - domain-specific code that still belongs in the
# central library, not duplicated per-program). fastprot_mpi has its own
# separate, untouched copies of these files (out of this plan's scope,
# deferred - see that plan doc's "Deferred: fastprot_mpi" section).
set(FASTPHYLO_PROTEIN_DIR ${CMAKE_SOURCE_DIR}/src/fastphylo/protein)
set(FASTPHYLO_PROTEIN_SRCS
${FASTPHYLO_PROTEIN_DIR}/ProtDistCalc.cpp
${FASTPHYLO_PROTEIN_DIR}/ModelMatrix.cpp
${FASTPHYLO_PROTEIN_DIR}/ExpectedDistance.cpp
${FASTPHYLO_PROTEIN_DIR}/Matrix.cpp
${FASTPHYLO_PROTEIN_DIR}/MaximumLikelihood.cpp
${FASTPHYLO_PROTEIN_DIR}/ProtSeqUtils.cpp
${FASTPHYLO_PROTEIN_DIR}/ProtSeqCode.cpp
${FASTPHYLO_PROTEIN_DIR}/ProtSeqCompare.cpp
)

# fastprot is the third app migrated off gengetopt onto CLI11
# (gengetopt_migration_plan.md, Phase D) - no generated cmdline.c to
# list here anymore.
set(FASTPROT_SRCS apps/fastprot/main.cpp
apps/fastprot/FastaInputStream.cpp
apps/fastprot/PhylipMaInputStream.cpp
)

if(WITH_LIBXML)
 set(FASTPROT_XML_SRCS apps/fastprot/XmlInputStream.cpp)
endif(WITH_LIBXML)

if(BUILD_WITH_MPI)
  # fastprot_mpi is the fourth app migrated off gengetopt onto CLI11
  # (gengetopt_migration_plan.md, Phase D) - no generated cmdline.c to
  # list here anymore.
  set(FASTPROT_MPI_SRCS apps/fastprot_mpi/main.cpp
  apps/fastprot_mpi/FastaInputStream.cpp
  apps/fastprot_mpi/DataOutputStream.cpp
  apps/fastprot_mpi/XmlOutputStream.cpp
  apps/fastprot_mpi/PhylipMaInputStream.cpp
  apps/fastprot_mpi/ProtDistCalc.cpp
  apps/fastprot_mpi/ModelMatrix.cpp
  apps/fastprot_mpi/ExpectedDistance.cpp
  apps/fastprot_mpi/Matrix.cpp
  apps/fastprot_mpi/MaximumLikelihood.cpp
  apps/fastprot_mpi/ProtSeqUtils.cpp
  )
  # fastprot_mpi is out of this plan's scope (project_layout_plan.md:
  # deferred, to be folded into fastprot on a future request) and was
  # never build-tested here (no MPI in this environment - configuring
  # with -DBUILD_WITH_MPI=ON already fails before reaching this point).
  # It used to reuse apps/fastprot/PhylipDmOutputStream.cpp and
  # BinaryDmOutputStream.cpp directly in this source list; both moved
  # into the shared fastphylo library as part of Phase C's io/
  # consolidation, so those two lines were removed above rather than
  # left pointing at deleted files. fastprot_mpi links fastphylo (below)
  # and so gets io::PhylipDmOutputStream/BinaryDmOutputStream from
  # there now - unverified whether fastprot_mpi's own (separate,
  # unmerged) DataOutputStream.hpp is still compatible with them; needs
  # checking when fastprot_mpi work is actually picked up.
  if(WITH_LIBXML)
    set(FASTPROT_MPI_XML_SRCS apps/fastprot_mpi/XmlInputStream.cpp)
  endif(WITH_LIBXML)
endif(BUILD_WITH_MPI)


foreach(i ${FASTPHYLO_SPECIAL_SRCS})
  set_source_files_properties(${i} PROPERTIES COMPILE_FLAGS "-DNDEBUG -fno-branch-count-reg -fomit-frame-pointer")
endforeach(i)

# Modernization Phase 5 (modernization_plan.md): target_include_directories
# replaces the old global include_directories() calls (every translation
# unit in the project used to get every directory below, regardless of
# whether it needed it). fastphylo's own directories are PUBLIC because
# fastdist/fnj/fastprot/fastprot_mpi #include core-library headers by
# bare filename (e.g. "SequenceTree.hpp", "NeighborJoining.hpp") -
# confirmed by grep before this rewrite, not assumed - so consumers
# linking fastphylo need these on their include path too, which PUBLIC
# gives them transitively. (DNA_b128/ and sequence_likelihood/ - once
# also part of this bare-filename include path - are gone: their
# surviving content moved into include/fastphylo/dna/ in
# project_layout_plan.md, and the rest was confirmed-dead and deleted
# in lint_plan.md's Phase 0.)
# Matrix.cpp/MaximumLikelihood.cpp use Eigen::Matrix<double,20,20> for
# the ML hot path's dense 20x20 operations (planning/
# fastprot_ml_speedup_implementation_plan.md) - header-only, so unlike
# BLAS/LAPACK (extern "C" declarations, no header needed until final
# link) this needs to be on fastphylo's own include path at compile
# time, not deferred to fastprot's target further down.
# No version pin (was "3.4 REQUIRED", the version available when this
# dependency was added): Homebrew's eigen formula has since moved its
# default/stable version to 5.0.1, and find_package(Eigen3 3.4)'s
# SameMajorVersion compatibility check rejects that outright (CI
# failure, 2026-08-11 - the first time this dependency was ever
# actually exercised on a real macOS runner, since pushes to this
# branch never triggered the workflow, only this PR did). Verified
# directly (full rebuild + ctest 6/6) against the Eigen 5.0.1 bottle:
# the small, stable subset used here (fixed-size Matrix, EigenSolver,
# SelfAdjointEigenSolver, basic array ops) is unaffected by whatever
# changed 3.4->5.0. No genuine technical floor is needed - these APIs
# have been stable since Eigen 3.0 (2010) - so accepting whatever
# Eigen3 a platform's package manager currently provides (3.4.x via
# Ubuntu's libeigen3-dev/vcpkg, 5.0.1 via current Homebrew) is more
# robust than re-pinning to today's latest and hitting this again next
# time a package manager moves on.
find_package(Eigen3 REQUIRED)

add_library(fastphylo STATIC ${FASTPHYLO_SRCS} ${FASTPHYLO_SPECIAL_SRCS} ${FASTPHYLO_IO_SRCS} ${FASTPHYLO_IO_XML_SRCS} ${FASTPHYLO_PROTEIN_SRCS})
# Modernization Phase B (project_layout_plan.md): the public API lives
# under include/fastphylo/{core,dna,io,protein}/ now, consumed as
# #include "fastphylo/<bucket>/X.hpp" - this is the one include root
# that needs to be on the path for that, replacing the four
# scattered src/c++-relative directories this used to list.
target_include_directories(fastphylo PUBLIC
  ${CMAKE_CURRENT_BINARY_DIR}
  ${CMAKE_SOURCE_DIR}/include
  ${FASTPHYLO_EXTRA_INCLUDE_DIRS}
)
# fastphylo's own XmlInputStream.cpp (when WITH_LIBXML) needs libxml2's
# headers to compile; PUBLIC so consumers linking fastphylo (which
# already separately link ${FASTPHYLO_XML_LIBS} themselves) get it
# transitively too rather than needing to know to add it twice.
target_link_libraries(fastphylo PUBLIC ${FASTPHYLO_XML_LIBS})
# MatrixExpm's private members (Matrix.hpp) are Eigen types now - the
# public Matrix class and MatrixExpm's own public methods are
# unchanged (still Matrix in, Matrix out), but Matrix.hpp itself
# #includes <Eigen/Dense> to declare that private state, so anything
# that #includes Matrix.hpp needs Eigen's include path too - PUBLIC,
# same reasoning as the XML libs above.
target_link_libraries(fastphylo PUBLIC Eigen3::Eigen)

find_package(Threads REQUIRED)

# libm is a separate library needing explicit linking on Unix (sin/cos/
# sqrt/etc historically lived outside libc there); on Windows every
# math function is already part of the standard CRT and there is no
# "m.lib" to link against at all - target_link_libraries(... m ...)
# there fails at link time with "cannot open input file 'm.lib'".
if(UNIX)
  set(FASTPHYLO_MATH_LIB m)
endif()

add_executable(fastdist ${FASTDIST_SRCS} ${FASTDIST_XML_SRCS})
target_include_directories(fastdist PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastdist
  ${CMAKE_CURRENT_SOURCE_DIR}/apps
  ${CMAKE_SOURCE_DIR}/third_party
)
target_link_libraries(fastdist PRIVATE fastphylo ${FASTPHYLO_MATH_LIB} ${FASTPHYLO_XML_LIBS})

add_executable(fnj ${FNJ_SRCS} ${FNJ_XML_SRCS})
target_include_directories(fnj PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/apps/fnj
  ${CMAKE_CURRENT_SOURCE_DIR}/apps
  ${CMAKE_SOURCE_DIR}/third_party
)
target_link_libraries(fnj PRIVATE fastphylo ${FASTPHYLO_MATH_LIB} ${FASTPHYLO_XML_LIBS})


find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_executable(fastprot ${FASTPROT_SRCS} ${FASTPROT_XML_SRCS})
target_include_directories(fastprot PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot
  ${CMAKE_CURRENT_SOURCE_DIR}/apps
  ${CMAKE_SOURCE_DIR}/third_party
)
target_link_libraries(fastprot PRIVATE fastphylo ${FASTPHYLO_MATH_LIB} ${FASTPHYLO_XML_LIBS} BLAS::BLAS LAPACK::LAPACK Threads::Threads)


if(BUILD_WITH_MPI)
  add_executable(fastprot_mpi ${FASTPROT_MPI_SRCS} ${FASTPROT_MPI_XML_SRCS})
  target_include_directories(fastprot_mpi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot_mpi
    ${CMAKE_CURRENT_SOURCE_DIR}/apps
    ${CMAKE_SOURCE_DIR}/third_party
  )
  target_link_libraries(fastprot_mpi PRIVATE fastphylo ${FASTPHYLO_MATH_LIB} ${FASTPHYLO_XML_LIBS} BLAS::BLAS LAPACK::LAPACK Threads::Threads MPI::MPI_CXX)
endif(BUILD_WITH_MPI)


install(TARGETS fastdist DESTINATION bin)
install(TARGETS fnj DESTINATION bin)
install(TARGETS fastprot DESTINATION bin)
if(BUILD_WITH_MPI)
  install(TARGETS fastprot_mpi DESTINATION bin)
endif(BUILD_WITH_MPI)

# Correctness tests for the speed2026a protein-code speedup (plan.md).
# tests/ (top-level, Layout Phase E - was src/c++/code_tests/) as a
# whole isn't wired into this build (see phase0_audit.md); these two
# are new and specific to this project, so they're wired in directly
# rather than expanding scope to the rest of tests/.
enable_testing()

set(FASTPHYLO_TESTS_DIR ${CMAKE_SOURCE_DIR}/tests)

# ProtSeqCompare.cpp needs FASTPHYLO_EXTRA_INCLUDE_DIRS (simde) directly:
# these three targets compile it as a plain extra source rather than
# linking the fastphylo library target, so they don't inherit it via
# fastphylo's PUBLIC target_include_directories. Likewise need
# ${CMAKE_SOURCE_DIR}/include directly (Layout Phase C2) now that
# ProtSeqCode.hpp/ProtSeqCompare.hpp live under include/fastphylo/protein/.
add_executable(ProtSeqCode_test ${FASTPHYLO_TESTS_DIR}/ProtSeqCode_test.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCode.cpp)
target_include_directories(ProtSeqCode_test PRIVATE ${CMAKE_SOURCE_DIR}/include ${FASTPHYLO_EXTRA_INCLUDE_DIRS})
add_test(NAME ProtSeqCode_test COMMAND ProtSeqCode_test)

add_executable(ProtSeqCompare_test ${FASTPHYLO_TESTS_DIR}/ProtSeqCompare_test.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCode.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCompare.cpp)
target_include_directories(ProtSeqCompare_test PRIVATE ${CMAKE_SOURCE_DIR}/include ${FASTPHYLO_EXTRA_INCLUDE_DIRS})
add_test(NAME ProtSeqCompare_test COMMAND ProtSeqCompare_test)

# phylip_reader_consolidation_plan.md, Phase C: SequenceTree::
# mapSequencesOntoTree(istream&) has no live caller anywhere in the
# codebase - its only two callers (Simulator.cpp and
# generating_block_transfer_data/Generate_sequence_data.cpp) were
# themselves unbuilt, hardcoded-external-tool wrapper code and were
# deleted outright, so this is the only way to verify this code path
# at all. Links the real fastphylo library target directly (unlike
# ProtSeqCode_test/ProtSeqCompare_test's narrow hand-picked source
# lists above) since SequenceTree.cpp pulls in most of the core
# library.
add_executable(SequenceTree_PhylipReader_test ${FASTPHYLO_TESTS_DIR}/SequenceTree_PhylipReader_test.cpp)
target_link_libraries(SequenceTree_PhylipReader_test PRIVATE fastphylo)
add_test(NAME SequenceTree_PhylipReader_test COMMAND SequenceTree_PhylipReader_test)

# fnj_binary_input_gap: BinaryInputStream::readDM(StrDblMatrix&, ...) was
# an unimplemented stub. BinaryInputStream.cpp/.hpp live under
# apps/fnj/ (not the shared fastphylo library, same as
# PhylipDmInputStream/XmlInputStream), so this test compiles that one
# source file directly - same pattern as ProtSeqCode_test/
# ProtSeqCompare_test above - alongside linking fastphylo for
# StrDblMatrix/BinaryDmOutputStream.
add_executable(BinaryDmIO_test ${FASTPHYLO_TESTS_DIR}/BinaryDmIO_test.cpp apps/fnj/BinaryInputStream.cpp)
target_include_directories(BinaryDmIO_test PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/apps/fnj)
target_link_libraries(BinaryDmIO_test PRIVATE fastphylo)
add_test(NAME BinaryDmIO_test COMMAND BinaryDmIO_test)

# likelihood_calc()'s safeguarded Newton-Raphson search, checked
# against real biological data (examples/globin_family.fasta) across
# all 5 ML models - see planning/fastprot_ml_speedup_investigation_
# plan.md's "Q3 results" for the bug this guards against (the previous
# finite-difference solver could return unverified, badly wrong
# answers on real, as opposed to narrow-divergence-synthetic, data).
# Links fastphylo directly (needs ModelMatrix/Matrix/MaximumLikelihood/
# ProtSeqCode together) rather than hand-picking sources.
add_executable(MaximumLikelihood_test ${FASTPHYLO_TESTS_DIR}/MaximumLikelihood_test.cpp)
target_link_libraries(MaximumLikelihood_test PRIVATE fastphylo BLAS::BLAS LAPACK::LAPACK)
target_compile_definitions(MaximumLikelihood_test PRIVATE GLOBIN_FIXTURE_PATH="${CMAKE_SOURCE_DIR}/examples/globin_family.fasta")
add_test(NAME MaximumLikelihood_test COMMAND MaximumLikelihood_test)

# Model data sanity check (row-sums-to-zero, frequencies-sum-to-one) -
# added after finding JTT's Q(0,0) was hardcoded 100x too large, a
# data-entry bug MaximumLikelihood_test's solver-self-consistency check
# couldn't catch (see ModelMatrix_test.cpp's header comment).
add_executable(ModelMatrix_test ${FASTPHYLO_TESTS_DIR}/ModelMatrix_test.cpp)
target_link_libraries(ModelMatrix_test PRIVATE fastphylo BLAS::BLAS LAPACK::LAPACK)
add_test(NAME ModelMatrix_test COMMAND ModelMatrix_test)

# Phase 5 performance benchmark (plan.md). Not an add_test() - it prints
# timing data, not a pass/fail result; see benchmarks/run_benchmarks.sh
# for the reproducible one-command entry point.
add_executable(bench_primitives ../../benchmarks/bench_primitives.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCode.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCompare.cpp)
target_include_directories(bench_primitives PRIVATE ${CMAKE_SOURCE_DIR}/include ${FASTPHYLO_EXTRA_INCLUDE_DIRS})
