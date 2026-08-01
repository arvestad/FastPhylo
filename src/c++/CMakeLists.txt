find_program(GENGETOPT gengetopt)
if(NOT GENGETOPT)
 set(GENGETOPTDEP "${CMAKE_CURRENT_BINARY_DIR}/external/install/bin/gengetopt")
if(WIN32)
 set(GENGETOPTDEP "${GENGETOPTDEP}.exe")
endif(WIN32)
 set(GENGETOPT "${GENGETOPTDEP}")
endif(NOT GENGETOPT)

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


if(STATIC)
 set(LIBXML2_INCLUDE_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2")
 set(LIBXML2_LIBRARIES "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libxml2.a")

# The next line seems to be need for macosx
 set(ZLIB_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a")
endif(STATIC)

set(GENGETOPT_VERSION 2.22.2)
set(LIBXML2_VERSION 2.7.2)
set(ZLIB_VERSION 1.2.3)

set(DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/downloads" INTERNAL)
set(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/build INTERNAL)
set(INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/install INTERNAL)


file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/download)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/build)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/install)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt)
file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt)
if(BUILD_WITH_MPI)
  file(MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt)
endif()

set(ARCH_FLAGS "-msse2")

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_CURRENT_BINARY_DIR}/script.nsi" DEPENDS fastdist fnj fastprot fastprot_mpi)

  add_custom_target(win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe")

  # STATIC's self-built gengetopt/libxml2/zlib toolchain: 2009-2011-era
  # library versions fetched from FTP mirrors (ftp.gnu.org, xmlsoft.org)
  # of unverified current availability. Not exercised by this project's
  # regular dev/CI builds (STATIC defaults OFF) and not touched by this
  # modernization pass beyond command-case consistency - flagged as a
  # candidate for its own future audit rather than silently assumed
  # working.
  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/gengetopt-${GENGETOPT_VERSION}.tar.gz" COMMAND wget --directory-prefix=external/download "ftp://ftp.gnu.org/gnu/gengetopt/gengetopt-${GENGETOPT_VERSION}.tar.gz")

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/libxml2-${LIBXML2_VERSION}.tar.gz" COMMAND wget --directory-prefix=external/download "ftp://xmlsoft.org/libxml2/libxml2-${LIBXML2_VERSION}.tar.gz")

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/zlib-${ZLIB_VERSION}.tar.gz" COMMAND wget --directory-prefix=external/download "http://www.gzip.org/zlib/zlib-${ZLIB_VERSION}.tar.gz")

  add_custom_command(OUTPUT "${GENGETOPTDEP}" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_gengetopt.sh" ${GENGETOPT_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install" DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/gengetopt-${GENGETOPT_VERSION}.tar.gz" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build")

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libxml2.a" "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_libxml2.sh" ${LIBXML2_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install" DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/libxml2-${LIBXML2_VERSION}.tar.gz" "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build")

  add_custom_command(OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_zlib.sh" ${ZLIB_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install" DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/zlib-${ZLIB_VERSION}.tar.gz" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build")


add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt/fastdist_gengetopt.c ${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt/fastdist_gengetopt.h
                     COMMAND ${GENGETOPT} -i ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastdist/gengetopt/cmdline.ggo --output-dir=${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastdist/gengetopt/cmdline.ggo ${GENGETOPTDEP}
                     VERBATIM)


add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt/fnj_gengetopt.c ${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt/fnj_gengetopt.h
                     COMMAND ${GENGETOPT} -i ${CMAKE_CURRENT_SOURCE_DIR}/apps/fnj/gengetopt/cmdline.ggo --output-dir=${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/apps/fnj/gengetopt/cmdline.ggo ${GENGETOPTDEP}
                     VERBATIM)

add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt/fastprot_gengetopt.c ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt/fastprot_gengetopt.h
                     COMMAND ${GENGETOPT} -i ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot/gengetopt/cmdline.ggo --output-dir=${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot/gengetopt/cmdline.ggo ${GENGETOPTDEP}
                     VERBATIM)

if(BUILD_WITH_MPI)
  add_custom_command(OUTPUT ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt/fastprot_mpi_gengetopt.c ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt/fastprot_mpi_gengetopt.h
                     COMMAND ${GENGETOPT} -i ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot_mpi/gengetopt/cmdline.ggo --output-dir=${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt
                     DEPENDS ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot_mpi/gengetopt/cmdline.ggo ${GENGETOPTDEP}
                     VERBATIM)
endif()

# Modernization Phase 5: STATIC's self-built libxml2 bypasses
# find_package's system detection by design (LIBXML2_INCLUDE_DIR/
# LIBXML2_LIBRARIES set above), so it can't use the LibXml2::LibXml2
# imported target directly - wrap it in a small INTERFACE library
# instead, so every consumer below links one FASTPHYLO_XML_LIBS
# regardless of STATIC. Preserves the pre-existing STATIC-vs-find_package
# ordering/behavior exactly; not attempting to fix whatever interaction
# there may be between the two (out of this pass's scope, STATIC is
# untested here).
if(WITH_LIBXML)
  if(STATIC)
    add_library(fastphylo_libxml2 INTERFACE)
    target_include_directories(fastphylo_libxml2 INTERFACE ${LIBXML2_INCLUDE_DIR})
    target_link_libraries(fastphylo_libxml2 INTERFACE ${LIBXML2_LIBRARIES})
    set(FASTPHYLO_XML_LIBS fastphylo_libxml2)
  else()
    find_package(LibXml2)
    if(TARGET LibXml2::LibXml2)
      set(FASTPHYLO_XML_LIBS LibXml2::LibXml2)
    endif()
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
${FASTPHYLO_CORE_DIR}/Exception.cpp
${FASTPHYLO_CORE_DIR}/InitAndPrintOn_utils.cpp
${FASTPHYLO_CORE_DIR}/Object.cpp
${FASTPHYLO_CORE_DIR}/Sequence.cpp
${FASTPHYLO_CORE_DIR}/SequenceTree.cpp
${FASTPHYLO_CORE_DIR}/SequenceTree_MostParsimonious.cpp
${FASTPHYLO_CORE_DIR}/Simulator.cpp
${FASTPHYLO_CORE_DIR}/arg_utils_ext.cpp
${FASTPHYLO_CORE_DIR}/file_utils.cpp
${FASTPHYLO_CORE_DIR}/stl_utils.cpp
${FASTPHYLO_DNA_DIR}/DNA_b128_String.cpp
${FASTPHYLO_DNA_DIR}/Sequences2DistanceMatrix.cpp
${FASTPHYLO_DNA_DIR}/LeastSquaresFit.cpp
${FASTPHYLO_DNA_DIR}/NeighborJoining.cpp
${FASTPHYLO_DNA_DIR}/Kimura2parameter.cpp
${FASTPHYLO_DNA_DIR}/TamuraNei.cpp
${FASTPHYLO_DNA_DIR}/ambiguity_nucleotide.cpp
${FASTPHYLO_DNA_DIR}/dna_pairwise_sequence_likelihood.cpp
${FASTPHYLO_DNA_DIR}/string_compare.cpp
${FASTPHYLO_CORE_DIR}/DistanceMatrix.cpp
${FASTPHYLO_CORE_DIR}/FloatDistanceMatrix.cpp
${FASTPHYLO_CORE_DIR}/DistanceRow.cpp
${FASTPHYLO_CORE_DIR}/arg_utils.c
${FASTPHYLO_CORE_DIR}/std_c_utils.c
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

if(STATIC AND WITH_LIBXML)
 set(LIBXML_EXTRA_COMPILE_FLAGS " -DLIBXML_STATIC ")
endif(STATIC AND WITH_LIBXML)

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

set(FASTDIST_SRCS apps/fastdist/main.cpp
apps/fastdist/PhylipMaInputStream.cpp
apps/fastdist/FastaInputStream.cpp
${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt/fastdist_gengetopt.c
)

if(WITH_LIBXML)
 if(STATIC)
  set(FASTDIST_XML_SRCS apps/fastdist/XmlInputStream.cpp "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h")
 else(STATIC)
  set(FASTDIST_XML_SRCS apps/fastdist/XmlInputStream.cpp)
 endif(STATIC)
endif(WITH_LIBXML)


set(FNJ_SRCS apps/fnj/main.cpp
apps/fnj/DataInputStream.cpp
apps/fnj/DataOutputStream.cpp
apps/fnj/TreeTextOutputStream.cpp
apps/fnj/XmlOutputStream.cpp
apps/fnj/PhylipDmInputStream.cpp
apps/fnj/BinaryInputStream.cpp
${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt/fnj_gengetopt.c
)

if(WITH_LIBXML)
 if(STATIC)
  set(FNJ_XML_SRCS apps/fnj/XmlInputStream.cpp "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h")
 else(STATIC)
  set(FNJ_XML_SRCS apps/fnj/XmlInputStream.cpp)
 endif(STATIC)
endif(WITH_LIBXML)

set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt/fnj_gengetopt.c PROPERTIES GENERATED true)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt/fnj_gengetopt.h PROPERTIES GENERATED true)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt/fastdist_gengetopt.c PROPERTIES GENERATED true)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt/fastdist_gengetopt.h PROPERTIES GENERATED true)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt/fastprot_gengetopt.c PROPERTIES GENERATED true)
set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt/fastprot_gengetopt.h PROPERTIES GENERATED true)
if(BUILD_WITH_MPI)
  set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt/fastprot_mpi_gengetopt.c PROPERTIES GENERATED true)
  set_source_files_properties(${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt/fastprot_mpi_gengetopt.h PROPERTIES GENERATED true)
endif()


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

set(FASTPROT_SRCS apps/fastprot/main.cpp
apps/fastprot/FastaInputStream.cpp
apps/fastprot/PhylipMaInputStream.cpp
${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt/fastprot_gengetopt.c
)

if(WITH_LIBXML)
 if(STATIC)
  set(FASTPROT_XML_SRCS apps/fastprot/XmlInputStream.cpp "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h")
 else(STATIC)
  set(FASTPROT_XML_SRCS apps/fastprot/XmlInputStream.cpp)
 endif(STATIC)
endif(WITH_LIBXML)

if(BUILD_WITH_MPI)
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
  ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt/fastprot_mpi_gengetopt.c
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
    if(STATIC)
      set(FASTPROT_MPI_XML_SRCS apps/fastprot_mpi/XmlInputStream.cpp "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h")
    else(STATIC)
      set(FASTPROT_MPI_XML_SRCS apps/fastprot_mpi/XmlInputStream.cpp)
    endif(STATIC)
  endif(WITH_LIBXML)
endif(BUILD_WITH_MPI)


foreach(i ${FASTPHYLO_SPECIAL_SRCS})
  set_source_files_properties(${i} PROPERTIES COMPILE_FLAGS "-DNDEBUG -fno-branch-count-reg -fomit-frame-pointer")
endforeach(i)

# Modernization Phase 5 (modernization_plan.md): target_include_directories
# replaces the old global include_directories() calls (every translation
# unit in the project used to get every directory below, regardless of
# whether it needed it). fastphylo's own directories are PUBLIC because
# fastdist/fnj/fastprot/fastprot_mpi #include core-library and DNA_b128/
# distance_methods/sequence_likelihood headers by bare filename (e.g.
# "SequenceTree.hpp", "NeighborJoining.hpp") - confirmed by grep before
# this rewrite, not assumed - so consumers linking fastphylo need these
# on their include path too, which PUBLIC gives them transitively.
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

find_package(Threads REQUIRED)

add_executable(fastdist ${FASTDIST_SRCS} ${FASTDIST_XML_SRCS})
target_include_directories(fastdist PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastdist
  ${CMAKE_CURRENT_BINARY_DIR}/apps/fastdist/gengetopt
)
target_link_libraries(fastdist PRIVATE fastphylo m ${FASTPHYLO_XML_LIBS} ${ZLIB_LIBRARY})

if(STATIC)
 set_property(TARGET fastdist PROPERTY LINK_SEARCH_END_STATIC ON)
endif(STATIC)

add_executable(fnj ${FNJ_SRCS} ${FNJ_XML_SRCS})
target_include_directories(fnj PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/apps/fnj
  ${CMAKE_CURRENT_BINARY_DIR}/apps/fnj/gengetopt
)
target_link_libraries(fnj PRIVATE fastphylo m ${FASTPHYLO_XML_LIBS} ${ZLIB_LIBRARY})

if(STATIC)
 set_property(TARGET fnj PROPERTY LINK_SEARCH_END_STATIC ON)
endif(STATIC)


find_package(BLAS REQUIRED)
find_package(LAPACK REQUIRED)
add_executable(fastprot ${FASTPROT_SRCS} ${FASTPROT_XML_SRCS})
target_include_directories(fastprot PRIVATE
  ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot
  ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot/gengetopt
)
target_link_libraries(fastprot PRIVATE fastphylo m ${FASTPHYLO_XML_LIBS} BLAS::BLAS LAPACK::LAPACK ${ZLIB_LIBRARY} Threads::Threads)

if(STATIC)
 set_property(TARGET fastprot PROPERTY LINK_SEARCH_END_STATIC ON)
endif(STATIC)


if(BUILD_WITH_MPI)
  add_executable(fastprot_mpi ${FASTPROT_MPI_SRCS} ${FASTPROT_MPI_XML_SRCS})
  target_include_directories(fastprot_mpi PRIVATE
    ${CMAKE_CURRENT_SOURCE_DIR}/apps/fastprot_mpi
    ${CMAKE_CURRENT_BINARY_DIR}/apps/fastprot_mpi/gengetopt
  )
  target_link_libraries(fastprot_mpi PRIVATE fastphylo m ${FASTPHYLO_XML_LIBS} BLAS::BLAS LAPACK::LAPACK ${ZLIB_LIBRARY} Threads::Threads MPI::MPI_CXX)

  if(STATIC)
    set_property(TARGET fastprot_mpi PROPERTY LINK_SEARCH_END_STATIC ON)
  endif(STATIC)
endif(BUILD_WITH_MPI)


install(TARGETS fastdist DESTINATION bin)
install(TARGETS fnj DESTINATION bin)
install(TARGETS fastprot DESTINATION bin)
if(BUILD_WITH_MPI)
  install(TARGETS fastprot_mpi DESTINATION bin)
endif(BUILD_WITH_MPI)

# Correctness tests for the speed2026a protein-code speedup (plan.md).
# code_tests/ as a whole isn't wired into this build (see
# phase0_audit.md); these two are new and specific to this project, so
# they're wired in directly rather than expanding scope to the rest of
# code_tests/.
enable_testing()

# ProtSeqCompare.cpp needs FASTPHYLO_EXTRA_INCLUDE_DIRS (simde) directly:
# these three targets compile it as a plain extra source rather than
# linking the fastphylo library target, so they don't inherit it via
# fastphylo's PUBLIC target_include_directories. Likewise need
# ${CMAKE_SOURCE_DIR}/include directly (Layout Phase C2) now that
# ProtSeqCode.hpp/ProtSeqCompare.hpp live under include/fastphylo/protein/.
add_executable(ProtSeqCode_test code_tests/ProtSeqCode_test.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCode.cpp)
target_include_directories(ProtSeqCode_test PRIVATE ${CMAKE_SOURCE_DIR}/include ${FASTPHYLO_EXTRA_INCLUDE_DIRS})
add_test(NAME ProtSeqCode_test COMMAND ProtSeqCode_test)

add_executable(ProtSeqCompare_test code_tests/ProtSeqCompare_test.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCode.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCompare.cpp)
target_include_directories(ProtSeqCompare_test PRIVATE ${CMAKE_SOURCE_DIR}/include ${FASTPHYLO_EXTRA_INCLUDE_DIRS})
add_test(NAME ProtSeqCompare_test COMMAND ProtSeqCompare_test)

# Phase 5 performance benchmark (plan.md). Not an add_test() - it prints
# timing data, not a pass/fail result; see benchmarks/run_benchmarks.sh
# for the reproducible one-command entry point.
add_executable(bench_primitives ../../benchmarks/bench_primitives.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCode.cpp ${FASTPHYLO_PROTEIN_DIR}/ProtSeqCompare.cpp)
target_include_directories(bench_primitives PRIVATE ${CMAKE_SOURCE_DIR}/include ${FASTPHYLO_EXTRA_INCLUDE_DIRS})
