# Modernization Phase 0 (modernization_plan.md): was VERSION 3.2
# (~2015). Phase 5 bumped 3.16 -> 3.18 (2020): still a conservative,
# widely-available floor (current Linux distros/macOS/CI images), and
# the version FindBLAS/FindLAPACK started reliably providing the
# BLAS::BLAS/LAPACK::LAPACK imported targets src/c++/CMakeLists.txt
# now relies on (FindLibXml2's LibXml2::LibXml2 has been available
# since 3.12, so it isn't the constraint here).
cmake_minimum_required(VERSION 3.18)
project(fastphylo)

# C++17 floor project-wide (modernization_plan.md). Portable across
# compilers (GCC/Clang/MSVC), unlike the old hand-appended "--std=c++11"
# GCC/Clang-only flag string this replaces (see src/c++/CMakeLists.txt).
set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

# CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR: this repo's own
# CMake floor is 3.18, no PROJECT_IS_TOP_LEVEL (added 3.21) - this is
# the portable equivalent, used here and again further down to guard
# enable_testing(), the docs subdirectory, and the CPack block: all
# only make sense when FastPhylo itself is the project being
# configured, not when it's pulled in as a subdirectory dependency
# (fastphylo-py's planned git-submodule consumption of just the
# fastphylo library target - see planning/ for that integration
# plan). PACKAGE_VERSION derivation just below stays unconditional -
# config.h.cmake's substitution (src/c++/CMakeLists.txt) needs it
# either way, since that's part of the fastphylo library target
# itself, not just the apps/packaging.
if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)
  enable_testing() # code_tests/ProtSeqCode_test and ProtSeqCompare_test as of speed2026a
endif()

# github_actions_release_builds_plan.md's Phase A/E: the tag drives the
# version, not the other way round - no checked-in number to keep in
# sync with a release tag by hand, so no way for the two to drift
# apart. Three ways PACKAGE_VERSION can end up set, in priority order:
#
# 1. Passed explicitly via -DPACKAGE_VERSION=X.Y.Z - what the
#    tag-triggered release workflow does, reading the pushed git tag
#    (v1.2.0 -> 1.2.0) and overriding whatever's checked in here.
# 2. Derived live via `git describe --tags --always --dirty` when a
#    .git directory is present (every ordinary dev/CI build - push/PR
#    CI, a local `cmake && make`) - informative strings like
#    "1.0.3-42-gabc1234", not a single frozen number nobody updates.
# 3. The static fallback below, only reached when neither applies -
#    e.g. a `make package_source` tarball extracted and built later,
#    with no git history at all to describe.
if(NOT DEFINED PACKAGE_VERSION)
  find_package(Git QUIET)
  if(GIT_FOUND AND EXISTS "${PROJECT_SOURCE_DIR}/.git")
    execute_process(
      COMMAND ${GIT_EXECUTABLE} describe --tags --always --dirty
      WORKING_DIRECTORY ${PROJECT_SOURCE_DIR}
      OUTPUT_VARIABLE PACKAGE_VERSION
      OUTPUT_STRIP_TRAILING_WHITESPACE
      RESULT_VARIABLE GIT_DESCRIBE_RESULT
      ERROR_QUIET
    )
    if(NOT GIT_DESCRIBE_RESULT EQUAL 0)
      unset(PACKAGE_VERSION)
    endif()
  endif()
endif()
if(NOT PACKAGE_VERSION)
  set(PACKAGE_VERSION "1.0.3") # static fallback - no git info available
endif()
message(STATUS "fastphylo PACKAGE_VERSION: ${PACKAGE_VERSION}")

option(BUILD_WITH_MPI "Build the MPI version of fastprot" OFF)
option(BUILD_DOCBOOK "Build Docbook html documentation" OFF)
option(WITH_LIBXML "Build xml functionality. This adds a dependeny for libxml" ON)
# The STATIC option (self-built gengetopt/libxml2/zlib from 2009-2011-era
# FTP-hosted tarballs, for portable win32/macOS binary packaging) was
# removed in lint_plan.md's Phase 6 - never exercised by this project's
# regular dev/CI builds, and superseded by a from-scratch, modern
# cross-platform release-build plan (see project memory/a future plan
# doc) rather than kept around as dead weight in the meantime.


add_subdirectory(src/c++)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_CURRENT_SOURCE_DIR)

add_subdirectory(docs)



# After running cmake, we can generate source tarballs ( in tar.gz and zip ) by running
#    make package_source
#
# Use
#    make package
# to get static binaries in a .tar.gz file!


#set(CPACK_GENERATOR RPM;DEB)

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxml2 (>= 2.6.16 )")


# Set directly rather than via CPACK_PACKAGE_VERSION_MAJOR/MINOR/PATCH
# - PACKAGE_VERSION (above) isn't always a clean MAJOR.MINOR.PATCH
# triplet (a git-describe dev build looks like "1.0.3-42-gabc1234").
set(CPACK_PACKAGE_VERSION ${PACKAGE_VERSION})
set(CPACK_SOURCE_GENERATOR TGZ;ZIP)
set(CPACK_PACKAGE_CONTACT "arve@su.se")
set(CPACK_PACKAGE_DESCRIPTION_SUMMARY "Fastphylo is software project containing the implementations of the algorithms \"Fast Computation of Distance Estimators\", \"Fast Neighbor Joining\", and more.")
set(CPACK_PACKAGE_VENDOR "Isaac Elias, Jensl Lagergren, Erik Sj�lund, Lars Arvestad, Mehmood Alam Khan, et al")
set(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README.md")
set(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
set(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
set(CPACK_RPM_PACKAGE_LICENSE "MIT License")
set(CPACK_RPM_PACKAGE_RELEASE "1")


set(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
include(CPack)

endif()
