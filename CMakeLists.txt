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

enable_testing() # code_tests/ProtSeqCode_test and ProtSeqCompare_test as of speed2026a

set(CPACK_PACKAGE_VERSION_MAJOR "1")
set(CPACK_PACKAGE_VERSION_MINOR "0")
set(CPACK_PACKAGE_VERSION_PATCH "3")
set(PACKAGE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})

option(BUILD_WITH_MPI "Build the MPI version of fastprot" OFF)
option(BUILD_DOCBOOK "Build Docbook html documentation" OFF)
option(WITH_LIBXML "Build xml functionality. This adds a dependeny for libxml" ON)
option(STATIC "Build static binaries. This is mainly intended for the win32-package and the MACOSX-bundle. This is will download gengetopt, lixml2, zlib, and build them. The libxml2-library (and therefore indirectly zlib) will be statically linked in" OFF)


add_subdirectory(src/c++)
add_subdirectory(src/docbook)



# After running cmake, we can generate source tarballs ( in tar.gz and zip ) by running
#    make package_source
#
# Use
#    make package
# to get static binaries in a .tar.gz file!


#set(CPACK_GENERATOR RPM;DEB)

set(CPACK_DEBIAN_PACKAGE_DEPENDS "libxml2 (>= 2.6.16 )")


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
