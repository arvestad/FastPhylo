cmake_minimum_required(VERSION 2.6.0)
project(fastphylo)

ENABLE_TESTING()

SET(CPACK_PACKAGE_VERSION_MAJOR "1")
SET(CPACK_PACKAGE_VERSION_MINOR "0")
SET(CPACK_PACKAGE_VERSION_PATCH "2")
set(PACKAGE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})

# OpenMPI support
FIND_PACKAGE(MPI REQUIRED)
INCLUDE(CMakeForceCompiler)
#mehmood's changes here
IF(CMAKE_HOST_APPLE)
	CMAKE_FORCE_CXX_COMPILER(mpic++ "MPI C++ Compiler for Apple")
ELSEIF(CMAKE_HOST_UNIX)
	CMAKE_FORCE_CXX_COMPILER(mpiCC "MPI C++ Compiler for Linux")
ENDIF()


# Ignore annoying deprecation warnings for now
#SET(CMAKE_CXX_FLAGS "-Wno-deprecated")

OPTION(BUILD_DOCBOOK "Build Docbook html documentation" OFF )
OPTION(WITH_LIBXML "Build xml functionality. This adds a dependeny for libxml" ON )
OPTION(STATIC "Build static binaries. This is mainly intended for the win32-package and the MACOSX-bundle. This is will download gengetopt, lixml2, zlib, and build them. The libxml2-library (and therefore indirectly zlib) will be statically linked in" OFF )


add_subdirectory(src/c++)
add_subdirectory(src/docbook)

# after running cmake, we can generate source tarballs ( in tar.gz and zip ) by running
# make package_source
#

#SET(CPACK_GENERATOR RPM;DEB)

SET(CPACK_DEBIAN_PACKAGE_DEPENDS "libxml2 (>= 2.6.16 )")


SET(CPACK_SOURCE_GENERATOR TGZ;ZIP)
SET(CPACK_PACKAGE_CONTACT "arve@su.se") 
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "FastPhylo is software project containing the implementations of the algorithms \"Fast Computation of Distance Estimators\", \"Fast Neighbor Joining\", and more.")
SET(CPACK_PACKAGE_VENDOR "Isaac Elias, Jensl Lagergren, Erik Sjölund, Lars Arvestad, Mehmood Alam Khan, et al")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
SET(CPACK_RPM_PACKAGE_LICENSE "MIT License")
SET(CPACK_RPM_PACKAGE_RELEASE "1")


SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
INCLUDE(CPack)
