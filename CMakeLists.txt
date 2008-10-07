cmake_minimum_required(VERSION 2.6.0)
project(fastphylo)

ENABLE_TESTING()

SET(CPACK_PACKAGE_VERSION_MAJOR "0")
SET(CPACK_PACKAGE_VERSION_MINOR "9")
SET(CPACK_PACKAGE_VERSION_PATCH "0")
set(PACKAGE_VERSION ${CPACK_PACKAGE_VERSION_MAJOR}.${CPACK_PACKAGE_VERSION_MINOR}.${CPACK_PACKAGE_VERSION_PATCH})


OPTION(BUILD_DOCBOOK "build Docbook html documentation" OFF )
OPTION(WITH_LIBXML "build xml functionality. This adds a dependeny for libxml" ON )


add_subdirectory(src/c++)
add_subdirectory(src/docbook)

# after running cmake, we can generate source tarballs ( in tar.gz and zip ) by running
# make package_source
#

SET(CPACK_GENERATOR RPM;DEB)
SET(CPACK_SOURCE_GENERATOR TGZ;ZIP)
# SET(CPACK_RPM_PACKAGE_ARCHITECTURE "i386")
SET(CPACK_PACKAGE_CONTACT "erik.sjolund@gmail.com") 
SET(CPACK_PACKAGE_DESCRIPTION_SUMMARY "fastphylo is software project containing the implementations of the algorithms \"Fast Computation of Distance Estimators\" and \"Fast Neighbor Joining\".")
SET(CPACK_PACKAGE_VENDOR "Isaac Elias, Jensl Lagergren, Erik Sjölund")
SET(CPACK_PACKAGE_DESCRIPTION_FILE "${CMAKE_CURRENT_SOURCE_DIR}/README")
SET(CPACK_RESOURCE_FILE_LICENSE "${CMAKE_CURRENT_SOURCE_DIR}/COPYING")
SET(CPACK_PACKAGE_INSTALL_DIRECTORY "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
SET(CPACK_SOURCE_PACKAGE_FILE_NAME "${CMAKE_PROJECT_NAME}-${PACKAGE_VERSION}")
INCLUDE(CPack)
