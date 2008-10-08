FIND_PROGRAM(GENGETOPT gengetopt)
IF(NOT GENGETOPT)
 SET(GENGETOPTDEP "${CMAKE_CURRENT_BINARY_DIR}/external/install/bin/gengetopt"  )
IF(WIN32)
 SET(GENGETOPTDEP "${GENGETOPTDEP}.exe"  )
ENDIF(WIN32)
 SET(GENGETOPT "${GENGETOPTDEP}" )
ENDIF(NOT GENGETOPT)

IF(STATIC)
 SET(LIBXML2_INCLUDE_DIR  "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2" )
 SET(LIBXML2_LIBRARIES  "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libxml2.a" )

# The next line seems to be need for macosx
 SET(ZLIB_LIBRARY "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a" )
ENDIF(STATIC)

#message(  FATAL_ERROR  "OU ${GENGETOPT} aa ${CMAKE_CURRENT_BINARY_DIR} bb" )

SET(GENGETOPT_VERSION 2.22 )
SET(LIBXML2_VERSION 2.7.2 )
SET(ZLIB_VERSION 1.2.3 )

SET(DOWNLOAD_DIR "${CMAKE_CURRENT_BINARY_DIR}/external/downloads" INTERNAL )
SET(BUILD_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/build INTERNAL )
SET(INSTALL_DIR ${CMAKE_CURRENT_BINARY_DIR}/external/install INTERNAL )


file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/download )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/build )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/external/install )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt )
file( MAKE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt )


  ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe" COMMAND makensis "${CMAKE_CURRENT_BINARY_DIR}/script.nsi" DEPENDS fastdist fnj )

  ADD_CUSTOM_TARGET( win32installer DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/fastphylo-${PACKAGE_VERSION}-win32.exe" )

  ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/gengetopt-${GENGETOPT_VERSION}.tar.gz" COMMAND  wget --directory-prefix=external/download "ftp://ftp.gnu.org/gnu/gengetopt/gengetopt-${GENGETOPT_VERSION}.tar.gz" )

  ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/libxml2-${LIBXML2_VERSION}.tar.gz" COMMAND  wget --directory-prefix=external/download  "ftp://xmlsoft.org/libxml2/libxml2-${LIBXML2_VERSION}.tar.gz" )

  ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/download/zlib-${ZLIB_VERSION}.tar.gz" COMMAND  wget --directory-prefix=external/download  "http://www.gzip.org/zlib/zlib-${ZLIB_VERSION}.tar.gz" )

  ADD_CUSTOM_COMMAND( OUTPUT "${GENGETOPTDEP}" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_gengetopt.sh" ${GENGETOPT_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install"  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/gengetopt-${GENGETOPT_VERSION}.tar.gz" WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build" )

  ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libxml2.a" "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_libxml2.sh" ${LIBXML2_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install"  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/libxml2-${LIBXML2_VERSION}.tar.gz" "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a"  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build" )

  ADD_CUSTOM_COMMAND( OUTPUT "${CMAKE_CURRENT_BINARY_DIR}/external/install/lib/libz.a" COMMAND sh "${CMAKE_CURRENT_SOURCE_DIR}/build_zlib.sh" ${ZLIB_VERSION} "${CMAKE_CURRENT_BINARY_DIR}/external/install"  DEPENDS "${CMAKE_CURRENT_BINARY_DIR}/external/download/zlib-${ZLIB_VERSION}.tar.gz"  WORKING_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/external/build" )


add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.c   ${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.h
                     COMMAND ${GENGETOPT} -i  ${CMAKE_CURRENT_SOURCE_DIR}/programs/fastdist/gengetopt/cmdline.ggo  --output-dir=${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt
                     DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/programs/fastdist/gengetopt/cmdline.ggo  ${GENGETOPTDEP}
                     VERBATIM)


add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.c   ${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.h
                     COMMAND ${GENGETOPT} -i  ${CMAKE_CURRENT_SOURCE_DIR}/programs/fnj/gengetopt/cmdline.ggo  --output-dir=${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt
                     DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/programs/fnj/gengetopt/cmdline.ggo  ${GENGETOPTDEP}
                     VERBATIM)


INCLUDE (CMakeDetermineCXXCompiler)
IF(NOT CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE(SEND_ERROR "Could not find gcc")
ENDIF(NOT CMAKE_COMPILER_IS_GNUCXX)


IF(WITH_LIBXML)
  FIND_PACKAGE(LibXml2)
ENDIF(WITH_LIBXML)


CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/config.h.cmake ${CMAKE_CURRENT_BINARY_DIR}/config.h)

CONFIGURE_FILE(${CMAKE_CURRENT_SOURCE_DIR}/script.nsi.cmake ${CMAKE_CURRENT_BINARY_DIR}/script.nsi)

INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
#LINK_DIRECTORIES(${LIBXML2_LIBRARIES})


INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/distance_methods)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/DNA_b128)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/programs)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/programs/fastdist)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/programs/fnj)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt)
INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/sequence_likelihood)
#INCLUDE_DIRECTORIES(${CMAKE_CURRENT_SOURCE_DIR}/aml)

SET(FASTPHYLO_SRCS BitVector.cpp 
Exception.cpp
InitAndPrintOn_utils.cpp
Object.cpp
Sequence.cpp
SequenceTree.cpp
SequenceTree_MostParsimonious.cpp
Simulator.cpp
arg_utils_ext.cpp
file_utils.cpp
stl_utils.cpp
DNA_b128/DNA_b128_String.cpp
DNA_b128/Sequences2DistanceMatrix.cpp
#aml/AML_LeafLifting.cpp
#aml/AML_given_edge_probabilities.cpp
#aml/AML_local_improve.cpp
#aml/AML_star.cpp
#aml/Big_AML.cpp
distance_methods/LeastSquaresFit.cpp
distance_methods/NeighborJoining.cpp
sequence_likelihood/Kimura2parameter.cpp
sequence_likelihood/TamuraNei.cpp
sequence_likelihood/ambiguity_nucleotide.cpp
sequence_likelihood/dna_pairwise_sequence_likelihood.cpp
sequence_likelihood/string_compare.cpp 
DistanceMatrix.cpp
arg_utils.c
std_c_utils.c
xml_output_global.cpp
)

SET(FASTPHYLO_SPECIAL_SRCS DNA_b128/sse2_wrapper.c 
                          DNA_b128/computeTAMURANEIDistance_DNA_b128_String.cpp 
                          DNA_b128/computeDistance_DNA_b128_String.cpp) 

# Erik Sjolund removed -fforce-mem from COMMON_FLAGS because of the g++ warning
# cc1plus: warning: -f[no-]force-mem is nop and option will be removed in 4.2

IF(STATIC AND WITH_LIBXML)
 set( LIBXML_EXTRA_COMPILE_FLAGS " -DLIBXML_STATIC " )
ENDIF(STATIC AND WITH_LIBXML)


SET(COMMON_FLAGS  "-march=pentium4 -msse2 -falign-functions -falign-jumps -falign-labels -falign-loops -fcse-follow-jumps -fcse-skip-blocks -fdelete-null-pointer-checks -fexpensive-optimizations -fgcse -fgcse-lm -fgcse-sm -finline-functions -foptimize-sibling-calls -fpeephole2 -fregmove -frename-registers -freorder-blocks -freorder-functions -frerun-cse-after-loop -frerun-loop-opt -fstrength-reduce -fstrict-aliasing ${LIBXML_EXTRA_COMPILE_FLAGS} -O1" )

SET(FASTDIST_SRCS programs/fastdist/main.cpp
programs/fastdist/PhylipMaInputStream.cpp
programs/fastdist/FastaInputStream.cpp
programs/fastdist/DataOutputStream.cpp
programs/fastdist/XmlOutputStream.cpp
${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.c 
)

IF(WITH_LIBXML)
 IF(STATIC)
  SET(FASTDIST_XML_SRCS programs/fastdist/XmlInputStream.cpp  "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h" )
 ELSE(STATIC)
  SET(FASTDIST_XML_SRCS programs/fastdist/XmlInputStream.cpp )
 ENDIF(STATIC)
ENDIF(WITH_LIBXML)


SET(FNJ_SRCS programs/fnj/main.cpp
programs/fnj/DataInputStream.cpp
programs/fnj/DataOutputStream.cpp
programs/fnj/XmlOutputStream.cpp
${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.c 
)

IF(WITH_LIBXML)
 IF(STATIC)
  SET(FNJ_XML_SRCS programs/fnj/XmlInputStream.cpp  "${CMAKE_CURRENT_BINARY_DIR}/external/install/include/libxml2/libxml/xmlreader.h" )
 ELSE(STATIC)
  SET(FNJ_XML_SRCS programs/fnj/XmlInputStream.cpp )
 ENDIF(STATIC)
ENDIF(WITH_LIBXML)

SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.c  PROPERTIES GENERATED true)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.h  PROPERTIES GENERATED true)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.c  PROPERTIES GENERATED true)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.h  PROPERTIES GENERATED true)

FOREACH(i ${FASTPHYLO_SPECIAL_SRCS})
  SET_SOURCE_FILES_PROPERTIES(${i} PROPERTIES COMPILE_FLAGS "${COMMON_FLAGS} -DNDEBUG  -fforce-addr  -fno-branch-count-reg  -fomit-frame-pointer ")
ENDFOREACH(i)


FOREACH(i ${FASTPHYLO_SRCS} ${FASTDIST_SRCS} ${FNJ_SRCS}  ${FASTDIST_XML_SRCS} ${FNJ_XML_SRCS} )
  SET_SOURCE_FILES_PROPERTIES(${i} PROPERTIES COMPILE_FLAGS "${COMMON_FLAGS}  -fcaller-saves  -ffast-math  -fno-default-inline  -fprefetch-loop-arrays  -fsched-interblock  -fsched-spec  -fschedule-insns  -fschedule-insns2  -ftracer  -funroll-loops ")
ENDFOREACH(i)

ADD_LIBRARY(fastphylo STATIC ${FASTPHYLO_SRCS} ${FASTPHYLO_SPECIAL_SRCS})

ADD_EXECUTABLE(fastdist ${FASTDIST_SRCS}  ${FASTDIST_XML_SRCS} )
TARGET_LINK_LIBRARIES(fastdist m ${LIBXML2_LIBRARIES} ${ZLIB_LIBRARY} fastphylo )

IF(STATIC)
 set_property(TARGET fastdist PROPERTY LINK_SEARCH_END_STATIC ON )
ENDIF(STATIC)

ADD_EXECUTABLE(fnj ${FNJ_SRCS} ${FNJ_XML_SRCS} ) 
TARGET_LINK_LIBRARIES(fnj m  ${LIBXML2_LIBRARIES}  ${ZLIB_LIBRARY} fastphylo )

IF(STATIC)
 set_property(TARGET fnj PROPERTY LINK_SEARCH_END_STATIC ON )
ENDIF(STATIC)


INSTALL(TARGETS fastdist DESTINATION bin)
INSTALL(TARGETS fnj DESTINATION bin)



SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/${LIBXML2_INCLUDE_DIR}/libxml/xmlreader.h PROPERTIES GENERATED TRUE HEADER_FILE_ONLY TRUE )
