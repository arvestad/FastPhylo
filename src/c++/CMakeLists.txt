CMAKE_MINIMUM_REQUIRED(VERSION 2.4) 

INCLUDE (CMakeDetermineCXXCompiler)
IF(NOT CMAKE_COMPILER_IS_GNUCXX)
      MESSAGE(SEND_ERROR "Could not find gcc")
ENDIF(NOT CMAKE_COMPILER_IS_GNUCXX)

PROJECT(fastphylo)
FIND_PACKAGE(LibXml2)
FIND_PROGRAM(GENGETOPT gengetopt)


INCLUDE_DIRECTORIES(${LIBXML2_INCLUDE_DIR})
LINK_DIRECTORIES(${LIBXML2_LIBRARY})


INCLUDE_DIRECTORIES(${fastphylo_BINARY_DIR}/)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/distance_methods)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/DNA_b128)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/programs)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/programs/fastdist)
INCLUDE_DIRECTORIES(${fastphylo_BINARY_DIR}/programs/fastdist/gengetopt)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/programs/fnj)
INCLUDE_DIRECTORIES(${fastphylo_BINARY_DIR}/programs/fnj/gengetopt)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/sequence_likelihood)
INCLUDE_DIRECTORIES(${fastphylo_SOURCE_DIR}/aml)


ADD_SUBDIRECTORY(programs/fnj/gengetopt)
ADD_SUBDIRECTORY(programs/fastdist/gengetopt)

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
aml/AML_LeafLifting.cpp
aml/AML_given_edge_probabilities.cpp
aml/AML_local_improve.cpp
aml/AML_star.cpp
aml/Big_AML.cpp
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

SET(COMMON_FLAGS  "-march=pentium4 -msse2 -falign-functions -falign-jumps -falign-labels -falign-loops -fcse-follow-jumps -fcse-skip-blocks -fdelete-null-pointer-checks -fexpensive-optimizations -fgcse -fgcse-lm -fgcse-sm -finline-functions -foptimize-sibling-calls -fpeephole2 -fregmove -frename-registers -freorder-blocks -freorder-functions -frerun-cse-after-loop -frerun-loop-opt -fstrength-reduce -fstrict-aliasing -O1" )

SET(FASTDIST_SRCS programs/fastdist/main.cpp
programs/fastdist/PhylipMaInputStream.cpp
programs/fastdist/FastaInputStream.cpp
programs/fastdist/XmlInputStream.cpp
programs/fastdist/DataOutputStream.cpp
programs/fastdist/XmlOutputStream.cpp
${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.c 
)

SET(FNJ_SRCS programs/fnj/main.cpp
programs/fnj/DataInputStream.cpp
programs/fnj/XmlInputStream.cpp
programs/fnj/DataOutputStream.cpp
programs/fnj/XmlOutputStream.cpp
${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.c 
)

SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.c  PROPERTIES GENERATED true)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fnj/gengetopt/fnj_gengetopt.h  PROPERTIES GENERATED true)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.c  PROPERTIES GENERATED true)
SET_SOURCE_FILES_PROPERTIES(${CMAKE_CURRENT_BINARY_DIR}/programs/fastdist/gengetopt/fastdist_gengetopt.h  PROPERTIES GENERATED true)



FOREACH(i ${FASTPHYLO_SPECIAL_SRCS})
  SET_SOURCE_FILES_PROPERTIES(${i} PROPERTIES COMPILE_FLAGS "${COMMON_FLAGS} -DNDEBUG  -fforce-addr  -fno-branch-count-reg  -fomit-frame-pointer ")
ENDFOREACH(i)

FOREACH(i ${FASTPHYLO_SRCS} ${FASTDIST_SRCS} ${FNJ_SRCS} )
  SET_SOURCE_FILES_PROPERTIES(${i} PROPERTIES COMPILE_FLAGS "${COMMON_FLAGS}  -fcaller-saves  -ffast-math  -fno-default-inline  -fprefetch-loop-arrays  -fsched-interblock  -fsched-spec  -fschedule-insns  -fschedule-insns2  -ftracer  -funroll-loops  ")
ENDFOREACH(i)

ADD_EXECUTABLE(fastdist ${FASTDIST_SRCS} ${FASTPHYLO_SRCS} ${FASTPHYLO_SPECIAL_SRCS})
TARGET_LINK_LIBRARIES(fastdist m xml2)

ADD_EXECUTABLE(fnj ${FNJ_SRCS} ${FASTPHYLO_SRCS} ${FASTPHYLO_SPECIAL_SRCS}  ) 
TARGET_LINK_LIBRARIES(fnj m xml2)


ADD_DEPENDENCIES(fnj fnj_gengetopt_target)
ADD_DEPENDENCIES(fastdist fastdist_gengetopt_target)

INSTALL(TARGETS fastdist DESTINATION bin)
INSTALL(TARGETS fnj DESTINATION bin)