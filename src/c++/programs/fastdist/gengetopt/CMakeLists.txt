add_custom_target(fastdist_gengetopt_target DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/fastdist_gengetopt.c ${CMAKE_CURRENT_BINARY_DIR}/fastdist_gengetopt.h )
add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/fastdist_gengetopt.c   ${CMAKE_CURRENT_BINARY_DIR}/fastdist_gengetopt.h
                     COMMAND ${GENGETOPT} -i  ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.ggo  
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     DEPENDS   ${CMAKE_CURRENT_SOURCE_DIR}/cmdline.ggo
                     VERBATIM)

