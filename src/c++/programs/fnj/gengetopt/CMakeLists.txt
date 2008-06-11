add_custom_target(fnj_gengetopt_target DEPENDS ${CMAKE_CURRENT_BINARY_DIR}/fnj_gengetopt.c ${CMAKE_CURRENT_BINARY_DIR}/fnj_gengetopt.h )
add_custom_command(OUTPUT  ${CMAKE_CURRENT_BINARY_DIR}/fnj_gengetopt.c   ${CMAKE_CURRENT_BINARY_DIR}/fnj_gengetopt.h
                     COMMAND ${GENGETOPT} -i  ${CMAKE_CURRENT_SOURCE_DIR}/fnj.ggo  
                     WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}
                     VERBATIM)

