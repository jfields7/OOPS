cmake_minimum_required(VERSION 3.0)
project(MultiGrid)

# We only need to do any of this if test building is enabled.

if(NOT BUILD_TESTS)
  return()
endif()

# Generate the Parameters and ParamParser files.
set(PARAM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/waveparser.cpp)
set(PARAM_INC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/waveparameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/waveparser.h
   )
set(SETUP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scripts/wave.json)

add_custom_command(
  OUTPUT ${PARAM_INC}
         ${PARAM_SRC}
  DEPENDS ${SETUP_SRC}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/genParams.py ${SETUP_SRC}
  COMMENT "Generating custom Parameters files"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM USES_TERMINAL
)

set(MULTIGRID_INCLUDE_FILES
    include/firstorderwave.h
   )
set(MULTIGRID_SOURCE_FILES
    src/multiGridTest.cpp
    src/firstorderwave.cpp
   )

set(SOURCE_FILES ${MULTIGRID_INCLUDE_FILES} ${MULTIGRID_SOURCE_FILES} ${PARAM_INC} ${PARAM_SRC})
add_executable(MultiGridTest ${SOURCE_FILES})
target_include_directories(MultiGridTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(MultiGridTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(MultiGridTest oops ${EXTRA_LIBS})
