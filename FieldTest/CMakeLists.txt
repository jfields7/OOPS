cmake_minimum_required(VERSION 3.0)
project(FieldTest)

# We only need to do any of this if test building is enabled.

if(NOT BUILD_TESTS)
  return()
endif()

# Generate the Parameters and ParamParser files.
set(PARAM_SRC ${CMAKE_CURRENT_SOURCE_DIR}/src/fieldtestparser.cpp)
set(PARAM_INC
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fieldtestparameters.h
    ${CMAKE_CURRENT_SOURCE_DIR}/include/fieldtestparser.h
   )
set(SETUP_SRC ${CMAKE_CURRENT_SOURCE_DIR}/scripts/fieldtest.json)

add_custom_command(
  OUTPUT ${PARAM_INC}
         ${PARAM_SRC}
  DEPENDS ${SETUP_SRC}
  COMMAND ${PYTHON_EXECUTABLE} ${CMAKE_SOURCE_DIR}/scripts/genParams.py ${SETUP_SRC}
  COMMENT "Generating custom Parameters files"
  WORKING_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}
  VERBATIM USES_TERMINAL
)

set(FIELDTEST_INCLUDE_FILES
    include/constraintwave.h
   )
set(FIELDTEST_SOURCE_FILES
    src/main.cpp
    src/constraintwave.cpp
   )

set(SOURCE_FILES ${FIELDTEST_INCLUDE_FILES} ${FIELDTEST_SOURCE_FILES} ${PARAM_INC} ${PARAM_SRC})
add_executable(FieldTest ${SOURCE_FILES})
target_include_directories(FieldTest PRIVATE ${CMAKE_CURRENT_SOURCE_DIR}/include)
target_include_directories(FieldTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(FieldTest oops ${EXTRA_LIBS})
