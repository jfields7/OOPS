cmake_minimum_required(VERSION 3.0)
project(Tests)

# We only need to do any of this if test building is enabled.

if(NOT BUILD_TESTS)
  return()
endif()

add_executable(DomainTest domainTest.cpp)
target_include_directories(DomainTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(DomainTest oops ${EXTRA_LIBS} m)

add_executable(GridTest gridTest.cpp)
target_include_directories(GridTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(GridTest oops ${EXTRA_LIBS} m)

add_executable(ParamTest paramTest.cpp)
target_include_directories(ParamTest PRIVATE ${CMAKE_SOURCE_DIR}/include)
target_link_libraries(ParamTest oops ${EXTRA_LIBS} m)
