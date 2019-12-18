cmake_minimum_required(VERSION 3.0)

# set the project name
project(GridTest VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)


set(GRID_INCLUDE_FILES
    include/grid.h
    include/domain.h
    include/rk4.h
    include/solver.h
    include/output.h
    include/cubic.h
    include/solverdata.h
    include/ode.h
    include/parameters.h
   )

set(GRID_SOURCE_FILES
    src/grid.cpp
    src/domain.cpp
    src/rk4.cpp
    src/ode.cpp
    src/solverdata.cpp
    src/cubic.cpp
   )

set(TEST_SOURCE_FILES
    tests/gridTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(GridTest ${SOURCE_FILES})
target_include_directories(GridTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(TEST_SOURCE_FILES
    tests/domainTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(DomainTest ${SOURCE_FILES})
target_include_directories(DomainTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(TEST_SOURCE_FILES
    tests/rkTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(RKTest ${SOURCE_FILES})
target_include_directories(RKTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(TEST_SOURCE_FILES
    tests/waveTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(WaveTest ${SOURCE_FILES})
target_link_libraries(WaveTest bbhutil)
target_include_directories(WaveTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(TEST_SOURCE_FILES
    tests/heatTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(HeatTest ${SOURCE_FILES})
target_link_libraries(HeatTest bbhutil)
target_include_directories(HeatTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)

set(MULTIGRID_INCLUDE_FILES
    MultiGrid/include/waveparameters.h
    MultiGrid/include/firstorderwave.h
   )
set(TEST_SOURCE_FILES
    MultiGrid/src/multiGridTest.cpp
    MultiGrid/src/firstorderwave.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES} ${MULTIGRID_INCLUDE_FILES})
add_executable(MultiGridTest ${SOURCE_FILES})
target_link_libraries(MultiGridTest bbhutil)
target_include_directories(MultiGridTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include 
                           ${CMAKE_CURRENT_SOURCE_DIR}/MultiGrid/include)
