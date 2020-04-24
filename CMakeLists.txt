cmake_minimum_required(VERSION 3.0)

# set the project name
project(OOPS VERSION 1.0)

# Specify the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

option(USE_SDF "Enable exporting to SDF file" OFF)

# Configure a header file that will pass on our settings.
configure_file(${PROJECT_SOURCE_DIR}/oopsconfig.h.in ${PROJECT_SOURCE_DIR}/oopsconfig.h)

set(EXTRA_LIBS
   )

if(USE_SDF)
  list(APPEND EXTRA_LIBS bbhutil)
endif()

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
    include/interpolator.h
    include/cubicinterpolator.h
    include/polynomialinterpolator.h
   )

set(GRID_SOURCE_FILES
    src/grid.cpp
    src/domain.cpp
    src/rk4.cpp
    src/ode.cpp
    src/solverdata.cpp
    src/cubic.cpp
    src/interpolator.cpp
    src/cubicinterpolator.cpp
    src/polynomialinterpolator.cpp
   )

set(TEST_SOURCE_FILES
    tests/gridTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(GridTest ${SOURCE_FILES})
target_include_directories(GridTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR})

set(TEST_SOURCE_FILES
    tests/domainTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(DomainTest ${SOURCE_FILES})
target_include_directories(DomainTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR})

set(TEST_SOURCE_FILES
    tests/rkTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(RKTest ${SOURCE_FILES})
target_include_directories(RKTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR})

set(TEST_SOURCE_FILES
    tests/waveTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(WaveTest ${SOURCE_FILES})
target_link_libraries(WaveTest ${EXTRA_LIBS})
target_include_directories(WaveTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR})

set(TEST_SOURCE_FILES
    tests/heatTest.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES})
add_executable(HeatTest ${SOURCE_FILES})
target_link_libraries(HeatTest ${EXTRA_LIBS})
target_include_directories(HeatTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include ${PROJECT_SOURCE_DIR})

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
target_link_libraries(MultiGridTest ${EXTRA_LIBS})
target_include_directories(MultiGridTest PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include 
                           ${CMAKE_CURRENT_SOURCE_DIR}/MultiGrid/include ${PROJECT_SOURCE_DIR})

set(FLUID_INCLUDE_FILES
    fluid/include/reconstruction.h
    fluid/include/minmod.h
    fluid/include/norecon.h
   )
set(TEST_SOURCE_FILES
    fluid/src/reconstruction.cpp
    fluid/src/minmod.cpp
    fluid/src/fluidmain.cpp
    fluid/src/norecon.cpp
   )

set(SOURCE_FILES ${GRID_SOURCE_FILES} ${GRID_INCLUDE_FILES} ${TEST_SOURCE_FILES} ${FLUID_INCLUDE_FILES})
add_executable(Fluid ${SOURCE_FILES})
target_link_libraries(Fluid ${EXTRA_LIBS})
target_include_directories(Fluid PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include
                           ${CMAKE_CURRENT_SOURCE_DIR}/fluid/include ${PROJECT_SOURCE_DIR})
