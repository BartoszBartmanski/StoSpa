cmake_minimum_required(VERSION 3.5)

# Project name
project(StoSpa)

# Add options for warning flags
set(WARNING_FLAGS "")
option(USE_Wall "Compile with Wall flag" OFF)
if(USE_Wall)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Wall ")
endif()
option(USE_Wextra "Compile with Wextra flag" OFF)
if(USE_Wextra)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Wextra ")
endif()
option(USE_Wpedantic "Compile with Wpedantic flag" OFF)
if(USE_Wpedantic)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Wpedantic ")
endif()
option(USE_Werror "Compile with Werror flag" OFF)
if(USE_Werror)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Werror ")
endif()

# Projects need c++14 standard - for make_unique and make_shared
set(CMAKE_CXX_STANDARD 14)

# Add the O3 optimisation flag and any specified warning flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -O3 ${WARNING_FLAGS}")

# Check for git version and put it in BUILD_VERSION variable
if(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)
    find_package(Git)
    if(GIT_FOUND)
        EXECUTE_PROCESS(
                COMMAND ${GIT_EXECUTABLE} rev-parse HEAD
                WORKING_DIRECTORY "${CMAKE_CURRENT_SOURCE_DIR}"
                OUTPUT_VARIABLE "BUILD_VERSION"
                ERROR_QUIET
                OUTPUT_STRIP_TRAILING_WHITESPACE)
        message( STATUS "Git version: ${BUILD_VERSION}" )
    else(GIT_FOUND)
        set(BUILD_VERSION 0)
    endif(GIT_FOUND)
endif(EXISTS ${CMAKE_CURRENT_SOURCE_DIR}/.git)

# Generate the Version.hpp file that will contain the git version hash key
configure_file(${CMAKE_CURRENT_SOURCE_DIR}/Misc/Version.hpp.in ${CMAKE_CURRENT_SOURCE_DIR}/Misc/Version.hpp @ONLY)

# Add Docopt to be compiled
add_subdirectory(Docopt)

# Add boost
set(Boost_USE_STATIC_LIBS OFF)
set(Boost_USE_MULTITHREADED ON)
set(Boost_USE_STATIC_RUNTIME OFF)
find_package(Boost)
if(Boost_FOUND)
    include_directories(${Boost_INCLUDE_DIRS})
#    add_executable(progname file1.cxx file2.cxx)
#    target_link_libraries(progname ${Boost_LIBRARIES})
endif()

# Add threads library and set threads flag on
set(THREADS_PREFER_PTHREAD_FLAG ON)
find_package(Threads REQUIRED)

# Add all the source files
include_directories(Simulation Reactions Misc Plots Tests)
file(GLOB SOURCE_FILES_SIMULATION "${CMAKE_CURRENT_SOURCE_DIR}/Simulation/*pp")
file(GLOB SOURCE_FILES_REACTIONS "${CMAKE_CURRENT_SOURCE_DIR}/Reactions/*pp")
file(GLOB SOURCE_FILES_MISC "${CMAKE_CURRENT_SOURCE_DIR}/Misc/*pp")
file(GLOB SOURCE_FILES_UNIT_TESTS "${CMAKE_CURRENT_SOURCE_DIR}/Tests/*pp")

# Place all the source files into a single variable (for convinence)
set(SOURCE_FILES ${SOURCE_FILES_SIMULATION} ${SOURCE_FILES_REACTIONS} ${SOURCE_FILES_MISC})

# Add an executable for Benchmarking
add_executable(benchmarking ${CMAKE_CURRENT_SOURCE_DIR}/Benchmarking/benchmarking.cpp ${SOURCE_FILES})
target_link_libraries(benchmarking ${Boost_LIBRARIES})

# Add an executable for Unit Testing
add_executable(unittest ${SOURCE_FILES_UNIT_TESTS} ${SOURCE_FILES})
target_link_libraries(unittest ${Boost_LIBRARIES})

# Add an executable for running simulations and link docopt to it
add_executable(stospa ${CMAKE_CURRENT_SOURCE_DIR}/main.cpp ${SOURCE_FILES})
target_link_libraries(stospa docopt ${Boost_LIBRARIES})

# Generate a list of source files in Plots directory
file(GLOB PLOTS_SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/Plots/" "${CMAKE_CURRENT_SOURCE_DIR}/Plots/*")

# For each source file in Plots directory generate an executable and link docopt and threads libraries
foreach( plot_source ${PLOTS_SOURCES} )
    string( REPLACE ".cpp" "" plot_exec ${plot_source} )
    add_executable( ${plot_exec} ${CMAKE_CURRENT_SOURCE_DIR}/Plots/${plot_source} ${SOURCE_FILES})
    target_link_libraries(${plot_exec} docopt Threads::Threads ${Boost_LIBRARIES})
endforeach( plot_source ${PLOTS_SOURCES} )
