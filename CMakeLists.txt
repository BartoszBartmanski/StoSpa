cmake_minimum_required(VERSION 3.5)

# Project name
project(StoSpa)

# Set the default build type to be Release
if (NOT EXISTS ${CMAKE_BINARY_DIR}/CMakeCache.txt)
    if (NOT CMAKE_BUILD_TYPE)
        set(CMAKE_BUILD_TYPE "Release" CACHE STRING "" FORCE)
    endif()
endif()

# Check if the default save dir exists. If not then source files directory is used.
set(SAVE_DIR "/scratch/bartmanski/Data/" CACHE STRING "Directory where to save data.")
if(EXISTS ${SAVE_DIR})
    message(STATUS "Save dir: ${SAVE_DIR}")
else(EXISTS ${SAVE_DIR})
    message(STATUS "${SAVE_DIR} is not a valid path!")
    set(SAVE_DIR "${CMAKE_CURRENT_SOURCE_DIR}" CACHE STRING "Directory where to save data." FORCE)
    message(STATUS "Save dir: ${SAVE_DIR}")
endif(EXISTS ${SAVE_DIR})

# Add options for warning flags
set(WARNING_FLAGS "")
option(USE_WALL "Compile with Wall flag" OFF)
if(USE_WALL)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Wall ")
endif()
option(USE_WEXTRA "Compile with Wextra flag" OFF)
if(USE_WEXTRA)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Wextra ")
endif()
option(USE_WPEDANTIC "Compile with Wpedantic flag" OFF)
if(USE_WPEDANTIC)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Wpedantic ")
endif()
option(USE_WERROR "Compile with Werror flag" OFF)
if(USE_WERROR)
    set(WARNING_FLAGS "${WARNING_FLAGS}-Werror ")
endif()

# Projects need c++14 standard - for make_unique and make_shared
set(CMAKE_CXX_STANDARD 14)

# Add the O3 optimisation flag and any specified warning flags
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} ${WARNING_FLAGS}")

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
include_directories(Docopt)

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
include_directories(Simulator Reactions Misc Simulations Tests)
file(GLOB SOURCE_FILES_SIMULATION "${CMAKE_CURRENT_SOURCE_DIR}/Simulator/*pp")
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

# Generate a list of source files in Simulations directory
file(GLOB SIM_SOURCES RELATIVE "${CMAKE_CURRENT_SOURCE_DIR}/Simulations/" "${CMAKE_CURRENT_SOURCE_DIR}/Simulations/*cpp")

# For each source file in Simulations directory generate an executable and link docopt and threads libraries
foreach( sim_source ${SIM_SOURCES} )
    string( REPLACE ".cpp" "" sim_exec ${sim_source} )
    add_executable( ${sim_exec} ${CMAKE_CURRENT_SOURCE_DIR}/Simulations/${sim_source} ${SOURCE_FILES})
    target_link_libraries(${sim_exec} docopt Threads::Threads)
endforeach( sim_source ${SIM_SOURCES} )
