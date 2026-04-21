cmake_minimum_required(VERSION 3.10)
project(OptimizationProject)

include(FetchContent)

set(BUILD_SHARED_LIBS OFF CACHE BOOL "" FORCE)

FetchContent_Declare(
    nlopt
    GIT_REPOSITORY https://github.com/stevengj/nlopt.git
    GIT_TAG v2.7.1
)

set(NLOPT_PYTHON OFF CACHE BOOL "" FORCE)
set(NLOPT_OCTAVE OFF CACHE BOOL "" FORCE)
set(NLOPT_MATLAB OFF CACHE BOOL "" FORCE)
set(NLOPT_GUILE OFF CACHE BOOL "" FORCE)
set(NLOPT_TESTS OFF CACHE BOOL "" FORCE)

FetchContent_MakeAvailable(nlopt)

FetchContent_GetProperties(nlopt)

set(CMAKE_CXX_STANDARD 17)

set(SOURCE_FILES
    src/main.cpp
    src/peano_mapper.cpp
    src/json_logger.cpp
    src/IOptProblem.cpp
    src/IOptProblemFamily.cpp
    src/IGeneralOptProblem.cpp
    src/IGeneralOptProblemFamily.cpp
    src/IConstrainedOptProblem.cpp
    src/IConstrainedOptProblemFamily.cpp
    Hill/HillProblem.cpp
    Shekel/ShekelProblem.cpp
    Grishagin/GrishaginOption.cpp
    Grishagin/grishagin_function.cpp
    GKLS/GKLSProblem.cpp
    GKLS/gkls_random.cpp
)

add_executable(OptimizationProject ${SOURCE_FILES})

target_link_libraries(OptimizationProject PRIVATE nlopt)

target_include_directories(OptimizationProject PUBLIC
    "${CMAKE_CURRENT_SOURCE_DIR}/include"
    "${CMAKE_CURRENT_SOURCE_DIR}/Hill"
    "${CMAKE_CURRENT_SOURCE_DIR}/Shekel"
    "${CMAKE_CURRENT_SOURCE_DIR}/src"
    "${CMAKE_CURRENT_SOURCE_DIR}/GKLS"
    "${CMAKE_CURRENT_SOURCE_DIR}/Grishagin"
    "${nlopt_SOURCE_DIR}/src/api"
    "${nlopt_BINARY_DIR}"
)