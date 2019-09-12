cmake_minimum_required(VERSION 3.1 FATAL_ERROR)
project(ExpressiveMonteCalro NONE)

include(config/core.cmake)

include_directories(${PROJECT_SOURCE_DIR}/include)

# Google test 
set(GOOGLETEST_VERSION 1)
add_subdirectory(tool/googletest)
include_directories(tool/googletest/include)

enable_testing()
add_subdirectory(test)
add_subdirectory(example)
