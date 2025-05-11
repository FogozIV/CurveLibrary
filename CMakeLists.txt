cmake_minimum_required(VERSION 3.31)
project(CurveLibrary)

set(CMAKE_CXX_STANDARD 20)

file(GLOB_RECURSE SOURCES CONFIGURE_DEPENDS src/*.cpp)
file(GLOB_RECURSE HEADERS CONFIGURE_DEPENDS include/*.h)
add_library(CurveLibrary STATIC ${SOURCES} ${HEADERS})

target_include_directories(CurveLibrary PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/include)