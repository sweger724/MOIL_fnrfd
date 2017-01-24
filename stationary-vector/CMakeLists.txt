cmake_minimum_required(VERSION 2.8.9)

project(stationary-vector CXX)
enable_language(CXX)

find_package(PkgConfig)
pkg_check_modules(ARPACK arpack)

include_directories(${ARPACK_INCLUDE_DIRS})

set(CXX_WARN_FLAGS "-W -Wall -Wextra")
set(CMAKE_CXX_FLAGS_DEBUG "-O0")
set(CMAKE_CXX_FLAGS_RELEASE "-O3 -DNDEBUG")
set(CMAKE_CXX_FLAGS "-g -std=c++11 ${CXX_WARN_FLAGS} ${ARPACK_CFLAGS}")

add_executable(stationary-vector stationary-vector.cpp linear-algebra.cpp eigenpair.cpp)
target_link_libraries(stationary-vector ${ARPACK_LDFLAGS})
