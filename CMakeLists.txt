cmake_minimum_required(VERSION 3.15)
project(skewed_pa)

set(CMAKE_CXX_STANDARD 20)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -march=native -ftemplate-depth=100000 -fopenmp")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -W -Wall")

find_package(OpenMP REQUIRED)

add_subdirectory(extlib/networkit)
add_subdirectory(extlib/stxxl)

add_library(libskewedpa INTERFACE)
target_include_directories(libskewedpa INTERFACE include/)
target_link_libraries(libskewedpa INTERFACE tlx foxxll stxxl networkit OpenMP::OpenMP_CXX)

add_executable(benchmark-polynomial-emgenpa main_polynomial_emgenpa.cpp)
target_link_libraries(benchmark-polynomial-emgenpa libskewedpa ${STXXL_LIBRARIES})

add_executable(benchmark-networkit-ba main_networkit_ba.cpp)
target_link_libraries(benchmark-networkit-ba libskewedpa ${STXXL_LIBRARIES})