cmake_minimum_required(VERSION 3.20)
project (lab2_1)

set(CMAKE_CXX_STANDARD 14)

# Add MPI Package to Project
find_package(MPI REQUIRED)

add_executable(lab2_1 main.cpp)
# Add libraries for code completion and compiling
target_link_libraries(lab2_1 PUBLIC MPI::MPI_CXX)