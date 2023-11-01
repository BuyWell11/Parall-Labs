cmake_minimum_required(VERSION 3.26)
project(lab1)

set(CMAKE_CXX_STANDARD 17)

add_executable(lab1 main.cpp
        timer.h)
target_link_libraries(lab1 PRIVATE pthread)