cmake_minimum_required(VERSION 3.23)
project(LDPC C)

set(CMAKE_C_STANDARD 99)

add_executable(LDPC test_mem.c)
target_link_libraries(LDPC m)
