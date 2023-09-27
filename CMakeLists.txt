cmake_minimum_required(VERSION 3.23)
project(SymNMF C)

set(CMAKE_C_STANDARD 11)

include_directories(.)

add_executable(SymNMF
        symnmf.c
        symnmf.h
        symnmfmodule.c)
