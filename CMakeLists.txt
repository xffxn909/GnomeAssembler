cmake_minimum_required(VERSION 3.20)

project(GnomeAssembler
    VERSION 0.1.0
    DESCRIPTION "Windows-focused de novo virus genome assembler"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

add_executable(gnome_assembler
    src/main.cpp)

target_compile_features(gnome_assembler PRIVATE cxx_std_17)

if(MSVC)
    target_compile_options(gnome_assembler PRIVATE /W4 /permissive-)
else()
    target_compile_options(gnome_assembler PRIVATE -Wall -Wextra -Wpedantic)
endif()
