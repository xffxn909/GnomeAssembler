cmake_minimum_required(VERSION 3.20)

project(GnomeAssembler
    VERSION 0.2.0
    DESCRIPTION "Windows-focused de novo virus genome assembler"
    LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 17)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)

find_package(Qt6 6.2 REQUIRED COMPONENTS Widgets)
qt_standard_project_setup()

# ---------------------------------------------------------------------------
# Core library — Qt-independent pipeline logic
# ---------------------------------------------------------------------------
add_library(gnome_assembler_core
    src/cli/command_line.cpp
    src/cleaning/cleaner.cpp
    src/io/fastq_reader.cpp
    src/io/output_writer.cpp
    src/kmer/kmer_counter.cpp
    src/graph/debruijn_graph.cpp
    src/graph/graph_pipeline.cpp)

target_include_directories(gnome_assembler_core PUBLIC include)
target_compile_features(gnome_assembler_core PUBLIC cxx_std_17)

# ---------------------------------------------------------------------------
# Main executable
# ---------------------------------------------------------------------------
add_executable(gnome_assembler
    src/main.cpp
    src/gui/main_window.cpp
    include/gnome_assembler/main_window.h)

target_link_libraries(gnome_assembler PRIVATE gnome_assembler_core Qt6::Widgets)

# ---------------------------------------------------------------------------
# Compiler warnings
# ---------------------------------------------------------------------------
if(MSVC)
    target_compile_options(gnome_assembler_core PRIVATE /W4 /permissive-)
    target_compile_options(gnome_assembler PRIVATE /W4 /permissive-)
else()
    target_compile_options(gnome_assembler_core PRIVATE -Wall -Wextra -Wpedantic)
    target_compile_options(gnome_assembler PRIVATE -Wall -Wextra -Wpedantic)
endif()

# ---------------------------------------------------------------------------
# Tests
# ---------------------------------------------------------------------------
include(CTest)
if(BUILD_TESTING)
    # Milestone 1: cleaning tests
    add_executable(gnome_assembler_tests tests/unit/test_cleaning.cpp)
    target_link_libraries(gnome_assembler_tests PRIVATE gnome_assembler_core)
    target_compile_features(gnome_assembler_tests PRIVATE cxx_std_17)
    add_test(NAME cleaning COMMAND gnome_assembler_tests)

    # Milestone 1: GUI window smoke test
    add_executable(gnome_assembler_gui_tests
        tests/unit/test_main_window.cpp
        src/gui/main_window.cpp
        include/gnome_assembler/main_window.h)
    target_link_libraries(gnome_assembler_gui_tests PRIVATE gnome_assembler_core Qt6::Widgets)
    target_compile_features(gnome_assembler_gui_tests PRIVATE cxx_std_17)
    add_test(NAME gui_window COMMAND gnome_assembler_gui_tests)

    # Milestone 2: k-mer and graph tests
    add_executable(gnome_assembler_graph_tests tests/unit/test_graph.cpp)
    target_link_libraries(gnome_assembler_graph_tests PRIVATE gnome_assembler_core)
    target_compile_features(gnome_assembler_graph_tests PRIVATE cxx_std_17)
    add_test(NAME graph COMMAND gnome_assembler_graph_tests)
endif()
