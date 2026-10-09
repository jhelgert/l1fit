cmake_minimum_required(VERSION 3.26)

project(l1fit LANGUAGES CXX Fortran)

# The Fortran solver is compiled into a static library and linked into the Python extension module.
# ---------------------------------------------------------------------------------------------------
add_library(l1fit_fortran STATIC
    src/l1_precision.f90
    src/l1_calgo552.f90
    src/l1_c_api.f90
)
set_target_properties(l1fit_fortran PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran_modules"
)

# The same flags as the test adapters (tests/legacy_cl1.py): no fused multiply-add, so that the results
# are reproducible and bit-identical to the tested reference builds on every CPU.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    target_compile_options(l1fit_fortran PRIVATE -O2 -ffp-contract=off)
endif()

# Python extension module (nanobind, stable ABI for Python >= 3.12)
# ---------------------------------------------------------------------------------------------------
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core STABLE_ABI NB_STATIC src/_core.cpp)
target_link_libraries(_core PRIVATE l1fit_fortran)

# CMake adds the Fortran runtime library only when Fortran is the linker language. The module is linked
# by the C++ driver, so the runtime that the compiler reported is added here.
option(L1FIT_STATIC_FORTRAN_RUNTIME "Link the Fortran runtime library statically" OFF)
if(L1FIT_STATIC_FORTRAN_RUNTIME AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    foreach(runtime_library gfortran quadmath)
        execute_process(
            COMMAND "${CMAKE_Fortran_COMPILER}" "-print-file-name=lib${runtime_library}.a"
            OUTPUT_VARIABLE runtime_path OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT IS_ABSOLUTE "${runtime_path}")
            message(FATAL_ERROR "no static lib${runtime_library}.a found for ${CMAKE_Fortran_COMPILER}")
        endif()
        target_link_libraries(_core PRIVATE "${runtime_path}")
    endforeach()
else()
    target_link_directories(_core PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
    target_link_libraries(_core PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
endif()

install(TARGETS _core LIBRARY DESTINATION l1fit)
