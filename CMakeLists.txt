cmake_minimum_required(VERSION 3.26)

project(l1fit LANGUAGES CXX)

# On Windows the C++ is built with MSVC and the Fortran with MinGW-w64 gfortran. That works because the Fortran
# code needs no runtime library (see below): its object files are linked by the MSVC linker only. CMake gives
# the Fortran compiler check the MSVC linker flags (/machine:x64), which gfortran rejects, so they are
# emptied while the language is enabled.
if(WIN32 AND CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
    set(l1fit_saved_exe_linker_flags "${CMAKE_EXE_LINKER_FLAGS}")
    set(CMAKE_EXE_LINKER_FLAGS "")
    enable_language(Fortran)
    set(CMAKE_EXE_LINKER_FLAGS "${l1fit_saved_exe_linker_flags}")
else()
    enable_language(Fortran)
endif()

# The Fortran solver is compiled into object files that are added to the Python extension module.
# (An OBJECT library, not a STATIC one: linking a Fortran static library makes CMake add the Fortran
# runtime libraries -lgfortran -lquadmath ... to the link line of the module, whether they are needed or not.)
# ---------------------------------------------------------------------------------------------------
add_library(l1fit_fortran OBJECT
    src/l1_precision.f90
    src/l1_calgo552.f90
    src/l1_c_api.f90
)
set_target_properties(l1fit_fortran PROPERTIES
    POSITION_INDEPENDENT_CODE ON
    Fortran_MODULE_DIRECTORY "${CMAKE_BINARY_DIR}/fortran_modules"
)

# The same flags as the test adapters (tests/legacy_cl1.py): no fused multiply-add, so that the results
# are reproducible and bit-identical to the tested reference builds on every CPU. The release default of
# CMake (-O3) is replaced, not extended, so that the optimisation level does not depend on flag order.
if(CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    set(CMAKE_Fortran_FLAGS_RELEASE "-O2")
    target_compile_options(l1fit_fortran PRIVATE -ffp-contract=off)
endif()

# Which Fortran runtime library to link: NONE (default) | SHARED | STATIC. See the end of this file.
set(L1FIT_FORTRAN_RUNTIME "NONE" CACHE STRING "Fortran runtime library to link: NONE, SHARED or STATIC")
set_property(CACHE L1FIT_FORTRAN_RUNTIME PROPERTY STRINGS NONE SHARED STATIC)

# CMake adds the implicit Fortran runtime libraries (-lgfortran -lquadmath ...) to every target that
# contains Fortran objects, whether they are needed or not (CMake documents this for
# CMAKE_<LANG>_IMPLICIT_LINK_LIBRARIES). With NONE they are removed, and checked for below.
if(L1FIT_FORTRAN_RUNTIME STREQUAL "NONE")
    set(CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES "")
    set(CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES "")
endif()

# Python extension module (nanobind, stable ABI for Python >= 3.12)
# ---------------------------------------------------------------------------------------------------
find_package(Python 3.12 REQUIRED COMPONENTS Interpreter Development.Module ${SKBUILD_SABI_COMPONENT})
find_package(nanobind CONFIG REQUIRED)

nanobind_add_module(_core STABLE_ABI NB_STATIC src/_core.cpp $<TARGET_OBJECTS:l1fit_fortran>)

# The Fortran solver needs no Fortran runtime library: it only calls malloc, free and memset. Not linking
# the runtime keeps the wheels free of a compiler runtime to bundle (and portable to machines without
# one). The check below fails the build if a change ever introduces a runtime call (an array temporary,
# a formatted write, error stop, ...): such a module would fail to import where the runtime is missing.
# L1FIT_FORTRAN_RUNTIME = NONE (default) | SHARED | STATIC links the runtime anyway, and skips the check.

if(L1FIT_FORTRAN_RUNTIME STREQUAL "STATIC" AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    foreach(runtime_library gfortran quadmath)
        execute_process(
            COMMAND "${CMAKE_Fortran_COMPILER}" "-print-file-name=lib${runtime_library}.a"
            OUTPUT_VARIABLE runtime_path OUTPUT_STRIP_TRAILING_WHITESPACE)
        if(NOT IS_ABSOLUTE "${runtime_path}")
            message(FATAL_ERROR "no static lib${runtime_library}.a found for ${CMAKE_Fortran_COMPILER}")
        endif()
        target_link_libraries(_core PRIVATE "${runtime_path}")
    endforeach()
elseif(L1FIT_FORTRAN_RUNTIME STREQUAL "SHARED")
    target_link_directories(_core PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_DIRECTORIES})
    target_link_libraries(_core PRIVATE ${CMAKE_Fortran_IMPLICIT_LINK_LIBRARIES})
elseif(NOT L1FIT_FORTRAN_RUNTIME STREQUAL "NONE")
    message(FATAL_ERROR "L1FIT_FORTRAN_RUNTIME must be NONE, SHARED or STATIC, not '${L1FIT_FORTRAN_RUNTIME}'")
endif()

# Guard: with no runtime linked, the Fortran code must not refer to one. The object files are checked, not
# the final module: the linker strips unreachable code, so a reference that is dead today would not show
# there, but would after the next refactoring. Needs `nm` (GNU Fortran on Unix).
find_program(L1FIT_NM nm)
if(L1FIT_FORTRAN_RUNTIME STREQUAL "NONE" AND L1FIT_NM AND CMAKE_Fortran_COMPILER_ID STREQUAL "GNU")
    add_custom_command(TARGET _core POST_BUILD
        COMMAND "${CMAKE_COMMAND}"
            "-DNM=${L1FIT_NM}" "-DOBJECTS=$<JOIN:$<TARGET_OBJECTS:l1fit_fortran>,$<COMMA>>"
            -P "${CMAKE_SOURCE_DIR}/cmake/check_no_fortran_runtime.cmake"
        COMMENT "Checking that the Fortran code needs no Fortran runtime library"
        VERBATIM)
endif()

install(TARGETS _core LIBRARY DESTINATION l1fit)
