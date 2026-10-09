# Fails if the Fortran object files refer to symbols of the GNU Fortran runtime (libgfortran, libquadmath).
# Usage: cmake -DNM=<nm> "-DOBJECTS=<object|object|...>" -P check_no_fortran_runtime.cmake
string(REPLACE "|" ";" object_files "${OBJECTS}")
set(offenders "")
foreach(object_file IN LISTS object_files)
    execute_process(
        COMMAND "${NM}" -u "${object_file}"
        OUTPUT_VARIABLE undefined_symbols
        RESULT_VARIABLE nm_result)
    if(NOT nm_result EQUAL 0)
        message(WARNING "nm failed on ${object_file}: the Fortran runtime check was skipped for it")
        continue()
    endif()
    string(REGEX MATCHALL "[^\n]*(_gfortran_|_gfortrani_|quadmath)[^\n]*" found "${undefined_symbols}")
    foreach(symbol IN LISTS found)
        string(STRIP "${symbol}" symbol)
        list(APPEND offenders "${object_file}: ${symbol}")
    endforeach()
endforeach()
if(offenders)
    string(REPLACE ";" "\n  " offenders "${offenders}")
    message(FATAL_ERROR
        "The Fortran code refers to the Fortran runtime, but it is not linked (L1FIT_FORTRAN_RUNTIME=NONE):\n"
        "  ${offenders}\n"
        "The module would fail to import on a machine without libgfortran. Remove the construct that needs "
        "the runtime (formatted I/O, array temporaries, error stop, ...) or build with "
        "-DL1FIT_FORTRAN_RUNTIME=STATIC.")
endif()
