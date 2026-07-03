cmake_minimum_required(VERSION 3.20)

if(NOT DEFINED SRC_DIR OR SRC_DIR STREQUAL "")
  set(SRC_DIR "$ENV{SRC_DIR}")
endif()
if(SRC_DIR STREQUAL "")
  get_filename_component(SRC_DIR "${CMAKE_CURRENT_LIST_DIR}/.." ABSOLUTE)
endif()

if(NOT DEFINED BUILD_ROOT OR BUILD_ROOT STREQUAL "")
  if(DEFINED ENV{BUILD_ROOT} AND NOT "$ENV{BUILD_ROOT}" STREQUAL "")
    set(BUILD_ROOT "$ENV{BUILD_ROOT}")
  else()
    set(BUILD_ROOT "${SRC_DIR}/_build_matrix_cmake")
  endif()
endif()

if(NOT DEFINED LOG_DIR OR LOG_DIR STREQUAL "")
  set(LOG_DIR "${BUILD_ROOT}/logs")
endif()

set(CMAKE_CMD "$ENV{CMAKE}")
if(CMAKE_CMD STREQUAL "")
  set(CMAKE_CMD "cmake")
endif()

set(CTEST_CMD "$ENV{CTEST}")
if(CTEST_CMD STREQUAL "")
  set(CTEST_CMD "ctest")
endif()
set(CTEST_COMMAND "${CTEST_CMD}")

set(GENERATOR "$ENV{GENERATOR}")
set(KEEP_BUILD "$ENV{KEEP_BUILD}")
if(KEEP_BUILD STREQUAL "")
  set(KEEP_BUILD "0")
endif()

if(NOT DEFINED ENABLE_MPFR_ORACLE OR ENABLE_MPFR_ORACLE STREQUAL "")
  set(ENABLE_MPFR_ORACLE "$ENV{ENABLE_MPFR_ORACLE}")
endif()
if(ENABLE_MPFR_ORACLE STREQUAL "")
  set(ENABLE_MPFR_ORACLE "OFF")
endif()

set(QD_TEST_SEED "$ENV{QD_TEST_SEED}")
if(QD_TEST_SEED STREQUAL "")
  set(QD_TEST_SEED "11400714819323198485")
endif()

set(CTEST_SOURCE_DIRECTORY "${SRC_DIR}")
set(CTEST_BINARY_DIRECTORY "${BUILD_ROOT}")

set(JOBS "$ENV{JOBS}")
if(JOBS STREQUAL "")
  cmake_host_system_information(RESULT JOBS QUERY NUMBER_OF_LOGICAL_CORES)
  if(NOT JOBS OR JOBS STREQUAL "0")
    set(JOBS "4")
  endif()
endif()

string(TOLOWER "${CMAKE_HOST_SYSTEM_PROCESSOR}" QD3_HOST_PROCESSOR)
if(QD3_HOST_PROCESSOR STREQUAL "")
  cmake_host_system_information(RESULT QD3_HOST_PROCESSOR QUERY OS_PLATFORM)
  string(TOLOWER "${QD3_HOST_PROCESSOR}" QD3_HOST_PROCESSOR)
endif()

function(qd3_probe_x86_64_v3 out_var)
  set(${out_var} FALSE PARENT_SCOPE)

  set(probe_cxx "$ENV{CXX}")
  if(probe_cxx STREQUAL "")
    find_program(probe_cxx NAMES c++ g++ clang++)
  endif()
  if(NOT probe_cxx)
    message(STATUS "x86-64-v3 probe skipped: no C++ compiler found")
    return()
  endif()

  set(probe_dir "${BUILD_ROOT}/_host_probes/x86-64-v3")
  set(probe_src "${probe_dir}/x86_64_v3.cpp")
  set(probe_exe "${probe_dir}/x86_64_v3")
  if(WIN32)
    set(probe_exe "${probe_exe}.exe")
  endif()
  file(MAKE_DIRECTORY "${probe_dir}")
  file(WRITE "${probe_src}" [[
#include <immintrin.h>

int main() {
  alignas(32) double input[4] = {1.0, 2.0, 3.0, 4.0};
  alignas(32) double output[4] = {};
  __m256d x = _mm256_load_pd(input);
  __m256d y = _mm256_set1_pd(2.0);
  __m256d z = _mm256_set1_pd(1.0);
  __m256d r = _mm256_fmadd_pd(x, y, z);
  _mm256_store_pd(output, r);
  return output[0] == 3.0 && output[3] == 9.0 ? 0 : 1;
}
]])

  execute_process(
    COMMAND "${probe_cxx}" -std=c++11 -march=x86-64-v3 "${probe_src}" -o "${probe_exe}"
    RESULT_VARIABLE build_rc
    OUTPUT_VARIABLE build_out
    ERROR_VARIABLE build_err
  )
  if(NOT build_rc STREQUAL "0")
    message(STATUS "x86-64-v3 probe build failed with ${probe_cxx}; treating FMA as unavailable")
    return()
  endif()

  execute_process(
    COMMAND "${probe_exe}"
    RESULT_VARIABLE run_rc
    OUTPUT_VARIABLE run_out
    ERROR_VARIABLE run_err
  )
  if(run_rc STREQUAL "0")
    set(${out_var} TRUE PARENT_SCOPE)
  else()
    message(STATUS "x86-64-v3 probe did not run cleanly; treating FMA as unavailable")
  endif()
endfunction()

set(QD3_FMA_ENABLED_MODE auto)
set(QD3_FMA_ENABLED_ARCH generic)
set(QD3_FMA_ENABLED_REASON "generic host target")
if(QD3_HOST_PROCESSOR MATCHES "^(x86_64|amd64)$")
  qd3_probe_x86_64_v3(QD3_HOST_HAS_X86_64_V3)
  if(QD3_HOST_HAS_X86_64_V3)
    set(QD3_FMA_ENABLED_ARCH x86-64-v3)
    set(QD3_FMA_ENABLED_REASON "x86-64-v3/FMA probe passed")
  else()
    set(QD3_FMA_ENABLED_MODE no)
    set(QD3_FMA_ENABLED_REASON "x86 host below runnable x86-64-v3")
  endif()
endif()

set(CORNER_REPORT_DIR "${BUILD_ROOT}/reports")
file(MAKE_DIRECTORY "${BUILD_ROOT}" "${LOG_DIR}" "${CORNER_REPORT_DIR}")

set(SUMMARY_OK "${BUILD_ROOT}/summary.ok")
set(SUMMARY_NG "${BUILD_ROOT}/summary.ng")
set(SUMMARY_TSV "${BUILD_ROOT}/summary.tsv")
file(WRITE "${SUMMARY_OK}" "")
file(WRITE "${SUMMARY_NG}" "")
file(WRITE "${SUMMARY_TSV}" "tag\tresult\tlog\n")

set(fail_count 0)

function(qd3_run_process log_file)
  execute_process(
    COMMAND ${ARGN}
    RESULT_VARIABLE rc
    OUTPUT_VARIABLE out
    ERROR_VARIABLE err
  )
  string(REPLACE ";" " " cmd_text "${ARGN}")
  file(APPEND "${log_file}" "$ ${cmd_text}\n")
  if(NOT out STREQUAL "")
    file(APPEND "${log_file}" "${out}")
  endif()
  if(NOT err STREQUAL "")
    file(APPEND "${log_file}" "${err}")
  endif()
  set(QD3_LAST_RC "${rc}" PARENT_SCOPE)
endfunction()

function(qd3_tail_file path max_lines out_var)
  if(NOT EXISTS "${path}")
    set(${out_var} "" PARENT_SCOPE)
    return()
  endif()
  file(STRINGS "${path}" lines)
  list(LENGTH lines line_count)
  math(EXPR start "${line_count} - ${max_lines}")
  if(start LESS 0)
    set(start 0)
  endif()
  set(tail "")
  foreach(index RANGE ${start} ${line_count})
    if(index LESS line_count)
      list(GET lines ${index} line)
      string(APPEND tail "${line}\n")
    endif()
  endforeach()
  set(${out_var} "${tail}" PARENT_SCOPE)
endfunction()

function(qd3_run_config tag)
  set(build_dir "${BUILD_ROOT}/${tag}")
  set(log_file "${LOG_DIR}/${tag}.log")
  set(report_file "${CORNER_REPORT_DIR}/${tag}.rounding_corners.csv")

  if(NOT KEEP_BUILD STREQUAL "1")
    file(REMOVE_RECURSE "${build_dir}")
  endif()
  file(MAKE_DIRECTORY "${build_dir}")
  file(WRITE "${log_file}" "")

  message("=== ${tag} ===")
  message("build_dir: ${build_dir}")
  message("log_file : ${log_file}")

  set(generator_args "")
  if(NOT GENERATOR STREQUAL "")
    list(APPEND generator_args -G "${GENERATOR}")
  endif()

  set(qa_config_args "")
  if(ENABLE_MPFR_ORACLE)
    list(APPEND qa_config_args -DQD3_ENABLE_MPFR_TESTS=ON)
  endif()

  qd3_run_process("${log_file}"
    "${CMAKE_CMD}" -S "${SRC_DIR}" -B "${build_dir}" ${generator_args} ${qa_config_args} ${ARGN}
  )
  set(rc "${QD3_LAST_RC}")

  if(rc EQUAL 0)
    qd3_run_process("${log_file}"
      "${CMAKE_CMD}" --build "${build_dir}" -j "${JOBS}"
    )
    set(rc "${QD3_LAST_RC}")
  endif()

  if(rc EQUAL 0)
    qd3_run_process("${log_file}"
      "${CMAKE_CMD}" -E env "QD_TEST_SEED=${QD_TEST_SEED}"
      "${CTEST_CMD}" --test-dir "${build_dir}" --output-on-failure
    )
    set(rc "${QD3_LAST_RC}")
  endif()

  if(rc EQUAL 0 AND ENABLE_MPFR_ORACLE)
    qd3_run_process("${log_file}"
      "${build_dir}/tests/oracle/test_rounding" --all
      "--seed=${QD_TEST_SEED}"
      "--worst-report=${report_file}"
    )
    set(rc "${QD3_LAST_RC}")
  endif()

  if(rc EQUAL 0)
    message("PASS ${tag}")
    file(APPEND "${SUMMARY_OK}" "${tag}\n")
    file(APPEND "${SUMMARY_TSV}" "${tag}\tPASS\t${log_file}\n")
  else()
    message("FAIL ${tag}")
    file(APPEND "${SUMMARY_NG}" "${tag}\n")
    file(APPEND "${SUMMARY_TSV}" "${tag}\tFAIL\t${log_file}\n")
    qd3_tail_file("${log_file}" 80 tail)
    if(NOT tail STREQUAL "")
      message("${tail}")
    endif()
    math(EXPR new_fail_count "${fail_count} + 1")
    set(fail_count "${new_fail_count}" PARENT_SCOPE)
  endif()
  message("")
endfunction()

qd3_run_config(default)

foreach(ieee_add OFF ON)
  foreach(sloppy_mul OFF ON)
    foreach(sloppy_div OFF ON)
      foreach(fma disabled enabled)
        set(tag "ieee_add-${ieee_add}__sloppy_mul-${sloppy_mul}__sloppy_div-${sloppy_div}__fma-${fma}")
        if(fma STREQUAL "enabled")
          qd3_run_config("${tag}"
            -DQD_ENABLE_IEEE_ADD=${ieee_add}
            -DQD_ENABLE_SLOPPY_MUL=${sloppy_mul}
            -DQD_ENABLE_SLOPPY_DIV=${sloppy_div}
            -DQD_FMA=${QD3_FMA_ENABLED_MODE}
            -DQD_ARCH=${QD3_FMA_ENABLED_ARCH}
          )
        else()
          qd3_run_config("${tag}"
            -DQD_ENABLE_IEEE_ADD=${ieee_add}
            -DQD_ENABLE_SLOPPY_MUL=${sloppy_mul}
            -DQD_ENABLE_SLOPPY_DIV=${sloppy_div}
            -DQD_FMA=no
            -DQD_ARCH=generic
          )
        endif()
      endforeach()
    endforeach()
  endforeach()
endforeach()

message("========================================")
message("Finished.")
message("PASS list : ${SUMMARY_OK}")
message("FAIL list : ${SUMMARY_NG}")
message("TSV       : ${SUMMARY_TSV}")
message("Logs      : ${LOG_DIR}")
message("MPFR mode : ${ENABLE_MPFR_ORACLE}")
message("Seed      : ${QD_TEST_SEED}")
message("Host CPU  : ${QD3_HOST_PROCESSOR}")
message("FMA mode  : ${QD3_FMA_ENABLED_MODE}")
message("FMA arch  : ${QD3_FMA_ENABLED_ARCH}")
message("FMA reason: ${QD3_FMA_ENABLED_REASON}")
message("Failures  : ${fail_count}")

if(ENABLE_MPFR_ORACLE AND EXISTS "${CORNER_REPORT_DIR}")
  execute_process(
    COMMAND python3 "${SRC_DIR}/qa/compare_rounding_matrix.py" "${CORNER_REPORT_DIR}"
    RESULT_VARIABLE compare_rc
  )
  if(NOT compare_rc EQUAL 0)
    message(FATAL_ERROR "CMake rounding matrix comparison failed")
  endif()
endif()

if(NOT fail_count EQUAL 0)
  message(FATAL_ERROR "CMake matrix failed: ${fail_count} configuration(s)")
endif()
