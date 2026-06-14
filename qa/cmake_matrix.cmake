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

set(JOBS "$ENV{JOBS}")
if(JOBS STREQUAL "")
  cmake_host_system_information(RESULT JOBS QUERY NUMBER_OF_LOGICAL_CORES)
  if(NOT JOBS OR JOBS STREQUAL "0")
    set(JOBS "4")
  endif()
endif()

file(MAKE_DIRECTORY "${BUILD_ROOT}" "${LOG_DIR}")

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
            -DQD_FMA=auto
            -DQD_ARCH=x86-64-v3
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
message("Failures  : ${fail_count}")

if(NOT fail_count EQUAL 0)
  message(FATAL_ERROR "CMake matrix failed: ${fail_count} configuration(s)")
endif()
