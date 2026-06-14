include(CheckCXXSourceRuns)
include(CheckCXXCompilerFlag)

function(qd3_fma_candidates_for_auto out_var)
  string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" proc)
  if(proc MATCHES "^(powerpc|ppc)")
    set(candidates ibm gnu c99)
  elseif(proc MATCHES "^ia64")
    set(candidates ia64 gnu c99)
  else()
    set(candidates gnu c99)
  endif()
  set(${out_var} ${candidates} PARENT_SCOPE)
endfunction()

function(qd3_fma_source candidate out_var)
  if(candidate STREQUAL "ibm")
    set(src [[
#include <cmath>
#include <builtins.h>
int main() {
  double d = std::ldexp(1.0, -52);
  double x = __fmadd(1.0 + d, 1.0 - d, -1.0);
  double y = __fmsub(1.0 + d, 1.0 - d, 1.0);
  return (x == -d*d && y == -d*d) ? 0 : 1;
}
]])
  elseif(candidate STREQUAL "gnu")
    set(src [[
#include <cmath>
int main() {
  double d = std::ldexp(1.0, -52);
  return (__builtin_fma(1.0 + d, 1.0 - d, -1.0) == -d*d ? 0 : 1);
}
]])
  elseif(candidate STREQUAL "ia64")
    set(src [[
#include <cmath>
int main() {
  double d = std::ldexp(1.0, -52);
  return (_Asm_fma(2, 1.0 + d, 1.0 - d, -1.0) == -d*d ? 0 : 1);
}
]])
  elseif(candidate STREQUAL "c99")
    set(src [[
#include <cmath>
int main() {
  double d = std::ldexp(1.0, -52);
  return (fma(1.0 + d, 1.0 - d, -1.0) == -d*d ? 0 : 1);
}
]])
  elseif(candidate STREQUAL "compiler")
    set(src [[
#include <cmath>
int main() {
  double d = std::ldexp(1.0, -52);
  return (((1.0 + d) * (1.0 - d) - 1.0) == -d*d ? 0 : 1);
}
]])
  else()
    message(FATAL_ERROR "Unknown FMA candidate '${candidate}'")
  endif()
  set(${out_var} "${src}" PARENT_SCOPE)
endfunction()

function(qd3_try_fma_candidate candidate out_var)
  qd3_fma_source("${candidate}" src)
  string(TOUPPER "${candidate}" upper)
  set(var "QD3_FMA_${upper}_RUNS")

  set(saved_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
  check_cxx_source_runs("${src}" ${var})
  if(NOT ${var} AND candidate STREQUAL "c99")
    set(CMAKE_REQUIRED_LIBRARIES ${saved_required_libraries} m)
    set(var_m "QD3_FMA_${upper}_RUNS_WITH_M")
    check_cxx_source_runs("${src}" ${var_m})
    set(CMAKE_REQUIRED_LIBRARIES "${saved_required_libraries}")
    if(${var_m})
      set(${out_var} TRUE PARENT_SCOPE)
      return()
    endif()
  endif()
  set(CMAKE_REQUIRED_LIBRARIES "${saved_required_libraries}")

  if(${var})
    set(${out_var} TRUE PARENT_SCOPE)
  else()
    set(${out_var} FALSE PARENT_SCOPE)
  endif()
endfunction()

function(qd3_check_fma)
  set(QD_FMA_SELECTED "none" PARENT_SCOPE)
  set(QD_FMS_SELECTED "none" PARENT_SCOPE)
  set(QD_VACPP_BUILTINS_H OFF PARENT_SCOPE)
  set(QD_FMA_CONFIG_DEFINE "/* #undef QD_FMA */" PARENT_SCOPE)
  set(QD_FMS_CONFIG_DEFINE "/* #undef QD_FMS */" PARENT_SCOPE)

  if((QD_FMA_EXPR AND NOT QD_FMS_EXPR) OR (QD_FMS_EXPR AND NOT QD_FMA_EXPR))
    message(FATAL_ERROR "QD_FMA_EXPR and QD_FMS_EXPR must be set together")
  endif()

  if(QD_FMA_EXPR AND QD_FMS_EXPR)
    message(WARNING "Using manual QD_FMA_EXPR/QD_FMS_EXPR overrides; no run-time FMA correctness probe will be performed")
    set(QD_FMA_SELECTED "manual" PARENT_SCOPE)
    set(QD_FMS_SELECTED "manual" PARENT_SCOPE)
    set(QD_FMA_CONFIG_DEFINE "#define QD_FMA(x,y,z) ${QD_FMA_EXPR}" PARENT_SCOPE)
    set(QD_FMS_CONFIG_DEFINE "#define QD_FMS(x,y,z) ${QD_FMS_EXPR}" PARENT_SCOPE)
    return()
  endif()

  if(QD_FMA STREQUAL "no")
    return()
  endif()

  if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    if(QD_FMA STREQUAL "auto")
      message(WARNING "Cross-compiling without CMAKE_CROSSCOMPILING_EMULATOR: FMA auto-detection is disabled. Set QD_FMA_EXPR/QD_FMS_EXPR for a manual override.")
      return()
    endif()
    message(FATAL_ERROR "QD_FMA=${QD_FMA} requires run tests while cross-compiling. Set CMAKE_CROSSCOMPILING_EMULATOR or QD_FMA_EXPR/QD_FMS_EXPR.")
  endif()

  if(QD_FMA STREQUAL "auto")
    qd3_fma_candidates_for_auto(candidates)
    set(fatal_if_missing FALSE)
  elseif(QD_FMA STREQUAL "yes")
    set(candidates ibm gnu c99 compiler)
    set(fatal_if_missing TRUE)
  else()
    set(candidates "${QD_FMA}")
    set(fatal_if_missing TRUE)
  endif()

  foreach(candidate IN LISTS candidates)
    message(STATUS "Checking FMA candidate '${candidate}'")
    qd3_try_fma_candidate("${candidate}" works)
    if(works)
      if(candidate STREQUAL "ibm")
        set(fma_expr "__fmadd(x,y,z)")
        set(fms_expr "__fmsub(x,y,z)")
        set(QD_VACPP_BUILTINS_H ON PARENT_SCOPE)
      elseif(candidate STREQUAL "gnu")
        set(fma_expr "__builtin_fma(x,y,z)")
        set(fms_expr "__builtin_fma(x,y,-z)")
      elseif(candidate STREQUAL "ia64")
        set(fma_expr "_Asm_fma(2, x,y,z)")
        set(fms_expr "_Asm_fms(2, x,y,z)")
      elseif(candidate STREQUAL "c99")
        set(fma_expr "fma(x,y,z)")
        set(fms_expr "fma(x,y,-z)")
      elseif(candidate STREQUAL "compiler")
        set(fma_expr "((x)*(y) + (z))")
        set(fms_expr "((x)*(y) - (z))")
      endif()
      set(QD_FMA_SELECTED "${fma_expr}" PARENT_SCOPE)
      set(QD_FMS_SELECTED "${fms_expr}" PARENT_SCOPE)
      set(QD_FMA_CONFIG_DEFINE "#define QD_FMA(x,y,z) ${fma_expr}" PARENT_SCOPE)
      set(QD_FMS_CONFIG_DEFINE "#define QD_FMS(x,y,z) ${fms_expr}" PARENT_SCOPE)
      return()
    endif()
  endforeach()

  if(fatal_if_missing)
    message(FATAL_ERROR "Cannot find a working fused multiply-add/subtract implementation for QD_FMA=${QD_FMA}")
  endif()
  message(WARNING "No working fused multiply-add/subtract implementation detected; QD_FMA/QD_FMS will be undefined")
endfunction()

function(qd3_check_fp_contract)
  check_cxx_compiler_flag("-fp-model strict" QD3_CXX_ACCEPTS_FP_MODEL_STRICT)
  if(QD3_CXX_ACCEPTS_FP_MODEL_STRICT)
    set(QD_FP_CONTRACT_FLAG "-fp-model strict" PARENT_SCOPE)
    set(QD_FP_CONTRACT_COMPILE_OPTIONS "-fp-model;strict" PARENT_SCOPE)
    return()
  endif()

  check_cxx_compiler_flag("-ffp-contract=off" QD3_CXX_ACCEPTS_FFP_CONTRACT_OFF)
  if(QD3_CXX_ACCEPTS_FFP_CONTRACT_OFF)
    set(QD_FP_CONTRACT_FLAG "-ffp-contract=off" PARENT_SCOPE)
    set(QD_FP_CONTRACT_COMPILE_OPTIONS "-ffp-contract=off" PARENT_SCOPE)
    return()
  endif()

  message(FATAL_ERROR "cannot find a compiler flag to suppress implicit FP contraction; need either -fp-model strict or -ffp-contract=off")
endfunction()

function(qd3_check_arch)
  set(QD_ARCH_CXX_FLAG "" PARENT_SCOPE)
  set(QD_ARCH_COMPILE_OPTIONS "" PARENT_SCOPE)
  if(QD_ARCH STREQUAL "generic")
    return()
  endif()

  if(QD_ARCH STREQUAL "x86-64-v3")
    string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" proc)
    if(NOT proc MATCHES "^(x86_64|amd64)$")
      message(FATAL_ERROR "QD_ARCH=x86-64-v3 is only valid for x86_64/amd64 targets; CMAKE_SYSTEM_PROCESSOR='${CMAKE_SYSTEM_PROCESSOR}'")
    endif()
    check_cxx_compiler_flag("-march=x86-64-v3" QD3_CXX_ACCEPTS_MARCH_X86_64_V3)
    if(NOT QD3_CXX_ACCEPTS_MARCH_X86_64_V3)
      message(FATAL_ERROR "${CMAKE_CXX_COMPILER} does not accept -march=x86-64-v3")
    endif()
    set(QD_ARCH_CXX_FLAG "-march=x86-64-v3" PARENT_SCOPE)
    set(QD_ARCH_COMPILE_OPTIONS "-march=x86-64-v3" PARENT_SCOPE)
    return()
  endif()

  if(QD_ARCH STREQUAL "native")
    if(CMAKE_CROSSCOMPILING)
      message(FATAL_ERROR "QD_ARCH=native is not allowed when cross-compiling")
    endif()
    check_cxx_compiler_flag("-march=native" QD3_CXX_ACCEPTS_MARCH_NATIVE)
    if(NOT QD3_CXX_ACCEPTS_MARCH_NATIVE)
      message(FATAL_ERROR "${CMAKE_CXX_COMPILER} does not accept -march=native")
    endif()
    set(QD_ARCH_CXX_FLAG "-march=native" PARENT_SCOPE)
    set(QD_ARCH_COMPILE_OPTIONS "-march=native" PARENT_SCOPE)
    return()
  endif()

  message(FATAL_ERROR "Unhandled QD_ARCH='${QD_ARCH}'")
endfunction()
