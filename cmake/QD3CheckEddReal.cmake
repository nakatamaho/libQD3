include(CheckCXXSourceCompiles)

function(qd3_check_edd_real)
  set(QD_HAVE_EDD_REAL OFF PARENT_SCOPE)
  set(QD_EDD_STATUS "disabled" PARENT_SCOPE)
  set(QD_EDD_WORD_IS_FLOAT64X OFF PARENT_SCOPE)
  set(QD_EDD_WORD_IS_LONG_DOUBLE OFF PARENT_SCOPE)
  set(QD_EDD_FLT64X_MANT_DIG "" PARENT_SCOPE)
  set(QD_EDD_FLT64X_MAX_EXP "" PARENT_SCOPE)
  set(QD_EDD_WORD_MANT_DIG "" PARENT_SCOPE)
  set(QD_EDD_WORD_MAX_EXP "" PARENT_SCOPE)
  set(QD_EDD_SPLIT_BITS "" PARENT_SCOPE)
  set(QD_EDD_SPLIT_SCALE_BITS "" PARENT_SCOPE)

  if(QD_ENABLE_EDD_REAL STREQUAL "OFF")
    return()
  endif()

  set(src_float64x [[
#ifndef __GNUC__
#error edd_real requires GNU C++
#endif
#if !defined(__x86_64__) && !defined(__i386__)
#error edd_real Phase 1 only supports x86 GNU binary80 targets
#endif
#ifndef __FLT64X_MANT_DIG__
#error missing __FLT64X_MANT_DIG__
#endif
#ifndef __FLT64X_MAX_EXP__
#error missing __FLT64X_MAX_EXP__
#endif
#if __FLT64X_MANT_DIG__ != 64
#error edd_real Phase 1 requires binary80 _Float64x with 64-bit mantissa
#endif
static_assert(__FLT64X_MANT_DIG__ == 64, "binary80 _Float64x required");
int main() {
  _Float64x x = (_Float64x)1.0;
  _Float64x y = __builtin_sqrtf64x(x);
  _Float64x z = __builtin_ldexpf64x(y, 32);
  z += __builtin_expf64x(-__builtin_logf64x(x));
  int e = 0;
  z += __builtin_frexpf64x(z, &e);
  z += __builtin_atan2f64x(z, x);
  z += __builtin_log10f64x(x);
  z += __builtin_floorf64x(x);
  int cls = __builtin_isfinite(z);
  return cls ? 0 : 1;
}
]])

  set(saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
  set(saved_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${QD_FP_CONTRACT_FLAG} ${QD_ARCH_CXX_FLAG}")
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${QD_SYSTEM_LIBS})
  check_cxx_source_compiles("${src_float64x}" QD3_HAVE_GNU_FLOAT64X_BINARY80)
  set(CMAKE_REQUIRED_FLAGS "${saved_required_flags}")
  set(CMAKE_REQUIRED_LIBRARIES "${saved_required_libraries}")

  if(NOT QD3_HAVE_GNU_FLOAT64X_BINARY80)
    set(src_long_double [[
#ifndef __GNUC__
#error edd_real requires GNU C++
#endif
#if !defined(__x86_64__) && !defined(__i386__)
#error edd_real Phase 1 only supports x86 GNU binary80 targets
#endif
#if !defined(__LDBL_MANT_DIG__) || !defined(__LDBL_MAX_EXP__)
#error missing long double format macros
#endif
#if __LDBL_MANT_DIG__ != 64 || __LDBL_MAX_EXP__ != 16384
#error edd_real requires binary80 long double with 64-bit mantissa
#endif
#ifndef __FLT64X_MANT_DIG__
#error long double backend requires _Float64x interop for the existing c_edd ABI
#endif
#include <cmath>
#include <limits>
int main() {
  long double x = 1.0L;
  _Float64x fx = (_Float64x)x;
  x = (long double)fx;
  long double y = std::sqrt(x);
  long double z = std::ldexp(y, 32);
  z += std::exp(-std::log(x));
  int e = 0;
  z += std::frexp(z, &e);
  z += std::atan2(z, x);
  z += std::log10(x);
  z += std::floor(x);
  return std::isfinite(z) ? 0 : 1;
}
]])

    set(saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
    set(saved_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
    set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${QD_FP_CONTRACT_FLAG} ${QD_ARCH_CXX_FLAG}")
    set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${QD_SYSTEM_LIBS})
    check_cxx_source_compiles("${src_long_double}" QD3_HAVE_GNU_LONG_DOUBLE_BINARY80)
    set(CMAKE_REQUIRED_FLAGS "${saved_required_flags}")
    set(CMAKE_REQUIRED_LIBRARIES "${saved_required_libraries}")
  endif()

  if(QD3_HAVE_GNU_FLOAT64X_BINARY80)
    set(QD_HAVE_EDD_REAL ON PARENT_SCOPE)
    set(QD_EDD_STATUS "enabled (_Float64x binary80)" PARENT_SCOPE)
    set(QD_EDD_WORD_IS_FLOAT64X ON PARENT_SCOPE)
    set(QD_EDD_WORD_IS_LONG_DOUBLE OFF PARENT_SCOPE)
    set(QD_EDD_FLT64X_MANT_DIG 64 PARENT_SCOPE)
    set(QD_EDD_FLT64X_MAX_EXP 16384 PARENT_SCOPE)
    set(QD_EDD_WORD_MANT_DIG 64 PARENT_SCOPE)
    set(QD_EDD_WORD_MAX_EXP 16384 PARENT_SCOPE)
    set(QD_EDD_SPLIT_BITS 32 PARENT_SCOPE)
    set(QD_EDD_SPLIT_SCALE_BITS 33 PARENT_SCOPE)
  elseif(QD3_HAVE_GNU_LONG_DOUBLE_BINARY80)
    set(QD_HAVE_EDD_REAL ON PARENT_SCOPE)
    set(QD_EDD_STATUS "enabled (long double binary80)" PARENT_SCOPE)
    set(QD_EDD_WORD_IS_FLOAT64X OFF PARENT_SCOPE)
    set(QD_EDD_WORD_IS_LONG_DOUBLE ON PARENT_SCOPE)
    set(QD_EDD_FLT64X_MANT_DIG 64 PARENT_SCOPE)
    set(QD_EDD_FLT64X_MAX_EXP 16384 PARENT_SCOPE)
    set(QD_EDD_WORD_MANT_DIG 64 PARENT_SCOPE)
    set(QD_EDD_WORD_MAX_EXP 16384 PARENT_SCOPE)
    set(QD_EDD_SPLIT_BITS 32 PARENT_SCOPE)
    set(QD_EDD_SPLIT_SCALE_BITS 33 PARENT_SCOPE)
  elseif(QD_ENABLE_EDD_REAL STREQUAL "ON")
    message(FATAL_ERROR "QD_ENABLE_EDD_REAL=ON requested, but GNU C++ x86 binary80 support was not detected")
  else()
    message(WARNING "edd_real disabled: GNU C++ binary80 support was not detected")
  endif()
endfunction()
