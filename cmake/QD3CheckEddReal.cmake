include(CheckCXXSourceCompiles)

function(qd3_check_edd_real)
  set(QD_HAVE_EDD_REAL OFF PARENT_SCOPE)
  set(QD_EDD_STATUS "disabled" PARENT_SCOPE)
  set(QD_EDD_FLT64X_MANT_DIG "" PARENT_SCOPE)
  set(QD_EDD_FLT64X_MAX_EXP "" PARENT_SCOPE)
  set(QD_EDD_SPLIT_BITS "" PARENT_SCOPE)
  set(QD_EDD_SPLIT_SCALE_BITS "" PARENT_SCOPE)

  if(QD_ENABLE_EDD_REAL STREQUAL "OFF")
    return()
  endif()

  set(src [[
#ifndef __GNUC__
#error edd_real requires GNU C++
#endif
#if !defined(__x86_64__) && !defined(__i386__)
#error edd_real Phase 1 only supports x86 GNU _Float64x binary80 targets
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
  int cls = __builtin_isfinite(z);
  return cls ? 0 : 1;
}
]])

  set(saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${QD_FP_CONTRACT_FLAG} ${QD_ARCH_CXX_FLAG}")
  check_cxx_source_compiles("${src}" QD3_HAVE_GNU_FLOAT64X_BINARY80)
  set(CMAKE_REQUIRED_FLAGS "${saved_required_flags}")

  if(QD3_HAVE_GNU_FLOAT64X_BINARY80)
    set(QD_HAVE_EDD_REAL ON PARENT_SCOPE)
    set(QD_EDD_STATUS "enabled" PARENT_SCOPE)
    set(QD_EDD_FLT64X_MANT_DIG 64 PARENT_SCOPE)
    set(QD_EDD_FLT64X_MAX_EXP 16384 PARENT_SCOPE)
    set(QD_EDD_SPLIT_BITS 32 PARENT_SCOPE)
    set(QD_EDD_SPLIT_SCALE_BITS 33 PARENT_SCOPE)
  elseif(QD_ENABLE_EDD_REAL STREQUAL "ON")
    message(FATAL_ERROR "QD_ENABLE_EDD_REAL=ON requested, but GNU C++ x86 binary80 _Float64x support was not detected")
  else()
    message(WARNING "edd_real disabled: GNU C++ _Float64x binary80 support was not detected")
  endif()
endfunction()
