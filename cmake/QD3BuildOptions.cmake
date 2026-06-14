# CMake build options for libQD3.

function(qd3_normalize_bool_option var)
  if(${var})
    set(${var} ON PARENT_SCOPE)
  else()
    set(${var} OFF PARENT_SCOPE)
  endif()
endfunction()

function(qd3_validate_choice var allowed)
  set(value "${${var}}")
  list(FIND allowed "${value}" found)
  if(found EQUAL -1)
    string(REPLACE ";" ", " allowed_text "${allowed}")
    message(FATAL_ERROR "Invalid ${var}='${value}'. Allowed values: ${allowed_text}")
  endif()
  set(${var} "${value}" PARENT_SCOPE)
endfunction()

function(qd3_define_options)
  option(QD_ENABLE_INLINE "Inline commonly used functions" ON)
  option(QD_ENABLE_IEEE_ADD "Use IEEE-style addition error bound" OFF)
  option(QD_ENABLE_SLOPPY_MUL "Use faster sloppy multiplication" ON)
  option(QD_ENABLE_SLOPPY_DIV "Use faster sloppy division" ON)

  set(QD_FMA "auto" CACHE STRING "FMA mode: auto, yes, no, gnu, c99, ibm, ia64, compiler")
  set_property(CACHE QD_FMA PROPERTY STRINGS auto yes no gnu c99 ibm ia64 compiler)

  set(QD_ARCH "generic" CACHE STRING "Code generation target: generic, x86-64-v3, native")
  set_property(CACHE QD_ARCH PROPERTY STRINGS generic x86-64-v3 native)

  set(QD_BUILD_FORTRAN "AUTO" CACHE STRING "Build Fortran interfaces: AUTO, ON, OFF")
  set_property(CACHE QD_BUILD_FORTRAN PROPERTY STRINGS AUTO ON OFF)

  set(QD_ENABLE_EDD_REAL "AUTO" CACHE STRING "Build edd_real sources: AUTO, ON, OFF")
  set_property(CACHE QD_ENABLE_EDD_REAL PROPERTY STRINGS AUTO ON OFF)

  option(BUILD_SHARED_LIBS "Build shared libraries" ON)
  option(QD_BUILD_STATIC "Build static libraries" ON)
  option(QD_PROPAGATE_FP_CONTRACT_FLAG "Propagate the selected FP-contraction flag to consumers" OFF)

  set(QD_FMA_EXPR "" CACHE STRING "Manual QD_FMA(x,y,z) replacement expression")
  set(QD_FMS_EXPR "" CACHE STRING "Manual QD_FMS(x,y,z) replacement expression")
endfunction()

function(qd3_normalize_options)
  qd3_normalize_bool_option(QD_ENABLE_INLINE)
  qd3_normalize_bool_option(QD_ENABLE_IEEE_ADD)
  qd3_normalize_bool_option(QD_ENABLE_SLOPPY_MUL)
  qd3_normalize_bool_option(QD_ENABLE_SLOPPY_DIV)
  qd3_normalize_bool_option(BUILD_SHARED_LIBS)
  qd3_normalize_bool_option(QD_BUILD_STATIC)
  qd3_normalize_bool_option(QD_PROPAGATE_FP_CONTRACT_FLAG)

  string(TOLOWER "${QD_FMA}" QD_FMA)
  set(allowed_fma auto yes no gnu c99 ibm ia64 compiler)
  qd3_validate_choice(QD_FMA "${allowed_fma}")

  string(TOLOWER "${QD_ARCH}" QD_ARCH)
  set(allowed_arch generic x86-64-v3 native)
  qd3_validate_choice(QD_ARCH "${allowed_arch}")

  string(TOUPPER "${QD_BUILD_FORTRAN}" QD_BUILD_FORTRAN)
  set(allowed_tristate AUTO ON OFF)
  qd3_validate_choice(QD_BUILD_FORTRAN "${allowed_tristate}")

  string(TOUPPER "${QD_ENABLE_EDD_REAL}" QD_ENABLE_EDD_REAL)
  qd3_validate_choice(QD_ENABLE_EDD_REAL "${allowed_tristate}")

  if(NOT BUILD_SHARED_LIBS AND NOT QD_BUILD_STATIC)
    message(FATAL_ERROR "At least one of BUILD_SHARED_LIBS or QD_BUILD_STATIC must be ON")
  endif()

  set(QD_FMA "${QD_FMA}" PARENT_SCOPE)
  set(QD_ARCH "${QD_ARCH}" PARENT_SCOPE)
  set(QD_BUILD_FORTRAN "${QD_BUILD_FORTRAN}" PARENT_SCOPE)
  set(QD_ENABLE_EDD_REAL "${QD_ENABLE_EDD_REAL}" PARENT_SCOPE)
endfunction()
