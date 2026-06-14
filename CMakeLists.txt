cmake_minimum_required(VERSION 3.20)

if(CMAKE_SOURCE_DIR STREQUAL CMAKE_BINARY_DIR)
  message(FATAL_ERROR "In-source builds are not supported. Use: cmake -S . -B build")
endif()

project(qd3 VERSION 1.1.0 LANGUAGES CXX C)

list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

include(GNUInstallDirs)
include(CTest)
include(CheckIncludeFileCXX)
include(CheckLibraryExists)
include(CheckSymbolExists)

include(QD3BuildOptions)
include(QD3CheckFMA)
include(QD3CheckEddReal)
include(QD3CheckFortran)
include(QD3Version)

qd3_define_options()
qd3_normalize_options()
qd3_assert_version_sync()

# Preserve Autotools ordering: FMA run probes happen before adding the
# contraction-suppression flag and before QD_ARCH target flags are applied.
qd3_check_fma()
qd3_check_fp_contract()
qd3_check_arch()

check_include_file_cxx(ieeefp.h HAVE_IEEEFP_H)

set(X86 OFF)
string(TOLOWER "${CMAKE_SYSTEM_PROCESSOR}" QD_SYSTEM_PROCESSOR_LOWER)
if(QD_SYSTEM_PROCESSOR_LOWER MATCHES "^(i[3-6]86|x86_64|amd64|k[0-9].*|athlon.*)$")
  set(X86 ON)
endif()
if(X86)
  check_include_file_cxx(fpu_control.h HAVE_FPU_CONTROL_H)
endif()

check_symbol_exists(gettimeofday "sys/time.h" HAVE_GETTIMEOFDAY)
check_symbol_exists(clock_gettime "time.h" QD_CLOCK_GETTIME_IN_LIBC)
set(QD_SYSTEM_LIBS "")
if(QD_CLOCK_GETTIME_IN_LIBC)
  set(HAVE_CLOCK_GETTIME ON)
else()
  check_library_exists(rt clock_gettime "" QD_CLOCK_GETTIME_IN_RT)
  if(QD_CLOCK_GETTIME_IN_RT)
    set(HAVE_CLOCK_GETTIME ON)
    list(APPEND QD_SYSTEM_LIBS rt)
  endif()
endif()
check_library_exists(m sqrt "" QD_HAVE_LIBM)
if(QD_HAVE_LIBM)
  list(APPEND QD_SYSTEM_LIBS m)
endif()

qd3_check_edd_real()
qd3_check_fortran()

set(QD_INLINE ${QD_ENABLE_INLINE})
set(QD_IEEE_ADD ${QD_ENABLE_IEEE_ADD})
set(QD_SLOPPY_MUL ${QD_ENABLE_SLOPPY_MUL})
set(QD_SLOPPY_DIV ${QD_ENABLE_SLOPPY_DIV})
set(HAVE_FORTRAN ${QD_HAVE_FORTRAN})

configure_file(
  "${CMAKE_SOURCE_DIR}/config.h.cmake.in"
  "${CMAKE_BINARY_DIR}/config.h"
  @ONLY
)
configure_file(
  "${CMAKE_SOURCE_DIR}/include/qd/qd_config.h.cmake.in"
  "${CMAKE_BINARY_DIR}/include/qd/qd_config.h"
  @ONLY
)

if(QD_HAVE_FORTRAN)
  set(ETIME "${QD_FORTRAN_ETIME}")
  configure_file(
    "${CMAKE_SOURCE_DIR}/fortran/second.f.in"
    "${CMAKE_BINARY_DIR}/fortran/second.f"
    @ONLY
  )
endif()

set(QD_BASE_SOURCES
  src/c_dd.cpp
  src/c_qd.cpp
  src/c_td.cpp
  src/dd_real.cpp
  src/dd_const.cpp
  src/td_real.cpp
  src/td_const.cpp
  src/fpu.cpp
  src/qd_real.cpp
  src/qd_const.cpp
  src/util.cpp
  src/bits.cpp
  src/util.h
)
set(QD_SOURCES ${QD_BASE_SOURCES})
if(QD_HAVE_EDD_REAL)
  list(APPEND QD_SOURCES
    src/edd_real.cpp
    src/edd_const.cpp
    src/edd_trans.cpp
    src/c_edd.cpp
  )
endif()

set(QD_PUBLIC_HEADERS
  "${CMAKE_SOURCE_DIR}/include/qd/c_dd.h"
  "${CMAKE_SOURCE_DIR}/include/qd/c_qd.h"
  "${CMAKE_SOURCE_DIR}/include/qd/c_td.h"
  "${CMAKE_SOURCE_DIR}/include/qd/dd_real.h"
  "${CMAKE_SOURCE_DIR}/include/qd/dd_inline.h"
  "${CMAKE_SOURCE_DIR}/include/qd/c_edd.h"
  "${CMAKE_SOURCE_DIR}/include/qd/edd_real.h"
  "${CMAKE_SOURCE_DIR}/include/qd/edd_inline.h"
  "${CMAKE_SOURCE_DIR}/include/qd/fpu.h"
  "${CMAKE_SOURCE_DIR}/include/qd/inline.h"
  "${CMAKE_SOURCE_DIR}/include/qd/td_real.h"
  "${CMAKE_SOURCE_DIR}/include/qd/td_inline.h"
  "${CMAKE_SOURCE_DIR}/include/qd/qd_real.h"
  "${CMAKE_SOURCE_DIR}/include/qd/qd_inline.h"
  "${CMAKE_SOURCE_DIR}/include/qd/bits.h"
)

if(QD3_ENABLE_MPFR_TESTS)
  find_package(PkgConfig REQUIRED)
  pkg_check_modules(MPFR REQUIRED IMPORTED_TARGET mpfr)
endif()

set(QD3_QA_COMPILE_OPTIONS "")
set(QD3_QA_LINK_OPTIONS "")

if(QD3_ENABLE_ASAN OR QD3_ENABLE_UBSAN OR QD3_ENABLE_COVERAGE)
  if(NOT CMAKE_CXX_COMPILER_ID MATCHES "GNU|Clang|AppleClang|IntelLLVM")
    message(FATAL_ERROR "QD3 QA instrumentation requires a GCC/Clang-compatible compiler")
  endif()
endif()

if(QD3_ENABLE_ASAN)
  list(APPEND QD3_QA_COMPILE_OPTIONS -fsanitize=address)
  list(APPEND QD3_QA_LINK_OPTIONS -fsanitize=address)
endif()

if(QD3_ENABLE_UBSAN)
  list(APPEND QD3_QA_COMPILE_OPTIONS -fsanitize=undefined -fno-sanitize-recover=all)
  list(APPEND QD3_QA_LINK_OPTIONS -fsanitize=undefined -fno-sanitize-recover=all)
endif()

if(QD3_ENABLE_COVERAGE)
  if(NOT BUILD_TESTING)
    message(FATAL_ERROR "QD3_ENABLE_COVERAGE requires BUILD_TESTING=ON")
  endif()
  list(APPEND QD3_QA_COMPILE_OPTIONS --coverage -O0 -g)
  list(APPEND QD3_QA_LINK_OPTIONS --coverage)
endif()

function(qd3_add_qa_options target)
  if(NOT QD3_QA_COMPILE_OPTIONS AND NOT QD3_QA_LINK_OPTIONS)
    return()
  endif()

  foreach(opt IN LISTS QD3_QA_COMPILE_OPTIONS)
    target_compile_options(${target} PRIVATE
      "$<$<COMPILE_LANGUAGE:C>:${opt}>"
      "$<$<COMPILE_LANGUAGE:CXX>:${opt}>"
    )
  endforeach()

  get_target_property(qd3_target_type ${target} TYPE)
  if(NOT qd3_target_type STREQUAL "STATIC_LIBRARY")
    foreach(opt IN LISTS QD3_QA_LINK_OPTIONS)
      target_link_options(${target} PRIVATE "${opt}")
    endforeach()
  endif()
endfunction()

function(qd3_add_cxx_private_options target)
  foreach(opt IN LISTS QD_FP_CONTRACT_COMPILE_OPTIONS QD_ARCH_COMPILE_OPTIONS)
    if(opt)
      target_compile_options(${target} PRIVATE "$<$<COMPILE_LANGUAGE:CXX>:${opt}>")
    endif()
  endforeach()
endfunction()

function(qd3_add_public_fp_contract_option target)
  if(QD_PROPAGATE_FP_CONTRACT_FLAG)
    foreach(opt IN LISTS QD_FP_CONTRACT_COMPILE_OPTIONS)
      if(opt)
        target_compile_options(${target} INTERFACE "$<$<COMPILE_LANGUAGE:CXX>:${opt}>")
      endif()
    endforeach()
  endif()
endfunction()

function(qd3_configure_qd_library target output_name)
  target_include_directories(${target} BEFORE
    PUBLIC
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}>
      $<BUILD_INTERFACE:${CMAKE_BINARY_DIR}/include>
      $<BUILD_INTERFACE:${CMAKE_SOURCE_DIR}/include>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}>
  )
  target_compile_features(${target} PUBLIC cxx_std_11)
  qd3_add_cxx_private_options(${target})
  qd3_add_qa_options(${target})
  qd3_add_public_fp_contract_option(${target})
  target_link_libraries(${target} PUBLIC ${QD_SYSTEM_LIBS})
  set_target_properties(${target} PROPERTIES OUTPUT_NAME ${output_name})
endfunction()

function(qd3_configure_internal_target target)
  target_include_directories(${target} BEFORE PRIVATE
    "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include"
    "${CMAKE_SOURCE_DIR}/tests"
  )
  target_compile_features(${target} PRIVATE cxx_std_11)
  qd3_add_cxx_private_options(${target})
  qd3_add_qa_options(${target})
endfunction()

set(QD3_INSTALL_TARGETS "")

if(BUILD_SHARED_LIBS)
  add_library(qd SHARED ${QD_SOURCES})
  set_target_properties(qd PROPERTIES VERSION 0.0.0 SOVERSION 0)
else()
  add_library(qd STATIC ${QD_SOURCES})
endif()
qd3_configure_qd_library(qd qd)
add_library(QD3::qd ALIAS qd)
list(APPEND QD3_INSTALL_TARGETS qd)

if(BUILD_SHARED_LIBS AND QD_BUILD_STATIC)
  add_library(qd_static STATIC ${QD_SOURCES})
  qd3_configure_qd_library(qd_static qd)
  add_library(QD3::qd_static ALIAS qd_static)
  list(APPEND QD3_INSTALL_TARGETS qd_static)
endif()

if(QD_HAVE_FORTRAN)
  set(QD_FORTRAN_SOURCES
    fortran/ddext.f
    fortran/ddmod.f
    fortran/qdext.f
    fortran/qdmod.f
    fortran/tdext.f
    fortran/tdmod.f
    fortran/f_dd.cpp
    fortran/f_qd.cpp
    fortran/f_td.cpp
  )
  set(QD_FORTRAN_FREE_FORM_SOURCES
    fortran/ddext.f
    fortran/ddmod.f
    fortran/qdext.f
    fortran/qdmod.f
    fortran/tdext.f
    fortran/tdmod.f
  )
  set_source_files_properties(${QD_FORTRAN_FREE_FORM_SOURCES} PROPERTIES Fortran_FORMAT FREE)

  if(BUILD_SHARED_LIBS)
    add_library(qdmod SHARED ${QD_FORTRAN_SOURCES})
    set_target_properties(qdmod PROPERTIES VERSION 0.0.0 SOVERSION 0)
  else()
    add_library(qdmod STATIC ${QD_FORTRAN_SOURCES})
  endif()
  qd3_configure_qd_library(qdmod qdmod)
  target_include_directories(qdmod PUBLIC
    $<BUILD_INTERFACE:${QD_FORTRAN_MODULE_DIR}>
    $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qd>
  )
  set_target_properties(qdmod PROPERTIES Fortran_MODULE_DIRECTORY "${QD_FORTRAN_MODULE_DIR}")
  target_link_libraries(qdmod PUBLIC qd)
  add_library(QD3::qdmod ALIAS qdmod)
  list(APPEND QD3_INSTALL_TARGETS qdmod)

  if(BUILD_SHARED_LIBS AND QD_BUILD_STATIC)
    set(QD_FORTRAN_MODULE_DIR_STATIC "${CMAKE_BINARY_DIR}/fortran/modules_static")
    file(MAKE_DIRECTORY "${QD_FORTRAN_MODULE_DIR_STATIC}")
    add_library(qdmod_static STATIC ${QD_FORTRAN_SOURCES})
    qd3_configure_qd_library(qdmod_static qdmod)
    target_include_directories(qdmod_static PUBLIC
      $<BUILD_INTERFACE:${QD_FORTRAN_MODULE_DIR_STATIC}>
      $<INSTALL_INTERFACE:${CMAKE_INSTALL_INCLUDEDIR}/qd>
    )
    set_target_properties(qdmod_static PROPERTIES Fortran_MODULE_DIRECTORY "${QD_FORTRAN_MODULE_DIR_STATIC}")
    target_link_libraries(qdmod_static PUBLIC qd_static)
    add_library(QD3::qdmod_static ALIAS qdmod_static)
    list(APPEND QD3_INSTALL_TARGETS qdmod_static)
  endif()

  if(BUILD_SHARED_LIBS)
    add_library(qd_f_main SHARED fortran/main.cpp)
    set_target_properties(qd_f_main PROPERTIES VERSION 0.0.0 SOVERSION 0)
  else()
    add_library(qd_f_main STATIC fortran/main.cpp)
  endif()
  qd3_configure_qd_library(qd_f_main qd_f_main)
  target_link_libraries(qd_f_main PUBLIC qdmod qd)
  add_library(QD3::qd_f_main ALIAS qd_f_main)
  list(APPEND QD3_INSTALL_TARGETS qd_f_main)

  if(BUILD_SHARED_LIBS AND QD_BUILD_STATIC)
    add_library(qd_f_main_static STATIC fortran/main.cpp)
    qd3_configure_qd_library(qd_f_main_static qd_f_main)
    target_link_libraries(qd_f_main_static PUBLIC qdmod_static qd_static)
    add_library(QD3::qd_f_main_static ALIAS qd_f_main_static)
    list(APPEND QD3_INSTALL_TARGETS qd_f_main_static)
  endif()
endif()

if(BUILD_TESTING)
  configure_file(
    "${CMAKE_SOURCE_DIR}/tests/coeff.dat"
    "${CMAKE_BINARY_DIR}/tests/coeff.dat"
    COPYONLY
  )

  add_executable(qd_test tests/qd_test.cpp)
  qd3_configure_internal_target(qd_test)
  target_link_libraries(qd_test PRIVATE qd)
  add_test(NAME qd_test COMMAND qd_test)
  add_test(NAME qd_test_dd COMMAND qd_test -dd)
  add_test(NAME qd_test_td COMMAND qd_test -td)
  add_test(NAME qd_test_qd COMMAND qd_test -qd)
  if(QD_HAVE_EDD_REAL)
    add_test(NAME qd_test_edd COMMAND qd_test -edd)
  endif()

  add_executable(pslq_test tests/pslq_test.cpp tests/tictoc.cpp tests/tictoc.h tests/pslq.h)
  qd3_configure_internal_target(pslq_test)
  target_link_libraries(pslq_test PRIVATE qd)
  add_test(NAME pslq_test COMMAND pslq_test)

  add_executable(c_test tests/c_test.c)
  target_include_directories(c_test BEFORE PRIVATE
    "${CMAKE_BINARY_DIR}"
    "${CMAKE_BINARY_DIR}/include"
    "${CMAKE_SOURCE_DIR}/include"
  )
  qd3_add_qa_options(c_test)
  set_target_properties(c_test PROPERTIES LINKER_LANGUAGE CXX)
  target_link_libraries(c_test PRIVATE qd)
  add_test(NAME c_test COMMAND c_test)

  add_executable(huge tests/huge.cpp)
  qd3_configure_internal_target(huge)
  target_link_libraries(huge PRIVATE qd)
  add_test(NAME huge COMMAND huge)

  add_executable(qd_timer EXCLUDE_FROM_ALL tests/qd_timer.cpp tests/tictoc.cpp tests/tictoc.h)
  qd3_configure_internal_target(qd_timer)
  target_link_libraries(qd_timer PRIVATE qd)

  add_executable(quadt_test EXCLUDE_FROM_ALL tests/quadt_test.cpp tests/tictoc.cpp tests/tictoc.h tests/quadt.h)
  qd3_configure_internal_target(quadt_test)
  target_link_libraries(quadt_test PRIVATE qd)
  set_target_properties(quadt_test PROPERTIES RUNTIME_OUTPUT_DIRECTORY "${CMAKE_BINARY_DIR}/tests")

  set(QD_DEMO_DEPS qd_test pslq_test c_test huge qd_timer quadt_test)

  if(QD_HAVE_FORTRAN)
    set_source_files_properties(tests/f_test.f PROPERTIES Fortran_FORMAT FREE)
    add_executable(f_test tests/f_test.f)
    target_include_directories(f_test BEFORE PRIVATE
      "${QD_FORTRAN_MODULE_DIR}"
      "${CMAKE_BINARY_DIR}"
      "${CMAKE_BINARY_DIR}/include"
      "${CMAKE_SOURCE_DIR}/include"
    )
    qd3_add_qa_options(f_test)
    set_target_properties(f_test PROPERTIES LINKER_LANGUAGE CXX)
    target_link_libraries(f_test PRIVATE qdmod qd_f_main qd)
    add_test(NAME f_test COMMAND f_test)
    list(APPEND QD_DEMO_DEPS f_test)
  endif()

  add_custom_target(demo DEPENDS ${QD_DEMO_DEPS})
  add_custom_target(time
    COMMAND qd_timer
    DEPENDS qd_timer
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
  )
endif()

if(QD3_ENABLE_COVERAGE)
  find_program(LCOV_EXECUTABLE lcov REQUIRED)
  find_program(GENHTML_EXECUTABLE genhtml REQUIRED)

  add_custom_target(coverage
    COMMAND "${LCOV_EXECUTABLE}" --zerocounters --directory "${CMAKE_BINARY_DIR}"
    COMMAND "${CMAKE_CTEST_COMMAND}" --test-dir "${CMAKE_BINARY_DIR}" --output-on-failure
    COMMAND "${LCOV_EXECUTABLE}" --capture --directory "${CMAKE_BINARY_DIR}" --output-file "${CMAKE_BINARY_DIR}/coverage.info"
    COMMAND "${LCOV_EXECUTABLE}" --extract "${CMAKE_BINARY_DIR}/coverage.info"
            "${CMAKE_SOURCE_DIR}/src/*"
            "${CMAKE_SOURCE_DIR}/include/qd/*_inline.h"
            --output-file "${CMAKE_BINARY_DIR}/coverage.filtered.info"
            --ignore-errors unused
    COMMAND "${GENHTML_EXECUTABLE}" "${CMAKE_BINARY_DIR}/coverage.filtered.info" --output-directory "${CMAKE_BINARY_DIR}/coverage-html"
    DEPENDS ${QD_DEMO_DEPS}
    BYPRODUCTS
      "${CMAKE_BINARY_DIR}/coverage.info"
      "${CMAKE_BINARY_DIR}/coverage.filtered.info"
      "${CMAKE_BINARY_DIR}/coverage-html/index.html"
    WORKING_DIRECTORY "${CMAKE_BINARY_DIR}"
    USES_TERMINAL
    COMMENT "Running tests and generating lcov HTML coverage"
  )
endif()

find_program(QD_MAKE_PROGRAM NAMES gmake make REQUIRED)
add_custom_target(doc
  COMMAND ${QD_MAKE_PROGRAM} -C "${CMAKE_SOURCE_DIR}/docs" qd.pdf
  COMMAND ${QD_MAKE_PROGRAM} -C "${CMAKE_SOURCE_DIR}/docs" edd.pdf
  COMMAND ${QD_MAKE_PROGRAM} -C "${CMAKE_SOURCE_DIR}/docs" td.pdf
  VERBATIM
)

add_custom_target(dist-cmake
  COMMAND sh -c "git -C '${CMAKE_SOURCE_DIR}' archive --format=tar --prefix='qd3-${PROJECT_VERSION}/' HEAD | xz -c > '${CMAKE_BINARY_DIR}/qd3-${PROJECT_VERSION}.tar.xz'"
  VERBATIM
)

set(QD_PC_LIBS "")
set(QD_CONFIG_LIBS "")
foreach(lib IN LISTS QD_SYSTEM_LIBS)
  if(lib STREQUAL "m" OR lib STREQUAL "rt")
    string(APPEND QD_PC_LIBS " -l${lib}")
    string(APPEND QD_CONFIG_LIBS " -l${lib}")
  endif()
endforeach()
if(QD_HAVE_FORTRAN)
  set(QD_CONFIG_FC "${CMAKE_Fortran_COMPILER}")
  set(QD_CONFIG_FCFLAGS "")
  set(QD_CONFIG_FC_LIBS "")
  set(QD_CONFIG_MODFLAG "${QD_FORTRAN_MODULE_FLAG}")
else()
  set(QD_CONFIG_FC "")
  set(QD_CONFIG_FCFLAGS "")
  set(QD_CONFIG_FC_LIBS "")
  set(QD_CONFIG_MODFLAG "")
endif()
configure_file("${CMAKE_SOURCE_DIR}/cmake/qd.pc.cmake.in" "${CMAKE_BINARY_DIR}/qd.pc" @ONLY)
configure_file("${CMAKE_SOURCE_DIR}/cmake/qd-config.cmake.in" "${CMAKE_BINARY_DIR}/qd-config" @ONLY)
file(CHMOD "${CMAKE_BINARY_DIR}/qd-config"
  PERMISSIONS OWNER_READ OWNER_WRITE OWNER_EXECUTE GROUP_READ GROUP_EXECUTE WORLD_READ WORLD_EXECUTE
)

include(QD3Install)
qd3_install_targets()

message(STATUS "")
message(STATUS "QD3 configuration summary")
message(STATUS "  version                    = ${PROJECT_VERSION}")
message(STATUS "  generator                  = ${CMAKE_GENERATOR}")
message(STATUS "  build type                 = ${CMAKE_BUILD_TYPE}")
message(STATUS "  install prefix             = ${CMAKE_INSTALL_PREFIX}")
message(STATUS "  C++ compiler               = ${CMAKE_CXX_COMPILER}")
message(STATUS "  C compiler                 = ${CMAKE_C_COMPILER}")
if(QD_HAVE_FORTRAN)
  message(STATUS "  Fortran compiler           = ${CMAKE_Fortran_COMPILER}")
else()
  message(STATUS "  Fortran compiler           = disabled")
endif()
message(STATUS "  QD_ENABLE_INLINE           = ${QD_ENABLE_INLINE}")
message(STATUS "  QD_ENABLE_IEEE_ADD         = ${QD_ENABLE_IEEE_ADD}")
message(STATUS "  QD_ENABLE_SLOPPY_MUL       = ${QD_ENABLE_SLOPPY_MUL}")
message(STATUS "  QD_ENABLE_SLOPPY_DIV       = ${QD_ENABLE_SLOPPY_DIV}")
message(STATUS "  QD_FMA requested           = ${QD_FMA}")
message(STATUS "  QD_FMA selected            = ${QD_FMA_SELECTED}")
message(STATUS "  QD_FMS selected            = ${QD_FMS_SELECTED}")
message(STATUS "  FP contraction flag        = ${QD_FP_CONTRACT_FLAG}")
message(STATUS "  QD_ARCH                    = ${QD_ARCH} (${QD_ARCH_CXX_FLAG})")
message(STATUS "  edd_real status            = ${QD_EDD_STATUS}")
message(STATUS "  BUILD_SHARED_LIBS          = ${BUILD_SHARED_LIBS}")
message(STATUS "  QD_BUILD_STATIC            = ${QD_BUILD_STATIC}")
message(STATUS "  tests enabled              = ${BUILD_TESTING}")
message(STATUS "  MPFR oracle tests          = ${QD3_ENABLE_MPFR_TESTS}")
message(STATUS "  AddressSanitizer           = ${QD3_ENABLE_ASAN}")
message(STATUS "  UndefinedBehaviorSanitizer = ${QD3_ENABLE_UBSAN}")
message(STATUS "  coverage instrumentation   = ${QD3_ENABLE_COVERAGE}")
message(STATUS "")
