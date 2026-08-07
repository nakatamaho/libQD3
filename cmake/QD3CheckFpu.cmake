include(CheckCXXSourceRuns)

function(qd3_check_fpu_fix_80bit)
  set(QD_HAVE_FPU_FIX_80BIT OFF PARENT_SCOPE)
  set(QD_FPU_FIX_80BIT_STATUS "no-op (non-x86 target)" PARENT_SCOPE)

  if(NOT X86)
    return()
  endif()

  set(QD_FPU_FIX_80BIT_STATUS "not available" PARENT_SCOPE)

  if(CMAKE_CROSSCOMPILING AND NOT CMAKE_CROSSCOMPILING_EMULATOR)
    set(QD_FPU_FIX_80BIT_STATUS "not probed (cross-compiling)" PARENT_SCOPE)
    return()
  endif()

  set(src_fpu_80bit [[
#if !defined(__i386__) && !defined(__x86_64__)
int main() { return 1; }
#else
#if !defined(__LDBL_MANT_DIG__) || __LDBL_MANT_DIG__ != 64
int main() { return 1; }
#else
#if defined(__has_include)
#  if __has_include(<fpu_control.h>)
#    include <fpu_control.h>
#  endif
#endif
#include <cmath>

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) __asm__ __volatile__("fnstcw %0":"=m"(x))
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) __asm__ __volatile__("fldcw %0": :"m"(x))
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

int main() {
  volatile unsigned short cw = 0;
  volatile unsigned short new_cw = 0;
  _FPU_GETCW(cw);
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_EXTENDED;
  _FPU_SETCW(new_cw);

  volatile long double one = 1.0L;
  volatile long double inc = std::ldexp(1.0L, -63);
  volatile long double sum = one + inc;
  volatile long double diff = sum - one;

  _FPU_SETCW(cw);
  return diff == inc ? 0 : 1;
}
#endif
#endif
]])

  set(saved_required_flags "${CMAKE_REQUIRED_FLAGS}")
  set(saved_required_libraries "${CMAKE_REQUIRED_LIBRARIES}")
  set(CMAKE_REQUIRED_FLAGS "${CMAKE_REQUIRED_FLAGS} ${QD_FP_CONTRACT_FLAG} ${QD_ARCH_CXX_FLAG}")
  set(CMAKE_REQUIRED_LIBRARIES ${CMAKE_REQUIRED_LIBRARIES} ${QD_SYSTEM_LIBS})
  check_cxx_source_runs("${src_fpu_80bit}" QD3_FPU_FIX_80BIT_RUNS)
  set(CMAKE_REQUIRED_FLAGS "${saved_required_flags}")
  set(CMAKE_REQUIRED_LIBRARIES "${saved_required_libraries}")

  if(QD3_FPU_FIX_80BIT_RUNS)
    set(QD_HAVE_FPU_FIX_80BIT ON PARENT_SCOPE)
    set(QD_FPU_FIX_80BIT_STATUS "available" PARENT_SCOPE)
  endif()
endfunction()
