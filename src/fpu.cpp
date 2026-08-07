/*
 * src/fpu.cc
 *
 * This work was supported by the Director, Office of Science, Division
 * of Mathematical, Information, and Computational Sciences of the
 * U.S. Department of Energy under contract number DE-AC03-76SF00098.
 *
 * Copyright (c) 2000-2001
 *
 * Contains functions to set and restore the round-to-double flag in the
 * control word of a x86 FPU.
 */

#include "config.h"
#include <qd/fpu.h>

#ifdef X86
#ifdef _WIN32
#include <float.h>
#endif

#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
#ifdef HAVE_FPU_CONTROL_H
#include <fpu_control.h>
#endif

#ifndef _FPU_GETCW
#define _FPU_GETCW(x) asm volatile ("fnstcw %0":"=m" (x));
#endif

#ifndef _FPU_SETCW
#define _FPU_SETCW(x) asm volatile ("fldcw %0": :"m" (x));
#endif
#endif

#ifndef _FPU_EXTENDED
#define _FPU_EXTENDED 0x0300
#endif

#ifndef _FPU_DOUBLE
#define _FPU_DOUBLE 0x0200
#endif

#endif /* X86 */

extern "C" {

void fpu_fix_start(unsigned int *old_cw) {
#ifdef X86
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);

  if (old_cw) {
    *old_cw = cw;
  }
#elif defined(_WIN32)
#ifdef __BORLANDC__
  /* Win 32 Borland C */
  unsigned short cw = _control87(0, 0);
  _control87(0x0200, 0x0300);
  if (old_cw) {
    *old_cw = cw;
  }
#else
  /* Win 32 MSVC */
  unsigned int cw = _control87(0, 0);
  _control87(0x00010000, 0x00030000);
  if (old_cw) {
    *old_cw = cw;
  }
#endif
#else
  /* Linux */
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_DOUBLE;
  _FPU_SETCW(new_cw);

  if (old_cw) {
    *old_cw = cw;
  }
#endif
#endif
}

void fpu_fix_start_80bit(unsigned int *old_cw) {
#ifdef X86
#if defined(__GNUC__) && (defined(__i386__) || defined(__x86_64__))
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  if (old_cw) {
    *old_cw = cw;
  }
#ifdef QD_HAVE_FPU_FIX_80BIT
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_EXTENDED;
  _FPU_SETCW(new_cw);
#endif
#elif defined(_WIN32)
#ifdef __BORLANDC__
  /* Win 32 Borland C */
  unsigned short cw = _control87(0, 0);
  if (old_cw) {
    *old_cw = cw;
  }
#ifdef QD_HAVE_FPU_FIX_80BIT
  _control87(0x0300, 0x0300);
#endif
#else
  /* Win 32 MSVC */
  unsigned int cw = _control87(0, 0);
  if (old_cw) {
    *old_cw = cw;
  }
#ifdef QD_HAVE_FPU_FIX_80BIT
  _control87(0x00000000, 0x00030000);
#endif
#endif
#else
  /* Linux */
  volatile unsigned short cw, new_cw;
  _FPU_GETCW(cw);

  if (old_cw) {
    *old_cw = cw;
  }
#ifdef QD_HAVE_FPU_FIX_80BIT
  new_cw = (cw & ~_FPU_EXTENDED) | _FPU_EXTENDED;
  _FPU_SETCW(new_cw);
#endif
#endif
#else
  if (old_cw) {
    *old_cw = 0;
  }
#endif
}

void fpu_fix_end(unsigned int *old_cw) {
#ifdef X86
#ifdef _WIN32

#ifdef __BORLANDC__
  /* Win 32 Borland C */
  if (old_cw) {
    unsigned short cw = (unsigned short) *old_cw;
    _control87(cw, 0xFFFF);
  }
#else
  /* Win 32 MSVC */
  if (old_cw) {
    _control87(*old_cw, 0xFFFFFFFF);
  }
#endif

#else
  /* Linux */
  if (old_cw) {
    int cw;
    cw = *old_cw;
    _FPU_SETCW(cw);
  }
#endif
#endif
}

#ifdef HAVE_FORTRAN

#define f_fpu_fix_start FC_FUNC_(f_fpu_fix_start, F_FPU_FIX_START)
#define f_fpu_fix_end   FC_FUNC_(f_fpu_fix_end,   F_FPU_FIX_END)

void f_fpu_fix_start(unsigned int *old_cw) {
  fpu_fix_start(old_cw);
}

void f_fpu_fix_end(unsigned int *old_cw) {
  fpu_fix_end(old_cw);
}

#endif

}
 
