! module file describing the C++-to-Fortran interface found in f_ts.cpp.
!
! Copyright (c) 2026, Nakata Maho
! All rights reserved.
!
! Redistribution and use in source and binary forms, with or without
! modification, are permitted provided that the following conditions are met:
!
! 1. Redistributions of source code must retain the above copyright notice,
!    this list of conditions and the following disclaimer.
! 2. Redistributions in binary form must reproduce the above copyright notice,
!    this list of conditions and the following disclaimer in the documentation
!    and/or other materials provided with the distribution.
!
! THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
! AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
! IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
! ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
! LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
! CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
! SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
! INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
! CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
! ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
! POSSIBILITY OF SUCH DAMAGE.

module tsext
  implicit none

  interface
    pure subroutine f_ts_add(a, b, c)
      real*4, intent(in) :: a(3), b(3)
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_add_ts_d(a, b, c)
      real*4, intent(in) :: a(3), b
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_sub(a, b, c)
      real*4, intent(in) :: a(3), b(3)
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_sub_ts_d(a, b, c)
      real*4, intent(in) :: a(3), b
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_sub_d_ts(a, b, c)
      real*4, intent(in) :: a, b(3)
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_mul(a, b, c)
      real*4, intent(in) :: a(3), b(3)
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_mul_ts_d(a, b, c)
      real*4, intent(in) :: a(3), b
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_div(a, b, c)
      real*4, intent(in) :: a(3), b(3)
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_div_ts_d(a, b, c)
      real*4, intent(in) :: a(3), b
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_div_d_ts(a, b, c)
      real*4, intent(in) :: a, b(3)
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_sqrt(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_sqr(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_abs(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_npwr(a, n, b)
      real*4, intent(in) :: a(3)
      integer, intent(in) :: n
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_nroot(a, n, b)
      real*4, intent(in) :: a(3)
      integer, intent(in) :: n
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_nint(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_aint(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_floor(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_ceil(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_log(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_log10(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_exp(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_sin(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_cos(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_tan(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_asin(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_acos(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_atan(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_atan2(a, b, c)
      real*4, intent(in) :: a(3), b(3)
      real*4, intent(out) :: c(3)
    end subroutine

    pure subroutine f_ts_sinh(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_cosh(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_tanh(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_asinh(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_acosh(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_atanh(a, b)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: b(3)
    end subroutine

    pure subroutine f_ts_sincos(a, s, c)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: s(3), c(3)
    end subroutine

    pure subroutine f_ts_sincosh(a, s, c)
      real*4, intent(in) :: a(3)
      real*4, intent(out) :: s(3), c(3)
    end subroutine

    subroutine f_ts_swrite(a, prec, str, maxlen)
      real*4, intent(in) :: a(3)
      integer, intent(in) :: prec, maxlen
      character, intent(out) :: str(maxlen)
    end subroutine

    subroutine f_ts_rand(a)
      real*4, intent(out) :: a(3)
    end subroutine

    pure subroutine f_ts_comp(a, b, r)
      real*4, intent(in) :: a(3), b(3)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_ts_comp_ts_d(a, b, r)
      real*4, intent(in) :: a(3), b
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_ts_comp_d_ts(a, b, r)
      real*4, intent(in) :: a, b(3)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_ts_pi(a)
      real*4, intent(out) :: a(3)
    end subroutine

    pure subroutine f_ts_nan(a)
      real*4, intent(out) :: a(3)
    end subroutine
  end interface
end
