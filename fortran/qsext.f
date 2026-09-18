! module file describing the C++-to-Fortran interface found in f_qs.cpp.
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

module qsext
  implicit none

  interface
    pure subroutine f_qs_add(a, b, c)
      real*4, intent(in) :: a(4), b(4)
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_add_qs_d(a, b, c)
      real*4, intent(in) :: a(4), b
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_sub(a, b, c)
      real*4, intent(in) :: a(4), b(4)
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_sub_qs_d(a, b, c)
      real*4, intent(in) :: a(4), b
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_sub_d_qs(a, b, c)
      real*4, intent(in) :: a, b(4)
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_mul(a, b, c)
      real*4, intent(in) :: a(4), b(4)
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_mul_qs_d(a, b, c)
      real*4, intent(in) :: a(4), b
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_div(a, b, c)
      real*4, intent(in) :: a(4), b(4)
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_div_qs_d(a, b, c)
      real*4, intent(in) :: a(4), b
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_div_d_qs(a, b, c)
      real*4, intent(in) :: a, b(4)
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_sqrt(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_sqr(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_abs(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_npwr(a, n, b)
      real*4, intent(in) :: a(4)
      integer, intent(in) :: n
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_nroot(a, n, b)
      real*4, intent(in) :: a(4)
      integer, intent(in) :: n
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_nint(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_aint(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_floor(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_ceil(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_log(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_log10(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_exp(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_sin(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_cos(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_tan(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_asin(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_acos(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_atan(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_atan2(a, b, c)
      real*4, intent(in) :: a(4), b(4)
      real*4, intent(out) :: c(4)
    end subroutine

    pure subroutine f_qs_sinh(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_cosh(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_tanh(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_asinh(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_acosh(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_atanh(a, b)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: b(4)
    end subroutine

    pure subroutine f_qs_sincos(a, s, c)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: s(4), c(4)
    end subroutine

    pure subroutine f_qs_sincosh(a, s, c)
      real*4, intent(in) :: a(4)
      real*4, intent(out) :: s(4), c(4)
    end subroutine

    subroutine f_qs_swrite(a, prec, str, maxlen)
      real*4, intent(in) :: a(4)
      integer, intent(in) :: prec, maxlen
      character, intent(out) :: str(maxlen)
    end subroutine

    subroutine f_qs_rand(a)
      real*4, intent(out) :: a(4)
    end subroutine

    pure subroutine f_qs_comp(a, b, r)
      real*4, intent(in) :: a(4), b(4)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_qs_comp_qs_d(a, b, r)
      real*4, intent(in) :: a(4), b
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_qs_comp_d_qs(a, b, r)
      real*4, intent(in) :: a, b(4)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_qs_pi(a)
      real*4, intent(out) :: a(4)
    end subroutine

    pure subroutine f_qs_nan(a)
      real*4, intent(out) :: a(4)
    end subroutine
  end interface
end
