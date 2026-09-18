! module file describing the C++-to-Fortran interface found in f_ds.cpp.
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

module dsext
  implicit none

  interface
    pure subroutine f_ds_add(a, b, c)
      real*4, intent(in) :: a(2), b(2)
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_add_ds_d(a, b, c)
      real*4, intent(in) :: a(2), b
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_sub(a, b, c)
      real*4, intent(in) :: a(2), b(2)
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_sub_ds_d(a, b, c)
      real*4, intent(in) :: a(2), b
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_sub_d_ds(a, b, c)
      real*4, intent(in) :: a, b(2)
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_mul(a, b, c)
      real*4, intent(in) :: a(2), b(2)
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_mul_ds_d(a, b, c)
      real*4, intent(in) :: a(2), b
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_div(a, b, c)
      real*4, intent(in) :: a(2), b(2)
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_div_ds_d(a, b, c)
      real*4, intent(in) :: a(2), b
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_div_d_ds(a, b, c)
      real*4, intent(in) :: a, b(2)
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_sqrt(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_sqr(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_abs(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_npwr(a, n, b)
      real*4, intent(in) :: a(2)
      integer, intent(in) :: n
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_nroot(a, n, b)
      real*4, intent(in) :: a(2)
      integer, intent(in) :: n
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_nint(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_aint(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_floor(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_ceil(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_log(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_log10(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_exp(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_sin(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_cos(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_tan(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_asin(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_acos(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_atan(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_atan2(a, b, c)
      real*4, intent(in) :: a(2), b(2)
      real*4, intent(out) :: c(2)
    end subroutine

    pure subroutine f_ds_sinh(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_cosh(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_tanh(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_asinh(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_acosh(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_atanh(a, b)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: b(2)
    end subroutine

    pure subroutine f_ds_sincos(a, s, c)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: s(2), c(2)
    end subroutine

    pure subroutine f_ds_sincosh(a, s, c)
      real*4, intent(in) :: a(2)
      real*4, intent(out) :: s(2), c(2)
    end subroutine

    subroutine f_ds_swrite(a, prec, str, maxlen)
      real*4, intent(in) :: a(2)
      integer, intent(in) :: prec, maxlen
      character, intent(out) :: str(maxlen)
    end subroutine

    subroutine f_ds_rand(a)
      real*4, intent(out) :: a(2)
    end subroutine

    pure subroutine f_ds_comp(a, b, r)
      real*4, intent(in) :: a(2), b(2)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_ds_comp_ds_d(a, b, r)
      real*4, intent(in) :: a(2), b
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_ds_comp_d_ds(a, b, r)
      real*4, intent(in) :: a, b(2)
      integer, intent(out) :: r
    end subroutine

    pure subroutine f_ds_pi(a)
      real*4, intent(out) :: a(2)
    end subroutine

    pure subroutine f_ds_nan(a)
      real*4, intent(out) :: a(2)
    end subroutine
  end interface
end
