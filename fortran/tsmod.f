!  tsmod.f
!
!  Fortran-90 module file to use with triple-single numbers.
!
!  Copyright (c) 2026, Nakata Maho
!  All rights reserved.
!
!  Redistribution and use in source and binary forms, with or without
!  modification, are permitted provided that the following conditions are met:
!
!  1. Redistributions of source code must retain the above copyright notice,
!     this list of conditions and the following disclaimer.
!  2. Redistributions in binary form must reproduce the above copyright notice,
!     this list of conditions and the following disclaimer in the documentation
!     and/or other materials provided with the distribution.
!
!  THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS "AS IS"
!  AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT LIMITED TO, THE
!  IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR A PARTICULAR PURPOSE
!  ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT HOLDER OR CONTRIBUTORS BE
!  LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL, EXEMPLARY, OR
!  CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO, PROCUREMENT OF
!  SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR PROFITS; OR BUSINESS
!  INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF LIABILITY, WHETHER IN
!  CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING NEGLIGENCE OR OTHERWISE)
!  ARISING IN ANY WAY OUT OF THE USE OF THIS SOFTWARE, EVEN IF ADVISED OF THE
!  POSSIBILITY OF SUCH DAMAGE.

module tsmodule
  use ddmodule, only: dd_real, dd_complex
  use qdmodule, only: qd_real, qd_complex
  use tsext
  implicit none

  type ts_real
    sequence
    real*4 :: re(3)
  end type ts_real

  type ts_complex
    sequence
    real*4 :: cmp(6)
  end type ts_complex

  real*4 d_ts_eps
  parameter (d_ts_eps = 8.4703295e-22)

  type (ts_real) ts_one, ts_zero, ts_eps, ts_huge, ts_tiny
  parameter (ts_one = ts_real((/1.0e0, 0.0e0, 0.0e0/)))
  parameter (ts_zero = ts_real((/0.0e0, 0.0e0, 0.0e0/)))
  parameter (ts_eps = ts_real((/d_ts_eps, 0.0e0, 0.0e0/)))
  parameter (ts_huge = ts_real((/3.4028235e+38, 0.0e0, 0.0e0/)))
  parameter (ts_tiny = ts_real((/3.3087225e-24, 0.0e0, 0.0e0/)))

  interface assignment (=)
    module procedure assign_ts_str
    module procedure assign_ts
    module procedure assign_ts_d
    module procedure assign_d_ts
    module procedure assign_ts_i
    module procedure assign_i_ts
    module procedure assign_ts_dd
    module procedure assign_dd_ts
    module procedure assign_ts_qd
    module procedure assign_qd_ts
    module procedure assign_tsc
    module procedure assign_tsc_ts
    module procedure assign_ts_tsc
    module procedure assign_tsc_d
    module procedure assign_tsc_i
    module procedure assign_d_tsc
    module procedure assign_tsc_dc
    module procedure assign_dc_tsc
    module procedure assign_tsc_ddc
    module procedure assign_ddc_tsc
    module procedure assign_tsc_qdc
    module procedure assign_qdc_tsc
  end interface

  interface operator (+)
    module procedure add_ts
    module procedure add_ts_d
    module procedure add_d_ts
    module procedure add_ts_i
    module procedure add_i_ts
    module procedure add_tsc
    module procedure add_tsc_ts
    module procedure add_ts_tsc
    module procedure add_tsc_d
    module procedure add_d_tsc
  end interface

  interface operator (-)
    module procedure sub_ts
    module procedure sub_ts_d
    module procedure sub_d_ts
    module procedure neg_ts
    module procedure sub_tsc
    module procedure sub_tsc_ts
    module procedure sub_ts_tsc
    module procedure sub_tsc_d
    module procedure sub_d_tsc
    module procedure neg_tsc
  end interface

  interface operator (*)
    module procedure mul_ts
    module procedure mul_ts_d
    module procedure mul_d_ts
    module procedure mul_ts_i
    module procedure mul_i_ts
    module procedure mul_tsc
    module procedure mul_tsc_ts
    module procedure mul_ts_tsc
    module procedure mul_tsc_d
    module procedure mul_d_tsc
    module procedure mul_tsc_i
    module procedure mul_i_tsc
  end interface

  interface operator (/)
    module procedure div_ts
    module procedure div_ts_d
    module procedure div_d_ts
    module procedure div_ts_i
    module procedure div_i_ts
    module procedure div_tsc
    module procedure div_tsc_ts
    module procedure div_ts_tsc
    module procedure div_tsc_d
  end interface

  interface operator (**)
    module procedure pwr_ts
    module procedure pwr_ts_i
    module procedure pwr_d_ts
    module procedure pwr_tsc_i
  end interface

  interface tsreal
    module procedure to_ts_i
    module procedure to_ts_d
    module procedure to_ts_dd
    module procedure to_ts_qd
    module procedure to_ts_ts
    module procedure to_ts_str
  end interface

  interface ddreal
    module procedure to_dd_ts
  end interface

  interface ddcomplex
    module procedure to_ddc_tsc
  end interface

  interface qdreal
    module procedure to_qd_ts
  end interface

  interface qdcomplex
    module procedure to_qdc_tsc
  end interface

  interface real
    module procedure to_d_ts
    module procedure to_ts_tsc
  end interface

  interface dble
    module procedure to_d_ts
    module procedure to_d_tsc
  end interface

  interface tscomplex
    module procedure to_tsc_ts
    module procedure to_tsc_ts2
    module procedure to_tsc_d
    module procedure to_tsc_dc
    module procedure to_tsc_ddc
    module procedure to_tsc_qdc
  end interface

  interface cmplx
    module procedure to_dc_tsc
  end interface

  interface conjg
    module procedure tscconjg
  end interface

  interface int
    module procedure to_int_ts
  end interface

  interface sin
    module procedure tssin
  end interface
  interface cos
    module procedure tscos
  end interface
  interface tan
    module procedure tstan
  end interface
  interface sincos
    module procedure tssincos
  end interface

  interface asin
    module procedure tsasin
  end interface
  interface acos
    module procedure tsacos
  end interface
  interface atan
    module procedure tsatan
  end interface
  interface atan2
    module procedure tsatan2
  end interface

  interface exp
    module procedure tsexp
    module procedure tscexp
  end interface
  interface log
    module procedure tslog
    module procedure tsclog
  end interface
  interface log10
    module procedure tslog10
  end interface

  interface sqrt
    module procedure tssqrt
  end interface
  interface sqr
    module procedure tssqr
  end interface
  interface nroot
    module procedure tsnroot
  end interface

  interface sinh
    module procedure tssinh
  end interface
  interface cosh
    module procedure tscosh
  end interface
  interface tanh
    module procedure tstanh
  end interface
  interface sincosh
    module procedure tssincosh
  end interface

  interface asinh
    module procedure tsasinh
  end interface
  interface acosh
    module procedure tsacosh
  end interface
  interface atanh
    module procedure tsatanh
  end interface

  interface aint
    module procedure tsaint
  end interface

  interface floor
    module procedure tsfloor
  end interface

  interface ceil
    module procedure tsceil
  end interface

  interface nint
    module procedure tsnint
  end interface

  interface anint
    module procedure tsanint
  end interface

  interface abs
    module procedure tsabs
    module procedure tscabs
  end interface

  interface sign
    module procedure tssign
    module procedure tssign_ts_d
  end interface

  interface random_number
    module procedure tsrand
  end interface

  interface aimag
    module procedure ts_aimag
  end interface

  interface operator (==)
    module procedure eq_ts
    module procedure eq_ts_d
    module procedure eq_d_ts
    module procedure eq_ts_i
    module procedure eq_i_ts
    module procedure eq_tsc
    module procedure eq_tsc_ts
    module procedure eq_ts_tsc
  end interface

  interface operator (/=)
    module procedure ne_ts
    module procedure ne_ts_d
    module procedure ne_d_ts
    module procedure ne_ts_i
    module procedure ne_i_ts
    module procedure ne_tsc
    module procedure ne_tsc_ts
    module procedure ne_ts_tsc
  end interface

  interface operator (>)
    module procedure gt_ts
    module procedure gt_ts_d
    module procedure gt_d_ts
    module procedure gt_ts_i
    module procedure gt_i_ts
  end interface

  interface operator (<)
    module procedure lt_ts
    module procedure lt_ts_d
    module procedure lt_d_ts
    module procedure lt_ts_i
    module procedure lt_i_ts
  end interface

  interface operator (>=)
    module procedure ge_ts
    module procedure ge_ts_d
    module procedure ge_d_ts
    module procedure ge_ts_i
    module procedure ge_i_ts
  end interface

  interface operator (<=)
    module procedure le_ts
    module procedure le_ts_d
    module procedure le_d_ts
    module procedure le_ts_i
    module procedure le_i_ts
  end interface

  interface read_scalar
    module procedure tsinpq
    module procedure tscinpq
  end interface

  interface write_scalar
    module procedure tsoutq
    module procedure tscoutq
  end interface

  interface tsread
    module procedure tsinpq
  end interface

  interface tswrite
    module procedure tsoutq
  end interface

  interface tscread
    module procedure tscinpq
  end interface

  interface tscwrite
    module procedure tscoutq
  end interface

  interface tspi
    module procedure ts_pi
  end interface

  interface huge
    module procedure tshuge
  end interface

  interface safe_huge
    module procedure ts_safe_huge
  end interface

  interface tiny
    module procedure tstiny
  end interface

  interface epsilon
    module procedure tsepsilon
  end interface

  interface radix
    module procedure ts_radix
  end interface

  interface digits
    module procedure ts_digits
  end interface

  interface maxexponent
    module procedure ts_max_expn
  end interface

  interface minexponent
    module procedure ts_min_expn
  end interface

  interface precision
    module procedure ts_precision
  end interface

  interface range
    module procedure ts_range
  end interface

  interface nan
    module procedure ts_nan
  end interface

  interface min
    module procedure tsmin
    module procedure tsmin2
  end interface

  interface max
    module procedure tsmax
    module procedure tsmax2
  end interface

  interface mod
    module procedure tsmod
  end interface

contains

  subroutine assign_ts_str(a, s)
    type (ts_real), intent(inout) :: a
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call tsinpc(t, a%re)
  end subroutine assign_ts_str

  elemental subroutine assign_ts(a, b)
    type (ts_real), intent(inout) :: a
    type (ts_real), intent(in) :: b
    a%re = b%re
  end subroutine assign_ts

  elemental subroutine assign_ts_d(a, d)
    type (ts_real), intent(inout) :: a
    real*4, intent(in) :: d
    a%re(1) = d
    a%re(2:3) = 0.0e0
  end subroutine assign_ts_d

  elemental subroutine assign_d_ts(d, a)
    real*4, intent(inout) :: d
    type (ts_real), intent(in) :: a
    d = a%re(1)
  end subroutine assign_d_ts

  elemental subroutine assign_ts_i(a, i)
    type (ts_real), intent(inout) :: a
    integer, intent(in) :: i
    a%re(1) = i
    a%re(2:3) = 0.0e0
  end subroutine assign_ts_i

  elemental subroutine assign_i_ts(i, a)
    integer, intent(inout) :: i
    type (ts_real), intent(in) :: a
    i = a%re(1)
  end subroutine assign_i_ts

  elemental subroutine assign_ts_dd(ts, dd)
    type (ts_real), intent(inout) :: ts
    type (dd_real), intent(in) :: dd
    ts%re(1:2) = dd%re
    ts%re(3) = 0.e0
  end subroutine assign_ts_dd

  elemental subroutine assign_dd_ts(dd, ts)
    type (dd_real), intent(inout) :: dd
    type (ts_real), intent(in) :: ts
    dd%re = ts%re(1:2)
  end subroutine assign_dd_ts

  elemental subroutine assign_ts_qd(ts, qd)
    type (ts_real), intent(inout) :: ts
    type (qd_real), intent(in) :: qd
    ts%re = qd%re(1:3)
  end subroutine assign_ts_qd

  elemental subroutine assign_qd_ts(qd, ts)
    type (qd_real), intent(inout) :: qd
    type (ts_real), intent(in) :: ts
    qd%re(1:3) = ts%re
    qd%re(4) = 0.e0
  end subroutine assign_qd_ts

  elemental subroutine assign_tsc(a, b)
    type (ts_complex), intent(inout) :: a
    type (ts_complex), intent(in) :: b
    a%cmp = b%cmp
  end subroutine assign_tsc

  elemental subroutine assign_tsc_ts(tsc, ts)
    type (ts_complex), intent(inout) :: tsc
    type (ts_real), intent(in) :: ts
    tsc%cmp(1:3) = ts%re
    tsc%cmp(4:6) = 0.e0
  end subroutine assign_tsc_ts

  elemental subroutine assign_ts_tsc(ts, tsc)
    type (ts_real), intent(inout) :: ts
    type (ts_complex), intent(in) :: tsc
    ts%re = tsc%cmp(1:3)
  end subroutine assign_ts_tsc

  elemental subroutine assign_tsc_d(tsc, d)
    type (ts_complex), intent(inout) :: tsc
    real*4, intent(in) :: d
    tsc%cmp(1) = d
    tsc%cmp(2:6) = 0.e0
  end subroutine assign_tsc_d

  elemental subroutine assign_tsc_i(tsc, i)
    type (ts_complex), intent(inout) :: tsc
    integer, intent(in) :: i
    tsc%cmp(1) = i
    tsc%cmp(2:6) = 0.e0
  end subroutine assign_tsc_i

  elemental subroutine assign_d_tsc(d, tsc)
    real*4, intent(inout) :: d
    type (ts_complex), intent(in) :: tsc
    d = tsc%cmp(1)
  end subroutine assign_d_tsc

  elemental subroutine assign_tsc_dc(tsc, dc)
    type (ts_complex), intent(inout) :: tsc
    complex(kind(0.e0)), intent(in) :: dc
    tsc%cmp(1) = dble(dc)
    tsc%cmp(2:3) = 0.e0
    tsc%cmp(4) = aimag(dc)
    tsc%cmp(5:6) = 0.e0
  end subroutine assign_tsc_dc

  elemental subroutine assign_dc_tsc(dc, tsc)
    complex(kind(0.e0)), intent(inout) :: dc
    type (ts_complex), intent(in) :: tsc
    dc = cmplx(tsc%cmp(1), tsc%cmp(4), kind(0.e0))
  end subroutine assign_dc_tsc

  elemental subroutine assign_tsc_ddc(tsc, ddc)
    type (ts_complex), intent(inout) :: tsc
    type (dd_complex), intent(in) :: ddc
    tsc%cmp(1:2) = ddc%cmp(1:2)
    tsc%cmp(3) = 0.e0
    tsc%cmp(4:5) = ddc%cmp(3:4)
    tsc%cmp(6) = 0.e0
  end subroutine assign_tsc_ddc

  elemental subroutine assign_ddc_tsc(ddc, tsc)
    type (dd_complex), intent(inout) :: ddc
    type (ts_complex), intent(in) :: tsc
    ddc%cmp(1:2) = tsc%cmp(1:2)
    ddc%cmp(3:4) = tsc%cmp(4:5)
  end subroutine assign_ddc_tsc

  elemental subroutine assign_tsc_qdc(tsc, qdc)
    type (ts_complex), intent(inout) :: tsc
    type (qd_complex), intent(in) :: qdc
    tsc%cmp(1:3) = qdc%cmp(1:3)
    tsc%cmp(4:6) = qdc%cmp(5:7)
  end subroutine assign_tsc_qdc

  elemental subroutine assign_qdc_tsc(qdc, tsc)
    type (qd_complex), intent(inout) :: qdc
    type (ts_complex), intent(in) :: tsc
    qdc%cmp(1:3) = tsc%cmp(1:3)
    qdc%cmp(4) = 0.e0
    qdc%cmp(5:7) = tsc%cmp(4:6)
    qdc%cmp(8) = 0.e0
  end subroutine assign_qdc_tsc

  elemental type (ts_real) function to_ts_i(ia)
    integer, intent(in) :: ia
    to_ts_i%re(1) = ia
    to_ts_i%re(2:3) = 0.e0
  end function to_ts_i

  elemental type (ts_real) function to_ts_d(a)
    real*4, intent(in) :: a
    to_ts_d%re(1) = a
    to_ts_d%re(2:3) = 0.e0
  end function to_ts_d

  elemental type (ts_real) function to_ts_dd(dd)
    type (dd_real), intent(in) :: dd
    to_ts_dd%re(1:2) = dd%re
    to_ts_dd%re(3) = 0.e0
  end function to_ts_dd

  elemental type (ts_real) function to_ts_qd(qd)
    type (qd_real), intent(in) :: qd
    to_ts_qd%re = qd%re(1:3)
  end function to_ts_qd

  elemental type (ts_real) function to_ts_ts(a)
    type (ts_real), intent(in) :: a
    to_ts_ts%re = a%re
  end function to_ts_ts

  type (ts_real) function to_ts_str(s)
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call tsinpc(t, to_ts_str%re)
  end function to_ts_str

  elemental type (dd_real) function to_dd_ts(ts)
    type (ts_real), intent(in) :: ts
    to_dd_ts%re = ts%re(1:2)
  end function to_dd_ts

  elemental type (qd_real) function to_qd_ts(ts)
    type (ts_real), intent(in) :: ts
    to_qd_ts%re(1:3) = ts%re
    to_qd_ts%re(4) = 0.e0
  end function to_qd_ts

  elemental type (dd_complex) function to_ddc_tsc(tsc)
    type (ts_complex), intent(in) :: tsc
    to_ddc_tsc%cmp(1:2) = tsc%cmp(1:2)
    to_ddc_tsc%cmp(3:4) = tsc%cmp(4:5)
  end function to_ddc_tsc

  elemental type (qd_complex) function to_qdc_tsc(tsc)
    type (ts_complex), intent(in) :: tsc
    to_qdc_tsc%cmp(1:3) = tsc%cmp(1:3)
    to_qdc_tsc%cmp(4) = 0.e0
    to_qdc_tsc%cmp(5:7) = tsc%cmp(4:6)
    to_qdc_tsc%cmp(8) = 0.e0
  end function to_qdc_tsc

  elemental real*4 function to_d_ts(a)
    type (ts_real), intent(in) :: a
    to_d_ts = a%re(1)
  end function to_d_ts

  elemental type (ts_real) function to_ts_tsc(tsc)
    type (ts_complex), intent(in) :: tsc
    to_ts_tsc%re = tsc%cmp(1:3)
  end function to_ts_tsc

  elemental type (ts_complex) function to_tsc_ts(ts)
    type (ts_real), intent(in) :: ts
    to_tsc_ts%cmp(1:3) = ts%re
    to_tsc_ts%cmp(4:6) = 0.e0
  end function to_tsc_ts

  elemental type (ts_complex) function to_tsc_ts2(x, y)
    type (ts_real), intent(in) :: x, y
    to_tsc_ts2%cmp(1:3) = x%re
    to_tsc_ts2%cmp(4:6) = y%re
  end function to_tsc_ts2

  elemental type (ts_complex) function to_tsc_d(d)
    real*4, intent(in) :: d
    to_tsc_d%cmp(1) = d
    to_tsc_d%cmp(2:6) = 0.e0
  end function to_tsc_d

  elemental complex(kind(0.e0)) function to_dc_tsc(tsc)
    type (ts_complex), intent(in) :: tsc
    to_dc_tsc = cmplx(tsc%cmp(1), tsc%cmp(4), kind(0.e0))
  end function to_dc_tsc

  elemental type (ts_complex) function to_tsc_dc(dc)
    complex(kind(0.e0)), intent(in) :: dc
    to_tsc_dc%cmp(1) = dble(dc)
    to_tsc_dc%cmp(2:3) = 0.e0
    to_tsc_dc%cmp(4) = aimag(dc)
    to_tsc_dc%cmp(5:6) = 0.e0
  end function to_tsc_dc

  elemental type (ts_complex) function to_tsc_ddc(ddc)
    type (dd_complex), intent(in) :: ddc
    to_tsc_ddc%cmp(1:2) = ddc%cmp(1:2)
    to_tsc_ddc%cmp(3) = 0.e0
    to_tsc_ddc%cmp(4:5) = ddc%cmp(3:4)
    to_tsc_ddc%cmp(6) = 0.e0
  end function to_tsc_ddc

  elemental type (ts_complex) function to_tsc_qdc(qdc)
    type (qd_complex), intent(in) :: qdc
    to_tsc_qdc%cmp(1:3) = qdc%cmp(1:3)
    to_tsc_qdc%cmp(4:6) = qdc%cmp(5:7)
  end function to_tsc_qdc

  elemental real*4 function to_d_tsc(tsc)
    type (ts_complex), intent(in) :: tsc
    to_d_tsc = tsc%cmp(1)
  end function to_d_tsc

  elemental integer function to_int_ts(a)
    type (ts_real), intent(in) :: a
    to_int_ts = a%re(1)
  end function to_int_ts

  elemental type (ts_complex) function tscconjg(tsc)
    type (ts_complex), intent(in) :: tsc
    tscconjg%cmp(1:3) = tsc%cmp(1:3)
    tscconjg%cmp(4:6) = -tsc%cmp(4:6)
  end function tscconjg

  elemental type (ts_real) function add_ts(a, b)
    type (ts_real), intent(in) :: a, b
    call f_ts_add(a%re, b%re, add_ts%re)
  end function add_ts

  elemental type (ts_real) function add_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ts_add_ts_d(a%re, b, add_ts_d%re)
  end function add_ts_d

  elemental type (ts_real) function add_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    add_d_ts = add_ts_d(b, a)
  end function add_d_ts

  elemental type (ts_real) function add_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    call f_ts_add_ts_d(a%re, ts_int_to_float(b), add_ts_i%re)
  end function add_ts_i

  elemental type (ts_real) function add_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    add_i_ts = add_ts_i(b, a)
  end function add_i_ts

  elemental type (ts_complex) function add_tsc(a, b)
    type (ts_complex), intent(in) :: a, b
    call f_ts_add(a%cmp(1:3), b%cmp(1:3), add_tsc%cmp(1:3))
    call f_ts_add(a%cmp(4:6), b%cmp(4:6), add_tsc%cmp(4:6))
  end function add_tsc

  elemental type (ts_complex) function add_tsc_ts(a, b)
    type (ts_complex), intent(in) :: a
    type (ts_real), intent(in) :: b
    call f_ts_add(a%cmp(1:3), b%re, add_tsc_ts%cmp(1:3))
    add_tsc_ts%cmp(4:6) = a%cmp(4:6)
  end function add_tsc_ts

  elemental type (ts_complex) function add_ts_tsc(a, b)
    type (ts_real), intent(in) :: a
    type (ts_complex), intent(in) :: b
    add_ts_tsc = add_tsc_ts(b, a)
  end function add_ts_tsc

  elemental type (ts_complex) function add_tsc_d(a, b)
    type (ts_complex), intent(in) :: a
    real*4, intent(in) :: b
    type (ts_real) :: tsb
    tsb%re(1) = b
    tsb%re(2:3) = 0.e0
    call f_ts_add(a%cmp(1:3), tsb%re, add_tsc_d%cmp(1:3))
    add_tsc_d%cmp(4:6) = a%cmp(4:6)
  end function add_tsc_d

  elemental type (ts_complex) function add_d_tsc(a, b)
    real*4, intent(in) :: a
    type (ts_complex), intent(in) :: b
    add_d_tsc = add_tsc_d(b, a)
  end function add_d_tsc

  elemental type (ts_real) function sub_ts(a, b)
    type (ts_real), intent(in) :: a, b
    call f_ts_sub(a%re, b%re, sub_ts%re)
  end function sub_ts

  elemental type (ts_real) function sub_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ts_sub_ts_d(a%re, b, sub_ts_d%re)
  end function sub_ts_d

  elemental type (ts_real) function sub_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    call f_ts_sub_d_ts(a, b%re, sub_d_ts%re)
  end function sub_d_ts

  elemental type (ts_real) function neg_ts(a)
    type (ts_real), intent(in) :: a
    neg_ts%re = -a%re
  end function neg_ts

  elemental type (ts_complex) function sub_tsc(a, b)
    type (ts_complex), intent(in) :: a, b
    call f_ts_sub(a%cmp(1:3), b%cmp(1:3), sub_tsc%cmp(1:3))
    call f_ts_sub(a%cmp(4:6), b%cmp(4:6), sub_tsc%cmp(4:6))
  end function sub_tsc

  elemental type (ts_complex) function sub_tsc_ts(a, b)
    type (ts_complex), intent(in) :: a
    type (ts_real), intent(in) :: b
    call f_ts_sub(a%cmp(1:3), b%re, sub_tsc_ts%cmp(1:3))
    sub_tsc_ts%cmp(4:6) = a%cmp(4:6)
  end function sub_tsc_ts

  elemental type (ts_complex) function sub_ts_tsc(a, b)
    type (ts_real), intent(in) :: a
    type (ts_complex), intent(in) :: b
    call f_ts_sub(a%re, b%cmp(1:3), sub_ts_tsc%cmp(1:3))
    sub_ts_tsc%cmp(4:6) = -b%cmp(4:6)
  end function sub_ts_tsc

  elemental type (ts_complex) function sub_tsc_d(a, b)
    type (ts_complex), intent(in) :: a
    real*4, intent(in) :: b
    type (ts_real) :: tsb
    tsb%re(1) = b
    tsb%re(2:3) = 0.e0
    call f_ts_sub(a%cmp(1:3), tsb%re, sub_tsc_d%cmp(1:3))
    sub_tsc_d%cmp(4:6) = a%cmp(4:6)
  end function sub_tsc_d

  elemental type (ts_complex) function sub_d_tsc(a, b)
    real*4, intent(in) :: a
    type (ts_complex), intent(in) :: b
    type (ts_real) :: tsa
    tsa%re(1) = a
    tsa%re(2:3) = 0.e0
    call f_ts_sub(tsa%re, b%cmp(1:3), sub_d_tsc%cmp(1:3))
    sub_d_tsc%cmp(4:6) = -b%cmp(4:6)
  end function sub_d_tsc

  elemental type (ts_complex) function neg_tsc(a)
    type (ts_complex), intent(in) :: a
    neg_tsc%cmp = -a%cmp
  end function neg_tsc

  elemental type (ts_real) function mul_ts(a, b)
    type (ts_real), intent(in) :: a, b
    call f_ts_mul(a%re, b%re, mul_ts%re)
  end function mul_ts

  elemental type (ts_real) function mul_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ts_mul_ts_d(a%re, b, mul_ts_d%re)
  end function mul_ts_d

  elemental type (ts_real) function mul_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    mul_d_ts = mul_ts_d(b, a)
  end function mul_d_ts

  elemental type (ts_real) function mul_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    call f_ts_mul_ts_d(a%re, ts_int_to_float(b), mul_ts_i%re)
  end function mul_ts_i

  elemental type (ts_real) function mul_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    mul_i_ts = mul_ts_i(b, a)
  end function mul_i_ts

  elemental type (ts_complex) function mul_tsc(a, b)
    type (ts_complex), intent(in) :: a, b
    type (ts_real) :: t1, t2
    call f_ts_mul(a%cmp(1:3), b%cmp(1:3), t1%re)
    call f_ts_mul(a%cmp(4:6), b%cmp(4:6), t2%re)
    call f_ts_sub(t1%re, t2%re, mul_tsc%cmp(1:3))
    call f_ts_mul(a%cmp(1:3), b%cmp(4:6), t1%re)
    call f_ts_mul(a%cmp(4:6), b%cmp(1:3), t2%re)
    call f_ts_add(t1%re, t2%re, mul_tsc%cmp(4:6))
  end function mul_tsc

  elemental type (ts_complex) function mul_tsc_ts(a, b)
    type (ts_complex), intent(in) :: a
    type (ts_real), intent(in) :: b
    call f_ts_mul(a%cmp(1:3), b%re, mul_tsc_ts%cmp(1:3))
    call f_ts_mul(a%cmp(4:6), b%re, mul_tsc_ts%cmp(4:6))
  end function mul_tsc_ts

  elemental type (ts_complex) function mul_ts_tsc(a, b)
    type (ts_real), intent(in) :: a
    type (ts_complex), intent(in) :: b
    call f_ts_mul(a%re, b%cmp(1:3), mul_ts_tsc%cmp(1:3))
    call f_ts_mul(a%re, b%cmp(4:6), mul_ts_tsc%cmp(4:6))
  end function mul_ts_tsc

  elemental type (ts_complex) function mul_tsc_d(a, b)
    type (ts_complex), intent(in) :: a
    real*4, intent(in) :: b
    call f_ts_mul_ts_d(a%cmp(1:3), b, mul_tsc_d%cmp(1:3))
    call f_ts_mul_ts_d(a%cmp(4:6), b, mul_tsc_d%cmp(4:6))
  end function mul_tsc_d

  elemental type (ts_complex) function mul_d_tsc(a, b)
    real*4, intent(in) :: a
    type (ts_complex), intent(in) :: b
    call f_ts_mul_ts_d(b%cmp(1:3), a, mul_d_tsc%cmp(1:3))
    call f_ts_mul_ts_d(b%cmp(4:6), a, mul_d_tsc%cmp(4:6))
  end function mul_d_tsc

  elemental type (ts_complex) function mul_tsc_i(a, b)
    type (ts_complex), intent(in) :: a
    integer, intent(in) :: b
    call f_ts_mul_ts_d(a%cmp(1:3), ts_int_to_float(b), mul_tsc_i%cmp(1:3))
    call f_ts_mul_ts_d(a%cmp(4:6), ts_int_to_float(b), mul_tsc_i%cmp(4:6))
  end function mul_tsc_i

  elemental type (ts_complex) function mul_i_tsc(a, b)
    integer, intent(in) :: a
    type (ts_complex), intent(in) :: b
    call f_ts_mul_ts_d(b%cmp(1:3), ts_int_to_float(a), mul_i_tsc%cmp(1:3))
    call f_ts_mul_ts_d(b%cmp(4:6), ts_int_to_float(a), mul_i_tsc%cmp(4:6))
  end function mul_i_tsc

  elemental type (ts_real) function div_ts(a, b)
    type (ts_real), intent(in) :: a, b
    call f_ts_div(a%re, b%re, div_ts%re)
  end function div_ts

  elemental type (ts_real) function div_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ts_div_ts_d(a%re, b, div_ts_d%re)
  end function div_ts_d

  elemental type (ts_real) function div_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    call f_ts_div_d_ts(a, b%re, div_d_ts%re)
  end function div_d_ts

  elemental type (ts_real) function div_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    call f_ts_div_ts_d(a%re, ts_int_to_float(b), div_ts_i%re)
  end function div_ts_i

  elemental type (ts_real) function div_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    call f_ts_div_d_ts(ts_int_to_float(a), b%re, div_i_ts%re)
  end function div_i_ts

  elemental type (ts_complex) function div_tsc(a, b)
    type (ts_complex), intent(in) :: a, b
    type (ts_real) :: t1, t2, t3, t4, t5
    call f_ts_mul(a%cmp(1:3), b%cmp(1:3), t1%re)
    call f_ts_mul(a%cmp(4:6), b%cmp(4:6), t2%re)
    call f_ts_add(t1%re, t2%re, t3%re)
    call f_ts_mul(a%cmp(1:3), b%cmp(4:6), t1%re)
    call f_ts_mul(a%cmp(4:6), b%cmp(1:3), t2%re)
    call f_ts_sub(t2%re, t1%re, t4%re)
    call f_ts_mul(b%cmp(1:3), b%cmp(1:3), t1%re)
    call f_ts_mul(b%cmp(4:6), b%cmp(4:6), t2%re)
    call f_ts_add(t1%re, t2%re, t5%re)
    call f_ts_div(t3%re, t5%re, div_tsc%cmp(1:3))
    call f_ts_div(t4%re, t5%re, div_tsc%cmp(4:6))
  end function div_tsc

  elemental type (ts_complex) function div_tsc_ts(a, b)
    type (ts_complex), intent(in) :: a
    type (ts_real), intent(in) :: b
    call f_ts_div(a%cmp(1:3), b%re, div_tsc_ts%cmp(1:3))
    call f_ts_div(a%cmp(4:6), b%re, div_tsc_ts%cmp(4:6))
  end function div_tsc_ts

  elemental type (ts_complex) function div_ts_tsc(a, b)
    type (ts_real), intent(in) :: a
    type (ts_complex), intent(in) :: b
    type (ts_real) :: t1, t2, t3, t4, t5
    call f_ts_mul(a%re, b%cmp(1:3), t1%re)
    call f_ts_mul(a%re, b%cmp(4:6), t2%re)
    t2%re = -t2%re
    call f_ts_mul(b%cmp(1:3), b%cmp(1:3), t3%re)
    call f_ts_mul(b%cmp(4:6), b%cmp(4:6), t4%re)
    call f_ts_add(t3%re, t4%re, t5%re)
    call f_ts_div(t1%re, t5%re, div_ts_tsc%cmp(1:3))
    call f_ts_div(t2%re, t5%re, div_ts_tsc%cmp(4:6))
  end function div_ts_tsc

  elemental type (ts_complex) function div_tsc_d(a, b)
    type (ts_complex), intent(in) :: a
    real*4, intent(in) :: b
    call f_ts_div_ts_d(a%cmp(1:3), b, div_tsc_d%cmp(1:3))
    call f_ts_div_ts_d(a%cmp(4:6), b, div_tsc_d%cmp(4:6))
  end function div_tsc_d

  elemental type (ts_real) function pwr_ts(a, b)
    type (ts_real), intent(in) :: a, b
    type (ts_real) t1, t2
    call f_ts_log(a%re, t1%re)
    call f_ts_mul(t1%re, b%re, t2%re)
    call f_ts_exp(t2%re, pwr_ts%re)
  end function pwr_ts

  elemental type (ts_real) function pwr_ts_i(a, n)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: n
    call f_ts_npwr(a%re, n, pwr_ts_i%re)
  end function pwr_ts_i

  elemental type (ts_real) function pwr_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    type (ts_real) t1, t2, t3
    t1%re(1) = a
    t1%re(2:3) = 0.e0
    call f_ts_log(t1%re, t2%re)
    call f_ts_mul(t2%re, b%re, t3%re)
    call f_ts_exp(t3%re, pwr_d_ts%re)
  end function pwr_d_ts

  elemental type (ts_complex) function pwr_tsc_i(a, n)
    type (ts_complex), intent(in) :: a
    integer, intent(in) :: n
    integer i2, n1
    type (ts_real) t1, t2, t3
    type (ts_complex) c1, c2

    intrinsic :: iabs, ishft

    if (n == 0) then
      if (all(a%cmp == 0.e0)) then
        call f_ts_nan(pwr_tsc_i%cmp(1:3))
        call f_ts_nan(pwr_tsc_i%cmp(4:6))
        return
      endif
      pwr_tsc_i%cmp(1) = 1.e0
      pwr_tsc_i%cmp(2:6) = 0.e0
      return
    endif

    n1 = iabs(n)
    i2 = ishft(1, n1 - 1)
    c1%cmp(1) = 1.e0
    c1%cmp(2:6) = 0.e0

110 continue
    if (n1 >= i2) then
      c2 = a * c1
      c1 = c2
      n1 = n1 - i2
    endif
    i2 = i2 / 2
    if (i2 >= 1) then
      c2 = c1 * c1
      c1 = c2
      goto 110
    endif

    if (n > 0) then
      pwr_tsc_i = c1
    else
      c1%cmp(4:6) = -c1%cmp(4:6)
      call f_ts_mul(c1%cmp(1:3), c1%cmp(1:3), t1%re)
      call f_ts_mul(c1%cmp(4:6), c1%cmp(4:6), t2%re)
      call f_ts_add(t1%re, t2%re, t3%re)
      call f_ts_div(c1%cmp(1:3), t3%re, pwr_tsc_i%cmp(1:3))
      call f_ts_div(c1%cmp(4:6), t3%re, pwr_tsc_i%cmp(4:6))
    endif
  end function pwr_tsc_i

  elemental type (ts_real) function tssin(a)
    type (ts_real), intent(in) :: a
    call f_ts_sin(a%re, tssin%re)
  end function tssin

  elemental type (ts_real) function tscos(a)
    type (ts_real), intent(in) :: a
    call f_ts_cos(a%re, tscos%re)
  end function tscos

  elemental type (ts_real) function tstan(a)
    type (ts_real), intent(in) :: a
    call f_ts_tan(a%re, tstan%re)
  end function tstan

  subroutine tssincos(a, s, c)
    type (ts_real), intent(in) :: a
    type (ts_real), intent(out) :: s, c
    call f_ts_sincos(a%re, s%re, c%re)
  end subroutine tssincos

  elemental type (ts_real) function tsasin(a)
    type (ts_real), intent(in) :: a
    call f_ts_asin(a%re, tsasin%re)
  end function tsasin

  elemental type (ts_real) function tsacos(a)
    type (ts_real), intent(in) :: a
    call f_ts_acos(a%re, tsacos%re)
  end function tsacos

  elemental type (ts_real) function tsatan(a)
    type (ts_real), intent(in) :: a
    call f_ts_atan(a%re, tsatan%re)
  end function tsatan

  elemental type (ts_real) function tsatan2(a, b)
    type (ts_real), intent(in) :: a, b
    call f_ts_atan2(a%re, b%re, tsatan2%re)
  end function tsatan2

  elemental type (ts_real) function tsexp(a)
    type (ts_real), intent(in) :: a
    call f_ts_exp(a%re, tsexp%re)
  end function tsexp

  elemental type (ts_complex) function tscexp(a)
    type (ts_complex), intent(in) :: a
    type (ts_real) :: t1, t2, t3
    call f_ts_exp(a%cmp(1:3), t1%re)
    call f_ts_sincos(a%cmp(4:6), t3%re, t2%re)
    call f_ts_mul(t1%re, t2%re, tscexp%cmp(1:3))
    call f_ts_mul(t1%re, t3%re, tscexp%cmp(4:6))
  end function tscexp

  elemental type (ts_real) function tslog(a)
    type (ts_real), intent(in) :: a
    call f_ts_log(a%re, tslog%re)
  end function tslog

  elemental type (ts_complex) function tsclog(a)
    type (ts_complex), intent(in) :: a
    type (ts_real) :: t1, t2, t3
    call f_ts_mul(a%cmp(1:3), a%cmp(1:3), t1%re)
    call f_ts_mul(a%cmp(4:6), a%cmp(4:6), t2%re)
    call f_ts_add(t1%re, t2%re, t3%re)
    call f_ts_log(t3%re, t1%re)
    tsclog%cmp(1:3) = 0.5e0 * t1%re
    call f_ts_atan2(a%cmp(4:6), a%cmp(1:3), tsclog%cmp(4:6))
  end function tsclog

  elemental type (ts_real) function tslog10(a)
    type (ts_real), intent(in) :: a
    call f_ts_log10(a%re, tslog10%re)
  end function tslog10

  elemental type (ts_real) function tssqrt(a)
    type (ts_real), intent(in) :: a
    call f_ts_sqrt(a%re, tssqrt%re)
  end function tssqrt

  elemental type (ts_real) function tssqr(a)
    type (ts_real), intent(in) :: a
    call f_ts_sqr(a%re, tssqr%re)
  end function tssqr

  elemental type (ts_real) function tsnroot(a, n)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: n
    call f_ts_nroot(a%re, n, tsnroot%re)
  end function tsnroot

  elemental type (ts_real) function tssinh(a)
    type (ts_real), intent(in) :: a
    call f_ts_sinh(a%re, tssinh%re)
  end function tssinh

  elemental type (ts_real) function tscosh(a)
    type (ts_real), intent(in) :: a
    call f_ts_cosh(a%re, tscosh%re)
  end function tscosh

  elemental type (ts_real) function tstanh(a)
    type (ts_real), intent(in) :: a
    call f_ts_tanh(a%re, tstanh%re)
  end function tstanh

  subroutine tssincosh(a, s, c)
    type (ts_real), intent(in) :: a
    type (ts_real), intent(out) :: s, c
    call f_ts_sincosh(a%re, s%re, c%re)
  end subroutine tssincosh

  elemental type (ts_real) function tsasinh(a)
    type (ts_real), intent(in) :: a
    call f_ts_asinh(a%re, tsasinh%re)
  end function tsasinh

  elemental type (ts_real) function tsacosh(a)
    type (ts_real), intent(in) :: a
    call f_ts_acosh(a%re, tsacosh%re)
  end function tsacosh

  elemental type (ts_real) function tsatanh(a)
    type (ts_real), intent(in) :: a
    call f_ts_atanh(a%re, tsatanh%re)
  end function tsatanh

  elemental type (ts_real) function tsaint(a)
    type (ts_real), intent(in) :: a
    call f_ts_aint(a%re, tsaint%re)
  end function tsaint

  elemental type (ts_real) function tsfloor(a)
    type (ts_real), intent(in) :: a
    call f_ts_floor(a%re, tsfloor%re)
  end function tsfloor

  elemental type (ts_real) function tsceil(a)
    type (ts_real), intent(in) :: a
    call f_ts_ceil(a%re, tsceil%re)
  end function tsceil

  elemental type (ts_real) function tsanint(a)
    type (ts_real), intent(in) :: a
    call f_ts_nint(a%re, tsanint%re)
  end function tsanint

  elemental integer function tsnint(a)
    type (ts_real), intent(in) :: a
    tsnint = to_int_ts(tsanint(a))
  end function tsnint

  elemental type (ts_real) function tsabs(a)
    type (ts_real), intent(in) :: a
    call f_ts_abs(a%re, tsabs%re)
  end function tsabs

  elemental type (ts_real) function tscabs(tsc)
    type (ts_complex), intent(in) :: tsc
    type (ts_real) :: t1, t2, t3
    call f_ts_mul(tsc%cmp(1:3), tsc%cmp(1:3), t1%re)
    call f_ts_mul(tsc%cmp(4:6), tsc%cmp(4:6), t2%re)
    call f_ts_add(t1%re, t2%re, t3%re)
    call f_ts_sqrt(t3%re, tscabs%re)
  end function tscabs

  elemental type (ts_real) function tssign(a, b) result (c)
    type (ts_real), intent(in) :: a, b
    if (b%re(1) .gt. 0.0e0) then
      if (a%re(1) .gt. 0.0e0) then
        c%re = a%re
      else
        c%re = -a%re
      end if
    else
      if (a%re(1) .gt. 0.0e0) then
        c%re = -a%re
      else
        c%re = a%re
      end if
    endif
  end function tssign

  elemental type (ts_real) function tssign_ts_d(a, b) result (c)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    if (b .gt. 0.0e0) then
      if (a%re(1) .gt. 0.0e0) then
        c%re = a%re
      else
        c%re = -a%re
      end if
    else
      if (a%re(1) .gt. 0.0e0) then
        c%re = -a%re
      else
        c%re = a%re
      end if
    endif
  end function tssign_ts_d

  subroutine tsrand(harvest)
    type (ts_real), intent(out) :: harvest
    call f_ts_rand(harvest%re)
  end subroutine tsrand

  elemental logical function eq_ts(a, b)
    type (ts_real), intent(in) :: a, b
    integer :: r
    call f_ts_comp(a%re, b%re, r)
    eq_ts = (r == 0)
  end function eq_ts

  elemental logical function eq_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ts_comp_ts_d(a%re, b, r)
    eq_ts_d = (r == 0)
  end function eq_ts_d

  elemental logical function eq_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    integer :: r
    call f_ts_comp_d_ts(a, b%re, r)
    eq_d_ts = (r == 0)
  end function eq_d_ts

  elemental logical function eq_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    eq_ts_i = eq_ts_d(a, ts_int_to_float(b))
  end function eq_ts_i

  elemental logical function eq_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    eq_i_ts = eq_d_ts(ts_int_to_float(a), b)
  end function eq_i_ts

  elemental logical function eq_tsc(a, b)
    type (ts_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_ts_comp(a%cmp(1:3), b%cmp(1:3), i1)
    call f_ts_comp(a%cmp(4:6), b%cmp(4:6), i2)
    eq_tsc = (i1 == 0 .and. i2 == 0)
  end function eq_tsc

  elemental logical function eq_tsc_ts(a, b)
    type (ts_complex), intent(in) :: a
    type (ts_real), intent(in) :: b
    integer :: i1
    call f_ts_comp(a%cmp(1:3), b%re, i1)
    eq_tsc_ts = (i1 == 0 .and. all(a%cmp(4:6) == 0.e0))
  end function eq_tsc_ts

  elemental logical function eq_ts_tsc(a, b)
    type (ts_real), intent(in) :: a
    type (ts_complex), intent(in) :: b
    integer :: i1
    call f_ts_comp(a%re, b%cmp(1:3), i1)
    eq_ts_tsc = (i1 == 0 .and. all(b%cmp(4:6) == 0.e0))
  end function eq_ts_tsc

  elemental logical function ne_ts(a, b)
    type (ts_real), intent(in) :: a, b
    ne_ts = .not. eq_ts(a, b)
  end function ne_ts

  elemental logical function ne_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    ne_ts_d = .not. eq_ts_d(a, b)
  end function ne_ts_d

  elemental logical function ne_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    ne_d_ts = .not. eq_d_ts(a, b)
  end function ne_d_ts

  elemental logical function ne_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    ne_ts_i = .not. eq_ts_i(a, b)
  end function ne_ts_i

  elemental logical function ne_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    ne_i_ts = .not. eq_i_ts(a, b)
  end function ne_i_ts

  elemental logical function ne_tsc(a, b)
    type (ts_complex), intent(in) :: a, b
    ne_tsc = .not. eq_tsc(a, b)
  end function ne_tsc

  elemental logical function ne_tsc_ts(a, b)
    type (ts_complex), intent(in) :: a
    type (ts_real), intent(in) :: b
    ne_tsc_ts = .not. eq_tsc_ts(a, b)
  end function ne_tsc_ts

  elemental logical function ne_ts_tsc(a, b)
    type (ts_real), intent(in) :: a
    type (ts_complex), intent(in) :: b
    ne_ts_tsc = .not. eq_ts_tsc(a, b)
  end function ne_ts_tsc

  elemental logical function gt_ts(a, b)
    type (ts_real), intent(in) :: a, b
    integer :: r
    call f_ts_comp(a%re, b%re, r)
    gt_ts = (r == 1)
  end function gt_ts

  elemental logical function gt_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ts_comp_ts_d(a%re, b, r)
    gt_ts_d = (r == 1)
  end function gt_ts_d

  elemental logical function gt_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    integer :: r
    call f_ts_comp_d_ts(a, b%re, r)
    gt_d_ts = (r == 1)
  end function gt_d_ts

  elemental logical function gt_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    gt_ts_i = gt_ts_d(a, ts_int_to_float(b))
  end function gt_ts_i

  elemental logical function gt_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    gt_i_ts = gt_d_ts(ts_int_to_float(a), b)
  end function gt_i_ts

  elemental logical function lt_ts(a, b)
    type (ts_real), intent(in) :: a, b
    integer :: r
    call f_ts_comp(a%re, b%re, r)
    lt_ts = (r == -1)
  end function lt_ts

  elemental logical function lt_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ts_comp_ts_d(a%re, b, r)
    lt_ts_d = (r == -1)
  end function lt_ts_d

  elemental logical function lt_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    integer :: r
    call f_ts_comp_d_ts(a, b%re, r)
    lt_d_ts = (r == -1)
  end function lt_d_ts

  elemental logical function lt_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    lt_ts_i = lt_ts_d(a, ts_int_to_float(b))
  end function lt_ts_i

  elemental logical function lt_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    lt_i_ts = lt_d_ts(ts_int_to_float(a), b)
  end function lt_i_ts

  elemental logical function ge_ts(a, b)
    type (ts_real), intent(in) :: a, b
    integer :: r
    call f_ts_comp(a%re, b%re, r)
    ge_ts = (r >= 0)
  end function ge_ts

  elemental logical function ge_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ts_comp_ts_d(a%re, b, r)
    ge_ts_d = (r >= 0)
  end function ge_ts_d

  elemental logical function ge_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    integer :: r
    call f_ts_comp_d_ts(a, b%re, r)
    ge_d_ts = (r >= 0)
  end function ge_d_ts

  elemental logical function ge_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    ge_ts_i = ge_ts_d(a, ts_int_to_float(b))
  end function ge_ts_i

  elemental logical function ge_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    ge_i_ts = ge_d_ts(ts_int_to_float(a), b)
  end function ge_i_ts

  elemental logical function le_ts(a, b)
    type (ts_real), intent(in) :: a, b
    integer :: r
    call f_ts_comp(a%re, b%re, r)
    le_ts = (r <= 0)
  end function le_ts

  elemental logical function le_ts_d(a, b)
    type (ts_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ts_comp_ts_d(a%re, b, r)
    le_ts_d = (r <= 0)
  end function le_ts_d

  elemental logical function le_d_ts(a, b)
    real*4, intent(in) :: a
    type (ts_real), intent(in) :: b
    integer :: r
    call f_ts_comp_d_ts(a, b%re, r)
    le_d_ts = (r <= 0)
  end function le_d_ts

  elemental logical function le_ts_i(a, b)
    type (ts_real), intent(in) :: a
    integer, intent(in) :: b
    le_ts_i = le_ts_d(a, ts_int_to_float(b))
  end function le_ts_i

  elemental logical function le_i_ts(a, b)
    integer, intent(in) :: a
    type (ts_real), intent(in) :: b
    le_i_ts = le_d_ts(ts_int_to_float(a), b)
  end function le_i_ts

  subroutine tsinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ts_real), intent(inout) :: q1
    type (ts_real), intent(inout), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call tsinp(u, q1%re)
    if (present(q2)) call tsinp(u, q2%re)
    if (present(q3)) call tsinp(u, q3%re)
    if (present(q4)) call tsinp(u, q4%re)
    if (present(q5)) call tsinp(u, q5%re)
    if (present(q6)) call tsinp(u, q6%re)
    if (present(q7)) call tsinp(u, q7%re)
    if (present(q8)) call tsinp(u, q8%re)
    if (present(q9)) call tsinp(u, q9%re)
  end subroutine tsinpq

  subroutine tscinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ts_complex), intent(inout) :: q1
    type (ts_complex), intent(inout), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call tsinp(u, q1%cmp(1:3))
    call tsinp(u, q1%cmp(4:6))
    if (present(q2)) then
      call tsinp(u, q2%cmp(1:3))
      call tsinp(u, q2%cmp(4:6))
    end if
    if (present(q3)) then
      call tsinp(u, q3%cmp(1:3))
      call tsinp(u, q3%cmp(4:6))
    end if
    if (present(q4)) then
      call tsinp(u, q4%cmp(1:3))
      call tsinp(u, q4%cmp(4:6))
    end if
    if (present(q5)) then
      call tsinp(u, q5%cmp(1:3))
      call tsinp(u, q5%cmp(4:6))
    end if
    if (present(q6)) then
      call tsinp(u, q6%cmp(1:3))
      call tsinp(u, q6%cmp(4:6))
    end if
    if (present(q7)) then
      call tsinp(u, q7%cmp(1:3))
      call tsinp(u, q7%cmp(4:6))
    end if
    if (present(q8)) then
      call tsinp(u, q8%cmp(1:3))
      call tsinp(u, q8%cmp(4:6))
    end if
    if (present(q9)) then
      call tsinp(u, q9%cmp(1:3))
      call tsinp(u, q9%cmp(4:6))
    end if
  end subroutine tscinpq

  subroutine tsoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ts_real), intent(in) :: q1
    type (ts_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call tsout(u, q1%re)
    if (present(q2)) call tsout(u, q2%re)
    if (present(q3)) call tsout(u, q3%re)
    if (present(q4)) call tsout(u, q4%re)
    if (present(q5)) call tsout(u, q5%re)
    if (present(q6)) call tsout(u, q6%re)
    if (present(q7)) call tsout(u, q7%re)
    if (present(q8)) call tsout(u, q8%re)
    if (present(q9)) call tsout(u, q9%re)
  end subroutine tsoutq

  subroutine tscoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ts_complex), intent(in) :: q1
    type (ts_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call tsout(u, q1%cmp(1:3))
    call tsout(u, q1%cmp(4:6))
    if (present(q2)) then
      call tsout(u, q2%cmp(1:3))
      call tsout(u, q2%cmp(4:6))
    end if
    if (present(q3)) then
      call tsout(u, q3%cmp(1:3))
      call tsout(u, q3%cmp(4:6))
    end if
    if (present(q4)) then
      call tsout(u, q4%cmp(1:3))
      call tsout(u, q4%cmp(4:6))
    end if
    if (present(q5)) then
      call tsout(u, q5%cmp(1:3))
      call tsout(u, q5%cmp(4:6))
    end if
    if (present(q6)) then
      call tsout(u, q6%cmp(1:3))
      call tsout(u, q6%cmp(4:6))
    end if
    if (present(q7)) then
      call tsout(u, q7%cmp(1:3))
      call tsout(u, q7%cmp(4:6))
    end if
    if (present(q8)) then
      call tsout(u, q8%cmp(1:3))
      call tsout(u, q8%cmp(4:6))
    end if
    if (present(q9)) then
      call tsout(u, q9%cmp(1:3))
      call tsout(u, q9%cmp(4:6))
    end if
  end subroutine tscoutq

  pure type (ts_real) function ts_pi()
    call f_ts_pi(ts_pi%re)
  end function ts_pi

  elemental type (ts_real) function tshuge(a)
    type (ts_real), intent(in) :: a
    tshuge = ts_huge
  end function tshuge

  elemental type (ts_real) function ts_safe_huge(a)
    type (ts_real), intent(in) :: a
    ts_safe_huge = ts_real((/3.4028235e+38, 0.0e0, 0.0e0/))
  end function ts_safe_huge

  elemental type (ts_real) function tstiny(a)
    type (ts_real), intent(in) :: a
    tstiny = ts_tiny
  end function tstiny

  elemental type (ts_real) function tsepsilon(a)
    type (ts_real), intent(in) :: a
    tsepsilon = ts_eps
  end function tsepsilon

  elemental integer function ts_radix(a)
    type (ts_real), intent(in) :: a
    ts_radix = 2
  end function ts_radix

  elemental integer function ts_digits(a)
    type (ts_real), intent(in) :: a
  ts_digits = 70
  end function ts_digits

  elemental integer function ts_max_expn(a)
    type (ts_real), intent(in) :: a
  ts_max_expn = 127
  end function ts_max_expn

  elemental integer function ts_min_expn(a)
    type (ts_real), intent(in) :: a
  ts_min_expn = -78
  end function ts_min_expn

  elemental integer function ts_precision(a)
    type (ts_real), intent(in) :: a
  ts_precision = 21
  end function ts_precision

  elemental integer function ts_range(a)
    type (ts_real), intent(in) :: a
  ts_range = 37
  end function ts_range

  elemental type (ts_real) function ts_nan(a)
    type (ts_real), intent(in) :: a
    call f_ts_nan(ts_nan%re)
  end function ts_nan

  elemental type (ts_real) function ts_aimag(a)
    type (ts_complex), intent(in) :: a
    ts_aimag%re = a%cmp(4:6)
  end function ts_aimag

  elemental type (ts_real) function tsmin2(a, b)
    type (ts_real), intent(in) :: a, b
    integer :: r
    call f_ts_comp(a%re, b%re, r)
    if (r <= 0) then
      tsmin2 = a
    else
      tsmin2 = b
    end if
  end function tsmin2

  elemental type (ts_real) function tsmin(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (ts_real), intent(in) :: a1, a2, a3
    type (ts_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    tsmin = tsmin2(tsmin2(a1, a2), a3)
    if (present(a4)) tsmin = tsmin2(tsmin, a4)
    if (present(a5)) tsmin = tsmin2(tsmin, a5)
    if (present(a6)) tsmin = tsmin2(tsmin, a6)
    if (present(a7)) tsmin = tsmin2(tsmin, a7)
    if (present(a8)) tsmin = tsmin2(tsmin, a8)
    if (present(a9)) tsmin = tsmin2(tsmin, a9)
  end function tsmin

  elemental type (ts_real) function tsmax2(a, b)
    type (ts_real), intent(in) :: a, b
    integer :: r
    call f_ts_comp(a%re, b%re, r)
    if (r >= 0) then
      tsmax2 = a
    else
      tsmax2 = b
    end if
  end function tsmax2

  elemental type (ts_real) function tsmax(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (ts_real), intent(in) :: a1, a2, a3
    type (ts_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    tsmax = tsmax2(tsmax2(a1, a2), a3)
    if (present(a4)) tsmax = tsmax2(tsmax, a4)
    if (present(a5)) tsmax = tsmax2(tsmax, a5)
    if (present(a6)) tsmax = tsmax2(tsmax, a6)
    if (present(a7)) tsmax = tsmax2(tsmax, a7)
    if (present(a8)) tsmax = tsmax2(tsmax, a8)
    if (present(a9)) tsmax = tsmax2(tsmax, a9)
  end function tsmax

  elemental type (ts_real) function tsmod(a, b)
    type (ts_real), intent(in) :: a, b
    type (ts_real) :: s1, s2
    call f_ts_div(a%re, b%re, s1%re)
    call f_ts_aint(s1%re, s2%re)
    call f_ts_mul(s2%re, b%re, s1%re)
    call f_ts_sub(a%re, s1%re, tsmod%re)
  end function tsmod

subroutine tsinp(iu, a)
  implicit none
  integer iu
  character*80 cs
  real*4 a(3)

  read (iu, '(a)', end = 100) cs
  call tsinpc(cs, a)
  goto 110

100 write (6, 1)
1 format ('*** tsinp: End-of-file encountered.')
  stop

110 return
end subroutine

subroutine tsinpc(a, b)
  implicit none
  integer i, id, ie, inz, ip, is, k, ln, lnn, beg
  parameter (ln = 80)
  real*4 bi
  character*80 a
  character*1 ai
  character*10 dig
  character*16 ca
  parameter (dig = '0123456789')
  real*4 b(3), f(3), s0(3), s1(3), s2(3)

  id = 0
  ip = -1
  is = 0
  inz = 0
  s1(1) = 0.e0
  s1(2) = 0.e0
  s1(3) = 0.e0

  beg = 0
  do i = 1, 80
    if (a(i:i) /= ' ') then
      beg = i
      goto 80
    end if
  end do

  goto 210
80 continue

  do i = beg, 80
    if (a(i:i) == ' ') then
      lnn = i - 1
      goto 90
    end if
  enddo

  lnn = 80
90 continue

  do i = beg, lnn
    ai = a(i:i)
    if (ai .eq. '.') then
      if (ip >= 0) goto 210
      ip = id
      inz = 1
    elseif (ai .eq. '+') then
      if (id .ne. 0 .or. ip >= 0 .or. is .ne. 0) goto 210
      is = 1
    elseif (ai .eq. '-') then
      if (id .ne. 0 .or. ip >= 0 .or. is .ne. 0) goto 210
      is = -1
    elseif (ai .eq. 'e' .or. ai .eq. 'E' .or. ai .eq. 'd' .or. ai .eq. 'D') then
      goto 100
    elseif (index(dig, ai) .eq. 0) then
      goto 210
    else
      bi = index(dig, ai) - 1
      if (inz > 0 .or. bi > 0.e0) then
        inz = 1
        id = id + 1
        call f_ts_mul_ts_d(s1, 10.e0, s0)
        f(1) = bi
        f(2) = 0.e0
        f(3) = 0.e0
        call f_ts_add(s0, f, s1)
      endif
    endif
  enddo

100 continue
  if (is .eq. -1) then
    s1(1) = -s1(1)
    s1(2) = -s1(2)
    s1(3) = -s1(3)
  endif
  k = i
  if (ip == -1) ip = id
  ie = 0
  is = 0
  ca = ' '

  do i = k + 1, lnn
    ai = a(i:i)
    if (ai .eq. ' ') then
    elseif (ai .eq. '+') then
      if (ie .ne. 0 .or. is .ne. 0) goto 210
      is = 1
    elseif (ai .eq. '-') then
      if (ie .ne. 0 .or. is .ne. 0) goto 210
      is = -1
    elseif (index(dig, ai) .eq. 0) then
      goto 210
    else
      ie = ie + 1
      if (ie .gt. 3) goto 210
      ca(ie:ie) = ai
    endif
  enddo

  ie = dddigin(ca, 4)
  if (is .eq. -1) ie = -ie
  ie = ie + ip - id
  s0(1) = 10.e0
  s0(2) = 0.e0
  s0(3) = 0.e0
  call f_ts_npwr(s0, ie, s2)
  call f_ts_mul(s1, s2, b)
  goto 220

210 write (6, 1) a
1 format ('*** tsinpc: Syntax error in literal string: ', a)
  stop

220 return
end subroutine

subroutine tsout(iu, a)
  implicit none
  integer iu
  character cs(57)
  real*4 a(3)

  call tsoutc(a, cs)
  write (iu, '(57a)') cs
end subroutine

subroutine tsoutc(a, b)
  implicit none
  real*4 a(3)
  character b(57)

  b(1) = ' '
  b(2) = ' '
  call f_ts_swrite(a, 47, b(3), 55)
end subroutine

  real*4 function dddigin(ca, n)
    implicit none
    real*4 d1
    character*(*), ca
    character*16 digits
    integer i, k, n
    parameter (digits = '0123456789')

    d1 = 0.e0
    do i = 1, n
      k = index(digits, ca(i:i)) - 1
      if (k < 0) then
        write (6, *) 'dddigin: non-digit in character string'
      elseif (k <= 9) then
        d1 = 10.e0 * d1 + k
      endif
    enddo
    dddigin = d1
  end function dddigin

elemental real*4 function ts_int_to_float(i)
  implicit none
  integer, intent(in) :: i
  intrinsic :: real
  ts_int_to_float = real(i, kind=4)
end function ts_int_to_float

end module tsmodule
