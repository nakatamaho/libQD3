!  dsmod.f
!
!  This work was supported by the Director, Office of Science, Division
!  of Mathematical, Information, and Computational Sciences of the
!  U.S. Department of Energy under contract number DE-AC03-76SF00098.
!
!  Copyright (c) 2000-2008
!
!  Fortran-90 module file to use with double-single numbers.
!
!  Yozo Hida
!  David H Bailey    2008-03-20

module dsmodule
  use ddmodule, only: dd_real, dd_complex
  use qdmodule, only: qd_real, qd_complex
  use dsext
  implicit none

  type ds_real
    sequence
    real*4 :: re(2)
  end type ds_real

  type ds_complex
    sequence
    real*4 :: cmp(4)
  end type ds_complex

  real*4 d_ds_eps
  parameter (d_ds_eps = 1.4210855e-14)

  type (ds_real) ds_one, ds_zero, ds_eps, ds_huge, ds_tiny
  parameter (ds_one = ds_real((/1.0e0, 0.0e0/)), &
             ds_zero = ds_real((/0.0e0, 0.0e0/)))
  parameter (ds_eps = ds_real((/d_ds_eps, 0.0e0/)))
  parameter (ds_huge = ds_real((/3.4028235e+38, 0.0e0/)))
  parameter (ds_tiny = ds_real((/1.9721523e-31, 0.0e0/)))


  interface assignment (=)
    module procedure assign_ds_str
    module procedure assign_ds
    module procedure assign_ds_d
    module procedure assign_d_ds
    module procedure assign_ds_i
    module procedure assign_i_ds
    module procedure assign_dsc
    module procedure assign_dsc_ds
    module procedure assign_ds_dsc
    module procedure assign_dsc_d
    module procedure assign_d_dsc
    module procedure assign_dsc_dc
    module procedure assign_dc_dsc
    module procedure assign_dsc_i
    module procedure assign_ds_dp
    module procedure assign_dp_ds
    module procedure assign_dsc_zc
    module procedure assign_ds_dd
    module procedure assign_dd_ds
    module procedure assign_dsc_ddc
    module procedure assign_ddc_dsc
    module procedure assign_ds_qd
    module procedure assign_qd_ds
    module procedure assign_dsc_qdc
    module procedure assign_qdc_dsc
  end interface

  interface operator (+)
    module procedure add_ds
    module procedure add_ds_d
    module procedure add_d_ds
    module procedure add_ds_i
    module procedure add_i_ds
    module procedure add_dsc
    module procedure add_dsc_ds
    module procedure add_ds_dsc
    module procedure add_dsc_d
    module procedure add_d_dsc
    module procedure add_dsc_i
    module procedure add_i_dsc
    module procedure add_ds_dp
    module procedure add_dp_ds
    module procedure add_dsc_dp
    module procedure add_dp_dsc
  end interface

  interface operator (-)
    module procedure sub_ds
    module procedure sub_ds_d
    module procedure sub_d_ds
    module procedure neg_ds
    module procedure sub_dsc
    module procedure sub_dsc_ds
    module procedure sub_ds_dsc
    module procedure sub_dsc_d
    module procedure sub_d_dsc
    module procedure neg_dsc
    module procedure sub_ds_i
    module procedure sub_i_ds
    module procedure sub_dsc_i
    module procedure sub_i_dsc
    module procedure sub_ds_dp
    module procedure sub_dp_ds
    module procedure sub_dsc_dp
    module procedure sub_dp_dsc
  end interface

  interface operator (*)
    module procedure mul_ds
    module procedure mul_ds_d
    module procedure mul_d_ds
    module procedure mul_ds_i
    module procedure mul_i_ds
    module procedure mul_dsc
    module procedure mul_dsc_ds
    module procedure mul_ds_dsc
    module procedure mul_dsc_d
    module procedure mul_d_dsc
    module procedure mul_dsc_i
    module procedure mul_i_dsc
    module procedure mul_ds_dp
    module procedure mul_dp_ds
    module procedure mul_dsc_dp
    module procedure mul_dp_dsc
  end interface

  interface operator (/)
    module procedure div_ds
    module procedure div_ds_d
    module procedure div_d_ds
    module procedure div_ds_i
    module procedure div_i_ds
    module procedure div_dsc
    module procedure div_dsc_ds
    module procedure div_ds_dsc
    module procedure div_dsc_d
    module procedure div_dsc_i
    module procedure div_i_dsc
    module procedure div_d_dsc
    module procedure div_ds_dp
    module procedure div_dp_ds
    module procedure div_dsc_dp
    module procedure div_dp_dsc
  end interface

  interface operator (**)
    module procedure pwr_ds
    module procedure pwr_ds_i
    module procedure pwr_d_ds
    module procedure pwr_dsc_i
    module procedure pwr_dsc
    module procedure pwr_dsc_ds
  end interface

  interface dsreal
    module procedure to_ds_i
    module procedure to_ds_d
    module procedure to_ds_ds
    module procedure to_ds_str
    module procedure to_ds_dsc
    module procedure to_ds_dp
    module procedure to_ds_dd
    module procedure to_ds_qd
  end interface

  interface dscomplex
     module procedure to_dsc_ds
     module procedure to_dsc_ds2
     module procedure to_dsc_d
     module procedure to_dsc_dc
    module procedure to_dsc_zc
    module procedure to_dsc_ddc
    module procedure to_dsc_qdc
  end interface

  interface real
    module procedure to_d_ds
    module procedure to_ds_dsc
  end interface

  interface int
    module procedure to_int_ds
  end interface

  interface sin
    module procedure dssin
    module procedure dsc_sin
  end interface
  interface cos
    module procedure dscos
    module procedure dsc_cos
  end interface
  interface tan
    module procedure dstan
    module procedure dsc_tan
  end interface
  interface sincos
    module procedure dssincos
  end interface

  interface asin
    module procedure dsasin
    module procedure dsc_asin
  end interface
  interface acos
    module procedure dsacos
    module procedure dsc_acos
  end interface
  interface atan
    module procedure dsatan
    module procedure dsc_atan
  end interface
  interface atan2
    module procedure dsatan2
  end interface

  interface exp
    module procedure dsexp
    module procedure dscexp
  end interface
  interface log
    module procedure dslog
    module procedure dsclog
  end interface
  interface log10
    module procedure dslog10
  end interface

  interface sqrt
    module procedure dssqrt
    module procedure dsc_sqrt
  end interface
  interface sqr
    module procedure dssqr
  end interface
  interface nroot
    module procedure dsnroot
  end interface

  interface sinh
    module procedure dssinh
    module procedure dsc_sinh
  end interface
  interface cosh
    module procedure dscosh
    module procedure dsc_cosh
  end interface
  interface tanh
    module procedure dstanh
    module procedure dsc_tanh
  end interface
  interface sincosh
    module procedure dssincosh
  end interface

  interface asinh
    module procedure dsasinh
    module procedure dsc_asinh
  end interface
  interface acosh
    module procedure dsacosh
    module procedure dsc_acosh
  end interface
  interface atanh
    module procedure dsatanh
    module procedure dsc_atanh
  end interface

  interface aint
    module procedure dsaint
  end interface

  interface anint
    module procedure dsanint
  end interface

  interface nint
    module procedure dsnint
  end interface

  interface abs
    module procedure dsabs
    module procedure dscabs
  end interface

  interface sign
    module procedure dssign
    module procedure dssign_ds_d
  end interface

  interface random_number
    module procedure dsrand
  end interface

  interface aimag
    module procedure ds_aimag
  end interface

  interface operator (==)
    module procedure eq_ds
    module procedure eq_ds_d
    module procedure eq_d_ds
    module procedure eq_ds_i
    module procedure eq_i_ds
    module procedure eq_dsc
    module procedure eq_dsc_ds
    module procedure eq_ds_dsc
    module procedure eq_dsc_d
    module procedure eq_d_dsc
    module procedure eq_ds_dp
    module procedure eq_dp_ds
    module procedure eq_dsc_dp
    module procedure eq_dp_dsc
  end interface

  interface operator (/=)
    module procedure ne_ds
    module procedure ne_ds_d
    module procedure ne_d_ds
    module procedure ne_ds_i
    module procedure ne_i_ds
    module procedure ne_dsc
    module procedure ne_dsc_ds
    module procedure ne_ds_dsc
    module procedure ne_dsc_d
    module procedure ne_d_dsc
    module procedure ne_ds_dp
    module procedure ne_dp_ds
    module procedure ne_dsc_dp
    module procedure ne_dp_dsc
  end interface

  interface operator (>)
    module procedure gt_ds
    module procedure gt_ds_d
    module procedure gt_d_ds
    module procedure gt_ds_i
    module procedure gt_i_ds
    module procedure gt_ds_dp
    module procedure gt_dp_ds
  end interface

  interface operator (<)
    module procedure lt_ds
    module procedure lt_ds_d
    module procedure lt_d_ds
    module procedure lt_ds_i
    module procedure lt_i_ds
    module procedure lt_ds_dp
    module procedure lt_dp_ds
  end interface

  interface operator (>=)
    module procedure ge_ds
    module procedure ge_ds_d
    module procedure ge_d_ds
    module procedure ge_ds_i
    module procedure ge_i_ds
    module procedure ge_ds_dp
    module procedure ge_dp_ds
  end interface

  interface operator (<=)
    module procedure le_ds
    module procedure le_ds_d
    module procedure le_d_ds
    module procedure le_ds_i
    module procedure le_i_ds
    module procedure le_ds_dp
    module procedure le_dp_ds
  end interface

  interface read_scalar
    module procedure dsinpq
    module procedure dscinpq
  end interface

  interface write_scalar
    module procedure dsoutq
    module procedure dscoutq
  end interface

  interface dsread
    module procedure dsinpq
  end interface

  interface dswrite
    module procedure dsoutq
  end interface

  interface dscread
    module procedure dscinpq
  end interface

  interface dscwrite
    module procedure dscoutq
  end interface

  interface dble
    module procedure to_dp_ds
    module procedure to_dp_dsc
  end interface

  interface cmplx
    module procedure to_dc_dsc
  end interface

  interface conjg
    module procedure dscconjg
  end interface

  interface min
    module procedure dsmin
    module procedure dsmin2
  end interface
  interface max
    module procedure dsmax
    module procedure dsmax2
  end interface
  interface mod
    module procedure dsmod
  end interface

  interface dspi
    module procedure ds_pi
  end interface

  interface huge
    module procedure dshuge
  end interface

  interface safe_huge
    module procedure ds_safe_huge
  end interface

  interface tiny
    module procedure dstiny
  end interface

  interface epsilon
    module procedure dsepsilon
  end interface

  interface radix
    module procedure ds_radix
  end interface

  interface digits
    module procedure ds_digits
  end interface

  interface maxexponent
    module procedure ds_max_expn
  end interface

  interface minexponent
    module procedure ds_min_expn
  end interface

  interface nan
    module procedure ds_nan
  end interface

  interface floor
    module procedure dsfloor
  end interface

  interface ceil
    module procedure dsceil
  end interface

  interface ceiling
    module procedure dsceil
  end interface

  interface precision
    module procedure ds_precision
  end interface

  interface range
    module procedure ds_range
  end interface

  interface hypot
    module procedure dshypot
  end interface

  interface modulo
    module procedure dsmodulo
  end interface

  interface dim
    module procedure dsdim
  end interface

  interface ddreal
    module procedure to_dd_ds
  end interface

  interface ddcomplex
    module procedure to_ddc_dsc
  end interface

  interface qdreal
    module procedure to_qd_ds
  end interface

  interface qdcomplex
    module procedure to_qdc_dsc
  end interface

contains

! Assignments
  subroutine assign_ds_str(a, s)
    type (ds_real), intent(inout) :: a
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call dsinpc (t, a%re)
  end subroutine assign_ds_str

  elemental subroutine assign_ds (a, b)
    type (ds_real), intent(inout) :: a
    type (ds_real), intent(in) :: b
    a%re = b%re
  end subroutine assign_ds


  elemental subroutine assign_ds_d(a, d)
    type (ds_real), intent(inout) :: a
    real*4, intent(in) :: d
    a%re(1) = d
    a%re(2) = 0.0e0
  end subroutine assign_ds_d

  elemental subroutine assign_d_ds(d, a)
    real*4, intent(inout) :: d
    type (ds_real), intent(in) :: a
    d = a%re(1)
  end subroutine assign_d_ds

  elemental subroutine assign_ds_i(a, i)
    type (ds_real), intent(inout) :: a
    integer, intent(in) :: i
    a%re = ds_int_limbs(i)
  end subroutine assign_ds_i

  elemental subroutine assign_i_ds(i, a)
    integer, intent(inout) :: i
    type (ds_real), intent(in) :: a
    i = to_int_ds(a)
  end subroutine assign_i_ds

  elemental subroutine assign_dsc (a, b)
    type (ds_complex), intent(inout) :: a
    type (ds_complex), intent(in) :: b
    a%cmp = b%cmp
  end subroutine assign_dsc

  elemental subroutine assign_dsc_ds (dsc, ds)
    type (ds_complex), intent (inout) :: dsc
    type (ds_real), intent(in) :: ds
    dsc%cmp(1:2) = ds%re
    dsc%cmp(3:4) = 0.e0
  end subroutine assign_dsc_ds

  elemental subroutine assign_ds_dsc (ds, dsc)
    type (ds_real), intent (inout) :: ds
    type (ds_complex), intent(in) :: dsc
    ds%re = dsc%cmp(1:2)
  end subroutine assign_ds_dsc

  elemental subroutine assign_dsc_d (dsc, d)
    type (ds_complex), intent (inout) :: dsc
    real*4, intent(in) :: d
    dsc%cmp(1) = d
    dsc%cmp(2:4) = 0.e0
  end subroutine assign_dsc_d

  elemental subroutine assign_dsc_i (dsc, i)
    type (ds_complex), intent (inout) :: dsc
    integer, intent(in) :: i
    dsc%cmp = 0.e0
    dsc%cmp(1:2) = ds_int_limbs(i)
  end subroutine assign_dsc_i

  elemental subroutine assign_d_dsc (d, dsc)
    real*4, intent(inout) :: d
    type (ds_complex), intent (in) :: dsc
    d = dsc%cmp(1)
  end subroutine assign_d_dsc

  elemental subroutine assign_dsc_dc (dsc, dc)
    type (ds_complex), intent (inout) :: dsc
    complex (kind (0.e0)), intent (in) :: dc
    dsc%cmp(1) = dble (dc)
    dsc%cmp(2) = 0.e0
    dsc%cmp(3) = aimag (dc)
    dsc%cmp(4) = 0.e0
  end subroutine assign_dsc_dc

  elemental subroutine assign_dc_dsc (dc, dsc)
    complex (kind (0.e0)), intent (inout) :: dc
    type (ds_complex), intent (in) :: dsc
    dc = cmplx (dsc%cmp(1), dsc%cmp(3), kind (0.e0))
  end subroutine assign_dc_dsc


! Conversions

  elemental type (ds_real) function to_ds_i(ia)
    integer, intent(in) :: ia
    to_ds_i%re = ds_int_limbs(ia)
  end function to_ds_i

  elemental type (ds_real) function to_ds_d(a)
    real*4, intent(in) :: a
    to_ds_d%re(1) = a
    to_ds_d%re(2) = 0.0e0
  end function to_ds_d

  elemental type (ds_real) function to_ds_ds(a)
    type (ds_real), intent(in) :: a
    to_ds_ds%re = a%re
  end function to_ds_ds

  elemental real*4 function to_d_ds(a)
    type (ds_real), intent(in) :: a
    to_d_ds = a%re(1)
  end function to_d_ds

  elemental integer function to_int_ds(a)
    type (ds_real), intent(in) :: a
    type (ds_real) :: t
    t = dsaint(a)
    to_int_ds = int(t%re(1)) + int(sum(t%re(2:)))
  end function to_int_ds

  type (ds_real) function to_ds_str(s)
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call dsinpc (t, to_ds_str%re)
  end function to_ds_str

  elemental type (ds_real) function to_ds_dsc(dsc)
    type (ds_complex), intent(in) :: dsc
    to_ds_dsc%re = dsc%cmp(1:2)
  end function to_ds_dsc

  elemental type (ds_complex) function to_dsc_ds(ds)
    type (ds_real), intent(in) :: ds
    to_dsc_ds%cmp(1:2) = ds%re
    to_dsc_ds%cmp(3:4) = 0.e0
  end function to_dsc_ds

  elemental type (ds_complex) function to_dsc_ds2(x, y)
    type (ds_real), intent(in) :: x, y
    to_dsc_ds2%cmp(1:2) = x%re
    to_dsc_ds2%cmp(3:4) = y%re
  end function to_dsc_ds2

  elemental type (ds_complex) function to_dsc_d(d)
    real*4, intent(in) :: d
    to_dsc_d%cmp(1) = d
    to_dsc_d%cmp(2:4) = 0.e0
  end function to_dsc_d

  elemental complex (kind (0.e0)) function to_dc_dsc (dsc)
    type (ds_complex), intent (in) :: dsc
    to_dc_dsc = cmplx (dsc%cmp(1), dsc%cmp(3), kind (0.e0))
  end function to_dc_dsc

  elemental type (ds_complex) function to_dsc_dc (dc)
    complex (kind (0.e0)), intent(in) :: dc
    to_dsc_dc%cmp(1) = dble (dc)
    to_dsc_dc%cmp(2) = 0.e0
    to_dsc_dc%cmp(3) = aimag (dc)
    to_dsc_dc%cmp(4) = 0.e0
  end function to_dsc_dc

!  Complex conjugation
  elemental type (ds_complex) function dscconjg (dsc)
    type (ds_complex), intent(in) :: dsc
    dscconjg%cmp(1:2) = dsc%cmp(1:2)
    dscconjg%cmp(3:4) = - dsc%cmp(3:4)
  end function dscconjg

! Adsitions
  elemental type (ds_real) function add_ds(a, b)
    type (ds_real), intent(in) :: a, b
    call f_ds_add(a%re, b%re, add_ds%re)
  end function add_ds

  elemental type (ds_real) function add_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ds_add_ds_d(a%re, b, add_ds_d%re)
  end function add_ds_d

  elemental type (ds_real) function add_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_add_ds_d(b%re, a, add_d_ds%re)
  end function add_d_ds

  elemental type (ds_real) function add_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_add(b%re, ds_int_limbs(a), add_i_ds%re)
  end function add_i_ds

  elemental type (ds_real) function add_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    call f_ds_add(a%re, ds_int_limbs(b), add_ds_i%re)
  end function add_ds_i

  elemental type (ds_complex) function add_dsc(a, b)
    type (ds_complex), intent(in) :: a, b
    call f_ds_add (a%cmp(1:2), b%cmp(1:2), add_dsc%cmp(1:2))
    call f_ds_add (a%cmp(3:4), b%cmp(3:4), add_dsc%cmp(3:4))
  end function add_dsc

  elemental type (ds_complex) function add_dsc_d(a, b)
    type (ds_complex), intent(in) :: a
    real*4, intent(in) :: b
    type (ds_real) :: dsb
    dsb%re(1) = b
    dsb%re(2) = 0.e0
    call f_ds_add (a%cmp(1:2), dsb%re, add_dsc_d%cmp(1:2))
    add_dsc_d%cmp(3:4) = a%cmp(3:4)
  end function add_dsc_d

  elemental type (ds_complex) function add_d_dsc(a, b)
    real*4, intent(in) :: a
    type (ds_complex), intent(in) :: b
    type (ds_real) dsa
    dsa%re(1) = a
    dsa%re(2) = 0.e0
    call f_ds_add (dsa%re, b%cmp(1:2), add_d_dsc%cmp(1:2))
    add_d_dsc%cmp(3:4) = b%cmp(3:4)
  end function add_d_dsc

  elemental type (ds_complex) function add_dsc_ds(a, b)
    type (ds_complex), intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_add (a%cmp(1:2), b%re, add_dsc_ds%cmp(1:2))
    add_dsc_ds%cmp(3:4) = a%cmp(3:4)
  end function add_dsc_ds

  elemental type (ds_complex) function add_ds_dsc(a, b)
    type (ds_real), intent(in) :: a
    type (ds_complex), intent(in) :: b
    call f_ds_add (a%re, b%cmp(1:2), add_ds_dsc%cmp(1:2))
    add_ds_dsc%cmp(3:4) = b%cmp(3:4)
  end function add_ds_dsc

! Subtractions
  elemental type (ds_real) function sub_ds(a, b)
    type (ds_real), intent(in) :: a, b
    call f_ds_sub(a%re, b%re, sub_ds%re)
  end function sub_ds

  elemental type (ds_real) function sub_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ds_sub_ds_d(a%re, b, sub_ds_d%re)
  end function sub_ds_d

  elemental type (ds_real) function sub_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_sub_d_ds(a, b%re, sub_d_ds%re)
  end function sub_d_ds

  elemental type (ds_complex) function sub_dsc(a, b)
    type (ds_complex), intent(in) :: a, b
    call f_ds_sub (a%cmp(1:2), b%cmp(1:2), sub_dsc%cmp(1:2))
    call f_ds_sub (a%cmp(3:4), b%cmp(3:4), sub_dsc%cmp(3:4))
  end function sub_dsc

  elemental type (ds_complex) function sub_dsc_d(a, b)
    type (ds_complex), intent(in) :: a
    real*4, intent(in) :: b
    type (ds_real) dsb
    dsb%re(1) = b
    dsb%re(2) = 0.e0
    call f_ds_sub (a%cmp(1:2), dsb%re, sub_dsc_d%cmp(1:2))
    sub_dsc_d%cmp(3:4) = a%cmp(3:4)
  end function sub_dsc_d

  elemental type (ds_complex) function sub_d_dsc(a, b)
    real*4, intent(in) :: a
    type (ds_complex), intent(in) :: b
    type (ds_real) dsa
    dsa%re(1) = a
    dsa%re(2) = 0.e0
    call f_ds_sub (dsa%re, b%cmp(1:2), sub_d_dsc%cmp(1:2))
    sub_d_dsc%cmp(3:4) = - b%cmp(3:4)
  end function sub_d_dsc

  elemental type (ds_complex) function sub_dsc_ds(a, b)
    type (ds_complex), intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_sub (a%cmp(1:2), b%re, sub_dsc_ds%cmp(1:2))
    sub_dsc_ds%cmp(3:4) = a%cmp(3:4)
  end function sub_dsc_ds

  elemental type (ds_complex) function sub_ds_dsc(a, b)
    type (ds_real), intent(in) :: a
    type (ds_complex), intent(in) :: b
    call f_ds_sub (a%re, b%cmp(1:2), sub_ds_dsc%cmp(1:2))
    sub_ds_dsc%cmp(3:4) = - b%cmp(3:4)
  end function sub_ds_dsc

! Unary Minus
  elemental type (ds_real) function neg_ds(a)
    type (ds_real), intent(in) :: a
    neg_ds%re = -a%re
  end function neg_ds

  elemental type (ds_complex) function neg_dsc(a)
    type (ds_complex), intent(in) :: a
    neg_dsc%cmp = - a%cmp
  end function neg_dsc

! Multiplications
  elemental type (ds_real) function mul_ds(a, b)
    type (ds_real), intent(in) :: a, b
    call f_ds_mul(a%re, b%re, mul_ds%re)
  end function mul_ds

  elemental type (ds_real) function mul_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ds_mul_ds_d(a%re, b, mul_ds_d%re)
  end function mul_ds_d

  elemental type (ds_real) function mul_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_mul_ds_d(b%re, a, mul_d_ds%re)
  end function mul_d_ds

  elemental type (ds_real) function mul_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    call f_ds_mul(a%re, ds_int_limbs(b), mul_ds_i%re)
  end function mul_ds_i

  elemental type (ds_real) function mul_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_mul(b%re, ds_int_limbs(a), mul_i_ds%re)
  end function mul_i_ds

  elemental type (ds_complex) function mul_dsc(a, b)
    type (ds_complex), intent(in) :: a, b
    type (ds_real) t1, t2
    call f_ds_mul (a%cmp(1:2), b%cmp(1:2), t1%re)
    call f_ds_mul (a%cmp(3:4), b%cmp(3:4), t2%re)
    call f_ds_sub (t1%re, t2%re, mul_dsc%cmp(1:2))
    call f_ds_mul (a%cmp(1:2), b%cmp(3:4), t1%re)
    call f_ds_mul (a%cmp(3:4), b%cmp(1:2), t2%re)
    call f_ds_add (t1%re, t2%re, mul_dsc%cmp(3:4))
  end function mul_dsc

  elemental type (ds_complex) function mul_dsc_d(a, b)
    type (ds_complex), intent(in) :: a
    real*4, intent(in) :: b
    call f_ds_mul_ds_d (a%cmp(1:2), b, mul_dsc_d%cmp(1:2))
    call f_ds_mul_ds_d (a%cmp(3:4), b, mul_dsc_d%cmp(3:4))
  end function mul_dsc_d

  elemental type (ds_complex) function mul_d_dsc(a, b)
    real*4, intent(in) :: a
    type (ds_complex), intent(in) :: b
    call f_ds_mul_ds_d (b%cmp(1:2), a, mul_d_dsc%cmp(1:2))
    call f_ds_mul_ds_d (b%cmp(3:4), a, mul_d_dsc%cmp(3:4))
  end function mul_d_dsc

  elemental type (ds_complex) function mul_dsc_i(a, b)
    type (ds_complex), intent(in) :: a
    integer, intent(in) :: b
    call f_ds_mul(a%cmp(1:2), ds_int_limbs(b), mul_dsc_i%cmp(1:2))
    call f_ds_mul(a%cmp(3:4), ds_int_limbs(b), mul_dsc_i%cmp(3:4))
  end function mul_dsc_i

  elemental type (ds_complex) function mul_i_dsc(a, b)
    integer, intent(in) :: a
    type (ds_complex), intent(in) :: b
    call f_ds_mul(b%cmp(1:2), ds_int_limbs(a), mul_i_dsc%cmp(1:2))
    call f_ds_mul(b%cmp(3:4), ds_int_limbs(a), mul_i_dsc%cmp(3:4))
  end function mul_i_dsc

  elemental type (ds_complex) function mul_dsc_ds(a, b)
    type (ds_complex), intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_mul (a%cmp(1:2), b%re, mul_dsc_ds%cmp(1:2))
    call f_ds_mul (a%cmp(3:4), b%re, mul_dsc_ds%cmp(3:4))
  end function mul_dsc_ds

  elemental type (ds_complex) function mul_ds_dsc(a, b)
    type (ds_real), intent(in) :: a
    type (ds_complex), intent(in) :: b
    call f_ds_mul (a%re, b%cmp(1:2), mul_ds_dsc%cmp(1:2))
    call f_ds_mul (a%re, b%cmp(3:4), mul_ds_dsc%cmp(3:4))
  end function mul_ds_dsc

! Divisions
  elemental type (ds_real) function div_ds(a, b)
    type (ds_real), intent(in) :: a, b
    call f_ds_div(a%re, b%re, div_ds%re)
  end function div_ds

  elemental type (ds_real) function div_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_ds_div_ds_d(a%re, b, div_ds_d%re)
  end function div_ds_d

  elemental type (ds_real) function div_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_div_d_ds(a, b%re, div_d_ds%re)
  end function div_d_ds

  elemental type (ds_real) function div_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    call f_ds_div(a%re, ds_int_limbs(b), div_ds_i%re)
  end function div_ds_i

  elemental type (ds_real) function div_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_div(ds_int_limbs(a), b%re, div_i_ds%re)
  end function div_i_ds

  elemental type (ds_complex) function div_dsc(a, b)
    type (ds_complex), intent(in) :: a, b
    call f_ds_cdiv(a%cmp, b%cmp, div_dsc%cmp)
  end function div_dsc

  elemental type (ds_complex) function div_dsc_d(a,b)
    type (ds_complex), intent(in) :: a
    real*4, intent(in) :: b
    call f_ds_div_ds_d(a%cmp(1:2), b, div_dsc_d%cmp(1:2))
    call f_ds_div_ds_d(a%cmp(3:4), b, div_dsc_d%cmp(3:4))
  end function div_dsc_d

  elemental type (ds_complex) function div_dsc_ds(a, b)
    type (ds_complex), intent(in) :: a
    type (ds_real), intent(in) :: b
    call f_ds_div (a%cmp(1:2), b%re, div_dsc_ds%cmp(1:2))
    call f_ds_div (a%cmp(3:4), b%re, div_dsc_ds%cmp(3:4))
  end function div_dsc_ds

  elemental type (ds_complex) function div_ds_dsc(a, b)
    type (ds_real), intent(in) :: a
    type (ds_complex), intent(in) :: b
    type (ds_complex) :: ac
    ac%cmp(1:2) = a%re
    ac%cmp(3:4) = 0
    call f_ds_cdiv(ac%cmp, b%cmp, div_ds_dsc%cmp)
  end function div_ds_dsc

! Power
  elemental type (ds_real) function pwr_ds (a, b)
    type (ds_real), intent(in) :: a, b
    type (ds_real) q1, q2
    call f_ds_log(a%re, q1%re)
    call f_ds_mul(q1%re, b%re, q2%re)
    call f_ds_exp(q2%re, pwr_ds%re)
  end function pwr_ds

  elemental type (ds_real) function pwr_ds_i(a, n)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: n
    call f_ds_npwr(a%re, n, pwr_ds_i%re)
  end function pwr_ds_i

  elemental type (ds_real) function pwr_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    type (ds_real) q1, q2, q3
    q1%re(1) = a
    q1%re(2) = 0.e0
    call f_ds_log(q1%re, q2%re)
    call f_ds_mul(q2%re, b%re, q3%re)
    call f_ds_exp(q3%re, pwr_d_ds%re)
  end function pwr_d_ds

  elemental type (ds_complex) function pwr_dsc_i(a, n)
    type (ds_complex), intent(in) :: a
    integer, intent(in) :: n
    integer i2, j, n1
    type (ds_real) t1, t2, t3
    type (ds_complex) c1, c2

    intrinsic :: iabs, ishft

    if (n == 0) then
      if (all(a%cmp == 0.e0)) then
        !write (6, *) 'pwr_dsc_i: a = 0 and n = 0'
        call f_ds_nan(pwr_dsc_i%cmp(1:2))
        call f_ds_nan(pwr_dsc_i%cmp(3:4))
        return
      endif
      pwr_dsc_i%cmp(1) = 1.e0
      pwr_dsc_i%cmp(2:4) = 0.e0
      return
    endif
    n1 = iabs (n)
    i2 = ishft(1, n1-1)

    c1%cmp(1) = 1.e0
    c1%cmp(2:4) = 0.e0

110 continue

    if (n1 >= i2) then
      call f_ds_mul (a%cmp(1:2), c1%cmp(1:2), t1%re)
      call f_ds_mul (a%cmp(3:4), c1%cmp(3:4), t2%re)
      call f_ds_sub (t1%re, t2%re, c2%cmp(1:2))
      call f_ds_mul (a%cmp(1:2), c1%cmp(3:4), t1%re)
      call f_ds_mul (a%cmp(3:4), c1%cmp(1:2), t2%re)
      call f_ds_add (t1%re, t2%re, c2%cmp(3:4))
      do j = 1, 4
        c1%cmp(j) = c2%cmp(j)
      enddo
      n1 = n1 - i2
    endif
    i2 = i2 / 2
    if (i2 >= 1) then
      call f_ds_mul (c1%cmp(1:2), c1%cmp(1:2), t1%re)
      call f_ds_mul (c1%cmp(3:4), c1%cmp(3:4), t2%re)
      call f_ds_sub (t1%re, t2%re, c2%cmp(1:2))
      call f_ds_mul (c1%cmp(1:2), c1%cmp(3:4), t1%re)
      c2%cmp(3:4) = 2.e0 * t1%re
      c1%cmp = c2%cmp
      goto 110
    endif

    if (n > 0) then
      pwr_dsc_i%cmp = c1%cmp
    else
      c1%cmp(3:4) = - c1%cmp(3:4)
      call f_ds_mul (c1%cmp(1:2), c1%cmp(1:2), t1%re)
      call f_ds_mul (c1%cmp(3:4), c1%cmp(3:4), t2%re)
      call f_ds_add (t1%re, t2%re, t3%re)
      call f_ds_div (c1%cmp(1:2), t3%re, pwr_dsc_i%cmp(1:2))
      call f_ds_div (c1%cmp(3:4), t3%re, pwr_dsc_i%cmp(3:4))
    endif

    return
  end function pwr_dsc_i

! Trigonometric Functions
  elemental type (ds_real) function dssin(a)
    type (ds_real), intent(in) :: a
    call f_ds_sin(a%re, dssin%re)
  end function dssin

  elemental type (ds_real) function dscos(a)
    type (ds_real), intent(in) :: a
    call f_ds_cos(a%re, dscos%re)
  end function dscos

  elemental type (ds_real) function dstan(a)
    type (ds_real), intent(in) :: a
    call f_ds_tan(a%re, dstan%re)
  end function dstan

  elemental subroutine dssincos(a, s, c)
    type (ds_real), intent(in) :: a
    type (ds_real), intent(out) :: s, c
    call f_ds_sincos(a%re, s%re, c%re)
  end subroutine dssincos


! Inverse Trigonometric Functions
  elemental type (ds_real) function dsasin(a)
    type (ds_real), intent(in) :: a
    call f_ds_asin(a%re, dsasin%re)
  end function dsasin

  elemental type (ds_real) function dsacos(a)
    type (ds_real), intent(in) :: a
    call f_ds_acos(a%re, dsacos%re)
  end function dsacos

  elemental type (ds_real) function dsatan(a)
    type (ds_real), intent(in) :: a
    call f_ds_atan(a%re, dsatan%re)
  end function dsatan

  elemental type (ds_real) function dsatan2(a, b)
    type (ds_real), intent(in) :: a, b
    call f_ds_atan2(a%re, b%re, dsatan2%re)
  end function dsatan2

! Exponential and Logarithms
  elemental type (ds_real) function dsexp(a)
    type (ds_real), intent(in) :: a
    call f_ds_exp(a%re, dsexp%re)
  end function dsexp

  elemental type (ds_complex) function dscexp (a)
    type (ds_complex), intent(in) :: a
    type (ds_real) t1, t2, t3
    call f_ds_exp (a%cmp(1:2), t1%re)
    call f_ds_sincos (a%cmp(3:4), t3%re, t2%re)
    call f_ds_mul (t1%re, t2%re, dscexp%cmp(1:2))
    call f_ds_mul (t1%re, t3%re, dscexp%cmp(3:4))
  end function dscexp

  elemental type (ds_real) function dslog(a)
    type (ds_real), intent(in) :: a
    call f_ds_log(a%re, dslog%re)
  end function dslog

  elemental type (ds_complex) function dsclog (a)
    type (ds_complex), intent(in) :: a
    call f_ds_clog(a%cmp, dsclog%cmp)
  end function dsclog


  elemental type (ds_real) function dslog10(a)
    type (ds_real), intent(in) :: a
    call f_ds_log10(a%re, dslog10%re)
  end function dslog10

! SQRT, etc.
  elemental type (ds_real) function dssqrt(a)
    type (ds_real), intent(in) :: a
    call f_ds_sqrt(a%re, dssqrt%re)
  end function dssqrt

  elemental type (ds_real) function dssqr(a)
    type (ds_real), intent(in) :: a
    call f_ds_sqr(a%re, dssqr%re)
  end function dssqr

  elemental type (ds_real) function dsnroot(a, n)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: n
    call f_ds_nroot(a%re, n, dsnroot%re)
  end function dsnroot


! Hyperbolic Functions
  elemental type (ds_real) function dssinh(a)
    type (ds_real), intent(in) :: a
    call f_ds_sinh(a%re, dssinh%re)
  end function dssinh

  elemental type (ds_real) function dscosh(a)
    type (ds_real), intent(in) :: a
    call f_ds_cosh(a%re, dscosh%re)
  end function dscosh

  elemental type (ds_real) function dstanh(a)
    type (ds_real), intent(in) :: a
    call f_ds_tanh(a%re, dstanh%re)
  end function dstanh

  elemental subroutine dssincosh(a, s, c)
    type (ds_real), intent(in) :: a
    type (ds_real), intent(out) :: s, c
    call f_ds_sincosh(a%re, s%re, c%re)
  end subroutine dssincosh

! Inverse Hyperbolic Functions
  elemental type (ds_real) function dsasinh(a)
    type (ds_real), intent(in) :: a
    call f_ds_asinh(a%re, dsasinh%re)
  end function dsasinh

  elemental type (ds_real) function dsacosh(a)
    type (ds_real), intent(in) :: a
    call f_ds_acosh(a%re, dsacosh%re)
  end function dsacosh

  elemental type (ds_real) function dsatanh(a)
    type (ds_real), intent(in) :: a
    call f_ds_atanh(a%re, dsatanh%re)
  end function dsatanh


! Rounding
  elemental type (ds_real) function dsaint(a)
    type (ds_real), intent(in) :: a
    call f_ds_aint(a%re, dsaint%re)
  end function dsaint

  elemental type (ds_real) function dsanint(a)
    type (ds_real), intent(in) :: a
    call f_ds_nint(a%re, dsanint%re)
  end function dsanint

  elemental integer function dsnint(a)
    type (ds_real), intent(in) :: a
    dsnint = to_int_ds(dsanint(a))
  end function dsnint


! Random Number Generator
  subroutine dsrand(harvest)
    type (ds_real), intent(out) :: harvest
    call f_ds_rand(harvest%re)
  end subroutine dsrand


! Equality
  elemental logical function eq_ds(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    if (r == 0) then
      eq_ds = .true.
    else
      eq_ds = .false.
    end if
  end function eq_ds

  elemental logical function eq_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(a%re, b, r)
    if (r == 0) then
      eq_ds_d = .true.
    else
      eq_ds_d = .false.
    end if
  end function eq_ds_d

  elemental logical function eq_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(b%re, a, r)
    if (r == 0) then
      eq_d_ds = .true.
    else
      eq_d_ds = .false.
    end if
  end function eq_d_ds

  elemental logical function eq_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    eq_ds_i = eq_ds(a, ds_real(ds_int_limbs(b)))
  end function eq_ds_i

  elemental logical function eq_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    eq_i_ds = eq_ds(ds_real(ds_int_limbs(a)), b)
  end function eq_i_ds

  elemental logical function eq_dsc (a, b)
    type (ds_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_ds_comp (a%cmp(1:2), b%cmp(1:2), i1)
    call f_ds_comp (a%cmp(3:4), b%cmp(3:4), i2)
    if (i1 == 0 .and. i2 == 0) then
      eq_dsc = .true.
    else
      eq_dsc = .false.
    endif
  end function eq_dsc

  elemental logical function eq_dsc_ds (a, b)
    type (ds_complex), intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: i1
    call f_ds_comp (a%cmp(1:2), b%re, i1)
    if (i1 == 0 .and. all(a%cmp(3:4) == 0.e0)) then
      eq_dsc_ds = .true.
    else
      eq_dsc_ds = .false.
    endif
  end function eq_dsc_ds

  elemental logical function eq_ds_dsc (a, b)
    type (ds_real), intent(in) :: a
    type (ds_complex), intent(in) :: b
    integer :: i1
    call f_ds_comp (a%re, b%cmp(1:2), i1)
    if (i1 == 0 .and. all(b%cmp(3:4) == 0.e0)) then
      eq_ds_dsc = .true.
    else
      eq_ds_dsc = .false.
    endif
  end function eq_ds_dsc


! Non-Equality
  elemental logical function ne_ds(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    if (r == 0) then
      ne_ds = .false.
    else
      ne_ds = .true.
    end if
  end function ne_ds

  elemental logical function ne_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(a%re, b, r)
    if (r == 0) then
      ne_ds_d = .false.
    else
      ne_ds_d = .true.
    end if
  end function ne_ds_d

  elemental logical function ne_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(b%re, a, r)
    if (r == 0) then
      ne_d_ds = .false.
    else
      ne_d_ds = .true.
    end if
  end function ne_d_ds

  elemental logical function ne_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    ne_ds_i = ne_ds(a, ds_real(ds_int_limbs(b)))
  end function ne_ds_i

  elemental logical function ne_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    ne_i_ds = ne_ds(ds_real(ds_int_limbs(a)), b)
  end function ne_i_ds

  elemental logical function ne_dsc (a, b)
    type (ds_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_ds_comp (a%cmp(1:2), b%cmp(1:2), i1)
    call f_ds_comp (a%cmp(3:4), b%cmp(3:4), i2)
    if (i1 /= 0 .or. i2 /= 0) then
      ne_dsc = .true.
    else
      ne_dsc = .false.
    endif
  end function ne_dsc

  elemental logical function ne_dsc_ds (a, b)
    type (ds_complex), intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: i1
    call f_ds_comp (a%cmp(1:2), b%re, i1)
    if (i1 /= 0 .or. any(a%cmp(3:4) /= 0.e0)) then
      ne_dsc_ds = .true.
    else
      ne_dsc_ds = .false.
    endif
  end function ne_dsc_ds

  elemental logical function ne_ds_dsc (a, b)
    type (ds_real), intent(in) :: a
    type (ds_complex), intent(in) :: b
    integer :: i1
    call f_ds_comp (a%re, b%cmp(1:2), i1)
    if (i1 /= 0 .or. any(b%cmp(3:4) /= 0.e0)) then
      ne_ds_dsc = .true.
    else
      ne_ds_dsc = .false.
    endif
  end function ne_ds_dsc


! Greater-Than
  elemental logical function gt_ds(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    if (r == 1) then
      gt_ds = .true.
    else
      gt_ds = .false.
    end if
  end function gt_ds

  elemental logical function gt_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(a%re, b, r)
    if (r == 1) then
      gt_ds_d = .true.
    else
      gt_ds_d = .false.
    end if
  end function gt_ds_d

  elemental logical function gt_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(b%re, a, r)
    if (r == -1) then
      gt_d_ds = .true.
    else
      gt_d_ds = .false.
    end if
  end function gt_d_ds

  elemental logical function gt_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    gt_ds_i = gt_ds(a, ds_real(ds_int_limbs(b)))
  end function gt_ds_i

  elemental logical function gt_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    gt_i_ds = gt_ds(ds_real(ds_int_limbs(a)), b)
  end function gt_i_ds


! Less-Than
  elemental logical function lt_ds(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    if (r == -1) then
      lt_ds = .true.
    else
      lt_ds = .false.
    end if
  end function lt_ds

  elemental logical function lt_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(a%re, b, r)
    if (r == -1) then
      lt_ds_d = .true.
    else
      lt_ds_d = .false.
    end if
  end function lt_ds_d

  elemental logical function lt_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(b%re, a, r)
    if (r == 1) then
      lt_d_ds = .true.
    else
      lt_d_ds = .false.
    end if
  end function lt_d_ds

  elemental logical function lt_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    lt_ds_i = lt_ds(a, ds_real(ds_int_limbs(b)))
  end function lt_ds_i

  elemental logical function lt_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    lt_i_ds = lt_ds(ds_real(ds_int_limbs(a)), b)
  end function lt_i_ds

! Greater-Than-Or-Equal-To
  elemental logical function ge_ds(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    if (r == 0 .or. r == 1) then
      ge_ds = .true.
    else
      ge_ds = .false.
    end if
  end function ge_ds

  elemental logical function ge_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(a%re, b, r)
    if (r == 0 .or. r == 1) then
      ge_ds_d = .true.
    else
      ge_ds_d = .false.
    end if
  end function ge_ds_d

  elemental logical function ge_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(b%re, a, r)
    if (r == 0 .or. r == -1) then
      ge_d_ds = .true.
    else
      ge_d_ds = .false.
    end if
  end function ge_d_ds

  elemental logical function ge_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    ge_ds_i = ge_ds(a, ds_real(ds_int_limbs(b)))
  end function ge_ds_i

  elemental logical function ge_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    ge_i_ds = ge_ds(ds_real(ds_int_limbs(a)), b)
  end function ge_i_ds

! Less-Than-Or-Equal-To
  elemental logical function le_ds(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    if (r == 0 .or. r == -1) then
      le_ds = .true.
    else
      le_ds = .false.
    end if
  end function le_ds

  elemental logical function le_ds_d(a, b)
    type (ds_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(a%re, b, r)
    if (r == 0 .or. r == -1) then
      le_ds_d = .true.
    else
      le_ds_d = .false.
    end if
  end function le_ds_d

  elemental logical function le_d_ds(a, b)
    real*4, intent(in) :: a
    type (ds_real), intent(in) :: b
    integer :: r
    call f_ds_comp_ds_d(b%re, a, r)
    if (r == 0 .or. r == 1) then
      le_d_ds = .true.
    else
      le_d_ds = .false.
    end if
  end function le_d_ds

  elemental logical function le_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    le_ds_i = le_ds(a, ds_real(ds_int_limbs(b)))
  end function le_ds_i

  elemental logical function le_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    le_i_ds = le_ds(ds_real(ds_int_limbs(a)), b)
  end function le_i_ds

! Absolute Value
  elemental type (ds_real) function dsabs(a)
    type (ds_real), intent(in) :: a
    call f_ds_abs(a%re, dsabs%re)
  end function dsabs

  elemental type (ds_real) function dscabs (dsc)
    type (ds_complex), intent(in) :: dsc
    call f_ds_cabs(dsc%cmp, dscabs%re)
  end function dscabs

! Sign transfer
  elemental type (ds_real) function dssign(a, b) result (c)
    type (ds_real), intent(in) :: a, b
    if ((sign(1.0e0, a%re(1)) .lt. 0) .eqv. (sign(1.0e0, b%re(1)) .lt. 0)) then
      c%re = a%re
    else
      c%re = -a%re
    end if
  end function dssign

  elemental type (ds_real) function dssign_ds_d(a, b) result (c)
    type (ds_real), intent(in) :: a
    real*4, intent(in) ::  b
    if ((sign(1.0e0, a%re(1)) .lt. 0) .eqv. (sign(1.0e0, b) .lt. 0)) then
      c%re = a%re
    else
      c%re = -a%re
    end if
  end function dssign_ds_d

! Input
  subroutine dsinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ds_real), intent(in) :: q1
    type (ds_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9
!    CHARACTER (LEN=72) :: str

    call dsinp (u, q1%re)

    if (present(q2)) then
      call dsinp (u, q2%re)
    end if

    if (present(q3)) then
      call dsinp (u, q3%re)
    end if

    if (present(q4)) then
      call dsinp (u, q4%re)
    end if

    if (present(q5)) then
      call dsinp (u, q5%re)
    end if

    if (present(q6)) then
      call dsinp (u, q6%re)
    end if

    if (present(q7)) then
      call dsinp (u, q7%re)
    end if

    if (present(q8)) then
      call dsinp (u, q8%re)
    end if

    if (present(q9)) then
      call dsinp (u, q9%re)
   end if

  end subroutine dsinpq

  subroutine dscinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ds_complex), intent(in) :: q1
    type (ds_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call dsinp (u, q1%cmp(1:2))
    call dsinp (u, q1%cmp(3:4))

    if (present(q2)) then
      call dsinp (u, q2%cmp(1:2))
      call dsinp (u, q2%cmp(3:4))
    end if

    if (present(q3)) then
      call dsinp (u, q3%cmp(1:2))
      call dsinp (u, q3%cmp(3:4))
    end if

    if (present(q4)) then
      call dsinp (u, q4%cmp(1:2))
      call dsinp (u, q4%cmp(3:4))
    end if

    if (present(q5)) then
      call dsinp (u, q5%cmp(1:2))
      call dsinp (u, q5%cmp(3:4))
    end if

    if (present(q6)) then
      call dsinp (u, q6%cmp(1:2))
      call dsinp (u, q6%cmp(3:4))
    end if

    if (present(q7)) then
      call dsinp (u, q7%cmp(1:2))
      call dsinp (u, q7%cmp(3:4))
    end if

    if (present(q8)) then
      call dsinp (u, q8%cmp(1:2))
      call dsinp (u, q8%cmp(3:4))
    end if

    if (present(q9)) then
      call dsinp (u, q9%cmp(1:2))
      call dsinp (u, q9%cmp(3:4))
    end if

  end subroutine dscinpq


! Output
  subroutine dsoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ds_real), intent(in) :: q1
    type (ds_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9
!    CHARACTER (LEN=72) :: str

    call dsout (u, q1%re)

    if (present(q2)) then
      call dsout (u, q2%re)
    end if

    if (present(q3)) then
      call dsout (u, q3%re)
    end if

    if (present(q4)) then
      call dsout (u, q4%re)
    end if

    if (present(q5)) then
      call dsout (u, q5%re)
    end if

    if (present(q6)) then
      call dsout (u, q6%re)
    end if

    if (present(q7)) then
      call dsout (u, q7%re)
    end if

    if (present(q8)) then
      call dsout (u, q8%re)
    end if

    if (present(q9)) then
      call dsout (u, q9%re)
   end if

  end subroutine dsoutq

  subroutine dscoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (ds_complex), intent(in) :: q1
    type (ds_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call dsout (u, q1%cmp(1:2))
    call dsout (u, q1%cmp(3:4))

    if (present(q2)) then
      call dsout (u, q2%cmp(1:2))
      call dsout (u, q2%cmp(3:4))
    end if

    if (present(q3)) then
      call dsout (u, q3%cmp(1:2))
      call dsout (u, q3%cmp(3:4))
    end if

    if (present(q4)) then
      call dsout (u, q4%cmp(1:2))
      call dsout (u, q4%cmp(3:4))
    end if

    if (present(q5)) then
      call dsout (u, q5%cmp(1:2))
      call dsout (u, q5%cmp(3:4))
    end if

    if (present(q6)) then
      call dsout (u, q6%cmp(1:2))
      call dsout (u, q6%cmp(3:4))
    end if

    if (present(q7)) then
      call dsout (u, q7%cmp(1:2))
      call dsout (u, q7%cmp(3:4))
    end if

    if (present(q8)) then
      call dsout (u, q8%cmp(1:2))
      call dsout (u, q8%cmp(3:4))
    end if

    if (present(q9)) then
      call dsout (u, q9%cmp(1:2))
      call dsout (u, q9%cmp(3:4))
    end if

  end subroutine dscoutq

  elemental type (ds_real) function dsmin2(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    ! Unordered (NaN) operands keep the pre-1.6 result: treat as equal.
    if (r == 2) r = 0
    if (r == 1) then
      dsmin2 = b
    else
      dsmin2 = a
    end if
  end function dsmin2

  elemental type (ds_real) function dsmin(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (ds_real), intent(in) :: a1, a2, a3
    type (ds_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    dsmin = dsmin2(dsmin2(a1, a2), a3)
    if (present(a4)) dsmin = dsmin2(dsmin, a4)
    if (present(a5)) dsmin = dsmin2(dsmin, a5)
    if (present(a6)) dsmin = dsmin2(dsmin, a6)
    if (present(a7)) dsmin = dsmin2(dsmin, a7)
    if (present(a8)) dsmin = dsmin2(dsmin, a8)
    if (present(a9)) dsmin = dsmin2(dsmin, a9)
  end function dsmin

  elemental type (ds_real) function dsmax2(a, b)
    type (ds_real), intent(in) :: a, b
    integer :: r
    call f_ds_comp(a%re, b%re, r)
    ! Unordered (NaN) operands keep the pre-1.6 result: treat as equal.
    if (r == 2) r = 0
    if (r == -1) then
      dsmax2 = b
    else
      dsmax2 = a
    end if
  end function dsmax2

  elemental type (ds_real) function dsmax(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (ds_real), intent(in) :: a1, a2, a3
    type (ds_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    dsmax = dsmax2(dsmax2(a1, a2), a3)
    if (present(a4)) dsmax = dsmax2(dsmax, a4)
    if (present(a5)) dsmax = dsmax2(dsmax, a5)
    if (present(a6)) dsmax = dsmax2(dsmax, a6)
    if (present(a7)) dsmax = dsmax2(dsmax, a7)
    if (present(a8)) dsmax = dsmax2(dsmax, a8)
    if (present(a9)) dsmax = dsmax2(dsmax, a9)
  end function dsmax

  elemental type (ds_real) function dsmod (a, b)
    type (ds_real), intent(in) :: a, b
    type (ds_real) :: s1, s2
    call f_ds_div (a%re, b%re, s1%re)
    call f_ds_aint(s1%re, s2%re)
    call f_ds_mul (s2%re, b%re, s1%re)
    call f_ds_sub (a%re, s1%re, dsmod%re)
  end function dsmod

  pure type (ds_real) function ds_pi()
    call f_ds_pi(ds_pi%re)
  end function ds_pi

subroutine dsinp (iu, a)

!   This routine reads the DS number A from logical unit IU.  The input
!   value must be placed on a single line of not more than 80 characters.

implicit none
integer iu, ln
parameter (ln = 80)
character*80 cs
real*4 a(2)

read (iu, '(a)', end = 100) cs
call dsinpc (cs, a)
goto 110

100 write (6, 1)
1  format ('*** dsinp: End-of-file encountered.')
! call dsabrt
stop

110 return
end subroutine

subroutine dsinpc (a, b)

!   Converts the CHARACTER*80 array A into the DS number B using the C++
!   decimal reader, so Fortran and C++ parse literals identically.  Fortran
!   'd'/'D' exponent markers are accepted.

implicit none
character*80 a
real*4 b(2)
character*80 t
integer i, n, ierr

t = adjustl(a)
n = len_trim(t)
do i = 1, n
  if (t(i:i) == 'd' .or. t(i:i) == 'D') t(i:i) = 'e'
end do
ierr = -1
if (n > 0) call f_ds_read(t, n, b, ierr)
if (ierr /= 0) then
  write (6, 1) a
1 format ('*** dsinpc: Syntax error in literal string: ', a)
  stop
end if
end subroutine

subroutine dsout (iu, a)

!   This routine writes the DD number A on logical unit iu using a standard
!   E format, with lines 40 characters long.

implicit none
integer iu, ln
parameter (ln = 40)
character cs(40)
real*4 a(2)

call dsoutc (a, cs)
write (iu, '(40a)') cs

return
end subroutine

subroutine dsoutc (a, b)
  implicit none
  real*4 a(2)
  character b(40)

  b(1) = ' '
  b(2) = ' '
  call f_ds_swrite(a, 31, b(3), 38)
end subroutine

  character*16 function dsdigout (a, n)
    implicit none
    real*4 a, d1, d2
    character*16 ca, digits
    parameter (digits = '0123456789')
    integer i, is, k, n

    intrinsic :: abs, aint, sign

    ca = ' '
    is = sign (1.e0, a)
    d1 = abs (a)

    do i = n, 1, -1
      d2 = aint (d1 / 10.e0)
      k = 1.e0 + (d1 - 10.e0 * d2)
      d1 = d2
      ca(i:i) = digits(k:k)
      if (d1 == 0.e0) goto 100
    enddo

    i = 0

100 continue

    if (is < 0 .and. i > 1) then
      ca(i-1:i-1) = '-'
    elseif (i == 0 .or. is < 0 .and. i == 1) then
      ca = '****************'
    endif

    dsdigout = ca
    return
  end function

elemental type (ds_real) function dshuge(a)
  type (ds_real), intent(in) :: a
  dshuge = ds_huge
end function dshuge

elemental type (ds_real) function ds_safe_huge(a)
  type (ds_real), intent(in) :: a
  ds_safe_huge = ds_real((/3.4028235e+38, 0.0e0/))
end function ds_safe_huge

elemental type (ds_real) function dstiny(a)
  type (ds_real), intent(in) :: a
  dstiny = ds_tiny
end function dstiny

elemental type (ds_real) function dsepsilon(a)
  type (ds_real), intent(in) :: a
  dsepsilon = ds_eps
end function dsepsilon

elemental integer function ds_radix(a)
  type (ds_real), intent(in) :: a
  ds_radix = 2
end function ds_radix

elemental integer function ds_digits(a)
  type (ds_real), intent(in) :: a
  ds_digits = 46
end function ds_digits

elemental integer function ds_max_expn(a)
  type (ds_real), intent(in) :: a
  ds_max_expn = 127
end function ds_max_expn

elemental integer function ds_min_expn(a)
  type (ds_real), intent(in) :: a
  ds_min_expn = -102
end function ds_min_expn

elemental integer function ds_precision(a)
  type (ds_real), intent(in) :: a
  ds_precision = 13
end function ds_precision

elemental integer function ds_range(a)
  type (ds_real), intent(in) :: a
  ds_range = 37
end function ds_range

elemental type (ds_real) function ds_nan(a)
  type (ds_real), intent(in) :: a
  call f_ds_nan(ds_nan%re)
end function ds_nan

elemental type (ds_real) function ds_aimag(a)
  type (ds_complex), intent(in) :: a
  ds_aimag%re = a%cmp(3:4)
end function

! Exact conversion of a default integer to 2 binary32 limbs: the leading
! limb is the rounded value and the remainder (at most 2**7) is exact.
pure function ds_int_limbs(i) result(r)
  implicit none
  integer, intent(in) :: i
  real*4 :: r(2)
  r = 0.0e0
  r(1) = real(i, kind=4)
  r(2) = real(int(i, kind=8) - int(r(1), kind=8), kind=4)
end function ds_int_limbs

! Mixed-mode operators added for interface completeness.
  elemental type (ds_real) function sub_ds_i(a, b)
    type (ds_real), intent(in) :: a
    integer, intent(in) :: b
    sub_ds_i = sub_ds(a, to_ds_i(b))
  end function sub_ds_i

  elemental type (ds_real) function sub_i_ds(a, b)
    integer, intent(in) :: a
    type (ds_real), intent(in) :: b
    sub_i_ds = sub_ds(to_ds_i(a), b)
  end function sub_i_ds

  elemental type (ds_complex) function add_dsc_i(a, b)
    type (ds_complex), intent(in) :: a
    integer, intent(in) :: b
    add_dsc_i = add_dsc_ds(a, to_ds_i(b))
  end function add_dsc_i

  elemental type (ds_complex) function add_i_dsc(a, b)
    integer, intent(in) :: a
    type (ds_complex), intent(in) :: b
    add_i_dsc = add_ds_dsc(to_ds_i(a), b)
  end function add_i_dsc

  elemental type (ds_complex) function sub_dsc_i(a, b)
    type (ds_complex), intent(in) :: a
    integer, intent(in) :: b
    sub_dsc_i = sub_dsc_ds(a, to_ds_i(b))
  end function sub_dsc_i

  elemental type (ds_complex) function sub_i_dsc(a, b)
    integer, intent(in) :: a
    type (ds_complex), intent(in) :: b
    sub_i_dsc = sub_ds_dsc(to_ds_i(a), b)
  end function sub_i_dsc

  elemental type (ds_complex) function div_dsc_i(a, b)
    type (ds_complex), intent(in) :: a
    integer, intent(in) :: b
    div_dsc_i = div_dsc_ds(a, to_ds_i(b))
  end function div_dsc_i

  elemental type (ds_complex) function div_i_dsc(a, b)
    integer, intent(in) :: a
    type (ds_complex), intent(in) :: b
    div_i_dsc = div_ds_dsc(to_ds_i(a), b)
  end function div_i_dsc

  elemental type (ds_complex) function div_d_dsc(a, b)
    real*4, intent(in) :: a
    type (ds_complex), intent(in) :: b
    div_d_dsc = div_ds_dsc(to_ds_d(a), b)
  end function div_d_dsc

  elemental logical function eq_dsc_d(a, b)
    type (ds_complex), intent(in) :: a
    real*4, intent(in) :: b
    eq_dsc_d = eq_dsc_ds(a, to_ds_d(b))
  end function eq_dsc_d

  elemental logical function eq_d_dsc(a, b)
    real*4, intent(in) :: a
    type (ds_complex), intent(in) :: b
    eq_d_dsc = eq_ds_dsc(to_ds_d(a), b)
  end function eq_d_dsc

  elemental logical function ne_dsc_d(a, b)
    type (ds_complex), intent(in) :: a
    real*4, intent(in) :: b
    ne_dsc_d = ne_dsc_ds(a, to_ds_d(b))
  end function ne_dsc_d

  elemental logical function ne_d_dsc(a, b)
    real*4, intent(in) :: a
    type (ds_complex), intent(in) :: b
    ne_d_dsc = ne_ds_dsc(to_ds_d(a), b)
  end function ne_d_dsc

! binary64 interoperability for the binary32-based type.
  elemental type (ds_real) function to_ds_dp(a)
    real*8, intent(in) :: a
    call f_ds_from_double(a, to_ds_dp%re)
  end function to_ds_dp

  elemental real*8 function to_dp_ds(a)
    type (ds_real), intent(in) :: a
    call f_ds_to_double(a%re, to_dp_ds)
  end function to_dp_ds

  elemental real*8 function to_dp_dsc(a)
    type (ds_complex), intent(in) :: a
    to_dp_dsc = to_dp_ds(to_ds_dsc(a))
  end function to_dp_dsc

  elemental type (ds_complex) function to_dsc_zc(a)
    complex*16, intent(in) :: a
    to_dsc_zc = to_dsc_ds2(to_ds_dp(dble(a)), to_ds_dp(aimag(a)))
  end function to_dsc_zc

  elemental subroutine assign_ds_dp(a, b)
    type (ds_real), intent(inout) :: a
    real*8, intent(in) :: b
    a = to_ds_dp(b)
  end subroutine assign_ds_dp

  elemental subroutine assign_dp_ds(a, b)
    real*8, intent(inout) :: a
    type (ds_real), intent(in) :: b
    a = to_dp_ds(b)
  end subroutine assign_dp_ds

  elemental subroutine assign_dsc_zc(a, b)
    type (ds_complex), intent(inout) :: a
    complex*16, intent(in) :: b
    a = to_dsc_zc(b)
  end subroutine assign_dsc_zc

  elemental type (ds_real) function add_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    add_ds_dp = add_ds(a, to_ds_dp(b))
  end function add_ds_dp

  elemental type (ds_real) function add_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    add_dp_ds = add_ds(to_ds_dp(a), b)
  end function add_dp_ds

  elemental type (ds_complex) function add_dsc_dp(a, b)
    type (ds_complex), intent(in) :: a
    real*8, intent(in) :: b
    add_dsc_dp = add_dsc_ds(a, to_ds_dp(b))
  end function add_dsc_dp

  elemental type (ds_complex) function add_dp_dsc(a, b)
    real*8, intent(in) :: a
    type (ds_complex), intent(in) :: b
    add_dp_dsc = add_ds_dsc(to_ds_dp(a), b)
  end function add_dp_dsc

  elemental type (ds_real) function sub_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    sub_ds_dp = sub_ds(a, to_ds_dp(b))
  end function sub_ds_dp

  elemental type (ds_real) function sub_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    sub_dp_ds = sub_ds(to_ds_dp(a), b)
  end function sub_dp_ds

  elemental type (ds_complex) function sub_dsc_dp(a, b)
    type (ds_complex), intent(in) :: a
    real*8, intent(in) :: b
    sub_dsc_dp = sub_dsc_ds(a, to_ds_dp(b))
  end function sub_dsc_dp

  elemental type (ds_complex) function sub_dp_dsc(a, b)
    real*8, intent(in) :: a
    type (ds_complex), intent(in) :: b
    sub_dp_dsc = sub_ds_dsc(to_ds_dp(a), b)
  end function sub_dp_dsc

  elemental type (ds_real) function mul_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    mul_ds_dp = mul_ds(a, to_ds_dp(b))
  end function mul_ds_dp

  elemental type (ds_real) function mul_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    mul_dp_ds = mul_ds(to_ds_dp(a), b)
  end function mul_dp_ds

  elemental type (ds_complex) function mul_dsc_dp(a, b)
    type (ds_complex), intent(in) :: a
    real*8, intent(in) :: b
    mul_dsc_dp = mul_dsc_ds(a, to_ds_dp(b))
  end function mul_dsc_dp

  elemental type (ds_complex) function mul_dp_dsc(a, b)
    real*8, intent(in) :: a
    type (ds_complex), intent(in) :: b
    mul_dp_dsc = mul_ds_dsc(to_ds_dp(a), b)
  end function mul_dp_dsc

  elemental type (ds_real) function div_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    div_ds_dp = div_ds(a, to_ds_dp(b))
  end function div_ds_dp

  elemental type (ds_real) function div_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    div_dp_ds = div_ds(to_ds_dp(a), b)
  end function div_dp_ds

  elemental type (ds_complex) function div_dsc_dp(a, b)
    type (ds_complex), intent(in) :: a
    real*8, intent(in) :: b
    div_dsc_dp = div_dsc_ds(a, to_ds_dp(b))
  end function div_dsc_dp

  elemental type (ds_complex) function div_dp_dsc(a, b)
    real*8, intent(in) :: a
    type (ds_complex), intent(in) :: b
    div_dp_dsc = div_ds_dsc(to_ds_dp(a), b)
  end function div_dp_dsc

  elemental logical function eq_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    eq_ds_dp = eq_ds(a, to_ds_dp(b))
  end function eq_ds_dp

  elemental logical function eq_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    eq_dp_ds = eq_ds(to_ds_dp(a), b)
  end function eq_dp_ds

  elemental logical function ne_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    ne_ds_dp = ne_ds(a, to_ds_dp(b))
  end function ne_ds_dp

  elemental logical function ne_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    ne_dp_ds = ne_ds(to_ds_dp(a), b)
  end function ne_dp_ds

  elemental logical function gt_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    gt_ds_dp = gt_ds(a, to_ds_dp(b))
  end function gt_ds_dp

  elemental logical function gt_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    gt_dp_ds = gt_ds(to_ds_dp(a), b)
  end function gt_dp_ds

  elemental logical function lt_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    lt_ds_dp = lt_ds(a, to_ds_dp(b))
  end function lt_ds_dp

  elemental logical function lt_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    lt_dp_ds = lt_ds(to_ds_dp(a), b)
  end function lt_dp_ds

  elemental logical function ge_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    ge_ds_dp = ge_ds(a, to_ds_dp(b))
  end function ge_ds_dp

  elemental logical function ge_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    ge_dp_ds = ge_ds(to_ds_dp(a), b)
  end function ge_dp_ds

  elemental logical function le_ds_dp(a, b)
    type (ds_real), intent(in) :: a
    real*8, intent(in) :: b
    le_ds_dp = le_ds(a, to_ds_dp(b))
  end function le_ds_dp

  elemental logical function le_dp_ds(a, b)
    real*8, intent(in) :: a
    type (ds_real), intent(in) :: b
    le_dp_ds = le_ds(to_ds_dp(a), b)
  end function le_dp_ds

  elemental logical function eq_dsc_dp(a, b)
    type (ds_complex), intent(in) :: a
    real*8, intent(in) :: b
    eq_dsc_dp = eq_dsc_ds(a, to_ds_dp(b))
  end function eq_dsc_dp

  elemental logical function eq_dp_dsc(a, b)
    real*8, intent(in) :: a
    type (ds_complex), intent(in) :: b
    eq_dp_dsc = eq_ds_dsc(to_ds_dp(a), b)
  end function eq_dp_dsc

  elemental logical function ne_dsc_dp(a, b)
    type (ds_complex), intent(in) :: a
    real*8, intent(in) :: b
    ne_dsc_dp = ne_dsc_ds(a, to_ds_dp(b))
  end function ne_dsc_dp

  elemental logical function ne_dp_dsc(a, b)
    real*8, intent(in) :: a
    type (ds_complex), intent(in) :: b
    ne_dp_dsc = ne_ds_dsc(to_ds_dp(a), b)
  end function ne_dp_dsc

! Complex elementary functions, evaluated by the C++ complex implementation.
  elemental type (ds_complex) function dsc_sqrt(z)
    type (ds_complex), intent(in) :: z
    call f_ds_csqrt(z%cmp, dsc_sqrt%cmp)
  end function dsc_sqrt

  elemental type (ds_complex) function dsc_sin(z)
    type (ds_complex), intent(in) :: z
    call f_ds_csin(z%cmp, dsc_sin%cmp)
  end function dsc_sin

  elemental type (ds_complex) function dsc_cos(z)
    type (ds_complex), intent(in) :: z
    call f_ds_ccos(z%cmp, dsc_cos%cmp)
  end function dsc_cos

  elemental type (ds_complex) function dsc_tan(z)
    type (ds_complex), intent(in) :: z
    call f_ds_ctan(z%cmp, dsc_tan%cmp)
  end function dsc_tan

  elemental type (ds_complex) function dsc_sinh(z)
    type (ds_complex), intent(in) :: z
    call f_ds_csinh(z%cmp, dsc_sinh%cmp)
  end function dsc_sinh

  elemental type (ds_complex) function dsc_cosh(z)
    type (ds_complex), intent(in) :: z
    call f_ds_ccosh(z%cmp, dsc_cosh%cmp)
  end function dsc_cosh

  elemental type (ds_complex) function dsc_tanh(z)
    type (ds_complex), intent(in) :: z
    call f_ds_ctanh(z%cmp, dsc_tanh%cmp)
  end function dsc_tanh

  elemental type (ds_complex) function dsc_asin(z)
    type (ds_complex), intent(in) :: z
    call f_ds_casin(z%cmp, dsc_asin%cmp)
  end function dsc_asin

  elemental type (ds_complex) function dsc_acos(z)
    type (ds_complex), intent(in) :: z
    call f_ds_cacos(z%cmp, dsc_acos%cmp)
  end function dsc_acos

  elemental type (ds_complex) function dsc_atan(z)
    type (ds_complex), intent(in) :: z
    call f_ds_catan(z%cmp, dsc_atan%cmp)
  end function dsc_atan

  elemental type (ds_complex) function dsc_asinh(z)
    type (ds_complex), intent(in) :: z
    call f_ds_casinh(z%cmp, dsc_asinh%cmp)
  end function dsc_asinh

  elemental type (ds_complex) function dsc_acosh(z)
    type (ds_complex), intent(in) :: z
    call f_ds_cacosh(z%cmp, dsc_acosh%cmp)
  end function dsc_acosh

  elemental type (ds_complex) function dsc_atanh(z)
    type (ds_complex), intent(in) :: z
    call f_ds_catanh(z%cmp, dsc_atanh%cmp)
  end function dsc_atanh

  elemental type (ds_complex) function pwr_dsc(z, w)
    type (ds_complex), intent(in) :: z, w
    call f_ds_cpow(z%cmp, w%cmp, pwr_dsc%cmp)
  end function pwr_dsc

  elemental type (ds_complex) function pwr_dsc_ds(z, x)
    type (ds_complex), intent(in) :: z
    type (ds_real), intent(in) :: x
    call f_ds_cpow_r(z%cmp, x%re, pwr_dsc_ds%cmp)
  end function pwr_dsc_ds

! Rounding and model intrinsics added for interface completeness.
  elemental type (ds_real) function dsfloor(a)
    type (ds_real), intent(in) :: a
    call f_ds_floor(a%re, dsfloor%re)
  end function dsfloor

  elemental type (ds_real) function dsceil(a)
    type (ds_real), intent(in) :: a
    call f_ds_ceil(a%re, dsceil%re)
  end function dsceil

  ! sqrt(a**2 + b**2) without intermediate overflow or underflow.
  elemental type (ds_real) function dshypot(a, b)
    type (ds_real), intent(in) :: a, b
    type (ds_real) :: x, y, m
    x = abs(a)
    y = abs(b)
    m = max(x, y)
    if (m == 0) then
      dshypot = m
    else
      dshypot = m * sqrt((x / m) ** 2 + (y / m) ** 2)
    end if
  end function dshypot

  ! a - floor(a / p) * p: the result has the sign of p.
  elemental type (ds_real) function dsmodulo(a, p)
    type (ds_real), intent(in) :: a, p
    dsmodulo = a - dsfloor(a / p) * p
  end function dsmodulo

  ! max(a - b, 0)
  elemental type (ds_real) function dsdim(a, b)
    type (ds_real), intent(in) :: a, b
    if (a > b) then
      dsdim = a - b
    else
      dsdim = a - a
    end if
  end function dsdim

! Conversions to and from the binary64-based dd/qd types.
  elemental type (ds_real) function to_ds_dd(dd)
    type (dd_real), intent(in) :: dd
    call f_ds_from_dd(dd%re, to_ds_dd%re)
  end function to_ds_dd

  elemental type (dd_real) function to_dd_ds(ds)
    type (ds_real), intent(in) :: ds
    call f_ds_to_dd(ds%re, to_dd_ds%re)
  end function to_dd_ds

  elemental subroutine assign_ds_dd(ds, dd)
    type (ds_real), intent(inout) :: ds
    type (dd_real), intent(in) :: dd
    call f_ds_from_dd(dd%re, ds%re)
  end subroutine assign_ds_dd

  elemental subroutine assign_dd_ds(dd, ds)
    type (dd_real), intent(inout) :: dd
    type (ds_real), intent(in) :: ds
    call f_ds_to_dd(ds%re, dd%re)
  end subroutine assign_dd_ds

  elemental type (ds_complex) function to_dsc_ddc(ddc)
    type (dd_complex), intent(in) :: ddc
    call f_ds_from_dd(ddc%cmp(1:2), to_dsc_ddc%cmp(1:2))
    call f_ds_from_dd(ddc%cmp(3:4), to_dsc_ddc%cmp(3:4))
  end function to_dsc_ddc

  elemental type (dd_complex) function to_ddc_dsc(dsc)
    type (ds_complex), intent(in) :: dsc
    call f_ds_to_dd(dsc%cmp(1:2), to_ddc_dsc%cmp(1:2))
    call f_ds_to_dd(dsc%cmp(3:4), to_ddc_dsc%cmp(3:4))
  end function to_ddc_dsc

  elemental subroutine assign_dsc_ddc(dsc, ddc)
    type (ds_complex), intent(inout) :: dsc
    type (dd_complex), intent(in) :: ddc
    dsc = to_dsc_ddc(ddc)
  end subroutine assign_dsc_ddc

  elemental subroutine assign_ddc_dsc(ddc, dsc)
    type (dd_complex), intent(inout) :: ddc
    type (ds_complex), intent(in) :: dsc
    ddc = to_ddc_dsc(dsc)
  end subroutine assign_ddc_dsc

  elemental type (ds_real) function to_ds_qd(qd)
    type (qd_real), intent(in) :: qd
    call f_ds_from_qd(qd%re, to_ds_qd%re)
  end function to_ds_qd

  elemental type (qd_real) function to_qd_ds(ds)
    type (ds_real), intent(in) :: ds
    call f_ds_to_qd(ds%re, to_qd_ds%re)
  end function to_qd_ds

  elemental subroutine assign_ds_qd(ds, qd)
    type (ds_real), intent(inout) :: ds
    type (qd_real), intent(in) :: qd
    call f_ds_from_qd(qd%re, ds%re)
  end subroutine assign_ds_qd

  elemental subroutine assign_qd_ds(qd, ds)
    type (qd_real), intent(inout) :: qd
    type (ds_real), intent(in) :: ds
    call f_ds_to_qd(ds%re, qd%re)
  end subroutine assign_qd_ds

  elemental type (ds_complex) function to_dsc_qdc(qdc)
    type (qd_complex), intent(in) :: qdc
    call f_ds_from_qd(qdc%cmp(1:4), to_dsc_qdc%cmp(1:2))
    call f_ds_from_qd(qdc%cmp(5:8), to_dsc_qdc%cmp(3:4))
  end function to_dsc_qdc

  elemental type (qd_complex) function to_qdc_dsc(dsc)
    type (ds_complex), intent(in) :: dsc
    call f_ds_to_qd(dsc%cmp(1:2), to_qdc_dsc%cmp(1:4))
    call f_ds_to_qd(dsc%cmp(3:4), to_qdc_dsc%cmp(5:8))
  end function to_qdc_dsc

  elemental subroutine assign_dsc_qdc(dsc, qdc)
    type (ds_complex), intent(inout) :: dsc
    type (qd_complex), intent(in) :: qdc
    dsc = to_dsc_qdc(qdc)
  end subroutine assign_dsc_qdc

  elemental subroutine assign_qdc_dsc(qdc, dsc)
    type (qd_complex), intent(inout) :: qdc
    type (ds_complex), intent(in) :: dsc
    qdc = to_qdc_dsc(dsc)
  end subroutine assign_qdc_dsc

end module dsmodule
