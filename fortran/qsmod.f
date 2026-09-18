!  qsmod.f
!
!  This work was supported by the Director, Office of Science, Division
!  of Mathematical, Information, and Computational Sciences of the
!  U.S. Department of Energy under contract number DE-AC03-76SF00098.
!
!  Copyright (c) 2000-2008
!
!  Fortran-90 module file to use with quad-single numbers.
!
!  Yozo Hida
!  David H Bailey    2008-02-20

module qsmodule
  use ddmodule
  use qsext
  implicit none

  type qs_real
    sequence
    real*4 :: re(4)
  end type qs_real

  type qs_complex
    sequence
    real*4 :: cmp(8)
  end type qs_complex

  real*4 d_qs_eps
  parameter (d_qs_eps = 6.3108872e-29)

  type (qs_real) qs_one, qs_zero, qs_eps, qs_huge, qs_tiny
  parameter (qs_one = qs_real((/1.0e0, 0.0e0, 0.0e0, 0.0e0/)))
  parameter (qs_zero = qs_real((/0.0e0, 0.0e0, 0.0e0, 0.0e0/)))
  parameter (qs_eps = qs_real((/d_qs_eps, 0.0e0, 0.0e0, 0.0e0/)))
  parameter (qs_huge = qs_real((/3.4028235e+38, 0.0e0, 0.0e0, 0.0e0/)))
  parameter (qs_tiny = qs_real((/5.5511151e-17, 0.0e0, 0.0e0, 0.0e0/)))

  interface assignment (=)
    module procedure assign_qs_str
    module procedure assign_qs
    module procedure assign_qs_d
    module procedure assign_d_qs
    module procedure assign_dd_qs
    module procedure assign_qs_dd
    module procedure assign_qs_i
    module procedure assign_i_qs
    module procedure assign_qsc
    module procedure assign_qsc_qs
    module procedure assign_qs_qsc
    module procedure assign_qsc_d
    module procedure assign_qsc_i
    module procedure assign_d_qsc
    module procedure assign_qsc_dc
    module procedure assign_dc_qsc
  end interface

  interface operator (+)
    module procedure add_qs
    module procedure add_qs_d
    module procedure add_d_qs
    module procedure add_qs_i
    module procedure add_i_qs
    module procedure add_qsc
    module procedure add_qsc_qs
    module procedure add_qs_qsc
    module procedure add_qsc_d
    module procedure add_d_qsc
  end interface

  interface operator (-)
    module procedure sub_qs
    module procedure sub_qs_d
    module procedure sub_d_qs
    module procedure neg_qs
    module procedure sub_qsc
    module procedure sub_qsc_qs
    module procedure sub_qs_qsc
    module procedure sub_qsc_d
    module procedure sub_d_qsc
    module procedure neg_qsc
  end interface

  interface operator (*)
    module procedure mul_qs
    module procedure mul_qs_d
    module procedure mul_d_qs
    module procedure mul_qs_i
    module procedure mul_i_qs
    module procedure mul_qsc
    module procedure mul_qsc_qs
    module procedure mul_qs_qsc
    module procedure mul_qsc_d
    module procedure mul_d_qsc
    module procedure mul_i_qsc
    module procedure mul_qsc_i
  end interface

  interface operator (/)
    module procedure div_qs
    module procedure div_qs_d
    module procedure div_d_qs
    module procedure div_qs_i
    module procedure div_i_qs
    module procedure div_qsc
    module procedure div_qsc_qs
    module procedure div_qs_qsc
    module procedure div_qsc_d
  end interface

  interface operator (**)
    module procedure pwr_qs
    module procedure pwr_qs_i
    module procedure pwr_d_qs
    module procedure pwr_qsc_i
  end interface

  interface qsreal
    module procedure to_qs_i
    module procedure to_qs_d
    module procedure to_qs_dd
    module procedure to_qs_qs
    module procedure to_qs_str
    module procedure to_qs_qsc
  end interface

  interface ddreal
    module procedure to_dd_qs
  end interface

  interface real
    module procedure to_d_qs
    module procedure to_qs_qsc
  end interface

  interface qscomplex
    module procedure to_qsc_qs
    module procedure to_qsc_qs2
    module procedure to_qsc_d
    module procedure to_qsc_dc
  end interface

  interface int
    module procedure to_int_qs
  end interface

  interface sin
    module procedure qssin
  end interface
  interface cos
    module procedure qscos
  end interface
  interface tan
    module procedure qstan
  end interface
  interface sincos
    module procedure qssincos
  end interface

  interface asin
    module procedure qsasin
  end interface
  interface acos
    module procedure qsacos
  end interface
  interface atan
    module procedure qsatan
  end interface
  interface atan2
    module procedure qsatan2
  end interface

  interface exp
    module procedure qsexp
    module procedure qscexp
  end interface
  interface log
    module procedure qslog
    module procedure qsclog
  end interface
  interface log10
    module procedure qslog10
  end interface

  interface sqrt
    module procedure qssqrt
  end interface
  interface sqr
    module procedure qssqr
  end interface
  interface nroot
    module procedure qsnroot
  end interface

  interface sinh
    module procedure qssinh
  end interface
  interface cosh
    module procedure qscosh
  end interface
  interface tanh
    module procedure qstanh
  end interface
  interface sincosh
    module procedure qssincosh
  end interface

  interface asinh
    module procedure qsasinh
  end interface
  interface acosh
    module procedure qsacosh
  end interface
  interface atanh
    module procedure qsatanh
  end interface

  interface aint
    module procedure qsaint
  end interface

  interface nint
    module procedure qsnint
  end interface

  interface anint
    module procedure qsanint
  end interface

  interface abs
    module procedure qsabs
    module procedure qscabs
  end interface

  interface sign
    module procedure qssign
    module procedure qssign_dd_d
  end interface

  interface random_number
    module procedure qsrand
  end interface

  interface aimag
    module procedure qs_aimag
  end interface

  interface operator (==)
    module procedure eq_qs
    module procedure eq_qs_d
    module procedure eq_d_qs
    module procedure eq_qs_i
    module procedure eq_i_qs
    module procedure eq_qsc
    module procedure eq_qsc_qs
    module procedure eq_qs_qsc
  end interface

  interface operator (/=)
    module procedure ne_qs
    module procedure ne_qs_d
    module procedure ne_d_qs
    module procedure ne_qs_i
    module procedure ne_i_qs
    module procedure ne_qsc
    module procedure ne_qsc_qs
    module procedure ne_qs_qsc
  end interface

  interface operator (>)
    module procedure gt_qs
    module procedure gt_qs_d
    module procedure gt_d_qs
    module procedure gt_qs_i
    module procedure gt_i_qs
  end interface

  interface operator (<)
    module procedure lt_qs
    module procedure lt_qs_d
    module procedure lt_d_qs
    module procedure lt_qs_i
    module procedure lt_i_qs
  end interface

  interface operator (>=)
    module procedure ge_qs
    module procedure ge_qs_d
    module procedure ge_d_qs
    module procedure ge_qs_i
    module procedure ge_i_qs
  end interface

  interface operator (<=)
    module procedure le_qs
    module procedure le_qs_d
    module procedure le_d_qs
    module procedure le_qs_i
    module procedure le_i_qs
  end interface

  interface read_scalar
    module procedure qsinpq
    module procedure qscinpq
  end interface

  interface write_scalar
    module procedure qsoutq
    module procedure qscoutq
  end interface

  interface qsread
    module procedure qsinpq
  end interface

  interface qswrite
    module procedure qsoutq
  end interface

  interface qscread
    module procedure qscinpq
  end interface

  interface qscwrite
    module procedure qscoutq
  end interface

  interface dble
    module procedure to_d_qs
    module procedure to_d_qsc
  end interface

  interface cmplx
    module procedure to_dc_qsc
  end interface

  interface conjg
    module procedure qscconjg
  end interface

  interface min
    module procedure qsmin
    module procedure qsmin2
  end interface
  interface max
    module procedure qsmax
    module procedure qsmax2
  end interface
  interface mod
     module procedure qsmod
  end interface

  interface qspi
    module procedure qs_pi
  end interface

  interface huge
    module procedure qshuge
  end interface

  interface safe_huge
    module procedure qs_safe_huge
  end interface

  interface tiny
    module procedure qstiny
  end interface

  interface epsilon
    module procedure qsepsilon
  end interface

  interface radix
    module procedure qs_radix
  end interface

  interface digits
    module procedure qs_digits
  end interface

  interface maxexponent
    module procedure qs_max_expn
  end interface

  interface minexponent
    module procedure qs_min_expn
  end interface

  interface precision
    module procedure qs_precision
  end interface

  interface range
    module procedure qs_range
  end interface

  interface nan
    module procedure qs_nan
  end interface

contains

! Assignments
  subroutine assign_qs_str(a, s)
    type (qs_real), intent(inout) :: a
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call qsinpc (t, a%re)
  end subroutine assign_qs_str

  elemental subroutine assign_qs (a, b)
    type (qs_real), intent(inout) :: a
    type (qs_real), intent(in) :: b
    a%re = b%re
  end subroutine assign_qs

  elemental subroutine assign_qs_d(a, d)
    type (qs_real), intent(inout) :: a
    real*4, intent(in) :: d
    a%re(1) = d
    a%re(2:4) = 0.0e0
  end subroutine assign_qs_d

  elemental subroutine assign_d_qs(d, a)
    real*4, intent(inout) :: d
    type (qs_real), intent(in) :: a
    d = a%re(1)
  end subroutine assign_d_qs

  elemental subroutine assign_qs_i(a, i)
    type (qs_real), intent(inout) :: a
    integer, intent(in) :: i
    a%re(1) = i
    a%re(2:4) = 0.0e0
  end subroutine assign_qs_i

  elemental subroutine assign_i_qs(i, a)
    integer, intent(inout) :: i
    type (qs_real), intent(in) :: a
    i = a%re(1)
  end subroutine assign_i_qs

  elemental subroutine assign_dd_qs(dd, qs)
    type (dd_real), intent(inout) :: dd
    type (qs_real), intent(in) :: qs
    dd%re(1:2) = qs%re(1:2)
  end subroutine assign_dd_qs

  elemental subroutine assign_qs_dd(qs, dd)
    type (qs_real), intent(inout) :: qs
    type (dd_real), intent(in) :: dd
    qs%re(1:2) = dd%re
    qs%re(3:4) = 0.e0
  end subroutine assign_qs_dd

  elemental subroutine assign_qsc (a, b)
    type (qs_complex), intent(inout) :: a
    type (qs_complex), intent(in) :: b
    a%cmp = b%cmp
  end subroutine assign_qsc

  elemental subroutine assign_qsc_qs (qsc, qs)
    type (qs_complex), intent (inout) :: qsc
    type (qs_real), intent(in) :: qs
    qsc%cmp(1:4) = qs%re
    qsc%cmp(5:8) = 0.e0
  end subroutine assign_qsc_qs

  elemental subroutine assign_qs_qsc (qs, qsc)
    type (qs_real), intent (inout) :: qs
    type (qs_complex), intent(in) :: qsc
    qs%re = qsc%cmp(1:4)
  end subroutine assign_qs_qsc

  elemental subroutine assign_qsc_d (qsc, d)
    type (qs_complex), intent (inout) :: qsc
    real*4, intent(in) :: d
    qsc%cmp(1) = d
    qsc%cmp(2:8) = 0.e0
  end subroutine assign_qsc_d

  elemental subroutine assign_qsc_i (qsc, i)
    type (qs_complex), intent (inout) :: qsc
    integer, intent(in) :: i
    qsc%cmp(1) = i
    qsc%cmp(2:8) = 0.e0
  end subroutine assign_qsc_i

  elemental subroutine assign_d_qsc (d, qsc)
    real*4, intent(inout) :: d
    type (qs_complex), intent (in) :: qsc
    d = qsc%cmp(1)
  end subroutine assign_d_qsc

  elemental subroutine assign_qsc_dc (qsc, dc)
    type (qs_complex), intent (inout) :: qsc
    complex (kind (0.e0)), intent (in) :: dc
    qsc%cmp(1) = dble (dc)
    qsc%cmp(2:4) = 0.e0
    qsc%cmp(5) = aimag (dc)
    qsc%cmp(6:8) = 0.e0
  end subroutine assign_qsc_dc

  elemental subroutine assign_dc_qsc (dc, qsc)
    complex (kind (0.e0)), intent (inout) :: dc
    type (qs_complex), intent (in) :: qsc
    dc = cmplx (qsc%cmp(1), qsc%cmp(5), kind (0.e0))
  end subroutine assign_dc_qsc

! Conversions

  elemental type (qs_real) function to_qs_i(ia)
    integer, intent(in) :: ia
    to_qs_i%re(1) = ia
    to_qs_i%re(2:4) = 0.e0
  end function to_qs_i

  elemental type (qs_real) function to_qs_d(d)
    real*4, intent(in) :: d
    to_qs_d%re(1) = d
    to_qs_d%re(2:4) = 0.0e0
  end function to_qs_d

  elemental real*4 function to_d_qs(qs)
    type (qs_real), intent(in) :: qs
    to_d_qs = qs%re(1)
  end function to_d_qs

  elemental integer function to_int_qs(a)
    type (qs_real), intent(in) :: a
    to_int_qs = a%re(1)
  end function to_int_qs

  elemental type (qs_real) function to_qs_dd (dd)
    type (dd_real), intent(in) :: dd
    to_qs_dd%re(1:2) = dd%re
    to_qs_dd%re(3:4) = 0.e0
  end function to_qs_dd

  elemental type (qs_real) function to_qs_qs (qs)
    type (qs_real), intent(in) :: qs
    to_qs_qs%re = qs%re
  end function to_qs_qs

  elemental type (dd_real) function to_dd_qs (qs)
    type (qs_real), intent(in) :: qs
    to_dd_qs%re = qs%re(1:2)
  end function to_dd_qs

  type (qs_real) function to_qs_str(s)
    character (len=*), intent(in) :: s
    character*80 t
    t = s
    call qsinpc (t, to_qs_str%re)
  end function to_qs_str

  elemental type (qs_real) function to_qs_qsc(qsc)
    type (qs_complex), intent(in) :: qsc
    to_qs_qsc%re = qsc%cmp(1:4)
  end function to_qs_qsc

  elemental type (qs_complex) function to_qsc_qs(qs)
    type (qs_real), intent(in) :: qs
    to_qsc_qs%cmp(1:4) = qs%re
    to_qsc_qs%cmp(5:8) = 0.e0
  end function to_qsc_qs

  elemental type (qs_complex) function to_qsc_qs2(x, y)
    type (qs_real), intent(in) :: x, y
    to_qsc_qs2%cmp(1:4) = x%re
    to_qsc_qs2%cmp(5:8) = y%re
  end function to_qsc_qs2

  elemental type (qs_complex) function to_qsc_d(d)
    real*4, intent(in) :: d
    to_qsc_d%cmp(1) = d
    to_qsc_d%cmp(2:8) = 0.e0
  end function to_qsc_d

  elemental complex (kind (0.e0)) function to_dc_qsc (qsc)
    type (qs_complex), intent (in) :: qsc
    to_dc_qsc = cmplx (qsc%cmp(1), qsc%cmp(5), kind (0.e0))
  end function to_dc_qsc

  elemental type (qs_complex) function to_qsc_dc (dc)
    complex (kind (0.e0)), intent(in) :: dc
    to_qsc_dc%cmp(1) = dble (dc)
    to_qsc_dc%cmp(2:4) = 0.e0
    to_qsc_dc%cmp(5) = aimag (dc)
    to_qsc_dc%cmp(6:8) = 0.e0
  end function to_qsc_dc

  elemental real*4 function to_d_qsc(qsc)
    type (qs_complex), intent(in) :: qsc
    to_d_qsc = qsc%cmp(1)
  end function to_d_qsc

!  Complex conjugation

  elemental type (qs_complex) function qscconjg (qsc)
    type (qs_complex), intent(in) :: qsc
    qscconjg%cmp(1:4) = qsc%cmp(1:4)
    qscconjg%cmp(5:8) = -qsc%cmp(5:8)
  end function qscconjg

! Additions
  elemental type (qs_real) function add_qs(a, b)
    type (qs_real), intent(in) :: a, b
    call f_qs_add(a%re, b%re, add_qs%re)
  end function add_qs

  elemental type (qs_real) function add_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_qs_add_qs_d(a%re, b, add_qs_d%re)
  end function add_qs_d

  elemental type (qs_real) function add_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    add_d_qs = add_qs_d(b, a)
  end function add_d_qs

  elemental type (qs_real) function add_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    call f_qs_add_qs_d(a%re, qs_int_to_float(b), add_qs_i%re)
  end function add_qs_i

  elemental type (qs_real) function add_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    add_i_qs = add_qs_i(b, a)
  end function add_i_qs

  elemental type (qs_complex) function add_qsc(a, b)
    type (qs_complex), intent(in) :: a, b
    call f_qs_add (a%cmp(1:4), b%cmp(1:4), add_qsc%cmp(1:4))
    call f_qs_add (a%cmp(5:8), b%cmp(5:8), add_qsc%cmp(5:8))
  end function add_qsc

  elemental type (qs_complex) function add_qsc_qs(a, b)
    type (qs_complex), intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_add (a%cmp(1:4), b%re, add_qsc_qs%cmp(1:4))
    add_qsc_qs%cmp(5:8) = a%cmp(5:8)
  end function add_qsc_qs

  elemental type (qs_complex) function add_qs_qsc(a, b)
    type (qs_real), intent(in) :: a
    type (qs_complex), intent(in) :: b
    add_qs_qsc = add_qsc_qs(b, a)
  end function add_qs_qsc

  elemental type (qs_complex) function add_qsc_d(a, b)
    type (qs_complex), intent(in) :: a
    real*4, intent(in) :: b
    type (qs_real) :: qsb
    qsb%re(1) = b
    qsb%re(2:4) = 0.e0
    call f_qs_add (a%cmp(1:4), qsb%re, add_qsc_d%cmp(1:4))
    add_qsc_d%cmp(5:8) = a%cmp(5:8)
  end function add_qsc_d

  elemental type (qs_complex) function add_d_qsc(a, b)
    real*4, intent(in) :: a
    type (qs_complex), intent(in) :: b
    add_d_qsc = add_qsc_d(b, a)
  end function add_d_qsc

! Subtractions
  elemental type (qs_real) function sub_qs(a, b)
    type (qs_real), intent(in) :: a, b
    call f_qs_sub(a%re, b%re, sub_qs%re)
  end function sub_qs

  elemental type (qs_real) function sub_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_qs_sub_qs_d(a%re, b, sub_qs_d%re)
  end function sub_qs_d

  elemental type (qs_real) function sub_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_sub_d_qs(a, b%re, sub_d_qs%re)
  end function sub_d_qs

  elemental type (qs_complex) function sub_qsc(a, b)
    type (qs_complex), intent(in) :: a, b
    call f_qs_sub (a%cmp(1:4), b%cmp(1:4), sub_qsc%cmp(1:4))
    call f_qs_sub (a%cmp(5:8), b%cmp(5:8), sub_qsc%cmp(5:8))
  end function sub_qsc

  elemental type (qs_complex) function sub_qsc_qs(a, b)
    type (qs_complex), intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_sub (a%cmp(1:4), b%re(1:4), sub_qsc_qs%cmp(1:4))
    sub_qsc_qs%cmp(5:8) = a%cmp(5:8)
  end function sub_qsc_qs

  elemental type (qs_complex) function sub_qs_qsc(a, b)
    type (qs_real), intent(in) :: a
    type (qs_complex), intent(in) :: b
    call f_qs_sub (a%re(1:4), b%cmp(1:4), sub_qs_qsc%cmp(1:4))
    sub_qs_qsc%cmp(5:8) = - b%cmp(5:8)
  end function sub_qs_qsc

  elemental type (qs_complex) function sub_qsc_d(a, b)
    type (qs_complex), intent(in) :: a
    real*4, intent(in) :: b
    type (qs_real) qsb
    qsb%re(1) = b
    qsb%re(2:4) = 0.e0
    call f_qs_sub (a%cmp(1:4), qsb%re, sub_qsc_d%cmp(1:4))
    sub_qsc_d%cmp(5:8) = a%cmp(5:8)
  end function sub_qsc_d

  elemental type (qs_complex) function sub_d_qsc(a, b)
    real*4, intent(in) :: a
    type (qs_complex), intent(in) :: b
    type (qs_real) qsa
    qsa%re(1) = a
    qsa%re(2:4) = 0.e0
    call f_qs_sub (qsa%re, b%cmp(1:4), sub_d_qsc%cmp(1:4))
    sub_d_qsc%cmp(5:8) = - b%cmp(5:8)
  end function sub_d_qsc

! Unary Minus
  elemental type (qs_real) function neg_qs(a)
    type (qs_real), intent(in) :: a
    neg_qs%re = -a%re
  end function neg_qs

  elemental type (qs_complex) function neg_qsc(a)
    type (qs_complex), intent(in) :: a
    neg_qsc%cmp = -a%cmp
  end function neg_qsc

! Multiplications
  elemental type (qs_real) function mul_qs(a, b)
    type (qs_real), intent(in) :: a, b
    call f_qs_mul(a%re, b%re, mul_qs%re)
  end function mul_qs

  elemental type (qs_real) function mul_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_qs_mul_qs_d(a%re, b, mul_qs_d%re)
  end function mul_qs_d

  elemental type (qs_real) function mul_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_mul_qs_d(b%re, a, mul_d_qs%re)
  end function mul_d_qs

  elemental type (qs_real) function mul_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    call f_qs_mul_qs_d(a%re, qs_int_to_float(b), mul_qs_i%re)
  end function mul_qs_i

  elemental type (qs_real) function mul_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_mul_qs_d(b%re, qs_int_to_float(a), mul_i_qs%re)
  end function mul_i_qs

  elemental type (qs_complex) function mul_qsc(a, b)
    type (qs_complex), intent(in) :: a, b
    type (qs_real) t1, t2
    call f_qs_mul (a%cmp(1:4), b%cmp(1:4), t1%re)
    call f_qs_mul (a%cmp(5:8), b%cmp(5:8), t2%re)
    call f_qs_sub (t1%re, t2%re, mul_qsc%cmp(1:4))
    call f_qs_mul (a%cmp(1:4), b%cmp(5:8), t1%re)
    call f_qs_mul (a%cmp(5:8), b%cmp(1:4), t2%re)
    call f_qs_add (t1%re, t2%re, mul_qsc%cmp(5:8))
  end function mul_qsc

  elemental type (qs_complex) function mul_qsc_qs(a, b)
    type (qs_complex), intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_mul (a%cmp(1:4), b%re, mul_qsc_qs%cmp(1:4))
    call f_qs_mul (a%cmp(5:8), b%re, mul_qsc_qs%cmp(5:8))
  end function mul_qsc_qs

  elemental type (qs_complex) function mul_qs_qsc(a, b)
    type (qs_real), intent(in) :: a
    type (qs_complex), intent(in) :: b
    call f_qs_mul (a%re, b%cmp(1:4), mul_qs_qsc%cmp(1:4))
    call f_qs_mul (a%re, b%cmp(5:8), mul_qs_qsc%cmp(5:8))
  end function mul_qs_qsc

  elemental type (qs_complex) function mul_qsc_d(a, b)
    type (qs_complex), intent(in) :: a
    real*4, intent(in) :: b
    call f_qs_mul_qs_d (a%cmp(1:4), b, mul_qsc_d%cmp(1:4))
    call f_qs_mul_qs_d (a%cmp(5:8), b, mul_qsc_d%cmp(5:8))
  end function mul_qsc_d

  elemental type (qs_complex) function mul_d_qsc(a, b)
    real*4, intent(in) :: a
    type (qs_complex), intent(in) :: b
    call f_qs_mul_qs_d (b%cmp(1:4), a, mul_d_qsc%cmp(1:4))
    call f_qs_mul_qs_d (b%cmp(5:8), a, mul_d_qsc%cmp(5:8))
  end function mul_d_qsc

  elemental type (qs_complex) function mul_qsc_i(a, b)
    type (qs_complex), intent(in) :: a
    integer, intent(in) :: b
    call f_qs_mul_qs_d (a%cmp(1:4), qs_int_to_float(b), mul_qsc_i%cmp(1:4))
    call f_qs_mul_qs_d (a%cmp(5:8), qs_int_to_float(b), mul_qsc_i%cmp(5:8))
  end function mul_qsc_i

  elemental type (qs_complex) function mul_i_qsc(a, b)
    integer, intent(in) :: a
    type (qs_complex), intent(in) :: b
    call f_qs_mul_qs_d (b%cmp(1:4), qs_int_to_float(a), mul_i_qsc%cmp(1:4))
    call f_qs_mul_qs_d (b%cmp(5:8), qs_int_to_float(a), mul_i_qsc%cmp(5:8))
  end function mul_i_qsc

! Divisions
  elemental type (qs_real) function div_qs(a, b)
    type (qs_real), intent(in) :: a, b
    call f_qs_div(a%re, b%re, div_qs%re)
  end function div_qs

  elemental type (qs_real) function div_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    call f_qs_div_qs_d(a%re, b, div_qs_d%re)
  end function div_qs_d

  elemental type (qs_real) function div_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_div_d_qs(a, b%re, div_d_qs%re)
  end function div_d_qs

  elemental type (qs_real) function div_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    call f_qs_div_qs_d(a%re, qs_int_to_float(b), div_qs_i%re)
  end function div_qs_i

  elemental type (qs_real) function div_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_div_d_qs(qs_int_to_float(a), b%re, div_i_qs%re)
  end function div_i_qs

  elemental type (qs_complex) function div_qsc(a, b)
    type (qs_complex), intent(in) :: a, b
    type (qs_real) t1, t2, t3, t4, t5
    call f_qs_mul (a%cmp(1:4), b%cmp(1:4), t1%re)
    call f_qs_mul (a%cmp(5:8), b%cmp(5:8), t2%re)
    call f_qs_add (t1%re, t2%re, t3%re)
    call f_qs_mul (a%cmp(1:4), b%cmp(5:8), t1%re)
    call f_qs_mul (a%cmp(5:8), b%cmp(1:4), t2%re)
    call f_qs_sub (t2%re, t1%re, t4%re)
    call f_qs_mul (b%cmp(1:4), b%cmp(1:4), t1%re)
    call f_qs_mul (b%cmp(5:8), b%cmp(5:8), t2%re)
    call f_qs_add (t1%re, t2%re, t5%re)
    call f_qs_div (t3%re, t5%re, div_qsc%cmp(1:4))
    call f_qs_div (t4%re, t5%re, div_qsc%cmp(5:8))
  end function div_qsc

  elemental type (qs_complex) function div_qsc_qs(a, b)
    type (qs_complex), intent(in) :: a
    type (qs_real), intent(in) :: b
    call f_qs_div (a%cmp(1:4), b%re, div_qsc_qs%cmp(1:4))
    call f_qs_div (a%cmp(5:8), b%re, div_qsc_qs%cmp(5:8))
  end function div_qsc_qs

  elemental type (qs_complex) function div_qs_qsc(a, b)
    type (qs_real), intent(in) :: a
    type (qs_complex), intent(in) :: b
    type (qs_real) t1, t2, t3, t4, t5
    call f_qs_mul (a%re, b%cmp(1:4), t1%re)
    call f_qs_mul (a%re, b%cmp(5:8), t2%re)
    t2%re = - t2%re
    call f_qs_mul (b%cmp(1:4), b%cmp(1:4), t3%re)
    call f_qs_mul (b%cmp(5:8), b%cmp(5:8), t4%re)
    call f_qs_add (t3%re, t4%re, t5%re)
    call f_qs_div (t1%re, t5%re, div_qs_qsc%cmp(1:4))
    call f_qs_div (t2%re, t5%re, div_qs_qsc%cmp(5:8))
  end function div_qs_qsc

  elemental type (qs_complex) function div_qsc_d(a,b)
    type (qs_complex), intent(in) :: a
    real*4, intent(in) :: b
    call f_qs_div_qs_d(a%cmp(1:4), b, div_qsc_d%cmp(1:4))
    call f_qs_div_qs_d(a%cmp(5:8), b, div_qsc_d%cmp(5:8))
  end function div_qsc_d

! Power
  elemental type (qs_real) function pwr_qs (a, b)
    type (qs_real), intent(in) :: a, b
    type (qs_real) q1, q2
    call f_qs_log(a%re, q1%re)
    call f_qs_mul(q1%re, b%re, q2%re)
    call f_qs_exp(q2%re, pwr_qs%re)
  end function pwr_qs

  elemental type (qs_real) function pwr_qs_i(a, n)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: n
    call f_qs_npwr(a%re, n, pwr_qs_i%re)
  end function pwr_qs_i

  elemental type (qs_real) function pwr_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    type (qs_real) q1, q2, q3
    q1%re(1) = a
    q1%re(2:4) = 0.e0
    call f_qs_log(q1%re, q2%re)
    call f_qs_mul(q2%re, b%re, q3%re)
    call f_qs_exp(q3%re, pwr_d_qs%re)
  end function pwr_d_qs

  elemental type (qs_complex) function pwr_qsc_i(a, n)
    type (qs_complex), intent(in) :: a
    integer, intent(in) :: n
    integer i2, j, n1
    type (qs_real) t1, t2, t3
    type (qs_complex) c1, c2

    intrinsic :: iabs, ishft

    if (n == 0) then
      if (all(a%cmp == 0.e0)) then
        !write (6, *) 'pwr_qsc_i: a = 0 and n = 0'
        call f_qs_nan(pwr_qsc_i%cmp(1:4))
        call f_qs_nan(pwr_qsc_i%cmp(5))
        return
      endif
      pwr_qsc_i%cmp(1) = 1.e0
      pwr_qsc_i%cmp(2:8) = 0.e0
      return
    endif
    n1 = iabs (n)
    i2 = ishft(1, n1-1)

    c1%cmp(1) = 1.e0
    c1%cmp(2:8) = 0.e0

110 continue

    if (n1 >= i2) then
      call f_qs_mul (a%cmp(1:4), c1%cmp(1:4), t1%re)
      call f_qs_mul (a%cmp(5:8), c1%cmp(5:8), t2%re)
      call f_qs_sub (t1%re, t2%re, c2%cmp(1:4))
      call f_qs_mul (a%cmp(1:4), c1%cmp(5:8), t1%re)
      call f_qs_mul (a%cmp(5:8), c1%cmp(1:4), t2%re)
      call f_qs_add (t1%re, t2%re, c2%cmp(5:8))
      c1%cmp = c2%cmp
      n1 = n1 - i2
    endif
    i2 = i2 / 2
    if (i2 >= 1) then
      call f_qs_mul (c1%cmp(1:4), c1%cmp(1:4), t1%re)
      call f_qs_mul (c1%cmp(5:8), c1%cmp(5:8), t2%re)
      call f_qs_sub (t1%re, t2%re, c2%cmp(1:4))
      call f_qs_mul (c1%cmp(1:4), c1%cmp(5:8), t1%re)
      c2%cmp(5:8) = 2.e0 * t1%re
      c1%cmp = c2%cmp
      goto 110
    endif

    if (n > 0) then
      pwr_qsc_i%cmp = c1%cmp
    else
      c1%cmp(5:8) = - c1%cmp(5:8)
      call f_qs_mul (c1%cmp(1:4), c1%cmp(1:4), t1%re)
      call f_qs_mul (c1%cmp(5:8), c1%cmp(5:8), t2%re)
      call f_qs_add (t1%re, t2%re, t3%re)
      call f_qs_div (c1%cmp(1:4), t3%re, pwr_qsc_i%cmp(1:4))
      call f_qs_div (c1%cmp(5:8), t3%re, pwr_qsc_i%cmp(5:8))
    endif

    return
  end function pwr_qsc_i


! Trigonometric Functions
  elemental type (qs_real) function qssin(a)
    type (qs_real), intent(in) :: a
    call f_qs_sin(a%re, qssin%re)
  end function qssin

  elemental type (qs_real) function qscos(a)
    type (qs_real), intent(in) :: a
    call f_qs_cos(a%re, qscos%re)
  end function qscos

  elemental type (qs_real) function qstan(a)
    type (qs_real), intent(in) :: a
    call f_qs_tan(a%re, qstan%re)
  end function qstan

  elemental subroutine qssincos(a, s, c)
    type (qs_real), intent(in) :: a
    type (qs_real), intent(out) :: s, c
    call f_qs_sincos(a%re, s%re, c%re)
  end subroutine qssincos


! Inverse Trigonometric Functions
  elemental type (qs_real) function qsasin(a)
    type (qs_real), intent(in) :: a
    call f_qs_asin(a%re, qsasin%re)
  end function qsasin

  elemental type (qs_real) function qsacos(a)
    type (qs_real), intent(in) :: a
    call f_qs_acos(a%re, qsacos%re)
  end function qsacos

  elemental type (qs_real) function qsatan(a)
    type (qs_real), intent(in) :: a
    call f_qs_atan(a%re, qsatan%re)
  end function qsatan

  elemental type (qs_real) function qsatan2(a, b)
    type (qs_real), intent(in) :: a, b
    call f_qs_atan2(a%re, b%re, qsatan2%re)
  end function qsatan2

! Exponential and Logarithms
  elemental type (qs_real) function qsexp(a)
    type (qs_real), intent(in) :: a
    call f_qs_exp(a%re, qsexp%re)
  end function qsexp

  elemental type (qs_complex) function qscexp (a)
    type (qs_complex), intent(in) :: a
    type (qs_real) t1, t2, t3
    call f_qs_exp (a%cmp(1:4), t1%re)
    call f_qs_sincos (a%cmp(5:8), t3%re, t2%re)
    call f_qs_mul (t1%re, t2%re, qscexp%cmp(1:4))
    call f_qs_mul (t1%re, t3%re, qscexp%cmp(5:8))
  end function qscexp

  elemental type (qs_real) function qslog(a)
    type (qs_real), intent(in) :: a
    call f_qs_log(a%re, qslog%re)
  end function qslog

  elemental type (qs_complex) function qsclog (a)
    type (qs_complex), intent(in) :: a
    type (qs_real) t1, t2, t3
    call f_qs_mul (a%cmp(1:4), a%cmp(1:4), t1%re)
    call f_qs_mul (a%cmp(5:8), a%cmp(5:8), t2%re)
    call f_qs_add (t1%re, t2%re, t3%re)
    call f_qs_log (t3%re, t1%re)
    qsclog%cmp(1:4) = 0.5e0 * t1%re
    call f_qs_atan2 (a%cmp(5:8), a%cmp(1:4), qsclog%cmp(5:8))
  end function qsclog

  elemental type (qs_real) function qslog10(a)
    type (qs_real), intent(in) :: a
    call f_qs_log10(a%re, qslog10%re)
  end function qslog10


! SQRT, etc.
  elemental type (qs_real) function qssqrt(a)
    type (qs_real), intent(in) :: a
    call f_qs_sqrt(a%re, qssqrt%re)
  end function qssqrt

  elemental type (qs_real) function qssqr(a)
    type (qs_real), intent(in) :: a
    call f_qs_sqr(a%re, qssqr%re)
  end function qssqr

  elemental type (qs_real) function qsnroot(a, n)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: n
    call f_qs_nroot(a%re, n, qsnroot%re)
  end function qsnroot


! Hyperbolic Functions
  elemental type (qs_real) function qssinh(a)
    type (qs_real), intent(in) :: a
    call f_qs_sinh(a%re, qssinh%re)
  end function qssinh

  elemental type (qs_real) function qscosh(a)
    type (qs_real), intent(in) :: a
    call f_qs_cosh(a%re, qscosh%re)
  end function qscosh

  elemental type (qs_real) function qstanh(a)
    type (qs_real), intent(in) :: a
    call f_qs_tanh(a%re, qstanh%re)
  end function qstanh

  elemental subroutine qssincosh(a, s, c)
    type (qs_real), intent(in) :: a
    type (qs_real), intent(out) :: s, c
    call f_qs_sincosh(a%re, s%re, c%re)
  end subroutine qssincosh


! Inverse Hyperbolic Functions
  elemental type (qs_real) function qsasinh(a)
    type (qs_real), intent(in) :: a
    call f_qs_asinh(a%re, qsasinh%re)
  end function qsasinh

  elemental type (qs_real) function qsacosh(a)
    type (qs_real), intent(in) :: a
    call f_qs_acosh(a%re, qsacosh%re)
  end function qsacosh

  elemental type (qs_real) function qsatanh(a)
    type (qs_real), intent(in) :: a
    call f_qs_atanh(a%re, qsatanh%re)
  end function qsatanh


! Rounding
  elemental type (qs_real) function qsaint(a)
    type (qs_real), intent(in) :: a
    call f_qs_aint(a%re, qsaint%re)
  end function qsaint

  elemental type (qs_real) function qsanint(a)
    type (qs_real), intent(in) :: a
    call f_qs_nint(a%re, qsanint%re)
  end function qsanint

  elemental integer function qsnint(a)
    type (qs_real), intent(in) :: a
    qsnint = to_int_qs(qsaint(a));
  end function qsnint


! Random Number Generator
  subroutine qsrand(harvest)
    type (qs_real), intent(out) :: harvest
    call f_qs_rand(harvest%re)
  end subroutine qsrand


! Equality
  elemental logical function eq_qs(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r == 0) then
      eq_qs = .true.
    else
      eq_qs = .false.
    end if
  end function eq_qs

  elemental logical function eq_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(a%re, b, r)
    if (r == 0) then
      eq_qs_d = .true.
    else
      eq_qs_d = .false.
    end if
  end function eq_qs_d

  elemental logical function eq_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(b%re, a, r)
    if (r == 0) then
      eq_d_qs = .true.
    else
      eq_d_qs = .false.
    end if
  end function eq_d_qs

  elemental logical function eq_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    eq_qs_i = eq_qs_d(a, qs_int_to_float(b))
  end function eq_qs_i

  elemental logical function eq_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    eq_i_qs = eq_d_qs(qs_int_to_float(a), b)
  end function eq_i_qs

  elemental logical function eq_qsc (a, b)
    type (qs_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_qs_comp (a%cmp(1:4), b%cmp(1:4), i1)
    call f_qs_comp (a%cmp(5:8), b%cmp(5:8), i2)
    if (i1 == 0 .and. i2 == 0) then
      eq_qsc = .true.
    else
      eq_qsc = .false.
    endif
  end function eq_qsc

  elemental logical function eq_qsc_qs (a, b)
    type (qs_complex), intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: i1
    call f_qs_comp (a%cmp(1:4), b%re, i1)
    if (i1 == 0 .and. all(a%cmp(5:8) == 0.e0)) then
      eq_qsc_qs = .true.
    else
      eq_qsc_qs = .false.
    endif
  end function eq_qsc_qs

  elemental logical function eq_qs_qsc (a, b)
    type (qs_real), intent(in) :: a
    type (qs_complex), intent(in) :: b
    integer :: i1
    call f_qs_comp (a%re, b%cmp(1:4), i1)
    if (i1 == 0 .and. all(b%cmp(5:8) == 0.e0)) then
      eq_qs_qsc = .true.
    else
      eq_qs_qsc = .false.
    endif
  end function eq_qs_qsc


! Non-Equality
  elemental logical function ne_qs(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r == 0) then
      ne_qs = .false.
    else
      ne_qs = .true.
    end if
  end function ne_qs

  elemental logical function ne_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(a%re, b, r)
    if (r == 0) then
      ne_qs_d = .false.
    else
      ne_qs_d = .true.
    end if
  end function ne_qs_d

  elemental logical function ne_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(b%re, a, r)
    if (r == 0) then
      ne_d_qs = .false.
    else
      ne_d_qs = .true.
    end if
  end function ne_d_qs

  elemental logical function ne_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    ne_qs_i = ne_qs_d(a, qs_int_to_float(b))
  end function ne_qs_i

  elemental logical function ne_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    ne_i_qs = ne_d_qs(qs_int_to_float(a), b)
  end function ne_i_qs

  elemental logical function ne_qsc (a, b)
    type (qs_complex), intent(in) :: a, b
    integer :: i1, i2
    call f_qs_comp (a%cmp(1:4), b%cmp(1:4), i1)
    call f_qs_comp (a%cmp(5:8), b%cmp(5:8), i2)
    if (i1 /= 0 .or. i2 /= 0) then
      ne_qsc = .true.
    else
      ne_qsc = .false.
    endif
  end function ne_qsc

  elemental logical function ne_qsc_qs (a, b)
    type (qs_complex), intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: i1
    call f_qs_comp (a%cmp(1:4), b%re, i1)
    if (i1 /= 0 .or. any(a%cmp(5:8) /= 0.e0)) then
      ne_qsc_qs = .true.
    else
      ne_qsc_qs = .false.
    endif
  end function ne_qsc_qs

  elemental logical function ne_qs_qsc (a, b)
    type (qs_real), intent(in) :: a
    type (qs_complex), intent(in) :: b
    integer :: i1
    call f_qs_comp (a%re, b%cmp(1:4), i1)
    if (i1 /= 0 .or. any(b%cmp(5:8) /= 0.e0)) then
      ne_qs_qsc = .true.
    else
      ne_qs_qsc = .false.
    endif
  end function ne_qs_qsc


! Greater-Than
  elemental logical function gt_qs(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r == 1) then
      gt_qs = .true.
    else
      gt_qs = .false.
    end if
  end function gt_qs

  elemental logical function gt_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(a%re, b, r)
    if (r == 1) then
      gt_qs_d = .true.
    else
      gt_qs_d = .false.
    end if
  end function gt_qs_d

  elemental logical function gt_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(b%re, a, r)
    if (r == -1) then
      gt_d_qs = .true.
    else
      gt_d_qs = .false.
    end if
  end function gt_d_qs

  elemental logical function gt_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    gt_qs_i = gt_qs_d(a, qs_int_to_float(b))
  end function gt_qs_i

  elemental logical function gt_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    gt_i_qs = gt_d_qs(qs_int_to_float(a), b)
  end function gt_i_qs

! Less-Than
  elemental logical function lt_qs(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r == -1) then
      lt_qs = .true.
    else
      lt_qs = .false.
    end if
  end function lt_qs

  elemental logical function lt_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(a%re, b, r)
    if (r == -1) then
      lt_qs_d = .true.
    else
      lt_qs_d = .false.
    end if
  end function lt_qs_d

  elemental logical function lt_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(b%re, a, r)
    if (r == 1) then
      lt_d_qs = .true.
    else
      lt_d_qs = .false.
    end if
  end function lt_d_qs

  elemental logical function lt_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    lt_qs_i = lt_qs_d(a, qs_int_to_float(b))
  end function lt_qs_i

  elemental logical function lt_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    lt_i_qs = lt_d_qs(qs_int_to_float(a), b)
  end function lt_i_qs

! Greater-Than-Or-Equal-To
  elemental logical function ge_qs(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r >= 0) then
      ge_qs = .true.
    else
      ge_qs = .false.
    end if
  end function ge_qs

  elemental logical function ge_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(a%re, b, r)
    if (r >= 0) then
      ge_qs_d = .true.
    else
      ge_qs_d = .false.
    end if
  end function ge_qs_d

  elemental logical function ge_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(b%re, a, r)
    if (r <= 0) then
      ge_d_qs = .true.
    else
      ge_d_qs = .false.
    end if
  end function ge_d_qs

  elemental logical function ge_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    ge_qs_i = ge_qs_d(a, qs_int_to_float(b))
  end function ge_qs_i

  elemental logical function ge_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    ge_i_qs = ge_d_qs(qs_int_to_float(a), b)
  end function ge_i_qs

! Less-Than-Or-Equal-To
  elemental logical function le_qs(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r <= 0) then
      le_qs = .true.
    else
      le_qs = .false.
    end if
  end function le_qs

  elemental logical function le_qs_d(a, b)
    type (qs_real), intent(in) :: a
    real*4, intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(a%re, b, r)
    if (r <= 0) then
      le_qs_d = .true.
    else
      le_qs_d = .false.
    end if
  end function le_qs_d

  elemental logical function le_d_qs(a, b)
    real*4, intent(in) :: a
    type (qs_real), intent(in) :: b
    integer :: r
    call f_qs_comp_qs_d(b%re, a, r)
    if (r >= 0) then
      le_d_qs = .true.
    else
      le_d_qs = .false.
    end if
  end function le_d_qs

  elemental logical function le_qs_i(a, b)
    type (qs_real), intent(in) :: a
    integer, intent(in) :: b
    le_qs_i = le_qs_d(a, qs_int_to_float(b))
  end function le_qs_i

  elemental logical function le_i_qs(a, b)
    integer, intent(in) :: a
    type (qs_real), intent(in) :: b
    le_i_qs = le_d_qs(qs_int_to_float(a), b)
  end function le_i_qs


! Absolute Value
  elemental type (qs_real) function qsabs(a)
    type (qs_real), intent(in) :: a
    call f_qs_abs(a%re, qsabs%re)
  end function qsabs

  elemental type (qs_real) function qscabs (qsc)
    type (qs_complex), intent(in) :: qsc
    type (qs_real) t1, t2, t3
    call f_qs_mul (qsc%cmp(1:4), qsc%cmp(1:4), t1%re)
    call f_qs_mul (qsc%cmp(5:8), qsc%cmp(5:8), t2%re)
    call f_qs_add (t1%re, t2%re, t3%re)
    call f_qs_sqrt (t3%re, qscabs%re)
  end function qscabs

! Sign transfer
  elemental type (qs_real) function qssign(a, b) result (c)
    type (qs_real), intent(in) :: a, b
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
  end function qssign

  elemental type (qs_real) function qssign_dd_d(a, b) result (c)
    type (qs_real), intent(in) :: a
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
  end function qssign_dd_d

! Input
  subroutine qsinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qs_real), intent(in) :: q1
    type (qs_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qsinp (u, q1%re)

    if (present(q2)) then
      call qsinp (u, q2%re)
    end if

    if (present(q3)) then
      call qsinp (u, q3%re)
    end if

    if (present(q4)) then
      call qsinp (u, q4%re)
    end if

    if (present(q5)) then
      call qsinp (u, q5%re)
    end if

    if (present(q6)) then
      call qsinp (u, q6%re)
    end if

    if (present(q7)) then
      call qsinp (u, q7%re)
    end if

    if (present(q8)) then
      call qsinp (u, q8%re)
    end if

    if (present(q9)) then
      call qsinp (u, q9%re)
    end if

  end subroutine qsinpq

  subroutine qscinpq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qs_complex), intent(in) :: q1
    type (qs_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qsinp (u, q1%cmp(1:4))
    call qsinp (u, q1%cmp(5:8))

    if (present(q2)) then
      call qsinp (u, q2%cmp(1:4))
      call qsinp (u, q2%cmp(5:8))
    end if

    if (present(q3)) then
      call qsinp (u, q3%cmp(1:4))
      call qsinp (u, q3%cmp(5:8))
    end if

    if (present(q4)) then
      call qsinp (u, q4%cmp(1:4))
      call qsinp (u, q4%cmp(5:8))
    end if

    if (present(q5)) then
      call qsinp (u, q5%cmp(1:4))
      call qsinp (u, q5%cmp(5:8))
    end if

    if (present(q6)) then
      call qsinp (u, q6%cmp(1:4))
      call qsinp (u, q6%cmp(5:8))
    end if

    if (present(q7)) then
      call qsinp (u, q7%cmp(1:4))
      call qsinp (u, q7%cmp(5:8))
    end if

    if (present(q8)) then
      call qsinp (u, q8%cmp(1:4))
      call qsinp (u, q8%cmp(5:8))
    end if

    if (present(q9)) then
      call qsinp (u, q9%cmp(1:4))
      call qsinp (u, q9%cmp(5:8))
    end if

  end subroutine qscinpq

! Output
  subroutine qsoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qs_real), intent(in) :: q1
    type (qs_real), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qsout (u, q1%re)

    if (present(q2)) then
      call qsout (u, q2%re)
    end if

    if (present(q3)) then
      call qsout (u, q3%re)
    end if

    if (present(q4)) then
      call qsout (u, q4%re)
    end if

    if (present(q5)) then
      call qsout (u, q5%re)
    end if

    if (present(q6)) then
      call qsout (u, q6%re)
    end if

    if (present(q7)) then
      call qsout (u, q7%re)
    end if

    if (present(q8)) then
      call qsout (u, q8%re)
    end if

    if (present(q9)) then
      call qsout (u, q9%re)
    end if

  end subroutine qsoutq

  subroutine qscoutq(u, q1, q2, q3, q4, q5, q6, q7, q8, q9)
    integer, intent(in) :: u
    type (qs_complex), intent(in) :: q1
    type (qs_complex), intent(in), optional :: q2, q3, q4, q5, q6, q7, q8, q9

    call qsout (u, q1%cmp(1:4))
    call qsout (u, q1%cmp(5:8))

    if (present(q2)) then
      call qsout (u, q2%cmp(1:4))
      call qsout (u, q2%cmp(5:8))
    end if

    if (present(q3)) then
      call qsout (u, q3%cmp(1:4))
      call qsout (u, q3%cmp(5:8))
    end if

    if (present(q4)) then
      call qsout (u, q4%cmp(1:4))
      call qsout (u, q4%cmp(5:8))
    end if

    if (present(q5)) then
      call qsout (u, q5%cmp(1:4))
      call qsout (u, q5%cmp(5:8))
    end if

    if (present(q6)) then
      call qsout (u, q6%cmp(1:4))
      call qsout (u, q6%cmp(5:8))
    end if

    if (present(q7)) then
      call qsout (u, q7%cmp(1:4))
      call qsout (u, q7%cmp(5:8))
    end if

    if (present(q8)) then
      call qsout (u, q8%cmp(1:4))
      call qsout (u, q8%cmp(5:8))
    end if

    if (present(q9)) then
      call qsout (u, q9%cmp(1:4))
      call qsout (u, q9%cmp(5:8))
    end if

  end subroutine qscoutq

  elemental real*4 function qs_to_d(a)
    type (qs_real), intent(in) :: a
    qs_to_d = a%re(1)
  end function qs_to_d

  elemental type (qs_real) function qsmin2(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r == 1) then
      qsmin2 = b
    else
      qsmin2 = a
    end if
  end function qsmin2

  elemental type (qs_real) function qsmin(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (qs_real), intent(in) :: a1, a2, a3
    type (qs_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    qsmin = qsmin2(qsmin2(a1, a2), a3)
    if (present(a4)) qsmin = qsmin2(qsmin, a4)
    if (present(a5)) qsmin = qsmin2(qsmin, a5)
    if (present(a6)) qsmin = qsmin2(qsmin, a6)
    if (present(a7)) qsmin = qsmin2(qsmin, a7)
    if (present(a8)) qsmin = qsmin2(qsmin, a8)
    if (present(a9)) qsmin = qsmin2(qsmin, a9)
  end function qsmin

  elemental type (qs_real) function qsmax2(a, b)
    type (qs_real), intent(in) :: a, b
    integer :: r
    call f_qs_comp(a%re, b%re, r)
    if (r == -1) then
      qsmax2 = b
    else
      qsmax2 = a
    end if
  end function qsmax2

  elemental type (qs_real) function qsmax(a1, a2, a3, a4, a5, a6, a7, a8, a9)
    type (qs_real), intent(in) :: a1, a2, a3
    type (qs_real), intent(in), optional :: a4, a5, a6, a7, a8, a9
    qsmax = qsmax2(qsmax2(a1, a2), a3)
    if (present(a4)) qsmax = qsmax2(qsmax, a4)
    if (present(a5)) qsmax = qsmax2(qsmax, a5)
    if (present(a6)) qsmax = qsmax2(qsmax, a6)
    if (present(a7)) qsmax = qsmax2(qsmax, a7)
    if (present(a8)) qsmax = qsmax2(qsmax, a8)
    if (present(a9)) qsmax = qsmax2(qsmax, a9)
  end function qsmax

  elemental type (qs_real) function qsmod (a, b)
    type (qs_real), intent(in) :: a, b
    type (qs_real) :: s1, s2
    call f_qs_div (a%re, b%re, s1%re)
    call f_qs_aint(s1%re, s2%re)
    call f_qs_mul (s2%re, b%re, s1%re)
    call f_qs_sub (a%re, s1%re, qsmod%re)
  end function qsmod

  pure type (qs_real) function qs_pi()
    call f_qs_pi(qs_pi%re)
  end function qs_pi

subroutine qsinp (iu, a)

!   This routine readd the DD number A from logical unit IU.  The input
!   value must be placed on a single line of not more than 80 characters.

implicit none
integer iu, ln
parameter (ln = 80)
character*80 cs
real*4 a(4)

read (iu, '(a)', end = 100) cs
call qsinpc (cs, a)
goto 110

100 write (6, 1)
1  format ('*** qsinp: End-of-file encountered.')
! call qsabrt
stop

110 return

end subroutine

subroutine qsinpc (a, b)

!   Converts the CHARACTER*80 array A into the DD number B.

implicit none
integer i, id, ie, inz, ip, is, k, ln, lnn, beg
parameter (ln = 80)
real*4 bi
character*80 a
character*1 ai
character*10 dig
character*16 ca
parameter (dig = '0123456789')
real*4 b(4), f(4), s0(4), s1(4), s2(4)

id = 0
ip = -1
is = 0
inz = 0
s1(1) = 0.e0
s1(2) = 0.e0
s1(3) = 0.e0
s1(4) = 0.e0

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
    lnn = i-1
    goto 90
  end if
enddo

lnn = 80
90 continue

!   Scan for digits, looking for the period also.

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
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  else
!    read (ai, '(f1.0)') bi
    bi = index (dig, ai) - 1
    if (inz > 0 .or. bi > 0.e0) then
      inz = 1
      id = id + 1
! call qsmuld (s1, 10.e0, s0)
      call f_qs_mul_qs_d (s1, 10.e0, s0)
      f(1) = bi
      f(2) = 0.e0
      f(3) = 0.e0
      f(4) = 0.e0
!    call qsdqc (bi, f)
!    call qsadd (s0, f, s1)
      call f_qs_add (s0, f, s1)
    endif
  endif
enddo

100   continue
if (is .eq. -1) then
  s1(1) = - s1(1)
  s1(2) = - s1(2)
  s1(3) = - s1(3)
  s1(4) = - s1(4)
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
  elseif (index (dig, ai) .eq. 0) then
    goto 210
  else
    ie = ie + 1
    if (ie .gt. 3) goto 210
    ca(ie:ie) = ai
  endif
enddo

! read (ca, '(i4)') ie
ie = dddigin (ca, 4)
if (is .eq. -1) ie = - ie
ie = ie + ip - id
s0(1) = 10.e0
s0(2) = 0.e0
s0(3) = 0.e0
s0(4) = 0.e0
! call qsnpwr (s0, ie, s2)
call f_qs_npwr (s0, ie, s2)
! call qsmul (s1, s2, b)
call f_qs_mul (s1, s2, b)
goto 220

210  write (6, 1) a
1 format ('*** qsinpc: Syntax error in literal string: ', a)
! call qsabrt
stop

220  return

end subroutine

subroutine qsout (iu, a)

!   This routine writes the QD number A on logical unit iu using a standard
!   E format, with lines 72 characters long.

implicit none
integer iu, ln
parameter (ln = 72)
character cs(72)
real*4 a(4)

call qsoutc (a, cs)
write (iu, '  (72a)') cs

return
end subroutine

subroutine qsoutc (a, b)
  implicit none
  real*4 a(4)
  character b(72)

  b(1) = ' '
  b(2) = ' '
  call f_qs_swrite(a, 62, b(3), 70)
end subroutine

elemental type (qs_real) function qshuge(a)
  type (qs_real), intent(in) :: a
  qshuge = qs_huge
end function qshuge

elemental type (qs_real) function qs_safe_huge(a)
  type (qs_real), intent(in) :: a
  qs_safe_huge = qs_real((/3.4028235e+38, 0.0e0, 0.0e0, 0.0e0/))
end function qs_safe_huge

elemental type (qs_real) function qstiny(a)
  type (qs_real), intent(in) :: a
  qstiny = qs_tiny
end function qstiny

elemental type (qs_real) function qsepsilon(a)
  type (qs_real), intent(in) :: a
  qsepsilon = qs_eps
end function qsepsilon

elemental integer function qs_radix(a)
  type (qs_real), intent(in) :: a
  qs_radix = 2
end function qs_radix

elemental integer function qs_digits(a)
  type (qs_real), intent(in) :: a
  qs_digits = 94
end function qs_digits

elemental integer function qs_max_expn(a)
  type (qs_real), intent(in) :: a
  qs_max_expn = 127
end function qs_max_expn

elemental integer function qs_min_expn(a)
  type (qs_real), intent(in) :: a
  qs_min_expn = -54
end function qs_min_expn

elemental integer function qs_precision(a)
  type (qs_real), intent(in) :: a
  qs_precision = 28
end function qs_precision

elemental integer function qs_range(a)
  type (qs_real), intent(in) :: a
  qs_range = 37
end function qs_range

elemental type (qs_real) function qs_nan(a)
  type (qs_real), intent(in) :: a
  call f_qs_nan(qs_nan%re)
end function qs_nan

elemental type (qs_real) function qs_aimag(a)
  type (qs_complex), intent(in) :: a
  qs_aimag%re = a%cmp(5:8)
end function

elemental real*4 function qs_int_to_float(i)
  implicit none
  integer, intent(in) :: i
  intrinsic :: real
  qs_int_to_float = real(i, kind=4)
end function qs_int_to_float

end module qsmodule
