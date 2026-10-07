! Type-generic body of the Fortran precision suite.
!
! The including subroutine must import one libQD3 Fortran module with the
! renames rt (real type), ct (complex type), mkr (real constructor),
! mkc (complex constructor) and pi_r (pi), use f_suite_support, and declare
! the character parameter tname.  Tolerances are multiples of epsilon(rt),
! so a check only passes at the precision of the type under test.

  type(rt) :: x, y, z, s, c, one
  type(ct) :: cx, cy, cz

  one = mkr(1)

  ! Arithmetic and parsing.
  x = mkr(1) / mkr(3)
  call report(near(x, mkr(c_third), 2), tname, 'div 1/3')
  call report(near(x * 3, one, 4), tname, 'mul_i (1/3)*3')
  call report(near(3 * x, one, 4), tname, 'mul_i 3*(1/3)')
  call report(near(x + x + x, one, 4), tname, 'add (1/3)*3')
  call report(near(one - x, mkr(2) / mkr(3), 4), tname, 'sub 1-1/3')
  call report(near(mkr(1) / 3, x, 2), tname, 'div_i 1/3')
  call report(near(1 / mkr(3), x, 2), tname, 'div_i_x 1/3')
  call report(near(mkr(c_sqrt2) * mkr(c_sqrt2), mkr(2), 4), tname, 'parse sqrt2')
  call report(near(mkr('0.1') * 10, one, 4), tname, 'parse 0.1')
  call report(mkr(2) + 1 == mkr(3), tname, 'add_i')
  call report(1 + mkr(2) == mkr(3), tname, 'add_i_x')
  call report(-x == mkr(0) - x, tname, 'neg')

  ! Powers and roots.
  call report(near(sqrt(mkr(2)), mkr(c_sqrt2), 4), tname, 'sqrt')
  call report(near(sqr(mkr(c_sqrt2)), mkr(2), 4), tname, 'sqr')
  call report(near(nroot(mkr(2), 3), mkr(c_cbrt2), 8), tname, 'nroot 3')
  call report(near(mkr(c_cbrt2) ** 3, mkr(2), 8), tname, 'pwr_i')
  call report(near(mkr('2.5') ** mkr('1.5'), mkr(c_pow25_15), 64), tname, 'pwr')

  ! Exponential and logarithm.
  call report(near(exp(one), mkr(c_e), 16), tname, 'exp')
  call report(near(log(mkr(2)), mkr(c_ln2), 16), tname, 'log')
  call report(near(log10(mkr(2)), mkr(c_log10_2), 16), tname, 'log10')

  ! Trigonometric and hyperbolic functions.
  x = mkr('0.5')
  call report(near(pi_r(), mkr(c_pi), 2), tname, 'pi')
  call report(near(sin(x), mkr(c_sin05), 16), tname, 'sin')
  call report(near(cos(x), mkr(c_cos05), 16), tname, 'cos')
  call report(near(tan(x), mkr(c_tan05), 16), tname, 'tan')
  call sincos(x, s, c)
  call report(near(s, mkr(c_sin05), 16) .and. near(c, mkr(c_cos05), 16), &
              tname, 'sincos')
  call report(near(asin(x), mkr(c_asin05), 16), tname, 'asin')
  call report(near(acos(x), mkr(c_acos05), 16), tname, 'acos')
  call report(near(atan(x), mkr(c_atan05), 16), tname, 'atan')
  call report(near(atan2(one, mkr(2)), mkr(c_atan2_1_2), 16), tname, 'atan2')
  call report(near(sinh(x), mkr(c_sinh05), 16), tname, 'sinh')
  call report(near(cosh(x), mkr(c_cosh05), 16), tname, 'cosh')
  call report(near(tanh(x), mkr(c_tanh05), 16), tname, 'tanh')
  call sincosh(x, s, c)
  call report(near(s, mkr(c_sinh05), 16) .and. near(c, mkr(c_cosh05), 16), &
              tname, 'sincosh')
  call report(near(asinh(x), mkr(c_asinh05), 16), tname, 'asinh')
  call report(near(acosh(mkr('1.5')), mkr(c_acosh15), 16), tname, 'acosh')
  call report(near(atanh(x), mkr(c_atanh05), 16), tname, 'atanh')

  ! Rounding, integer conversion and sign manipulation.
  x = mkr('2.75')
  call report(aint(x) == mkr(2) .and. aint(-x) == mkr(-2), tname, 'aint')
  call report(anint(x) == mkr(3) .and. anint(-x) == mkr(-3), tname, 'anint')
  call report(nint(x) == 3 .and. nint(-x) == -3, tname, 'nint')
  call report(nint(mkr('2.25')) == 2, tname, 'nint 2.25')
  call report(int(x) == 2 .and. int(-x) == -2, tname, 'int')
  y = mkr(3) - epsilon(x)
  call report(int(y) == 2, tname, 'int below integer')
  call report(int(-y) == -2, tname, 'int above negative integer')
  call report(near(mod(mkr('7.5'), mkr(2)), mkr('1.5'), 4), tname, 'mod')
  call report(sign(mkr(2), mkr(-1)) == mkr(-2), tname, 'sign negative')
  call report(sign(mkr(-2), mkr(1)) == mkr(2), tname, 'sign positive')
  call report(sign(mkr(2), mkr(0)) == mkr(2), tname, 'sign +0')
  call report(sign(mkr(-2), mkr(0)) == mkr(2), tname, 'sign -a +0')
  call report(abs(mkr(-2)) == mkr(2), tname, 'abs')
  call report(min(mkr(3), one, mkr(2)) == one, tname, 'min')
  call report(max(mkr(3), one, mkr(2)) == mkr(3), tname, 'max')

  ! Comparisons, including mixed integer forms.
  x = one + epsilon(one)
  call report(x > one .and. one < x .and. x /= one, tname, 'compare eps')
  call report(mkr(2) > 1 .and. 1 < mkr(2) .and. mkr(2) == 2 .and. &
              2 == mkr(2) .and. mkr(2) >= 2 .and. mkr(2) <= 2, &
              tname, 'compare_i')

  ! Metadata.
  call report(epsilon(x) > mkr(0) .and. one + epsilon(x) > one, tname, &
              'epsilon')
  call report(huge(x) > one .and. tiny(x) > mkr(0) .and. tiny(x) < one, &
              tname, 'huge/tiny')
  call report(radix(x) == 2 .and. digits(x) > 0, tname, 'radix/digits')
  call random_number(x)
  call report(x >= mkr(0) .and. x < one, tname, 'random_number')

  ! Complex arithmetic and functions.
  cx = mkc(one, mkr(2))
  cy = mkc(mkr(3), mkr(-1))
  cz = cx * conjg(cx)
  call report(near(real(cz), mkr(5), 4) .and. near(aimag(cz), mkr(0), 4), &
              tname, 'complex conjg product')
  cz = (cx * cy) / cy
  call report(near(real(cz), one, 8) .and. near(aimag(cz), mkr(2), 8), &
              tname, 'complex mul/div')
  cz = (cx + cy) - cy
  call report(real(cz) == one .and. aimag(cz) == mkr(2), tname, &
              'complex add/sub')
  call report(near(abs(cx) ** 2, mkr(5), 8), tname, 'complex abs')
  cz = exp(cx)
  call report(near(real(cz), mkr(c_cexp_re), 32) .and. &
              near(aimag(cz), mkr(c_cexp_im), 32), tname, 'complex exp')
  cz = log(cx)
  call report(near(real(cz), mkr(c_clog_re), 32) .and. &
              near(aimag(cz), mkr(c_clog_im), 32), tname, 'complex log')
  cz = cx ** 2
  call report(near(real(cz), mkr(-3), 8) .and. near(aimag(cz), mkr(4), 8), &
              tname, 'complex pwr_i')

contains

  ! True when got is within k ulps-of-epsilon of ref (relative above 1).
  logical function near(got, ref, k)
    type(rt), intent(in) :: got, ref
    integer, intent(in) :: k
    type(rt) :: scale
    scale = abs(ref)
    if (scale < 1) scale = mkr(1)
    near = abs(got - ref) <= k * epsilon(ref) * scale
  end function near
