! Type-generic body of the Fortran precision suite.
!
! The including subroutine must import one libQD3 Fortran module with the
! renames rt (real type), ct (complex type), mkr (real constructor),
! mkc (complex constructor) and pi_r (pi), use f_suite_support, and declare
! the character parameter tname; it must also provide extra_checks after the
! include (it may be empty).  Tolerances are multiples of epsilon(rt),
! so a check only passes at the precision of the type under test.

  type(rt) :: x, y, z, s, c, one
  type(ct) :: cx, cy, cz
  integer :: i
  logical :: ok

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
  call report(floor(x) == mkr(2) .and. floor(-x) == mkr(-3), tname, 'floor')
  call report(ceiling(x) == mkr(3) .and. ceiling(-x) == mkr(-2) .and. &
              ceil(x) == mkr(3), tname, 'ceiling')
  call report(floor(mkr(3) - epsilon(x)) == mkr(2), tname, 'floor below integer')
  call report(near(modulo(mkr('-7.5'), mkr(2)), mkr('0.5'), 4) .and. &
              near(modulo(mkr('7.5'), mkr(-2)), mkr('-0.5'), 4), tname, 'modulo')
  call report(dim(mkr(5), mkr(3)) == mkr(2) .and. dim(mkr(3), mkr(5)) == mkr(0), &
              tname, 'dim')
  call report(near(hypot(mkr(3), mkr(4)), mkr(5), 4), tname, 'hypot')
  call report(near(hypot(mkr(1), mkr(2)), sqrt(mkr(5)), 4), tname, 'hypot inexact')
  call report(hypot(huge(x) / 2, huge(x) / 2) > huge(x) / 2, tname, 'hypot no overflow')
  call report(precision(x) > 0 .and. range(x) > 0, tname, 'precision/range')
  call report(sign(mkr(2), mkr(-1)) == mkr(-2), tname, 'sign negative')
  call report(sign(mkr(-2), mkr(1)) == mkr(2), tname, 'sign positive')
  call report(sign(mkr(2), mkr(0)) == mkr(2), tname, 'sign +0')
  call report(sign(mkr(-2), mkr(0)) == mkr(2), tname, 'sign -a +0')
  call report(abs(mkr(-2)) == mkr(2), tname, 'abs')
  call report(min(mkr(3), one, mkr(2)) == one, tname, 'min')
  call report(max(mkr(3), one, mkr(2)) == mkr(3), tname, 'max')

  ! Mixed integer arithmetic, including integers wider than one limb.
  call report(mkr(5) - 2 == mkr(3) .and. 5 - mkr(2) == mkr(3), tname, 'sub_i')
  call report(near(one - 1 + x, x, 1), tname, 'sub_i keeps precision')
  call report(mkr(16777217) == mkr('16777217'), tname, 'to_real large int')
  call report(mkr('16777217') == 16777217 .and. 16777217 == mkr('16777217'), &
              tname, 'compare large int')
  call report(mkr(0) + 16777217 == mkr('16777217'), tname, 'add large int')
  call report(mkr('16777217') - 16777217 == mkr(0), tname, 'sub large int')
  call report(mkr(2147483647) == mkr('2147483647') .and. &
              int(mkr(2147483647)) == 2147483647, tname, 'int range limit')
  i = mkr(3) - epsilon(x)
  call report(i == 2, tname, 'assign integer = real')

  ! Truth table of every comparison operator, for real/real and both
  ! integer orders, over a < b, a == b and a > b.
  ok = .true.
  do i = 1, 3
    x = mkr(2)
    y = mkr(i + 0)
    ok = ok .and. ((x == y) .eqv. (2 == i)) .and. ((x /= y) .eqv. (2 /= i)) &
         .and. ((x < y) .eqv. (2 < i)) .and. ((x <= y) .eqv. (2 <= i)) &
         .and. ((x > y) .eqv. (2 > i)) .and. ((x >= y) .eqv. (2 >= i))
    ok = ok .and. ((x == i) .eqv. (2 == i)) .and. ((x /= i) .eqv. (2 /= i)) &
         .and. ((x < i) .eqv. (2 < i)) .and. ((x <= i) .eqv. (2 <= i)) &
         .and. ((x > i) .eqv. (2 > i)) .and. ((x >= i) .eqv. (2 >= i))
    ok = ok .and. ((i == x) .eqv. (i == 2)) .and. ((i /= x) .eqv. (i /= 2)) &
         .and. ((i < x) .eqv. (i < 2)) .and. ((i <= x) .eqv. (i <= 2)) &
         .and. ((i > x) .eqv. (i > 2)) .and. ((i >= x) .eqv. (i >= 2))
  end do
  call report(ok, tname, 'comparison truth table')

  ! Comparisons with NaN are unordered: only /= is true.  They used to
  ! report "equal", so ==, <= and >= were true and hid NaN results.
  y = nan(x)
  call report(.not. (y == y) .and. (y /= y) .and. .not. (y < one) .and. &
              .not. (y <= one) .and. .not. (y > one) .and. &
              .not. (y >= one) .and. .not. (y >= 1) .and. .not. (1 <= y), &
              tname, 'NaN comparisons unordered')

  ! Comparisons, including mixed integer forms.
  x = one + epsilon(one)
  call report(x > one .and. one < x .and. x /= one, tname, 'compare eps')
  call report(mkr(2) > 1 .and. 1 < mkr(2) .and. mkr(2) == 2 .and. &
              2 == mkr(2) .and. mkr(2) >= 2 .and. mkr(2) <= 2, &
              tname, 'compare_i')

  ! Metadata.
  call report(epsilon(x) > mkr(0) .and. one + epsilon(x) > one, tname, &
              'epsilon')
  ! epsilon() must match the C++ type: 2**-104/-157/-209 for dd/td/qd and
  ! 2**-46/-70/-94 for ds/ts/qs.
  select case (tname)
  case ('dd'); i = 104
  case ('td'); i = 157
  case ('qd'); i = 209
  case ('ds'); i = 46
  case ('ts'); i = 70
  case default; i = 94
  end select
  call report(epsilon(x) == mkr(2) ** (-i), tname, 'epsilon value')
  call report(huge(x) > one .and. tiny(x) > mkr(0) .and. tiny(x) < one, &
              tname, 'huge/tiny')
  call report(radix(x) == 2 .and. digits(x) > 0, tname, 'radix/digits')
  ! random_number must cover [0, 1) uniformly and follow qd_random_seed.
  y = mkr(0)
  z = one
  s = mkr(0)
  do i = 1, 2000
    call random_number(x)
    if (x < mkr(0) .or. x >= one) y = mkr(-1)
    if (x < z) z = x
    if (x > c .or. i == 1) c = x
    s = s + x
  end do
  call report(y == mkr(0), tname, 'random_number in [0, 1)')
  call report(z < mkr('0.01') .and. c > mkr('0.99'), tname, 'random_number range')
  call report(abs(s / 2000 - mkr('0.5')) < mkr('0.05'), tname, 'random_number mean')
  call qd_random_seed(2024)
  call random_number(x)
  call qd_random_seed(2024)
  call random_number(y)
  call report(x == y, tname, 'qd_random_seed reproducible')

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
  cz = cx + 1
  call report(real(cz) == mkr(2) .and. aimag(cz) == mkr(2), tname, 'complex add_i')
  cz = 1 - cx
  call report(real(cz) == mkr(0) .and. aimag(cz) == mkr(-2), tname, 'complex sub_i')
  cz = cx / 2
  call report(real(cz) == mkr('0.5') .and. aimag(cz) == one, tname, 'complex div_i')
  cz = 5 / cx
  call report(near(real(cz), one, 8) .and. near(aimag(cz), mkr(-2), 8), &
              tname, 'complex i/complex')
  call report(mkc(mkr(3), mkr(0)) == mkr(3) .and. cx /= one, tname, &
              'complex compare real')
  cz = cx ** 2
  call report(near(real(cz), mkr(-3), 8) .and. near(aimag(cz), mkr(4), 8), &
              tname, 'complex pwr_i')

  call extra_checks()

  ! Complex division, abs and log must not form re**2 + im**2, which
  ! overflows above sqrt(huge) and underflows below sqrt(tiny).
  x = sqrt(huge(x)) * 4
  y = sqrt(tiny(x)) / 4
  cx = mkc(x, x)
  cz = cx / cx
  call report(near(real(cz), one, 8) .and. abs(aimag(cz)) <= 8 * epsilon(x), &
              tname, 'complex div large')
  call report(near(abs(cx), x * mkr(c_sqrt2), 8), tname, 'complex abs large')
  ! The real log of ~1e154 itself is only good to ~13 eps (dd), so this
  ! uses the transcendental tolerance; the old code returned NaN/Inf here.
  cz = log(cx)
  call report(near(real(cz), log(x) + log(mkr(2)) / 2, 64) .and. &
              near(aimag(cz), pi_r() / 4, 8), tname, 'complex log large')
  cz = x / cx
  call report(near(real(cz), mkr('0.5'), 8) .and. near(aimag(cz), mkr('-0.5'), 8), &
              tname, 'real / complex large')
  cy = mkc(y, y)
  cz = cy / cy
  call report(near(real(cz), one, 8) .and. abs(aimag(cz)) <= 8 * epsilon(x), &
              tname, 'complex div small')
  call report(abs(abs(cy) - y * mkr(c_sqrt2)) <= 8 * epsilon(x) * y, tname, &
              'complex abs small')

  ! Complex elementary functions against MPC references for z = 1+2i.
  cx = mkc(one, mkr(2))
  call report(cnear(sqrt(cx), c_csqrt_re, c_csqrt_im, 8), tname, 'complex sqrt')
  call report(cnear(sin(cx), c_csin_re, c_csin_im, 64), tname, 'complex sin')
  call report(cnear(cos(cx), c_ccos_re, c_ccos_im, 64), tname, 'complex cos')
  call report(cnear(tan(cx), c_ctan_re, c_ctan_im, 64), tname, 'complex tan')
  call report(cnear(sinh(cx), c_csinh_re, c_csinh_im, 64), tname, 'complex sinh')
  call report(cnear(cosh(cx), c_ccosh_re, c_ccosh_im, 64), tname, 'complex cosh')
  call report(cnear(tanh(cx), c_ctanh_re, c_ctanh_im, 64), tname, 'complex tanh')
  call report(cnear(asin(cx), c_casin_re, c_casin_im, 64), tname, 'complex asin')
  call report(cnear(acos(cx), c_cacos_re, c_cacos_im, 64), tname, 'complex acos')
  call report(cnear(atan(cx), c_catan_re, c_catan_im, 64), tname, 'complex atan')
  call report(cnear(asinh(cx), c_casinh_re, c_casinh_im, 64), tname, 'complex asinh')
  call report(cnear(acosh(cx), c_cacosh_re, c_cacosh_im, 64), tname, 'complex acosh')
  call report(cnear(atanh(cx), c_catanh_re, c_catanh_im, 64), tname, 'complex atanh')
  call report(cnear(cx ** mkc(mkr('0.5'), mkr('-1.5')), c_cpowc_re, c_cpowc_im, 64), &
              tname, 'complex ** complex')
  call report(cnear(cx ** mkr('1.5'), c_cpowr_re, c_cpowr_im, 64), tname, &
              'complex ** real')
  cz = sqrt(mkc(mkr(-4), mkr(0)))
  call report(real(cz) == mkr(0) .and. aimag(cz) == mkr(2), tname, &
              'complex sqrt branch')

contains

  ! True when both parts of got are within k epsilon of the references.
  logical function cnear(got, re_ref, im_ref, k)
    type(ct), intent(in) :: got
    character(len=*), intent(in) :: re_ref, im_ref
    integer, intent(in) :: k
    cnear = near(real(got), mkr(re_ref), k) .and. near(aimag(got), mkr(im_ref), k)
  end function cnear

  ! True when got is within k ulps-of-epsilon of ref (relative above 1).
  logical function near(got, ref, k)
    type(rt), intent(in) :: got, ref
    integer, intent(in) :: k
    type(rt) :: scale
    scale = abs(ref)
    if (scale < 1) scale = mkr(1)
    near = abs(got - ref) <= k * epsilon(ref) * scale
  end function near
