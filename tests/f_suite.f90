! Fortran precision suite: runs tests/f_suite_body.inc against every libQD3
! Fortran module with epsilon-relative tolerances.

subroutine check_dd()
  use ddmodule, rt => dd_real, ct => dd_complex, mkr => ddreal, &
                mkc => ddcomplex, pi_r => ddpi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'dd'
  include 'f_suite_body.inc'

  subroutine extra_checks()
  end subroutine extra_checks
end subroutine check_dd

subroutine check_td()
  use tdmodule, rt => td_real, ct => td_complex, mkr => tdreal, &
                mkc => tdcomplex, pi_r => tdpi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'td'
  include 'f_suite_body.inc'

  subroutine extra_checks()
  end subroutine extra_checks
end subroutine check_td

subroutine check_qd()
  use qdmodule, rt => qd_real, ct => qd_complex, mkr => qdreal, &
                mkc => qdcomplex, pi_r => qdpi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'qd'
  include 'f_suite_body.inc'

  subroutine extra_checks()
  end subroutine extra_checks
end subroutine check_qd

subroutine check_ds()
  use dsmodule, rt => ds_real, ct => ds_complex, mkr => dsreal, &
                mkc => dscomplex, pi_r => dspi
  use ddmodule
  use qdmodule
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'ds'
  include 'f_suite_body.inc'

  subroutine extra_checks()
    call check_binary64()
    call check_dd_conversion()
    call check_qd_conversion()
  end subroutine extra_checks

  include 'f_suite_single.inc'

  ! Conversions between ds_real and dd_real must keep the precision of the
  ! narrower type rather than truncating limbs.
  subroutine check_dd_conversion()
    type(dd_real) :: w, ref
    type(rt) :: v
    real*8 :: tol
    type(dd_complex) :: wc
    type(ct) :: vc
    ref = ddreal(1) / ddreal(3)
    tol = 4 * max(dble(epsilon(v)), dble(epsilon(ref)))
    v = ref
    call report(near(v, mkr(c_third), 4), tname, 'ds = dd')
    call report(near(mkr(ref), mkr(c_third), 4), tname, 'dsreal(dd)')
    w = mkr(c_third)
    call report(abs(dble(w - ref)) <= tol, tname, 'dd = ds')
    call report(abs(dble(ddreal(mkr(c_third)) - ref)) <= tol, tname, 'ddreal(ds)')
    wc = ddcomplex(ref, -ref)
    vc = wc
    call report(near(real(vc), mkr(c_third), 4) .and. near(aimag(vc), -mkr(c_third), 4), tname, 'dsc = ddc')
    vc = mkc(mkr(c_third), -mkr(c_third))
    wc = vc
    call report(abs(dble(real(wc) - ref)) <= tol .and. abs(dble(aimag(wc) + ref)) <= tol, tname, 'ddc = dsc')
  end subroutine check_dd_conversion

  ! Conversions between ds_real and qd_real must keep the precision of the
  ! narrower type rather than truncating limbs.
  subroutine check_qd_conversion()
    type(qd_real) :: w, ref
    type(rt) :: v
    real*8 :: tol
    type(qd_complex) :: wc
    type(ct) :: vc
    ref = qdreal(1) / qdreal(3)
    tol = 4 * max(dble(epsilon(v)), dble(epsilon(ref)))
    v = ref
    call report(near(v, mkr(c_third), 4), tname, 'ds = qd')
    call report(near(mkr(ref), mkr(c_third), 4), tname, 'dsreal(qd)')
    w = mkr(c_third)
    call report(abs(dble(w - ref)) <= tol, tname, 'qd = ds')
    call report(abs(dble(qdreal(mkr(c_third)) - ref)) <= tol, tname, 'qdreal(ds)')
    wc = qdcomplex(ref, -ref)
    vc = wc
    call report(near(real(vc), mkr(c_third), 4) .and. near(aimag(vc), -mkr(c_third), 4), tname, 'dsc = qdc')
    vc = mkc(mkr(c_third), -mkr(c_third))
    wc = vc
    call report(abs(dble(real(wc) - ref)) <= tol .and. abs(dble(aimag(wc) + ref)) <= tol, tname, 'qdc = dsc')
  end subroutine check_qd_conversion
end subroutine check_ds

subroutine check_ts()
  use tsmodule, rt => ts_real, ct => ts_complex, mkr => tsreal, &
                mkc => tscomplex, pi_r => tspi
  use ddmodule
  use qdmodule
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'ts'
  include 'f_suite_body.inc'

  subroutine extra_checks()
    call check_binary64()
    call check_dd_conversion()
    call check_qd_conversion()
  end subroutine extra_checks

  include 'f_suite_single.inc'

  ! Conversions between ts_real and dd_real must keep the precision of the
  ! narrower type rather than truncating limbs.
  subroutine check_dd_conversion()
    type(dd_real) :: w, ref
    type(rt) :: v
    real*8 :: tol
    type(dd_complex) :: wc
    type(ct) :: vc
    ref = ddreal(1) / ddreal(3)
    tol = 4 * max(dble(epsilon(v)), dble(epsilon(ref)))
    v = ref
    call report(near(v, mkr(c_third), 4), tname, 'ts = dd')
    call report(near(mkr(ref), mkr(c_third), 4), tname, 'tsreal(dd)')
    w = mkr(c_third)
    call report(abs(dble(w - ref)) <= tol, tname, 'dd = ts')
    call report(abs(dble(ddreal(mkr(c_third)) - ref)) <= tol, tname, 'ddreal(ts)')
    wc = ddcomplex(ref, -ref)
    vc = wc
    call report(near(real(vc), mkr(c_third), 4) .and. near(aimag(vc), -mkr(c_third), 4), tname, 'tsc = ddc')
    vc = mkc(mkr(c_third), -mkr(c_third))
    wc = vc
    call report(abs(dble(real(wc) - ref)) <= tol .and. abs(dble(aimag(wc) + ref)) <= tol, tname, 'ddc = tsc')
  end subroutine check_dd_conversion

  ! Conversions between ts_real and qd_real must keep the precision of the
  ! narrower type rather than truncating limbs.
  subroutine check_qd_conversion()
    type(qd_real) :: w, ref
    type(rt) :: v
    real*8 :: tol
    type(qd_complex) :: wc
    type(ct) :: vc
    ref = qdreal(1) / qdreal(3)
    tol = 4 * max(dble(epsilon(v)), dble(epsilon(ref)))
    v = ref
    call report(near(v, mkr(c_third), 4), tname, 'ts = qd')
    call report(near(mkr(ref), mkr(c_third), 4), tname, 'tsreal(qd)')
    w = mkr(c_third)
    call report(abs(dble(w - ref)) <= tol, tname, 'qd = ts')
    call report(abs(dble(qdreal(mkr(c_third)) - ref)) <= tol, tname, 'qdreal(ts)')
    wc = qdcomplex(ref, -ref)
    vc = wc
    call report(near(real(vc), mkr(c_third), 4) .and. near(aimag(vc), -mkr(c_third), 4), tname, 'tsc = qdc')
    vc = mkc(mkr(c_third), -mkr(c_third))
    wc = vc
    call report(abs(dble(real(wc) - ref)) <= tol .and. abs(dble(aimag(wc) + ref)) <= tol, tname, 'qdc = tsc')
  end subroutine check_qd_conversion
end subroutine check_ts

subroutine check_qs()
  use qsmodule, rt => qs_real, ct => qs_complex, mkr => qsreal, &
                mkc => qscomplex, pi_r => qspi
  use ddmodule
  use qdmodule
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'qs'
  include 'f_suite_body.inc'

  subroutine extra_checks()
    call check_binary64()
    call check_dd_conversion()
    call check_qd_conversion()
  end subroutine extra_checks

  include 'f_suite_single.inc'

  ! Conversions between qs_real and dd_real must keep the precision of the
  ! narrower type rather than truncating limbs.
  subroutine check_dd_conversion()
    type(dd_real) :: w, ref
    type(rt) :: v
    real*8 :: tol
    type(dd_complex) :: wc
    type(ct) :: vc
    ref = ddreal(1) / ddreal(3)
    tol = 4 * max(dble(epsilon(v)), dble(epsilon(ref)))
    v = ref
    call report(near(v, mkr(c_third), 4), tname, 'qs = dd')
    call report(near(mkr(ref), mkr(c_third), 4), tname, 'qsreal(dd)')
    w = mkr(c_third)
    call report(abs(dble(w - ref)) <= tol, tname, 'dd = qs')
    call report(abs(dble(ddreal(mkr(c_third)) - ref)) <= tol, tname, 'ddreal(qs)')
    wc = ddcomplex(ref, -ref)
    vc = wc
    call report(near(real(vc), mkr(c_third), 4) .and. near(aimag(vc), -mkr(c_third), 4), tname, 'qsc = ddc')
    vc = mkc(mkr(c_third), -mkr(c_third))
    wc = vc
    call report(abs(dble(real(wc) - ref)) <= tol .and. abs(dble(aimag(wc) + ref)) <= tol, tname, 'ddc = qsc')
  end subroutine check_dd_conversion

  ! Conversions between qs_real and qd_real must keep the precision of the
  ! narrower type rather than truncating limbs.
  subroutine check_qd_conversion()
    type(qd_real) :: w, ref
    type(rt) :: v
    real*8 :: tol
    type(qd_complex) :: wc
    type(ct) :: vc
    ref = qdreal(1) / qdreal(3)
    tol = 4 * max(dble(epsilon(v)), dble(epsilon(ref)))
    v = ref
    call report(near(v, mkr(c_third), 4), tname, 'qs = qd')
    call report(near(mkr(ref), mkr(c_third), 4), tname, 'qsreal(qd)')
    w = mkr(c_third)
    call report(abs(dble(w - ref)) <= tol, tname, 'qd = qs')
    call report(abs(dble(qdreal(mkr(c_third)) - ref)) <= tol, tname, 'qdreal(qs)')
    wc = qdcomplex(ref, -ref)
    vc = wc
    call report(near(real(vc), mkr(c_third), 4) .and. near(aimag(vc), -mkr(c_third), 4), tname, 'qsc = qdc')
    vc = mkc(mkr(c_third), -mkr(c_third))
    wc = vc
    call report(abs(dble(real(wc) - ref)) <= tol .and. abs(dble(aimag(wc) + ref)) <= tol, tname, 'qdc = qsc')
  end subroutine check_qd_conversion
end subroutine check_qs

subroutine f_main
  use f_suite_support
  implicit none
  integer*4 old_cw

  call f_fpu_fix_start(old_cw)
  call check_dd()
  call check_td()
  call check_qd()
  call check_ds()
  call check_ts()
  call check_qs()
  call f_fpu_fix_end(old_cw)

  write (*, '(a,i0,a,i0,a)') 'f_suite: ', n_checks - n_failures, '/', &
    n_checks, ' checks passed'
  if (n_failures /= 0) stop 1
end subroutine f_main
