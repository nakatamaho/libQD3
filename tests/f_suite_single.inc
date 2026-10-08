! binary64 interoperability checks for the binary32-based modules
! (ds/ts/qs).  Included in the contains section of check_ds/check_ts/check_qs.

  subroutine check_binary64()
    type(rt) :: v
    type(ct) :: cv
    real*8 :: d, tol
    complex*16 :: zc

    ! ds_real is narrower than binary64, so allow its own epsilon as well.
    tol = max(2 * epsilon(1.0d0), dble(epsilon(v)))
    d = 1.0d0 / 3.0d0
    v = mkr(d)
    call report(abs(dble(v) - d) <= tol * d, tname, 'real*8 round trip')
    v = d
    call report(abs(dble(v) - d) <= tol * d, tname, 'assign from real*8')
    d = mkr(c_third)
    call report(abs(d - 1.0d0 / 3.0d0) <= tol * d, tname, 'assign to real*8')
    call report(kind(dble(v)) == kind(1.0d0), tname, 'dble returns real*8')
    call report(near(mkr(1) + 0.5d0, mkr('1.5'), 1) .and. &
                near(0.5d0 + mkr(1), mkr('1.5'), 1) .and. &
                near(mkr(1) - 0.5d0, mkr('0.5'), 1) .and. &
                near(2.0d0 - mkr(1), mkr(1), 1) .and. &
                near(mkr(3) * 0.5d0, mkr('1.5'), 1) .and. &
                near(0.5d0 * mkr(3), mkr('1.5'), 1) .and. &
                near(mkr(3) / 2.0d0, mkr('1.5'), 1) .and. &
                near(3.0d0 / mkr(2), mkr('1.5'), 1), tname, 'real*8 arithmetic')
    call report(mkr(2) > 1.0d0 .and. 1.0d0 < mkr(2) .and. mkr(2) == 2.0d0 &
                .and. 2.0d0 == mkr(2) .and. mkr(2) /= 1.0d0 .and. &
                mkr(2) >= 2.0d0 .and. mkr(2) <= 2.0d0, tname, 'real*8 compare')
    zc = cmplx(0.25d0, -0.5d0, kind(0.d0))
    cv = mkc(zc)
    call report(real(cv) == mkr('0.25') .and. aimag(cv) == mkr('-0.5'), &
                tname, 'complex from complex*16')
    cv = zc
    call report(real(cv) == mkr('0.25') .and. dble(cv) == 0.25d0, tname, &
                'complex assign/dble')
    cv = mkc(mkr(1), mkr(2)) * 2.0d0 + 1.0d0
    call report(cv == mkc(mkr(3), mkr(4)), tname, 'complex real*8 arithmetic')
  end subroutine check_binary64
