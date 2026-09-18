! program fortran_test
subroutine f_main
! A simple test of the fortran wrappers

  use ddmodule
  use dsmodule
  use qdmodule
  use qsmodule
  use tdmodule
  use tsmodule
  implicit none
  integer*4 old_cw
  integer i
  type (dd_complex) ddcx
  type (dd_real) ddx
  type (ds_complex) dscx, dscy, dscz
  type (ds_real) dsx, dsy, dsz
  type (qd_complex) qdcx
  type (qs_complex) qscx, qscy, qscz
  type (qs_real) qsx, qsy, qsz
  type (td_complex) tcx, tcy, tcz
  type (td_real) tx, ty
  type (ts_complex) tscx, tscy, tscz
  type (ts_real) tsx, tsy, tsz
  type (qd_real) x, y, z
  complex(kind(0.d0)) dc

  call f_fpu_fix_start (old_cw)

  ! Test for read/write
  z = "3.14159265358979323846264338327950288419716939937510582097494459230"
  call write_scalar(6, z)

  ! Test for atan/write
  do i=1,3
    x = qdreal(dble(i))
    call write_scalar(6, x)
    y = atan(x)
    call write_scalar(6, y)
  end do

  call write_scalar(6, nan(x))
  call write_scalar(6, qdcomplex(x, y))

  tx = tdreal("1.4142135623730950488016887242096980785696718753769")
  call write_scalar(6, tx)
  ddx = tx
  y = tx
  ty = ddx
  if (abs(real(ty) - 1.4142135623730951d0) .gt. 1.0d-30) stop 11
  ty = y
  if (abs(real(ty) - 1.4142135623730951d0) .gt. 1.0d-30) stop 12

  tx = exp(log(tdreal(2.0d0)))
  if (abs(tx - tdreal(2.0d0)) .gt. tdreal("1.0d-40")) stop 13

  ty = tdreal("1.25d0")
  ddx = ty
  y = ty
  if (abs(ddx%re(1) - 1.25d0) .gt. 1.0d-30) stop 14
  if (abs(y%re(1) - 1.25d0) .gt. 1.0d-30) stop 15

  call tdrand(tx)
  if (tx < tdreal(0.d0) .or. tx > tdreal(1.d0)) stop 16

  tx = tdreal("2.75d0")
  if (floor(tx) /= tdreal(2.d0)) stop 17
  if (ceil(tx) /= tdreal(3.d0)) stop 18
  if (aint(-tx) /= tdreal(-2.d0)) stop 19
  if (abs(mod(tdreal("7.5d0"), tdreal("2.d0")) - tdreal("1.5d0")) > tdreal("1.0d-28")) stop 20
  if (min(tdreal("3.d0"), tdreal("1.d0"), tdreal("2.d0")) /= tdreal("1.d0")) stop 21
  if (max(tdreal("3.d0"), tdreal("1.d0"), tdreal("2.d0")) /= tdreal("3.d0")) stop 22

  tcx = tdcomplex(tdreal("1.25d0"), tdreal("-0.5d0"))
  tcy = conjg(tcx)
  if (aimag(tcy) /= tdreal("0.5d0")) stop 23
  tcz = tcx * tcy
  if (abs(aimag(tcz)) > tdreal("1.0d-40")) stop 24
  if (abs(real(tcz) - (tdreal("1.25d0") * tdreal("1.25d0") + tdreal("0.5d0") * tdreal("0.5d0"))) > tdreal("1.0d-40")) stop 25
  dc = cmplx(1.0d0, -2.0d0, kind(0.d0))
  tcx = dc
  if (real(tcx) /= tdreal(1.0d0)) stop 26
  if (aimag(tcx) /= tdreal(-2.0d0)) stop 27

  tcx = tdcomplex(tdreal("1.0d0") + td_eps, tdreal("-2.0d0") + tdreal("3.0d-32"))
  ddcx = tcx
  tcy = ddcx
  if (abs(real(tcy) - tdreal(ddreal(real(tcx)))) > tdreal("1.0d-40")) stop 28
  if (abs(aimag(tcy) - tdreal(ddreal(aimag(tcx)))) > tdreal("1.0d-40")) stop 29
  ddcx = ddcomplex(tcx)
  tcy = tdcomplex(ddcx)
  if (abs(real(tcy) - tdreal(ddreal(real(tcx)))) > tdreal("1.0d-40")) stop 30
  if (abs(aimag(tcy) - tdreal(ddreal(aimag(tcx)))) > tdreal("1.0d-40")) stop 31

  qdcx = tcx
  tcz = qdcx
  if (abs(real(tcz) - real(tcx)) > tdreal("1.0d-40")) stop 32
  if (abs(aimag(tcz) - aimag(tcx)) > tdreal("1.0d-40")) stop 33
  qdcx = qdcomplex(tcx)
  tcz = tdcomplex(qdcx)
  if (abs(real(tcz) - real(tcx)) > tdreal("1.0d-40")) stop 34
  if (abs(aimag(tcz) - aimag(tcx)) > tdreal("1.0d-40")) stop 35

  ! Test the double-single, triple-single, and quad-single Fortran modules.
  dsx = dsreal("1.25")
  dsy = dsreal("2.0")
  dsz = dsx + dsy
  if (abs(dsz - dsreal("3.25")) > dsreal("1.0e-5")) stop 36
  dsz = dsx * dsy
  if (abs(dsz - dsreal("2.5")) > dsreal("1.0e-5")) stop 37
  dsz = dsz / dsy
  if (abs(dsz - dsx) > dsreal("1.0e-5")) stop 38
  dsz = exp(log(dsy))
  if (abs(dsz - dsy) > dsreal("1.0e-5")) stop 39
  if (digits(dsx) /= 46 .or. maxexponent(dsx) /= 127 .or. &
      minexponent(dsx) /= -102) stop 40
  dscx = dscomplex(dsx, dsy)
  dscy = conjg(dscx)
  dscz = dscx * dscy
  if (abs(aimag(dscz)) > dsreal("1.0e-5")) stop 41
  call dsrand(dsx)
  if (dsx < dsreal(0.0e0) .or. dsx > dsreal(1.0e0)) stop 42

  tsx = tsreal("1.25")
  tsy = tsreal("2.0")
  tsz = tsx + tsy
  if (abs(tsz - tsreal("3.25")) > tsreal("1.0e-5")) stop 43
  tsz = tsx * tsy
  if (abs(tsz - tsreal("2.5")) > tsreal("1.0e-5")) stop 44
  tsz = tsz / tsy
  if (abs(tsz - tsx) > tsreal("1.0e-5")) stop 45
  tsz = exp(log(tsy))
  if (abs(tsz - tsy) > tsreal("1.0e-5")) stop 46
  tsz = floor(tsreal("2.75"))
  if (abs(tsz - tsreal("2.0")) > tsreal("1.0e-5")) stop 47
  if (abs(ceil(tsreal("2.25")) - tsreal("3.0")) > tsreal("1.0e-5")) stop 48
  if (digits(tsx) /= 70 .or. maxexponent(tsx) /= 127 .or. &
      minexponent(tsx) /= -78) stop 49
  tscx = tscomplex(tsx, tsy)
  tscy = conjg(tscx)
  tscz = tscx * tscy
  if (abs(aimag(tscz)) > tsreal("1.0e-5")) stop 50
  call tsrand(tsx)
  if (tsx < tsreal(0.0e0) .or. tsx > tsreal(1.0e0)) stop 51

  qsx = qsreal("1.25")
  qsy = qsreal("2.0")
  qsz = qsx + qsy
  if (abs(qsz - qsreal("3.25")) > qsreal("1.0e-5")) stop 52
  qsz = qsx * qsy
  if (abs(qsz - qsreal("2.5")) > qsreal("1.0e-5")) stop 53
  qsz = qsz / qsy
  if (abs(qsz - qsx) > qsreal("1.0e-5")) stop 54
  qsz = exp(log(qsy))
  if (abs(qsz - qsy) > qsreal("1.0e-5")) stop 55
  if (digits(qsx) /= 94 .or. maxexponent(qsx) /= 127 .or. &
      minexponent(qsx) /= -54) stop 56
  qscx = qscomplex(qsx, qsy)
  qscy = conjg(qscx)
  qscz = qscx * qscy
  if (abs(aimag(qscz)) > qsreal("1.0e-5")) stop 57
  call qsrand(qsx)
  if (qsx < qsreal(0.0e0) .or. qsx > qsreal(1.0e0)) stop 58

  call f_fpu_fix_end (old_cw)
end
