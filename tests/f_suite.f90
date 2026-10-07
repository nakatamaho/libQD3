! Fortran precision suite: runs tests/f_suite_body.inc against every libQD3
! Fortran module with epsilon-relative tolerances.

subroutine check_dd()
  use ddmodule, rt => dd_real, ct => dd_complex, mkr => ddreal, &
                mkc => ddcomplex, pi_r => ddpi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'dd'
  include 'f_suite_body.inc'
end subroutine check_dd

subroutine check_td()
  use tdmodule, rt => td_real, ct => td_complex, mkr => tdreal, &
                mkc => tdcomplex, pi_r => tdpi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'td'
  include 'f_suite_body.inc'
end subroutine check_td

subroutine check_qd()
  use qdmodule, rt => qd_real, ct => qd_complex, mkr => qdreal, &
                mkc => qdcomplex, pi_r => qdpi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'qd'
  include 'f_suite_body.inc'
end subroutine check_qd

subroutine check_ds()
  use dsmodule, rt => ds_real, ct => ds_complex, mkr => dsreal, &
                mkc => dscomplex, pi_r => dspi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'ds'
  include 'f_suite_body.inc'
end subroutine check_ds

subroutine check_ts()
  use tsmodule, rt => ts_real, ct => ts_complex, mkr => tsreal, &
                mkc => tscomplex, pi_r => tspi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'ts'
  include 'f_suite_body.inc'
end subroutine check_ts

subroutine check_qs()
  use qsmodule, rt => qs_real, ct => qs_complex, mkr => qsreal, &
                mkc => qscomplex, pi_r => qspi
  use f_suite_support
  implicit none
  character(len=*), parameter :: tname = 'qs'
  include 'f_suite_body.inc'
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
