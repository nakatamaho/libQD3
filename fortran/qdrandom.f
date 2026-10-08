!  qdrandom.f
!
!  Seeding of the random number source shared by every libQD3 module's
!  random_number interface (and by the C and C++ random functions).
!  See include/qd/qd_random.h.

module qdrandom
  use, intrinsic :: iso_c_binding, only: c_int64_t
  implicit none
  private
  public :: qd_random_seed

  interface
    subroutine qd_srand(seed) bind(C, name='qd_srand')
      import :: c_int64_t
      integer(c_int64_t), value :: seed
    end subroutine qd_srand
  end interface

  interface qd_random_seed
    module procedure qd_random_seed_i4
    module procedure qd_random_seed_i8
  end interface

contains

  subroutine qd_random_seed_i4(seed)
    integer(kind=4), intent(in) :: seed
    call qd_srand(int(seed, c_int64_t))
  end subroutine qd_random_seed_i4

  subroutine qd_random_seed_i8(seed)
    integer(kind=8), intent(in) :: seed
    call qd_srand(int(seed, c_int64_t))
  end subroutine qd_random_seed_i8

end module qdrandom
