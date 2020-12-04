!> @brief Other misc utis
!> @file misc_utils.f90
!> @author Peter Somkuti
!>
!> @details

module misc_utils_mod

  ! System modules
  use :: OMP_LIB

  implicit none

contains


  !> @brief Returns cpu time
  function get_cpu_time() result(cpu_t)

    double precision :: cpu_t

#ifdef _OPENMP
    cpu_t = omp_get_wtime()
#else
    call cpu_time(cpu_t)
#endif

  end function get_cpu_time


end module misc_utils_mod
