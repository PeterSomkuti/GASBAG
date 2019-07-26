  !> @brief Aerosol module
  !> @file Aerosols.f90
  !> @author Peter Somkuti


module aerosols_mod

  ! User modules
  use file_utils_mod, only: read_mom_file, read_mie_file
  use control_mod, only: MCS, CS_aerosol
  ! Third-party modules
  use stringifor
  use logger_mod, only: logger => master_logger

  public :: read_aerosol_file

contains


  !> @brief Feed in the aerosol data, given a mom/mie file combination
  !> @param name Identifier
  !> @param aerosol Aerosol object
  !> @param mom_filename Location of mom file
  !> @param mie_filename Location of mie file
  subroutine read_aerosol_file(aerosol)

    implicit none
    type(CS_aerosol), intent(inout) :: aerosol
    logical :: success

    character(len=*), parameter :: fname = "create_aerosol"
    double precision, allocatable :: wavelengths_tmp(:)
    integer :: i
    character(len=999) :: tmp_str

    ! Attempt at read-in of mom file
    call logger%debug(fname, "Reading in mom file: " //  aerosol%mom_filename%chars())
    call read_mom_file(aerosol%mom_filename%chars(), &
         aerosol%wavelengths, aerosol%coef, success)

    ! Attempt at read-in of mie file
    call logger%debug(fname, "Reading in mie file: " //  aerosol%mie_filename%chars())
    call read_mie_file(aerosol%mie_filename%chars(), &
         wavelengths_tmp, aerosol%qext, &
         aerosol%qsca, aerosol%ssa, &
         aerosol%sigma_ext, aerosol%reff, success)

    ! Do some checks to confirm that the mie and mom file
    ! 'belong' together by checking if the wavelengths are the same.
    if (size(wavelengths_tmp) /= size(aerosol%wavelengths)) then
       call logger%fatal(fname, "Mie and Mom files have different numbers of wavelengths.")
       stop 1
    end if

    ! Check if the wavelength values in both mom and mie file correspond. There are obviously
    ! rounding errors, and probably depend on the code that wrote the files, so rather than
    ! terminating the program, just spit out a warning. The user should know whether this
    ! is acceptable or not.
    do i=1, size(aerosol%wavelengths)
       if (abs(aerosol%wavelengths(i) - wavelengths_tmp(i)) > 1.0d-3) then
          call logger%warning(fname, "Mie and Mom files have different wavelength values (>1nm).")
          write(tmp_str, '(A, ES10.6)') "Mie file: ", wavelengths_tmp(i)
          call logger%warning(fname, trim(tmp_str))
          write(tmp_str, '(A, ES10.6)') "Mom file: ", aerosol%wavelengths(i)
          call logger%warning(fname, trim(tmp_str))
       end if
    end do

  end subroutine read_aerosol_file


end module aerosols_mod
