!> @brief Solar module
!> @file solar_model.f90
!> @author Peter Somkuti
!>
!! This module contains routines related to solar irradiance,
!! or the loading of the measurement-based solar transmission
!! spectra.


module solar_model_mod

  ! User modules
  use math_utils_mod

  ! Third-party modules
  use logger_mod, only: logger => master_logger

  implicit none
  
  public :: read_toon_spectrum, calculate_solar_planck_function
  private
  
contains

  !> @brief Subroutine to load a text file based solar transmittance spectrum
  !> (like the Toon solar files)
  !> @param filename Path to the text file
  !> @param solar_spectrum Output array containing the transmittance spectrum
  !> @param wl_min Lower wavelength at which the spectrum will be truncated
  !> @param wl_max Higher wavelength at which the spectrum will be truncated  
  subroutine read_toon_spectrum(filename, solar_spectrum, wl_min, wl_max)

    implicit none
    character(len=*), intent(in) :: filename
    double precision, dimension(:,:), allocatable, intent(inout) :: solar_spectrum
    double precision, optional, intent(in) :: wl_min
    double precision, optional, intent(in) :: wl_max

    ! Function name for error logging
    character(len=*), parameter :: fname = "read_toon_spectrum"
    ! File unit integer
    integer :: funit
    ! Line count variable
    integer:: line_count
    ! Loop variable
    integer :: i
    ! Variable to check IO status
    integer :: io_status
    ! Does the file exist?
    logical :: exist_file
    ! Dumme for file reading
    double precision :: dummy
    ! Array which will contain the read-in data
    double precision, allocatable :: tmp_array(:,:)
    ! Indices at which the array will be truncated
    integer :: idx_min, idx_max

    call logger%debug(fname, "Reading Toon solar spectrum at: " // trim(filename))

    ! Check if the file exists, and abort if not
    inquire(file=filename, exist=exist_file)
    if (exist_file .eqv. .false.) then
       call logger%fatal(fname, "Toon spectrum file could not be found.")
       stop 1
    end if

    ! First count how many lines we have
    open(newunit=funit, file=filename, status='old', action='read')
    line_count = 0
    do
       read(funit, *, iostat=io_status) dummy
       if (io_status < 0) then
          exit
       end if
       line_count = line_count + 1
    end do
    close(funit)

    ! And then read it all in, properly
    allocate(solar_spectrum(line_count, 2))
    open(newunit=funit, file=filename, status='old', action='read')
    do i=1, line_count
       read(funit, *) solar_spectrum(i,1), solar_spectrum(i,2)
    end do
    close(funit)
    call logger%debug(fname, "Finished reading in Toon spectrum.")

    ! Now turn them into wavelengths (microns) and re-order the array, so that
    ! we have the spectrum in units of increasing wavelength
    solar_spectrum(:,:) = solar_spectrum(line_count:1:-1,:)
    solar_spectrum(:,1) = 1e4 / solar_spectrum(:,1)
    !solar_spectrum(:,1) = 1.0d-3 *  solar_spectrum(:,1)

    ! If the user wants to, we can cut off the fairly large solar spectrum
    ! array, depending on the wavelength limits supplied to the function

    if(present(wl_min) .and. present(wl_max)) then
       idx_min = -1
       idx_max = -1

       idx_min = minloc(abs(solar_spectrum(:,1) - wl_min), dim=1)
       idx_max = minloc(abs(solar_spectrum(:,1) - wl_max), dim=1)

       allocate(tmp_array(idx_max - idx_min + 1, 2))
       tmp_array(:,:) = solar_spectrum(idx_min:idx_max, :)

       deallocate(solar_spectrum)
       allocate(solar_spectrum(idx_max - idx_min + 1, 2))
       solar_spectrum(:,:) = tmp_array(:,:)
    end if


  end subroutine read_toon_spectrum

  !> @brief Calculates per-pixel solar continuum, meant to be multiplied by transmittance
  !>
  !> @param temp Temperature in Kelvin
  !> @param solar_distance Distance between Earth and the Sun in meters
  !> @param wavelength Array containing wavelengths in micrometers
  !> @param irradiance Where the solar continuum will go
  subroutine calculate_solar_planck_function(temp, solar_distance, wavelength, irradiance)
    implicit none
    double precision, intent(in) :: temp ![K]
    double precision, intent(in) :: solar_distance ! [m]
    double precision, dimension(:), intent(in) :: wavelength ! in microns
    double precision, dimension(:), allocatable, intent(out) :: irradiance

    ! Some constants that we need for the calculations
    double precision :: solar_solid_angle, solar_angular_radius

    double precision, dimension(:), allocatable :: denom

    allocate(irradiance(size(wavelength)))
    allocate(denom(size(wavelength)))

    ! Wavelength comes in as microns, so convert them to meters

    solar_angular_radius = atan(SOLAR_RADIUS / solar_distance)
    solar_solid_angle = PI * solar_angular_radius * solar_angular_radius

    ! This is all in OCO-like units of photons/s/m2/sr/um
    irradiance(:) = 2d18 * SPEED_OF_LIGHT / ((wavelength(:))**4)
    denom(:) = exp(1d6 * PLANCK_CONST * SPEED_OF_LIGHT / (wavelength(:) * BOLTZMANN_CONST * temp)) - 1.0d0
    irradiance(:) = solar_solid_angle * irradiance(:) / denom(:)

    ! OLD calculation
    !irradiance(:) = 2 * PLANCK_CONST * (SPEED_OF_LIGHT**2) / ((wavelength(:) * 1d-6) ** 5)
    !denom(:) = exp(PLANCK_CONST * SPEED_OF_LIGHT / (wavelength(:) * 1d-6 * BOLTZMANN_CONST * temp)) - 1.0d0
    !irradiance(:) = irradiance(:) / denom(:)
    ! And convert to ph/s/m2/sr/um
    !irradiance(:) = SUN_APP_SIZE * 1d-6 * irradiance(:) / (PLANCK_CONST * SPEED_OF_LIGHT / (wavelength(:) * 1d-6))
  end subroutine calculate_solar_planck_function


end module solar_model_mod
