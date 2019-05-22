!> @brief Module to house radiative transfer models
!> @file Radiance.f90
!> @author Peter Somkuti

module radiance_mod

  ! User modules
  use math_utils_mod

  implicit none

  public :: calculate_Beer_Lambert

contains

  !> @brief Beer-Lambert type RT
  !> @param wavelengths Per-pixel wavelength array
  !> @param mu0 cos(SZA)
  !> @param mu cos(VZA)
  !> @param albedo Per-pixel albedo value
  !> @param tau Per-pixel total column optical depth
  !> @param radiance Per-pixel TOA transmission spectrum
  subroutine calculate_Beer_Lambert(wavelengths, mu0, mu, albedo, tau, &
       radiance)

        implicit none
        double precision, intent(in) :: wavelengths(:)
        double precision, intent(in) :: mu0, mu, albedo(:)
        double precision, allocatable, intent(in) :: tau(:)
        double precision, intent(inout) :: radiance(:)

        ! First, calculate the radiance reflected JUST above the surface
        radiance(:) = albedo(:) / PI * mu0

        ! .. and if there are gases in the atmosphere modify that TOA reflected
        ! radiance by the optical depths, which pass through the atmosphere twice.
        if (allocated(tau)) then
           ! Have gas absorbers in the atmosphere?
           radiance(:) = radiance(:) * exp(-1.0d0/mu0 * tau(:)) * exp(-1.0d0/mu * tau(:))
        endif

    end subroutine



end module
