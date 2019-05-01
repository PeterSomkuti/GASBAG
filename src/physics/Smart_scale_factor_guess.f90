!> @brief Smart first guess for scale factors
!> @file Smart_scale_factor_guess.f90
!> @author Peter Somkuti
!>
!! Add text here

module smart_first_guess_mod

  use math_utils_mod

  implicit none

  public :: estimate_first_guess_scale_factor

contains

  !> @brief Make a smart guess about the scale factor first guess
  subroutine estimate_first_guess_scale_factor(wavelengths, radiance, &
       expected_wavelengths_in, expected_wavelengths_out, &
       expected_delta_tau, &
       SZA, VZA, scale_first_guess)

    implicit none

    double precision, intent(in) :: wavelengths(:)
    double precision, intent(in) :: radiance(:)
    double precision, intent(in) :: expected_wavelengths_in(:)
    double precision, intent(in) :: expected_wavelengths_out(:)
    double precision, intent(in) :: expected_delta_tau(:)
    double precision, intent(in) :: SZA
    double precision, intent(in) :: VZA
    double precision, intent(inout) :: scale_first_guess(:)

    integer :: i
    integer :: N_guess
    double precision :: mu_bar
    integer :: idx_wl_in, idx_wl_out
    double precision :: observed_delta_tau

    mu_bar = (cos(DEG2RAD * SZA) + cos(DEG2RAD * VZA)) / &
         (cos(DEG2RAD * SZA) * cos(DEG2RAD * VZA))

    N_guess = size(expected_delta_tau)

    do i=1, N_guess

       ! Find at which pixel our requested wavelength is
       idx_wl_in = minloc(abs(wavelengths(:) - expected_wavelengths_in(i)), dim=1)
       idx_wl_out = minloc(abs(wavelengths(:) - expected_wavelengths_out(i)), dim=1)

       observed_delta_tau = log(radiance(idx_wl_out) / radiance(idx_wl_in)) / mu_bar

       scale_first_guess(i) = observed_delta_tau / expected_delta_tau(i)

    end do


  end subroutine estimate_first_guess_scale_factor


end module smart_first_guess_mod
