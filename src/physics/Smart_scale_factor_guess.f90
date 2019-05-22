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
  !> @param wavelengths Low-res per-pixel wavelength array
  !> @param radiance Measured radiance for which the delta-tau is calculated
  !> @param expected_wavelengths_in Wavelegths of selected line cores
  !> @param expected_wavelengths_out Wavelenghts of selected continua
  !> @param expected_delta_tau The expected delta-tau between *_in and *_out
  !> @param SZA Solar zenith angle
  !> @param VZA Viewing zenith angle
  !> @param scale_first_guess Estimated first-guess for the scale factor
  !> @detail The first guess works in the following way. For a non-scattering scene,
  !> Beer-Lambert holds somewhat reasonably, so the radiance can be written like
  !> I = I_o * albedo * mu0 / PI * exp(-tau/mu0 - tau/mu). Let's now take two
  !> wavelengths, wl_in (inside line core) and wl_out (continuum), and calculate the
  !> ratio between them. If we consider albedo to be the same and create the ratio
  !> R = I(out)/I(in) we end up with R = exp[-(tau(out) + tau(in))] * mu_bar, where
  !> mu_bar = 1/mu0 + 1/mu. Rearranging gives: tau(in) - tau(out) = ln(R) / mu_bar,
  !> which is the difference in optical depth as a function of the ratio I(out)/I(in).
  !> Calculating this difference from the measured radiance and comparing it to an
  !> expected value (derived from clear-sky cases or simulations), allows to
  !> create a better first guess for the scale factor and thus resulting in
  !> less iterations for the retrieval.
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

    ! Loop variable
    integer :: i
    ! Number of wavelength pairs used for guesses
    integer :: N_guess
    ! mu_bar is 1/mu0 + 1/mu
    double precision :: mu_bar
    ! Needed to find the indices corresponding to the wavelengths in/out
    integer :: idx_wl_in, idx_wl_out
    ! The delta tau from the measurement, which is compared to the expectation
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
