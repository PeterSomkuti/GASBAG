  !> @brief Module to house the Beer-Lambert RT model
  !> @file Beer_Lambert.f90
  !> @author Peter Somkuti



module Beer_Lambert_mod

  ! User modules
  use math_utils_mod

  implicit none

  public :: calculate_BL_radiance, calculate_BL_psurf_jacobian, &
       calculate_BL_temp_jacobian, calculate_BL_albedo_jacobian, &
       calculate_BL_gas_subcolumn_jacobian


contains

  !> @brief Beer-Lambert type TOA radiances
  !> @param wavelengths Per-pixel wavelength array
  !> @param mu0 cos(SZA)
  !> @param mu cos(VZA)
  !> @param albedo Per-pixel albedo value
  !> @param tau Per-pixel total column optical depth
  !> @param radiance Per-pixel TOA transmission spectrum
  subroutine calculate_BL_radiance(wavelengths, mu0, mu, albedo, tau, &
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

    !> @brief Surface pressure Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param gas_tau_dpsurf dTau/dpsurf array (wavelength, layer, gas)
    !> @param mu0 1/cos(SZA)
    !> @param mu 1/cos(VZA)
    !> @param psurf_jacobian dI/dpsurf array (wavelength)
    subroutine calculate_BL_psurf_jacobian(TOA_radiance, gas_tau_dpsurf, mu0, mu, psurf_jacobian)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      double precision, intent(in) :: gas_tau_dpsurf(:,:,:)
      double precision, intent(in) :: mu0
      double precision, intent(in) :: mu
      double precision, intent(inout) :: psurf_jacobian(:)


      psurf_jacobian(:) = TOA_radiance(:) * (1.0d0 / mu0 + 1.0d0 / mu) &
           * (sum(sum(gas_tau_dpsurf, dim=2), dim=2))

    end subroutine calculate_BL_psurf_jacobian

    !> @brief Temperature offset Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param gas_tau Tau array (wavelength, layer, gas)
    !> @param gas_tau_dtemp Tau array (wavelength, layer, gas) with +1K perturbation
    !> @param mu0 1/cos(SZA)
    !> @param mu 1/cos(VZA)
    !> @param temp_jacobian dI/dTemp array (wavelength)
    subroutine calculate_BL_temp_jacobian(TOA_radiance, gas_tau, gas_tau_dtemp, mu0, mu, temp_jacobian)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      double precision, intent(in) :: gas_tau(:,:,:)
      double precision, intent(in) :: gas_tau_dtemp(:,:,:)
      double precision, intent(in) :: mu0
      double precision, intent(in) :: mu
      double precision, intent(inout) :: temp_jacobian(:)

      temp_jacobian(:) =  -TOA_radiance(:) * ((1.0d0 / mu0) + (1.0d0 / mu)) &
           * (sum(sum(gas_tau_dtemp, dim=2), dim=2) - sum(sum(gas_tau, dim=2), dim=2))

    end subroutine calculate_BL_temp_jacobian

    !> @brief Albedo Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param albedo_array Per-wavelength albedo value
    !> @param hires_grid Highres wavelength grid
    !> @param albedo_coeff Which albedo order?
    !> @param albedo_jacobian dI/dAlbedoCoeff
    subroutine calculate_BL_albedo_jacobian(TOA_radiance, albedo_array, hires_grid, &
         albedo_coeff, albedo_jacobian)
      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      double precision, intent(in) :: albedo_array(:)
      double precision, intent(in) :: hires_grid(:)
      integer, intent(in) :: albedo_coeff
      double precision, intent(inout) :: albedo_jacobian(:)

      albedo_jacobian(:) = TOA_radiance(:) / albedo_array(:) * &
           ((hires_grid(:) - hires_grid(1)) ** dble(albedo_coeff-1))

    end subroutine calculate_BL_albedo_jacobian

    !> @brief Gas sub-column Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param mu0 1/cos(SZA)
    !> @param mu 1/cos(VZA)
    !> @param gas_tau_subcolumn Array for gas OD subcolumns (wavelenght, layers)
    !> @param scale_factor the currently applied scale factor to this subcolumn
    !> @param gas_jacobian dI/dSubColumnScaleFactor
    subroutine calculate_BL_gas_subcolumn_jacobian(TOA_radiance, mu0, mu, gas_tau_subcolumn, &
         scale_factor, gas_jacobian)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      double precision, intent(in) :: mu0
      double precision, intent(in) :: mu
      double precision, intent(in) :: gas_tau_subcolumn(:,:)
      double precision, intent(in) :: scale_factor
      double precision, intent(inout) :: gas_jacobian(:)

      gas_jacobian(:) = -TOA_radiance(:) * ((1.0d0 / mu0) + (1.0d0 / mu)) &
           * sum(gas_tau_subcolumn(:,:), dim=2) / scale_factor

    end subroutine calculate_BL_gas_subcolumn_jacobian

end module Beer_Lambert_mod
