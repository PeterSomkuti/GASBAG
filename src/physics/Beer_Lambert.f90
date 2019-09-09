!> @brief Module to house the Beer-Lambert RT model
!> @file Beer_Lambert.f90
!> @author Peter Somkuti



module Beer_Lambert_mod

  ! User modules
  use physical_model_addon_mod
  use math_utils_mod

  implicit none

  public :: calculate_BL_radiance, calculate_BL_psurf_jacobian, &
       calculate_BL_temp_jacobian, calculate_BL_albedo_jacobian, &
       calculate_BL_gas_subcolumn_jacobian


contains

  !> @brief Beer-Lambert type TOA radiances
  !> @param scn Scene object
  !> @param radiance Per-pixel TOA transmission spectrum
  subroutine calculate_BL_radiance(scn, radiance)

        implicit none
        type(scene), intent(in) :: scn
        double precision, intent(inout) :: radiance(:)

        ! First, calculate the radiance reflected JUST above the surface
        radiance(:) = scn%op%albedo(:) / PI * scn%mu0

        ! .. and if there are gases in the atmosphere modify that TOA reflected
        ! radiance by the optical depths, which pass through the atmosphere twice.
        if (allocated(scn%op%total_tau)) then
           ! Have gas absorbers in the atmosphere?
           radiance(:) = radiance(:) * exp(-1.0d0/scn%mu0 * scn%op%total_tau(:)) &
                * exp(-1.0d0/scn%mu * scn%op%total_tau(:))
        endif

    end subroutine

    !> @brief Surface pressure Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param scn Scene object
    !> @param psurf_jacobian dI/dpsurf array (wavelength)
    subroutine calculate_BL_psurf_jacobian(TOA_radiance, scn, psurf_jacobian)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      type(scene), intent(in) :: scn
      double precision, intent(inout) :: psurf_jacobian(:)


      psurf_jacobian(:) = TOA_radiance(:) * (1.0d0 / scn%mu0 + 1.0d0 / scn%mu) &
           * (sum(sum(scn%op%gas_tau_dpsurf, dim=2), dim=2))

    end subroutine calculate_BL_psurf_jacobian

    !> @brief Temperature offset Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param scn Scene object
    !> @param temp_jacobian dI/dTemp array (wavelength)
    subroutine calculate_BL_temp_jacobian(TOA_radiance, scn, temp_jacobian)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      type(scene), intent(in) :: scn
      double precision, intent(inout) :: temp_jacobian(:)

      temp_jacobian(:) =  -TOA_radiance(:) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
           * (sum(sum(scn%op%gas_tau_dtemp, dim=2), dim=2) &
           - sum(sum(scn%op%gas_tau, dim=2), dim=2))

    end subroutine calculate_BL_temp_jacobian

    !> @brief Albedo Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param scn Scene object
    !> @param center_pixel Which pixel is the reference?
    !> @param albedo_coeff Which albedo order?
    !> @param albedo_jacobian dI/dAlbedoCoeff
    subroutine calculate_BL_albedo_jacobian(TOA_radiance, scn, &
         center_pixel, albedo_coeff, albedo_jacobian)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      type(scene), intent(in) :: scn
      integer, intent(in) :: center_pixel
      integer, intent(in) :: albedo_coeff
      double precision, intent(inout) :: albedo_jacobian(:)

      albedo_jacobian(:) = TOA_radiance(:) / scn%op%albedo(:) * &
           ((scn%op%wl(:) - scn%op%wl(center_pixel)) ** dble(albedo_coeff-1))

    end subroutine calculate_BL_albedo_jacobian

    !> @brief Gas sub-column Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param scn Scene object
    !> @param idx_start Gas sub-column start layer index
    !> @param idx_stop Gas sub-column end layer index
    !> @param idx_gas Gas index
    !> @param scale_factor the currently applied scale factor to this subcolumn
    !> @param gas_jacobian dI/dSubColumnScaleFactor
    subroutine calculate_BL_gas_subcolumn_jacobian(TOA_radiance, scn, &
         idx_start, idx_stop, idx_gas, scale_factor, gas_jacobian)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      type(scene), intent(in) :: scn
      integer, intent(in) :: idx_start
      integer, intent(in) :: idx_stop
      integer, intent(in) :: idx_gas
      double precision, intent(in) :: scale_factor
      double precision, intent(inout) :: gas_jacobian(:)

      gas_jacobian(:) = -TOA_radiance(:) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
           * sum(scn%op%gas_tau(:,idx_start:idx_stop,idx_gas), dim=2) / scale_factor

    end subroutine calculate_BL_gas_subcolumn_jacobian

end module Beer_Lambert_mod
