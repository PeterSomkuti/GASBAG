!> @brief Module to house the Beer-Lambert RT model
!> @file Beer_Lambert.f90
!> @author Peter Somkuti



module Beer_Lambert_mod

  ! User modules
  use physical_model_addon_mod
  use math_utils_mod
  use statevector_mod

  implicit none

  public :: calculate_BL_radiance, calculate_BL_psurf_jacobian, &
       calculate_BL_temp_jacobian, calculate_BL_albedo_jacobian, &
       calculate_BL_gas_subcolumn_jacobian, &
       calculate_BL_scale_AK_corr


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
           * (sum(sum(scn%op%gas_tau_dtemp, dim=2), dim=2))

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

    !> @brief Solar irradiance scaling  Jacobian for BL-type RT
    !> @param TOA_radiance The TOA radiance array WITHOUT ZLO and SIF
    !> @param scn Scene object
    !> @param center_pixel Which pixel is the reference?
    !> @param scale_coeff Which albedo order?
    !> @param scale_jacobian dI/dSolarIrradScaleCoeff
    subroutine calculate_BL_solar_irrad_scale_jacobian(solar_radiance_unscaled, &
         scn, center_pixel, scale_coeff, scale_jacobian)

      implicit none
      double precision, intent(in) :: solar_radiance_unscaled(:)
      type(scene), intent(in) :: scn
      integer, intent(in) :: center_pixel
      integer, intent(in) :: scale_coeff
      double precision, intent(inout) :: scale_jacobian(:)

      scale_jacobian(:) = solar_radiance_unscaled(:) * &
           ((scn%op%wl(:) - scn%op%wl(center_pixel)) ** dble(scale_coeff-1))
      
    end subroutine calculate_BL_solar_irrad_scale_jacobian

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


    subroutine calculate_BL_scale_AK_corr(TOA_radiance, scn, SV, &
         gain_matrix, &
         ILS_delta_lambda, ILS_relative_response, dispersion, &
         psurf, &
         num_active_levels, N_spec, &
         idx_start, idx_stop, &
         col_AK)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      type(scene), intent(in) :: scn
      type(statevector), intent(in) :: SV
      double precision, intent(in) :: gain_matrix(:,:)
      double precision, intent(in) :: ILS_delta_lambda(:,:)
      double precision, intent(in) :: ILS_relative_response(:,:)
      double precision, intent(in) :: dispersion(:)
      double precision, intent(in) :: psurf
      integer, intent(in) :: num_active_levels
      integer, intent(in) :: N_spec
      integer, intent(in) :: idx_start(:)
      integer, intent(in) :: idx_stop(:)
      double precision, intent(inout) :: col_AK(:,:)

      integer :: N_gas
      integer :: i_gas, i_SV
      integer :: i

      logical :: ILS_success

      double precision, allocatable :: dI_dVMR(:)
      double precision :: dI_dVMR_conv(N_spec, num_active_levels)
      double precision :: pwgts(num_active_levels)
      double precision :: prior_VMR(num_active_levels), this_VMR(num_active_levels)
      double precision :: s_bar(num_active_levels)

      double precision, allocatable :: AK_profile(:,:), AK_profile_total(:)
      double precision, allocatable :: tmp_v1(:), tmp_v2(:)

      N_gas = size(scn%op%gas_tau, 3)

      allocate(dI_dVMR(size(TOA_radiance)))
      allocate(tmp_v1(size(SV%svap)))
      allocate(tmp_v2(num_active_levels))
      allocate(AK_profile(SV%num_gas, num_active_levels))
      allocate(AK_profile_total(num_active_levels))

      ! Calculate this for every retrieved gas
      do i_gas=1, N_gas

         ! Is this gas found in ANY of the state vector elements?
         ! If not - skip this gas
         if (count(SV%gas_idx_lookup(:) == i_gas) == 0) cycle

         dI_dVMR(:) = 0.0d0
         dI_dVMR_conv(:,:) = 0.0d0
         pwgts(:) = 0.0d0
         prior_VMR(:) = 0.0d0
         this_VMR(:) = 0.0d0

         ! Partial derivatives w.r.t. a change in a level
         ! VMR can be computed outside of the SV-gas loop, since
         ! those values are independent of the sub-column setup

         do i=1, num_active_levels
            if (i == 1) then
               ! TOA layer
               dI_dVMR(:) = -TOA_radiance(:) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
                    * scn%op%gas_tau_dvmr(1,:,i,i_gas)
            else if (i == num_active_levels) then
               ! BOA layer
               dI_dVMR(:) = -TOA_radiance(:) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
                    * scn%op%gas_tau_dvmr(1,:,i-1,i_gas)
            else
               ! everything in between
               dI_dVMR(:) = -TOA_radiance(:) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) * &
                    (scn%op%gas_tau_dvmr(2,:,i-1,i_gas) + &
                    scn%op%gas_tau_dvmr(1,:,i,i_gas))
            end if

            call oco_type_convolution(scn%op%wl, &
                 dI_dVMR(:), &
                 ILS_delta_lambda(:,:), &
                 ILS_relative_response(:,:), &
                 dispersion(:), &
                 dI_dVMR_conv(:,i), &
                 ILS_success)
         end do


         ! Compute / adjust the VMRs given the prior scale factor (not necessarily 1.0)
         ! and the region of the column in which the sub-column is retrieved
         prior_VMR = scn%atm%gas_vmr(1:num_active_levels, i_gas)
         this_VMR = scn%atm%gas_vmr(1:num_active_levels, i_gas)

         do i_SV=1, SV%num_gas
            if (i_gas /= SV%gas_idx_lookup(i_SV)) cycle

            prior_VMR(idx_start(i_SV):idx_stop(i_SV)) = prior_VMR(idx_start(i_SV):idx_stop(i_SV)) &
                 * SV%svap(SV%idx_gas(i_SV, 1))
            this_VMR(idx_start(i_SV):idx_stop(i_SV)) = this_VMR(idx_start(i_SV):idx_stop(i_SV)) &
                 * SV%svsv(SV%idx_gas(i_SV, 1))
         end do

         ! Calculate pressure weights for the prior VMR / column
         call pressure_weighting_function( &
              num_active_levels, &
              scn%atm%p(1:num_active_levels), &
              psurf, &
              prior_VMR(1:num_active_levels), &
              pwgts(1:num_active_levels))

         AK_profile_total(:) = 0.0d0
         do i_SV=1, SV%num_gas
            if (i_gas /= SV%gas_idx_lookup(i_SV)) cycle

            s_bar(:) = 0.0d0
            s_bar(idx_start(i_SV):idx_stop(i_SV)) = 1.0d0

            tmp_v2(:) = 0.0d0
            do i=1, num_active_levels
               ! G dot dI/dVMR is a vector of length(SV)
               tmp_v1(:) = matmul(gain_matrix, dI_dVMR_conv(:, i))
               ! Now multiply by (h^T dot [prior_VMR * s_bar])
               AK_profile(i_SV, i) = tmp_v1(SV%idx_gas(i_SV, 1)) * dot_product(pwgts, prior_VMR * s_bar)
               AK_profile_total(i) = AK_profile_total(i) + AK_profile(i_SV, i)

            end do

         end do
         col_AK(i_gas, 1:num_active_levels) = AK_profile_total(:)

      end do


    end subroutine calculate_BL_scale_AK_corr


end module Beer_Lambert_mod
