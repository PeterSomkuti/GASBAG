!> @brief Module to house the Beer-Lambert RT model
!> @file Beer_Lambert.f90
!> @author Peter Somkuti
!> @detail
!> Contains all needed functions for an absorption-only atmosphere


module Beer_Lambert_mod

  ! User modules
  use physical_model_addon_mod
  use math_utils_mod
  use statevector_mod

  implicit none

  logical :: BL_SPHERICAL = .true.

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
        type(scene), intent(inout) :: scn
        double precision, intent(inout) :: radiance(:)

        integer :: Nlay, Nspec
        integer :: i, j

        double precision :: this_chapman_fac
        double precision, allocatable :: cum_slopd_solar(:), cum_slopd_viewing(:)
        double precision, allocatable :: slopd_solar(:), slopd_viewing(:)

        Nlay = scn%num_active_levels - 1
        Nspec = size(radiance)

        ! Step zero:
        ! If the user wants a spherical Earth,
        ! we must modify the plane-parallel optical
        ! depths do a slant-path formulation.

        ! The scene Chapman factors
        !   scn%atm%chapman_solar
        !   scn%atm%chapman_viewing
        ! have a (1/mu) factor
        ! built into them, so we can use the same formalism
        ! for both plane-parallel and curved atmospheres, without
        ! needing to change the equations.
        ! We just modify the optical depths.

        if (BL_SPHERICAL) then

           allocate(slopd_solar(Nlay)) ! this is 1/mu * tau
           allocate(cum_slopd_solar(Nlay)) ! this is cumulative 1/mu * tau
           allocate(slopd_viewing(Nlay)) ! this is 1/mu * tau
           allocate(cum_slopd_viewing(Nlay)) ! this is cumulative 1/mu * tau

           ! This quantity represents the ratio
           !
           ! slant optical depth (spherical)
           ! ===============================
           ! optical depth (plane-parallel)

           allocate(scn%op%sph_tau_factor_solar(Nspec, Nlay))
           allocate(scn%op%sph_tau_factor_viewing(Nspec, Nlay))

           do j = 1, size(radiance)

              ! This bit here is confusing - taken from MS3
              ! Calculates cumulative optical depth from layer i
              ! up to TOA?
              do i = Nlay, 1, -1
                 cum_slopd_solar(i) = sum( &
                      scn%op%layer_tau(j,i:1:-1) * &
                      scn%atm%chapman_solar(i,i:1:-1) &
                      )
                 cum_slopd_viewing(i) = sum( &
                      scn%op%layer_tau(j,i:1:-1) * &
                      scn%atm%chapman_viewing(i,i:1:-1) &
                      )
              end do

              ! Calculates the per-layer optical depth
              do i = Nlay, 2, -1
                 slopd_solar(i) = cum_slopd_solar(i) - cum_slopd_solar(i-1)
                 slopd_viewing(i) = cum_slopd_viewing(i) - cum_slopd_viewing(i-1)
              end do
              slopd_solar(1) = scn%op%layer_tau(j, 1) * scn%atm%chapman_solar(1, 1)
              slopd_viewing(1) = scn%op%layer_tau(j, 1) * scn%atm%chapman_viewing(1, 1)

              ! Factor needed to upscale plane parallel optical depths
              do i = 1, Nlay
                 scn%op%sph_tau_factor_solar(j,i) = slopd_solar(i) / (scn%op%layer_tau(j,i) / scn%mu0)
                 scn%op%sph_tau_factor_viewing(j,i) = slopd_viewing(i) / (scn%op%layer_tau(j,i) / scn%mu)
              end do

           end do

        end if



        ! First, calculate the radiance reflected JUST above the surface
        radiance(:) = scn%op%albedo(:) / PI * scn%mu0


        ! .. and if there are gases in the atmosphere modify that TOA reflected
        ! radiance by the optical depths, which pass through the atmosphere twice.
        if (allocated(scn%op%total_tau)) then
           ! Have gas absorbers in the atmosphere?

           if (BL_SPHERICAL) then
              do j = 1, size(radiance)

                 radiance(j) = radiance(j) * exp(&
                      -1.0d0/scn%mu0 &
                      * sum(scn%op%layer_tau(j,1:nlay) * scn%op%sph_tau_factor_solar(j,1:nlay)) &
                      ) &
                      * exp( &
                      -1.0d0 / scn%mu &
                      * sum(scn%op%layer_tau(j,1:nlay) * scn%op%sph_tau_factor_viewing(j,1:nlay)) &
                      )

                 !radiance(j) = radiance(j) * exp(&
                 !     -1.0d0/scn%mu0 &
                 !     * sum(scn%op%layer_tau(j,1:nlay) * scn%atm%ps_factors_solar(1:nlay)) &
                 !     ) &
                 !     * exp( &
                 !     -1.0d0/scn%mu &
                 !     * sum(scn%op%layer_tau(j,1:nlay) * scn%atm%ps_factors_viewing(1:nlay)) &
                 !     )
              end do
           else
              do j = 1, size(radiance)
                 radiance(j) = radiance(j) * exp(-1.0d0/scn%mu0 &
                      * sum(scn%op%layer_tau(j,1:scn%num_active_levels-1))) &
                      * exp(-1.0d0/scn%mu &
                      * sum(scn%op%layer_tau(j,1:scn%num_active_levels-1)) &
                      )
              end do
           end if
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

      integer :: j, nlay

      nlay = scn%num_active_levels - 1

      if (BL_SPHERICAL) then

         do j = 1, size(psurf_jacobian)

            psurf_jacobian(j) = TOA_radiance(j) * sum( &
                 ( &
                 scn%op%sph_tau_factor_solar(j, 1:nlay) / scn%mu0 &
                 + 1.0d0 / scn%mu &
                 ) * sum(scn%op%gas_tau_dpsurf(j,1:nlay,:), dim=2) &
                 )

            !psurf_jacobian(j) = TOA_radiance(j) * sum( &
            !     ( &
            !     scn%atm%ps_factors_solar(1:nlay) / scn%mu0 &
            !     + scn%atm%ps_factors_viewing(1:nlay) / scn%mu &
            !     ) * sum(scn%op%gas_tau_dpsurf(j,1:nlay,:), dim=2) &
            !     )
         end do
      else
         do j = 1, size(psurf_jacobian)
            psurf_jacobian(j) = TOA_radiance(j) * (1.0d0 / scn%mu0 + 1.0d0 / scn%mu) &
                 * (sum(sum(scn%op%gas_tau_dpsurf(j,1:nlay,:), dim=1), dim=1))
         end do
      end if

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

      integer :: j, nlay

      nlay = scn%num_active_levels - 1

      if (BL_SPHERICAL) then
         do j = 1, size(temp_jacobian)

            temp_jacobian(j) = -TOA_radiance(j) * sum( &
                 ( &
                 scn%op%sph_tau_factor_solar(j, 1:nlay) / scn%mu0 &
                 + scn%op%sph_tau_factor_viewing(j, 1:nlay) / scn%mu &
                 ) * sum(scn%op%gas_tau_dtemp(j,1:nlay,:), dim=2) &
                 )


            !temp_jacobian(j) = -TOA_radiance(j) * sum( &
            !     ( &
            !     scn%atm%ps_factors_solar(1:nlay) / scn%mu0 &
            !     + scn%atm%ps_factors_viewing(1:nlay) / scn%mu &
            !     ) * sum(scn%op%gas_tau_dtemp(j,1:nlay,:), dim=2) &
            !     )
         end do
      else
         do j = 1, size(temp_jacobian)
            temp_jacobian(j) =  -TOA_radiance(j) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
                 * (sum(sum(scn%op%gas_tau_dtemp(j,1:nlay,:), dim=1), dim=1))
         end do
      end if


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

      integer :: j, nlay
      nlay = scn%num_active_levels - 1


      if (BL_SPHERICAL) then
         do j = 1, size(gas_jacobian)

            gas_jacobian(j) = -TOA_radiance(j) * sum( &
                 ( &
                 scn%op%sph_tau_factor_solar(j, idx_start:idx_stop) / scn%mu0 &
                 + scn%op%sph_tau_factor_viewing(j, idx_start:idx_stop) / scn%mu &
                 ) * scn%op%gas_tau(j,idx_start:idx_stop,idx_gas) &
                 ) / scale_factor

            !gas_jacobian(j) = -TOA_radiance(j) * sum( &
            !     ( &
            !     scn%atm%ps_factors_solar(idx_start:idx_stop) / scn%mu0 &
            !     + scn%atm%ps_factors_viewing(idx_start:idx_stop) / scn%mu &
            !     ) * scn%op%gas_tau(j,idx_start:idx_stop,idx_gas) &
            !     ) / scale_factor
         end do
      else
         do j = 1, size(gas_jacobian)
            gas_jacobian(j) = -TOA_radiance(j) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
                 * sum(scn%op%gas_tau(j,idx_start:idx_stop,idx_gas), dim=1) / scale_factor
         end do
      end if


    end subroutine calculate_BL_gas_subcolumn_jacobian


    !> @brief Column averaging kernel for scale retrieval
    !> @param TOA_radiance top-of-atmosphere radiance (sans SIF and ZLO)
    !> @param scn Scene object
    !> @param SV Statevector object
    !> @param gain_matrix Gain matrix (see Rodgers)
    !> @param ILS_delta_lambda ILS delta lambda array (sample, delta lambda)
    !> @param ILS_relative_response ILS relative response array (sample, delta_lambda)
    !> @param dispersion Dispersion coefficients
    !> @param psurf Current surface pressure
    !> @param num_active_levels Number of active levels in the retrieval
    !> @param N_spec Number of spectral points
    !> @param idx_start Array of start indices for retrieved gas scale coefficients
    !> @param idx_stop Array of stop indices for retrieved gas scale coefficients
    !> @param col_AK Result - column averaging kernel for each gas
    subroutine calculate_BL_scale_AK_corr( &
         TOA_radiance, &
         scn, &
         SV, &
         gain_matrix, &
         ILS_delta_lambda, &
         ILS_relative_response, &
         dispersion, &
         psurf, &
         num_active_levels, &
         N_spec, &
         idx_start, &
         idx_stop, &
         col_AK)

      implicit none
      double precision, intent(in) :: TOA_radiance(:)
      type(scene), intent(in) :: scn
      type(statevector), intent(in) :: SV
      double precision, intent(in) :: gain_matrix(:,:)
      real, intent(in) :: ILS_delta_lambda(:,:)
      real, intent(in) :: ILS_relative_response(:,:)
      double precision, intent(in) :: dispersion(:)
      double precision, intent(in) :: psurf
      integer, intent(in) :: num_active_levels
      integer, intent(in) :: N_spec
      integer, intent(in) :: idx_start(:)
      integer, intent(in) :: idx_stop(:)
      double precision, intent(inout) :: col_AK(:,:)

      integer :: N_gas, nlay
      integer :: i_gas, i_SV
      integer :: i, j

      logical :: ILS_success

      double precision, allocatable :: dI_dVMR(:)
      double precision :: dI_dVMR_conv(N_spec, num_active_levels)
      double precision :: pwgts(num_active_levels)
      double precision :: prior_VMR(num_active_levels), this_VMR(num_active_levels)
      double precision :: s_bar(num_active_levels)

      double precision, allocatable :: AK_profile(:,:), AK_profile_total(:)
      double precision, allocatable :: tmp_v1(:), tmp_v2(:)

      integer :: idx1, idx2

      N_gas = size(scn%op%gas_tau, 3)
      nlay = scn%num_active_levels - 1

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


            if (BL_SPHERICAL) then

               if (i == 1) then
                  ! TOA layer
                  do j = 1, size(dI_dVMR) ! Loop over spectral indices
                     dI_dVMR(j) = -TOA_radiance(j) * ( &
                          (scn%op%sph_tau_factor_solar(j,i) / scn%mu0) &
                          + (scn%op%sph_tau_factor_viewing(j,i) / scn%mu) &
                          ) * scn%op%gas_tau_dvmr(1,j,i,i_gas)
                  end do

               else if (i == num_active_levels) then
                  ! BOA layer
                  do j = 1, size(dI_dVMR) ! Loop over spectral indices
                     dI_dVMR(j) = -TOA_radiance(j) * ( &
                          (scn%op%sph_tau_factor_solar(j,i-1) / scn%mu0) &
                          + (scn%op%sph_tau_factor_viewing(j,i-1) / scn%mu) &
                          ) * scn%op%gas_tau_dvmr(2,j,i-1,i_gas)
                  end do

               else
                  ! everything in between
                  do j = 1, size(dI_dVMR) ! Loop over spectral indices
                     dI_dVMR(j) = -TOA_radiance(j) * ( &
                          (scn%op%sph_tau_factor_solar(j,i) / scn%mu0) &
                          + (scn%op%sph_tau_factor_solar(j,i) / scn%mu) &
                          ) &
                          * (scn%op%gas_tau_dvmr(2,j,i-1,i_gas) + &
                          scn%op%gas_tau_dvmr(1,j,i,i_gas))
                  end do

               end if

            else

               if (i == 1) then
                  ! TOA layer
                  do j = 1, size(dI_dVMR)
                     dI_dVMR(j) = -TOA_radiance(j) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
                          * scn%op%gas_tau_dvmr(1,j,i,i_gas)
                  end do
               else if (i == num_active_levels) then
                  ! BOA layer
                  do j = 1, size(dI_dVMR)
                     dI_dVMR(j) = -TOA_radiance(j) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) &
                          * scn%op%gas_tau_dvmr(2,j,i-1,i_gas)
                  end do
               else
                  ! everything in between
                  do j = 1, size(dI_dVMR)
                     dI_dVMR(j) = -TOA_radiance(j) * ((1.0d0 / scn%mu0) + (1.0d0 / scn%mu)) * &
                          (scn%op%gas_tau_dvmr(2,j,i-1,i_gas) + &
                          scn%op%gas_tau_dvmr(1,j,i,i_gas))
                  end do
               end if

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

            idx1 = max(idx_start(i_SV), 1)
            idx2 = min(idx_stop(i_SV), num_active_levels)

            prior_VMR(idx1:idx2) = prior_VMR(idx1:idx2) &
                 * SV%svap(SV%idx_gas(i_SV, 1))
            this_VMR(idx1:idx2) = this_VMR(idx1:idx2) &
                 * SV%svsv(SV%idx_gas(i_SV, 1))
         end do

         AK_profile_total(:) = 0.0d0
         do i_SV=1, SV%num_gas
            if (i_gas /= SV%gas_idx_lookup(i_SV)) cycle

            idx1 = max(idx_start(i_SV), 1)
            idx2 = min(idx_stop(i_SV), num_active_levels)

            s_bar(:) = 0.0d0
            s_bar(idx1:idx2) = 1.0d0

            tmp_v2(:) = 0.0d0
            do i=1, num_active_levels
               ! G dot dI/dVMR is a vector of length(SV)
               tmp_v1(:) = matmul(gain_matrix, dI_dVMR_conv(:, i))
               ! Now multiply by (h^T dot [prior_VMR * s_bar])
               AK_profile(i_SV, i) = tmp_v1(SV%idx_gas(i_SV, 1)) * dot_product(scn%atm%pwgts, prior_VMR * s_bar)
               AK_profile_total(i) = AK_profile_total(i) + AK_profile(i_SV, i)

            end do

         end do

         col_AK(i_gas, 1:num_active_levels) = AK_profile_total(:)

      end do

    end subroutine calculate_BL_scale_AK_corr


end module Beer_Lambert_mod
