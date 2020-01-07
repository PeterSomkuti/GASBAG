

module physical_model_addon_mod

  ! User modules
  use math_utils_mod
  use statevector_mod
  use Rayleigh_mod
  use scene_mod
  use aerosols_mod

  ! Third-party modules
  use stringifor
  use mod_datetime
  use logger_mod, only: logger => master_logger

  ! System modules
  use ISO_FORTRAN_ENV
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

  !> This structure contains the result data that will be stored in the output HDF file.
  type result_container
     !> State vector names (SV number)
     type(string), allocatable :: sv_names(:)
     !> Retrieved state vector (SV number, footprint, frame)
     double precision, allocatable :: sv_retrieved(:,:,:)
     !> State vector prior (SV number, footprint, frame)
     double precision, allocatable :: sv_prior(:,:,:)
     !> State vector posterior uncertainty (SV number, footprint, frame)
     double precision, allocatable :: sv_uncertainty(:,:,:)
     !> Column-average dry air mixing ratio for retrieved gases (footprint, frame, gas_number)
     double precision, allocatable :: xgas(:,:,:)
     !> Column-average dry air mixing ratio for prior gases! (footprint, frame, gas_number)
     double precision, allocatable :: xgas_prior(:,:,:)
     !> Pressure levels (footprint, frame, pressure level)
     double precision, allocatable :: pressure_levels(:,:,:)
     !> Prior gas VMRs per level (footprint, frame, gas_number, pressure level)
     double precision, allocatable :: vmr_prior(:,:,:,:)
     !> Retrieved gas VMRs per level (footprint, frame, gas_number, pressure level)
     double precision, allocatable :: vmr_retrieved(:,:,:,:)
     !> Pressure weighting functions (footprint, frame, gas_number, pressure level)
     double precision, allocatable :: pwgts(:,:,:,:)
     !> Column averaging kernels (footprint, frame, gas_number, pressure level)
     double precision, allocatable :: col_ak(:,:,:,:)
     !> Final Chi2 (footprint, frame)
     double precision, allocatable :: chi2(:,:)
     !> Final residual RMS (footprint, frame)
     double precision, allocatable :: residual_rms(:,:)
     !> Final dsigma-squared
     double precision, allocatable :: dsigma_sq(:,:)
     !> Final number of iterations
     integer, allocatable :: num_iterations(:,:)
     !> Converged or not? (1=converged, 0=not converged, -1=not properly run)
     integer, allocatable :: converged(:,:)
     !> SNR estimate (mean of per-pixel SNR)
     double precision, allocatable :: SNR(:,:)
     !> SNR standard deviation (std of per-pixel SNR)
     double precision, allocatable :: SNR_std(:,:)
     !> Continuum level radiance estimate
     double precision, allocatable :: continuum(:,:)
     !> Single-sounding retrieval processing time
     double precision, allocatable :: processing_time(:,:)
     !> Number of moles of dry air per m2 for various sections
     !> of the model atmosphere - corresponding to retrieved
     !> gas scale factors.
     double precision, allocatable :: ndry(:,:,:)

  end type result_container


  public scene_altitude

contains


  !> @brief Calculates mid-layer pressures, taking into account surface pressure
  !> @param scn Scene object
  subroutine calculate_layer_pressure(scn)
    type(scene), intent(inout) :: scn

    ! Loop variable
    integer :: l

    ! Get rid of any existing mid-layer pressures
    if (allocated(scn%atm%p_layers)) deallocate(scn%atm%p_layers)
    ! Construct a new array
    allocate(scn%atm%p_layers(scn%num_active_levels - 1))

    scn%atm%p_layers(:) = -1.0d0

    ! Loop over all layers and compute mid-layer pressure
    do l = 1, scn%num_active_levels - 1
       if (l < scn%num_active_levels - 2) then
          scn%atm%p_layers(l) = 0.5d0 * (scn%atm%p(l) + scn%atm%p(l+1))
       else
          ! For bottom-most layer, use surface pressure instead
          ! of the last active pressure level
          scn%atm%p_layers(l) = 0.5d0 * (scn%atm%p(l) + scn%atm%psurf)
       end if
    end do


  end subroutine calculate_layer_pressure


  subroutine allocate_optical_properties(scn, N_hires, N_gases)

    type(scene), intent(inout) :: scn
    integer, intent(in) :: N_hires
    integer, intent(in) :: N_gases

    allocate(scn%op%wl(N_hires))
    scn%op%wl(:) = 0.0d0
    allocate(scn%op%albedo(N_hires))
    scn%op%albedo(:) = 0.0d0


    ! NOTE
    ! This explicitly does NOT include aerosols, as those are handled
    ! in the aerosol initialization routine.

    ! These arrays are only allocated if we have an atmosphere with at
    ! least one gas in it.
    if (N_gases > 0) then
       allocate(scn%op%gas_tau(N_hires, scn%num_levels-1, N_gases))
       scn%op%gas_tau(:,:,:) = 0.0d0
       allocate(scn%op%gas_tau_dpsurf(N_hires, scn%num_levels-1, N_gases))
       scn%op%gas_tau_dpsurf(:,:,:) = 0.0d0
       allocate(scn%op%gas_tau_dvmr(N_hires, scn%num_levels-1, N_gases, 2))
       scn%op%gas_tau_dvmr(:,:,:,:) = 0.0d0
       allocate(scn%op%gas_tau_dtemp(N_hires, scn%num_levels-1, N_gases))
       scn%op%gas_tau_dtemp(:,:,:) = 0.0d0
       allocate(scn%op%ray_tau(N_hires, scn%num_levels-1))
       scn%op%ray_tau(:,:) = 0.0d0
       allocate(scn%op%ray_depolf(N_hires))
       scn%op%ray_depolf(:) = 0.0d0
       allocate(scn%op%total_tau(N_hires))
       scn%op%total_tau(:) = 0.0d0
       allocate(scn%op%layer_tau(N_hires, scn%num_levels-1))
       scn%op%layer_tau(:,:) = 0.0d0
       allocate(scn%op%layer_omega(N_hires, scn%num_levels-1))
       scn%op%layer_omega(:,:) = 0.0d0
    end if

  end subroutine allocate_optical_properties


  subroutine destroy_optical_properties(scn)
    type(scene), intent(inout) :: scn

    deallocate(scn%op%wl)
    deallocate(scn%op%albedo)

    if (allocated(scn%op%gas_tau)) deallocate(scn%op%gas_tau)
    if (allocated(scn%op%gas_tau_dpsurf)) deallocate(scn%op%gas_tau_dpsurf)
    if (allocated(scn%op%gas_tau_dvmr)) deallocate(scn%op%gas_tau_dvmr)
    if (allocated(scn%op%gas_tau_dtemp)) deallocate(scn%op%gas_tau_dtemp)
    if (allocated(scn%op%ray_tau)) deallocate(scn%op%ray_tau)
    if (allocated(scn%op%ray_depolf)) deallocate(scn%op%ray_depolf)
    if (allocated(scn%op%total_tau)) deallocate(scn%op%total_tau)
    if (allocated(scn%op%layer_tau)) deallocate(scn%op%layer_tau)
    if (allocated(scn%op%layer_omega)) deallocate(scn%op%layer_omega)

  end subroutine destroy_optical_properties


  subroutine precompute_all_coef(scn, SV, n_stokes, n_derivs, constant_coef, coef, lcoef)

    type(scene), intent(in) :: scn
    type(statevector), intent(in) :: SV
    integer, intent(in) :: n_stokes
    integer, intent(in) :: n_derivs
    logical, intent(in) :: constant_coef
    double precision, allocatable, intent(inout) :: coef(:,:,:,:)
    double precision, allocatable, intent(inout) :: lcoef(:,:,:,:,:)


    character(len=*), parameter :: fname = "precompute_all_coef"
    double precision, allocatable :: ray_coef(:,:,:)

    double precision, allocatable :: coef_left(:,:,:)
    double precision, allocatable :: coef_right(:,:,:)
    double precision, allocatable :: lcoef_left(:,:,:,:)
    double precision, allocatable :: lcoef_right(:,:,:,:)

    integer :: aer_idx

    double precision, allocatable :: aer_sca_left(:,:), aer_sca_right(:,:)
    double precision, allocatable :: aer_ext_left(:,:), aer_ext_right(:,:)
    double precision :: denom
    double precision :: left_wl, right_wl
    double precision :: wl
    double precision :: fac
    integer :: n_mom
    integer :: n_coefs
    integer :: i, l, p, a
    integer :: n_layers

    call logger%debug(fname, "FUNCTION START")

    n_layers = scn%num_active_levels - 1

    ! Depending on whether we use polarization or not,
    ! we only need to do a certain number of phase matrix
    ! elements.

    n_mom = -1
    if (n_stokes == 1) then
       call logger%debug(fname, "Setting number of PFmom elements to 1")
       n_mom = 1
    else if (n_stokes == 3) then
       call logger%debug(fname, "Setting number of PFmom elements to 6")
       n_mom = 6
    else
       call logger%fatal(fname, "Number of stokes elements is neither 1 or 3!")
       stop 1
    end if

    allocate(ray_coef(3, n_mom, size(scn%op%wl)))

    if (scn%num_aerosols > 0) then
       call logger%debug(fname, "We have aerosols - allocating scattering optical depths.")
       left_wl = MCS%aerosol(1)%wavelengths(scn%op%aer_wl_idx_l(1))
       right_wl = MCS%aerosol(1)%wavelengths(scn%op%aer_wl_idx_r(1))

       allocate(aer_sca_left(n_layers, scn%num_aerosols))
       allocate(aer_sca_right(n_layers, scn%num_aerosols))

       allocate(aer_ext_left(n_layers, scn%num_aerosols))
       allocate(aer_ext_right(n_layers, scn%num_aerosols))

       aer_sca_left(:,:) = scn%op%aer_sca_tau(1,1:n_layers,:)
       aer_sca_right(:,:) = scn%op%aer_sca_tau(size(scn%op%wl),1:n_layers,:)

       aer_ext_left(:,:) = scn%op%aer_ext_tau(1,1:n_layers,:)
       aer_ext_right(:,:) = scn%op%aer_ext_tau(size(scn%op%wl),1:n_layers,:)

    else
       call logger%debug(fname, "No aerosols present. Just grabbing band edge wavelengths.")
       left_wl = scn%op%wl(1)
       right_wl = scn%op%wl(size(scn%op%wl))
    end if

    call logger%debug(fname, "Calculating coefficients at left edge")
    call compute_coef_at_wl(scn, SV, left_wl, n_mom, n_derivs, &
         scn%op%ray_tau(1,:), aer_sca_left, aer_ext_left, coef_left, lcoef_left)

    call logger%debug(fname, "Calculating coefficients at right edge")
    call compute_coef_at_wl(scn, SV, right_wl, n_mom, n_derivs, &
         scn%op%ray_tau(size(scn%op%wl),:), aer_sca_right, aer_ext_right, coef_right, lcoef_right)

    if (constant_coef) then
       n_coefs = 1
    else
       n_coefs = size(scn%op%wl)
    end if

    ! This array might turn out to be huge, and the code will crash
    ! if it is too large for the operating system to handle.

    allocate(coef(size(coef_left, 1), &
         size(coef_left, 2), size(coef_left, 3), n_coefs))
    allocate(lcoef(size(lcoef_left, 1), size(lcoef_left, 2), &
         size(lcoef_left, 3), size(lcoef_left, 4), n_coefs))

    call logger%debug(fname, "Setting coef and lcoef arrays to zero.")
    coef(:,:,:,:) = 0.0d0
    lcoef(:,:,:,:,:) = 0.0d0


    ! Assume the scattering properties from the middle of the reference wavelengths
    if (constant_coef) then
       call logger%debug(fname, "Using constant scattering coefficients for the band.")
       fac = 0.5
       coef(:,:,:,1) = (1.0d0 - fac) * coef_left(:,:,:) + fac * coef_right(:,:,:)
       lcoef(:,:,:,:,1) = (1.0d0 - fac) * lcoef_left(:,:,:,:) + fac * lcoef_right(:,:,:,:)
    else
       ! NOTE
       ! This here is a heinously slow operation

       call logger%debug(fname, "Creating big coef and lcoef arrays (SLOW).")

       do i = 1, size(scn%op%wl)
          wl = scn%op%wl(i)

          fac = (wl - left_wl)  / (right_wl - left_wl)

          coef(:,:,:,i) = (1.0d0 - fac) * coef_left(:,:,:) + fac * coef_right(:,:,:)
          lcoef(:,:,:,:,i) = (1.0d0 - fac) * lcoef_left(:,:,:,:) + fac * lcoef_right(:,:,:,:)

       end do
    end if

    call logger%debug(fname, "FUNCTION END")

  end subroutine precompute_all_coef


  subroutine compute_coef_at_wl(scn, SV, wl, n_mom, n_derivs, ray_tau, &
       aer_sca_tau, aer_ext_tau, coef, lcoef)

    type(scene), intent(in) :: scn
    type(statevector), intent(in) :: SV
    double precision, intent(in) :: wl
    integer, intent(in) :: n_mom ! 1 for scalar, 6 for vector
    integer, intent(in) :: n_derivs
    double precision, intent(in) :: ray_tau(:)
    double precision, intent(in) :: aer_sca_tau(:,:)
    double precision, intent(in) :: aer_ext_tau(:,:)
    double precision, allocatable, intent(inout) :: coef(:,:,:)
    double precision, allocatable, intent(inout) :: lcoef(:,:,:,:)

    integer :: n_aer
    integer :: n_layer
    integer :: n_pfmom
    integer :: a, l, p, i, k
    integer :: aer_idx
    integer :: aer_sv_idx
    integer :: l_aero_idx

    double precision, allocatable :: aer_fac(:)

    double precision :: fac
    double precision :: denom
    double precision, allocatable :: aer_height_sum(:,:)
    double precision :: this_aero_height
    double precision, allocatable :: ray_coef(:,:)
    double precision, allocatable :: aerpmom(:,:,:,:)

    n_layer = scn%num_active_levels - 1
    n_aer = scn%num_aerosols
    ! Set number of phase function coefficients to either 3 (Rayleigh only)
    ! or whatever number we need if we have aerosols
    n_pfmom = max(3, scn%max_pfmom)

    allocate(ray_coef(3, n_mom))
    allocate(aerpmom(n_aer, n_pfmom, n_mom, n_layer))
    allocate(coef(n_pfmom, n_mom, n_layer))
    allocate(lcoef(n_pfmom, n_mom, n_derivs, n_layer))

    coef(:,:,:) = 0.0d0
    lcoef(:,:,:,:) = 0.0d0

    call calculate_rayleigh_scatt_matrix(&
         calculate_rayleigh_depolf(wl), ray_coef(:,:))

    ! Remember, the total phasefunction coefficients are calculated as follows
    ! beta = (beta_aer * tau_aer_sca + beta_ray * tau_ray) / (tau_aer_sca + tau_ray)

    ! Loop over every type of aerosol used
    do l = 1, n_layer

       ! The denominator is the sum of aerosol and Rayleigh scattering optical depth
       denom = ray_tau(l)
       if (n_aer > 0) denom = denom + sum(aer_sca_tau(l, :))

       do a = 1, n_aer

          aer_idx = scn%op%aer_mcs_map(a)

          ! Calculate the "interpolation" factor needed for coef interpolation
          ! Using this here makes the subroutine a bit more general and useful
          fac = (wl - MCS%aerosol(aer_idx)%wavelengths(scn%op%aer_wl_idx_l(a))) &
               / (MCS%aerosol(aer_idx)%wavelengths(scn%op%aer_wl_idx_r(a)) &
               - MCS%aerosol(aer_idx)%wavelengths(scn%op%aer_wl_idx_l(a)) )

          ! This gives us the phase function moments from the mom file, but interpolated
          ! AT the wavelength requested by the user. If this is used to calculate the
          ! total moments at the wavelengths that are given by the mom file itself, "fac"
          ! should be either 0 or 1, and essentially just use the numbers as they are
          ! stored in the file.

          do p = 1, n_mom
             aerpmom(a,1:MCS%aerosol(aer_idx)%max_n_coef,p,l) = &
               (1.0d0 - fac) * MCS%aerosol(aer_idx)%coef(:,p,scn%op%aer_wl_idx_l(a)) &
               + fac * MCS%aerosol(aer_idx)%coef(:,p,scn%op%aer_wl_idx_r(a))

             ! Add aerosol contributions here
             coef(:, p, l) = coef(:, p, l) + aerpmom(a, :, p, l) * aer_sca_tau(l, a)
          end do

       end do
       
       ! And add Rayleigh contributions here (after partial aerosol sum)
       coef(1:3, :, l) = coef(1:3, :, l) + ray_coef(:, :) * ray_tau(l)
       ! and divide the entire layer-moments by the denominator
       coef(:, :, l) = coef(:, :, l) / denom

       ! Now that we have beta, we can 'simply' calculate the derivative inputs
       ! needed by the RT model(s).

       ! Aerosol AODs
       do i = 1, SV%num_aerosol_aod

          ! The position of the AOD derivative inputs must essentially
          ! match the SV structure, and we are generally using this ordering
          ! TODO: this is a bit hacky, would be nice to have some global
          ! dictionary where one can look these up

          l_aero_idx = SV%num_gas + SV%num_temp + SV%num_albedo + i
          ! Which aerosol belongs to SV index 'i'?
          aer_sv_idx = SV%aerosol_aod_idx_lookup(i)
          ! What is the corresponding aerosol in the MCS?
          aer_idx = scn%op%aer_mcs_map(aer_sv_idx)

          ! And calculate dBeta/dAOD
          lcoef(:,:,l_aero_idx,l) = aer_sca_tau(l, aer_sv_idx) / scn%op%reference_aod(aer_sv_idx) * &
               (aerpmom(aer_sv_idx,:,:,l) - coef(:,:,l)) / (aer_sca_tau(l, aer_sv_idx) + ray_tau(l))

       end do

       ! Aerosol heights
       do i = 1, SV%num_aerosol_height

          allocate(aer_fac(n_layer))

          ! The position of the AOD derivative inputs must essentially
          ! match the SV structure, and we are generally using this ordering
          ! TODO: this is a bit hacky, would be nice to have some global
          ! dictionary where one can look these up

          l_aero_idx = SV%num_gas + SV%num_temp + SV%num_albedo + SV%num_aerosol_aod + i
          ! Which aerosol belongs to SV index 'i'?
          aer_sv_idx = SV%aerosol_height_idx_lookup(i)
          ! What is the corresponding aerosol in the MCS?
          aer_idx = scn%op%aer_mcs_map(aer_sv_idx)

          this_aero_height = exp(SV%svsv(SV%idx_aerosol_height(i))) * scn%atm%psurf

          call calculate_aero_height_factors( &
               scn%atm%p_layers(1:n_layer), &
               this_aero_height, &
               MCS%aerosol(scn%op%aer_mcs_map(aer_sv_idx))%default_width, &
               aer_fac)

          ! And calculate dBeta/dAerosolHeight for layer l
          lcoef(:,:,l_aero_idx,l) = aer_sca_tau(l, aer_sv_idx) * aer_fac(l) * &
               (aerpmom(aer_sv_idx,:,:,l) - coef(:,:,l)) / (aer_sca_tau(l, aer_sv_idx) + ray_tau(l))

          deallocate(aer_fac)

       end do


    end do

  end subroutine compute_coef_at_wl



  subroutine calculate_active_levels(scn)

    type(scene), intent(inout) :: scn

    ! Value to reduce surface pressure if too close to layer boundary
    double precision, parameter :: psurf_bump = 1e-2
    character(len=99) :: tmp_str
    integer :: i

    ! Counting from TOA down to BOA, the layer for which
    ! psurf is located in, is the last layer - hence the
    ! lower boundary of that layer is considered the last
    ! level and thus sets the number of active levels.

    do i = 1, scn%num_levels - 1
       if ((scn%atm%psurf > scn%atm%p(i)) .and. &
            (scn%atm%psurf < scn%atm%p(i+1))) then

          scn%num_active_levels = i + 1
          exit
       end if

       ! Fringe case:
       ! In the event that the surface pressure is too close to the
       ! as a certain pressure level, it will cause problems later on.
       ! So if we find the surface pressure is too close to a pressure
       ! level, we reduce the surface pressure by a small amount to
       ! avoid the problem.

       if (abs(scn%atm%psurf - scn%atm%p(i+1)) < 1d-5) then

          write(tmp_str, '(A, F14.4, A, F14.4)') "Surface pressure bumped up: ", &
          scn%atm%psurf, " -> ", scn%atm%psurf - psurf_bump

          call logger%debug("calculate_active_levels", trim(tmp_str))

          scn%atm%psurf = scn%atm%psurf - psurf_bump
          scn%num_active_levels = i + 1

          exit
       end if

    end do

  end subroutine calculate_active_levels

  !> @brief Calculate altitude levels, gravity levels and ndry
  !> @detail This code is borrowed from MS3 (O'Brien, O'Dell et al.)
  !> @param scn Scene object
  subroutine scene_altitude(scn)
    type(scene), intent(inout) :: scn

    ! Local stuff
    double precision :: SH_layer, g_layer, p_layer, T_layer, Tv, dP
    double precision :: logratio, dz, constant
    ! Loop variable
    integer :: i

    scn%atm%altitude_layers(:) = 0.0d0
    ! Set altitudes to zero, and the lowest level to the altitude
    scn%atm%altitude_levels(:) = 0.0d0
    scn%atm%altitude_levels(scn%num_levels) = scn%alt
    scn%atm%grav(scn%num_levels) = jpl_gravity(scn%lat, scn%alt)
    scn%atm%ndry(:) = 0.0d0

    ! Loop through layers, starting with the bottom-most (surface) one
    do i = scn%num_levels - 1, 1, -1

       SH_layer = (scn%atm%sh(i) + scn%atm%sh(i+1)) * 0.5d0
       g_layer = jpl_gravity(scn%lat, scn%atm%altitude_levels(i+1))
       p_layer = (scn%atm%p(i) + scn%atm%p(i+1)) * 0.5d0
       dP = scn%atm%p(i+1) - scn%atm%p(i)
       T_layer = (scn%atm%T(i) + scn%atm%T(i+1)) * 0.5d0
       Tv = T_layer * (1.0d0 + SH_layer * (1.0d0 - EPSILON) / EPSILON)
       logratio = log(scn%atm%p(i+1) / scn%atm%p(i))
       dz = logratio * Tv * Rd / g_layer
       g_layer = jpl_gravity(scn%lat, scn%atm%altitude_levels(i+1) + 0.5d0 * dz)
       dz = logratio * Tv * Rd / g_layer
       constant = dP / (DRY_AIR_MASS * g_layer)

       ! Write the important stuff back into the scene object
       !-----------------------------------------------------
       scn%atm%ndry(i) = constant * (1.0d0 - SH_layer)
       scn%atm%altitude_levels(i) = scn%atm%altitude_levels(i+1) + dz
       scn%atm%grav(i) = jpl_gravity(scn%lat, scn%atm%altitude_levels(i))
       !-----------------------------------------------------

    end do

    ! Some calculations want the layer altitude, so might as well compute them
    ! here and store them.
    do i = 1, scn%num_levels - 1
       scn%atm%altitude_layers(i) = 0.5d0 * (scn%atm%altitude_levels(i) + scn%atm%altitude_levels(i+1))
    end do


  end subroutine scene_altitude

  !> @brief Calculate acceleration g(lat, altitude)
  !> @param gdlat Geodetic latitude [deg]
  !> @param altit Geometric altitude [m]
  pure function jpl_gravity(gdlat,altit) result(gravity)
    ! NOTE
    ! This code was taken from MS3 and just modified slightly for data
    ! types.
    !
    ! Computes the effective Earth gravity at a given latitude and altitude.
    ! This is the sum of the gravitational and centripital accelerations.
    ! These are based on equation I.2.4-(17) in US Standard Atmosphere 1962
    ! The Earth is assumed to be an oblate ellipsoid, with a ratio of the
    ! major to minor axes = sqrt(1+con) where con=.006738
    ! This eccentricity makes the Earth's gravititational field smaller at the
    ! poles and larger at the equator than if the Earth were a sphere of the
    ! same mass. It also makes the local mid-latitude gravity field not point
    ! toward the center of mass.
    !
    ! Input Parameters:
    !   gdlat       GeoDetric Latitude (degrees)
    !   altit       Geometric Altitude (meters) ! CHANGED CWO 6/22/2009
    !
    ! Output Parameter:
    !   gravity     Effective Gravitational Acceleration (m/s2)
    !
    ! Interestingly, since the centripital effect of the Earth's rotation
    ! (-ve at equator, 0 at poles) has almost the opposite shape to the
    ! second order gravitational field (+ve at equator, -ve at poles), their
    ! sum is almost constant so that the surface gravity can be approximated
    ! (to .07%) by the simple expression g = 0.99746*GM/radius**2, the latitude
    ! variation coming entirely from the variation of surface r with latitude.

    implicit none

    ! In/Out variables
    double precision, intent(in) :: gdlat   ! geodetic latitude [degrees]
    double precision, intent(in) :: altit   ! geometric altitude [meters]
    double precision :: gravity ! gravitational acceleration [m/s^2]

    ! Local Variables
    double precision :: radius        ! radial distance (metres)
    double precision :: gclat         ! geocentric latitude [radians]
    double precision :: ff, hh, ge       ! scratch variables

    ! Parameters
    ! Gravitational constant times Earth's Mass (m3/s2)
    double precision, parameter  :: gm = 3.9862216d14
    ! Earth's angular rotational velocity (radians/s)
    double precision, parameter  :: omega = 7.292116d-5
    ! (a/b)**2-1 where a & b are equatorial & polar radii
    double precision, parameter  :: con = 0.006738d0
    ! 2nd harmonic coefficient of Earth's gravity field
    double precision, parameter  :: shc = 1.6235d-3

    ! Convert from geodetic latitude (GDLAT) to geocentric latitude (GCLAT).
    gclat = atan(tan(DEG2RAD*gdlat)/(1+con))  ! radians
    ! On computers which crash at the poles try the following expression
    ! gclat=d2r*gdlat-con*sin(d2r*gdlat)*cos(d2r*gdlat)/(1+con*cos(d2r*gdlat)**2)
    radius = altit + EARTH_EQUATORIAL_RADIUS/sqrt(1.0d0+con*sin(gclat)**2)
    ff = (radius/ EARTH_EQUATORIAL_RADIUS)**2
    hh = radius*omega**2
    ge = gm / EARTH_EQUATORIAL_RADIUS**2                       ! = gravity at Re
    gravity = (ge*(1-shc*(3.0d0*sin(gclat)**2-1.0d0)/ff)/ff-hh*cos(gclat)**2) &
         *(1.0d0+0.5d0*(sin(gclat)*cos(gclat)*(hh/ge+2.0d0*shc/ff**2))**2)

  end function jpl_gravity


  !> @brief Reads in the atmosphere file, which contains column-based profiles.
  !>
  !> @param filename Path to the atmosphere file
  !> @param Array of 'strings' containing names of required gases
  !> @param Atmosphere object that will be populated by this function
  !>
  !> @detail We supply here a filename, a list of gas strings and the
  !> atmosphere-object. The function checks whether all required gases
  !> are present in the atmosphere file, and will throw a fatal error if
  !> that is not the case.
  subroutine read_atmosphere_file(filename, gas_strings, atm)

    implicit none
    character(len=*), intent(in) :: filename
    type(string), intent(in) :: gas_strings(:)
    type(atmosphere), intent(inout) :: atm

    ! Function name
    character(len=*), parameter :: fname = "read_atmosphere_file"
    ! File handler and IO stat variable
    integer :: funit, iostat
    ! Whether file exists
    logical :: file_exist
    ! Various counting variables to figure out where the
    ! contents of the atmosphere file start.
    integer :: line_count, nonempty_line_count, file_start, level_count
    ! Indices which tell us which columns are pressure and temperature
    integer :: idx_p, idx_t
    ! Character variables to store lines
    character(len=999) :: dummy, tmp_str
    ! This dummy is for reading in numerical values
    double precision :: dummy_dp
    ! This dummy is for reading in strings
    type(string) :: dummy_string
    ! Lines will be split using this dummy variable
    type(string), allocatable :: split_string(:)
    ! Various loop counters
    integer :: i, j, cnt
    ! Too keep track of whether we have found the required gases, and
    ! at which position/column they are.
    integer, allocatable :: this_gas_index(:)

    integer :: num_gases

    ! Check whether the file exists or not.
    inquire(file=filename, exist=file_exist)
    if (.not. file_exist) then
       call logger%fatal(fname, "Atmosphere file does not exist: " // filename)
       stop 1
    else
       call logger%debug(fname, "File does exist.")
    end if

    ! First pass: we scan the file to see how many levels our atmosphere has
    open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')
    rewind(unit=funit, iostat=iostat)

    line_count = 0
    level_count = 0
    nonempty_line_count = 0
    file_start = -1

    ! Loop through the file until we have reached the end
    do
       read(funit, '(A)', iostat=iostat) dummy

       if (iostat == iostat_end) then
          ! End of file?
          exit
       end if

       ! Keep track of how many lines we have traversed
       line_count = line_count + 1

       if (scan(dummy, "!#;") > 0) then
          ! Skip this line, as it's commented
          cycle
       else if (trim(dummy) == "") then
          ! Skip empty lines
          cycle
       end if

       ! And keep track of now many lines are non-empty
       nonempty_line_count = nonempty_line_count + 1

       ! First non-empty line should be saved here.
       ! file_start is where the real file contents begin
       if (file_start == -1) then
          file_start = line_count
       else
          if (nonempty_line_count > 1) then
             level_count = level_count + 1
          end if
       end if

    end do

    ! Go back to the top of the file.
    rewind(unit=funit, iostat=iostat)

    idx_p = -1
    idx_t = -1

    ! Reset the line counter since we're starting from the top again.
    line_count = 0
    do
       ! Read the line
       read(funit, '(A)', iostat=iostat) dummy
       line_count = line_count + 1

       ! .. and immediately skip until we are at the
       ! beginning of the contents of the file.
       if (line_count < file_start) cycle

       if (line_count == file_start) then
          ! This is the proper atmosphere header that should contain
          ! the information about the gases. So first, let's check if
          ! the numbers match

          ! Read the line into a string object and split it by whitespaces
          dummy_string = dummy
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now that we know both the number of levels and gases, we can allocate the
          ! arrays in the atmosphere structure.
          num_gases = 0
          idx_p = -1
          idx_t = -1
          do j=1, size(split_string)
             ! Skip temp or pressure - this requires the pressure and temperature
             ! columns to be explicitly labeled p/P and t/T
             if (split_string(j)%lower() == "p") then
                idx_p = j
                cycle
             end if
             if (split_string(j)%lower() == "t") then
                idx_t = j
                cycle
             end if
             num_gases = num_gases + 1
          end do

          ! Let the user know how many gases and levels we have found
          write(tmp_str, '(A,G0.1,A,A)') "There seem to be ", num_gases, " gases in ", filename
          call logger%info(fname, trim(tmp_str))
          write(tmp_str, '(A, G0.1)') "The number of atmospheric levels is: ", level_count
          call logger%info(fname, trim(tmp_str))

          ! If the atmosphere structure was already allocated, deallocate it first
          if (allocated(atm%p)) deallocate(atm%p)
          if (allocated(atm%T)) deallocate(atm%T)
          if (allocated(atm%sh)) deallocate(atm%sh)
          if (allocated(atm%gas_names)) deallocate(atm%gas_names)
          if (allocated(atm%gas_index)) deallocate(atm%gas_index)
          if (allocated(atm%gas_vmr)) deallocate(atm%gas_vmr)
          if (allocated(atm%altitude_levels)) deallocate(atm%altitude_levels)
          if (allocated(atm%altitude_layers)) deallocate(atm%altitude_layers)
          if (allocated(atm%grav)) deallocate(atm%grav)
          if (allocated(atm%ndry)) deallocate(atm%ndry)

          ! Allocate according to the file structure
          atm%num_levels = level_count
          atm%num_gases = num_gases

          allocate(atm%T(level_count))
          allocate(atm%p(level_count))
          allocate(atm%sh(level_count))
          allocate(atm%gas_names(num_gases))
          allocate(atm%gas_vmr(level_count, num_gases))
          allocate(atm%gas_index(num_gases))
          allocate(atm%altitude_levels(level_count))
          allocate(atm%altitude_layers(level_count - 1))
          allocate(atm%grav(level_count))
          allocate(atm%ndry(level_count))

          allocate(this_gas_index(size(gas_strings)))

          ! But we also want to know what gas index to use for storage
          do i=1, size(gas_strings)
             this_gas_index(i) = -1
             cnt = 1
             do j=1, size(split_string)
                ! Skip temp or pressure
                if (split_string(j)%lower() == "p") cycle
                if (split_string(j)%lower() == "t") cycle
                ! Check if gas description matches gases we know
                if (split_string(j) == gas_strings(i)) then
                   this_gas_index(i) = j
                   atm%gas_index(i) = j
                   atm%gas_names(i) = split_string(j)%lower()

                   write(tmp_str, '(A,A,A,G0.1)') "Index for atmosphere gas ", &
                        split_string(j)%chars(), ": ", j
                   call logger%debug(fname, trim(tmp_str))
                   exit
                end if
                cnt = cnt + 1
             end do
          end do

          ! And last check - do all required gases have a VMR column in the
          ! atmosphere file.

          do i=1, size(gas_strings)
             if (this_gas_index(i) == -1) then
                ! Uh-oh, gas that was speficied in the "gases" option of the window
                ! could not be found in the atmosphere! Hard exit.
                write(tmp_str, "(A,A)") "The following gas was not found in the atmosphere file: " &
                     // gas_strings(i)%chars()
                call logger%fatal(fname, trim(tmp_str))
                stop 1
             end if
          end do

       end if


       if (line_count > file_start) then
          ! Right after the header, we should have the data in rows. We
          ! use the string split option to split the row string into substrings,
          ! and then convert each into double precision values and feed them into
          ! the atmosphere structure - at the right position!

          dummy_string = dummy
          ! Need to deallocate the split_string object first
          if (allocated(split_string)) deallocate(split_string)
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now here we need to check again whether a certain line has more than
          ! num_gases+1 columns.
          if (size(split_string) /= (num_gases + 2)) then
             write(tmp_str, '(A, G0.1)') "Too many values in line ", line_count
             call logger%fatal(fname, trim(tmp_str))
             stop 1
          end if

          ! Get the pressure value
          tmp_str = split_string(idx_p)%chars()
          read(tmp_str, *) dummy_dp
          atm%p(line_count - file_start) = dummy_dp

          ! Get the temperature value
          tmp_str = split_string(idx_t)%chars()
          read(tmp_str, *) dummy_dp
          atm%T(line_count - file_start) = dummy_dp

          ! And the gas value(s) from the other column(s)
          do i=1, size(gas_strings)
             tmp_str = split_string(this_gas_index(i))%chars()
             read(tmp_str, *) dummy_dp
             atm%gas_vmr(line_count - file_start, i) = dummy_dp
          end do

       end if

       if (line_count == (file_start + level_count)) exit

    end do

    close(funit)

  end subroutine read_atmosphere_file

  !> @brief Creates / allocates the "results" container to hold all retrieval results
  !>
  !> @param results Result container
  !> @param num_frames Number of frames
  !> @param num_fp Number of footprints
  !> @param num_SV Number of state vector elements
  !> @param num_gas Number of retrieved (!) gases
  subroutine create_result_container(results, num_frames, num_fp, num_SV, num_gas, num_level)
    implicit none
    type(result_container), intent(inout) :: results
    integer, intent(in) :: num_frames, num_fp, num_SV, num_gas, num_level

    allocate(results%sv_names(num_SV))

    allocate(results%sv_retrieved(num_fp, num_frames, num_SV))
    allocate(results%sv_prior(num_fp, num_frames, num_SV))
    allocate(results%sv_uncertainty(num_fp, num_frames, num_SV))
    allocate(results%xgas(num_fp, num_frames, num_gas))
    allocate(results%xgas_prior(num_fp, num_frames, num_gas))
    allocate(results%pwgts(num_fp, num_frames, num_gas, num_level))
    allocate(results%col_ak(num_fp, num_frames, num_gas, num_level))
    allocate(results%pressure_levels(num_fp, num_frames, num_level))
    allocate(results%vmr_prior(num_fp, num_frames, num_gas, num_level))
    allocate(results%vmr_retrieved(num_fp, num_frames, num_gas, num_level))
    allocate(results%chi2(num_fp, num_frames))
    allocate(results%residual_rms(num_fp, num_frames))
    allocate(results%dsigma_sq(num_fp, num_frames))
    allocate(results%SNR(num_fp, num_frames))
    allocate(results%SNR_std(num_fp, num_frames))
    allocate(results%continuum(num_fp, num_frames))
    allocate(results%processing_time(num_fp, num_frames))

    allocate(results%num_iterations(num_fp, num_frames))
    allocate(results%converged(num_fp, num_frames))

    allocate(results%ndry(num_fp, num_frames, num_gas))

    results%sv_names = "NONE"

    ! This might cause problems for some, but I find it convenient
    ! to set 'unused' fields to NaNs. Remember this only works for
    ! reals/double precision, but not for integers.
    results%sv_retrieved = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%sv_uncertainty = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%xgas = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%xgas_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%pwgts = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%vmr_prior = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%vmr_retrieved = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%pressure_levels = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%col_ak = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%chi2 = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%residual_rms = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%dsigma_sq = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%SNR = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%SNR_std = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%continuum = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    results%processing_time = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

    results%num_iterations = -1
    results%converged = -1

  end subroutine create_result_container


  !> @brief Destroys the "results" container allocated in "create_result_container"
  !>
  !> @param results Result container
  subroutine destroy_result_container(results)
    implicit none
    type(result_container), intent(inout) :: results

    deallocate(results%sv_names)
    deallocate(results%sv_retrieved)
    deallocate(results%sv_prior)
    deallocate(results%sv_uncertainty)
    deallocate(results%xgas)
    deallocate(results%xgas_prior)
    deallocate(results%pwgts)
    deallocate(results%vmr_prior)
    deallocate(results%vmr_retrieved)
    deallocate(results%pressure_levels)
    deallocate(results%col_ak)
    deallocate(results%chi2)
    deallocate(results%residual_rms)
    deallocate(results%dsigma_sq)
    deallocate(results%SNR)
    deallocate(results%SNR_std)
    deallocate(results%continuum)
    deallocate(results%num_iterations)
    deallocate(results%converged)
    deallocate(results%ndry)
    deallocate(results%processing_time)

  end subroutine destroy_result_container


  !> @brief Creates human-readable names for state vector elements
  !>
  !> The idea is faily simple: we loop through all the state vector elements,
  !> and then check for each one if there is a corresponding SV\%idx_* associated with
  !> that element position. Based on that, we create a name for the state vector
  !> element, which usually has the parameter number (e.g. albedo order) baked in.
  !> @param results Result container
  !> @param SV State vector object
  !> @param i_win Retrieval window index for MCS
  subroutine assign_SV_names_to_result(results, SV, i_win)
    implicit none

    type(result_container), intent(inout) :: results
    type(statevector), intent(in) :: SV
    integer, intent(in) :: i_win

    type(string) :: lower_str
    character(len=999) :: tmp_str
    integer :: i,j,k

    i = 1
    do while (i <= size(SV%svsv))

       ! Albedo names
       do j=1, SV%num_albedo
          if (SV%idx_albedo(j) == i) then
             write(tmp_str, '(A,G0.1)') "albedo_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! SIF name (really only one at this point)
       do j=1, SV%num_sif
          if (SV%idx_sif(j) == i) then
             write(tmp_str, '(A)') "sif_radiance"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! ZLO name
       do j=1, SV%num_zlo
          if (SV%idx_zlo(j) == i) then
             write(tmp_str, '(A)') "zero_level_offset"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       do j=1, SV%num_temp
          if (SV%idx_temp(j) == i) then
             write(tmp_str, '(A)') "temperature_offset"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Solar shift name
       if (SV%idx_solar_shift(1) == i) then
          write(tmp_str, '(A,A)') "solar_shift"
          results%sv_names(i) = trim(tmp_str)
       end if

       ! Solar stretch name
       if (SV%idx_solar_stretch(1) == i) then
          write(tmp_str, '(A)') "solar_stretch"
          results%sv_names(i) = trim(tmp_str)
       end if

       ! Surface pressure name
       do j=1, SV%num_psurf
          if (SV%idx_psurf(j) == i) then
             write(tmp_str, '(A)') "surface_pressure"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Dispersion parameter names
       do j=1, SV%num_dispersion
          if (SV%idx_dispersion(j) == i) then
             write(tmp_str, '(A,G0.1)') "dispersion_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! ILS parameter names
       do j=1, SV%num_ils_stretch
          if (SV%idx_ils_stretch(j) == i) then
             write(tmp_str, '(A,G0.1)') "ils_stretch_order_", j-1
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Retrieved aerosol AOD names
       do j=1, SV%num_aerosol_aod
          if (SV%idx_aerosol_aod(j) == i) then
             lower_str = MCS%window(i_win)%aerosols(sv%aerosol_aod_idx_lookup(j))%lower()
             write(tmp_str, '(A,A)') lower_str%chars(), "_aod"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Retrieved aerosol AOD names
       do j=1, SV%num_aerosol_height
          if (SV%idx_aerosol_height(j) == i) then
             lower_str = MCS%window(i_win)%aerosols(sv%aerosol_height_idx_lookup(j))%lower()
             write(tmp_str, '(A,A)') lower_str%chars(), "_height"
             results%sv_names(i) = trim(tmp_str)
          end if
       end do

       ! Retrieved gas names (scalar retrieval only so far)
       k = 1
       do j=1, SV%num_gas
          ! Check if this SV element is a scalar retrieval
          if (SV%idx_gas(j,1) == i) then
             if (MCS%window(i_win)%gas_retrieve_scale(sv%gas_idx_lookup(j))) then
                if (MCS%window(i_win)%gas_retrieve_scale_start(sv%gas_idx_lookup(j), k) == -1.0) cycle

                lower_str = MCS%window(i_win)%gases(sv%gas_idx_lookup(j))%lower()
                write(tmp_str, '(A)') trim(lower_str%chars() // "_scale_")
                write(tmp_str, '(A, F4.2)') trim(tmp_str), &
                     SV%gas_retrieve_scale_start(j)
                write(tmp_str, '(A,A,F4.2)') trim(tmp_str), "_" , &
                     SV%gas_retrieve_scale_stop(j)
                results%sv_names(i) = trim(tmp_str)

                k = k + 1
             end if
          end if
       end do

       i = i+1
    end do

  end subroutine assign_SV_names_to_result

  !> @begin Wrapper to replace prior VMRs with special functions
  subroutine replace_prior_VMR(scn, prior_types)

    type(scene), intent(inout) :: scn
    type(string), intent(in) :: prior_types(:)

    ! Function name
    character(len=*), parameter :: fname = "replace_prior_VMR"
    character(len=999) :: tmp_str

    integer :: i

    ! The prior_types are in order of scn%atm%gas_vmr, so
    ! we can calculate a new prior using prior_type(i) and stick
    ! it into scn%atm%gas_vmr(:,i).

    do i=1, size(prior_types)



       ! Nothing to do if this string is empty
       if (prior_types(i) == "") cycle


       if (prior_types(i) == "SC4C2018") then
          ! Max Reuter / Oliver Schneising
          ! CH4 profiles as derived from a climatology file


       else
          ! If the prior type is not implemented, the user has made a mistake.
          ! Again, we are terminating here immediately, since falling back to
          ! some default behavior is not a good option..

          write(tmp_str, '(A,A)') "Sorry, the following prior VMR function " &
              // "is not implemented: ", prior_types(i)%chars()
          call logger%fatal(fname, trim(tmp_str))
          stop 1
       end if

    end do

  end subroutine replace_prior_VMR



end module physical_model_addon_mod
