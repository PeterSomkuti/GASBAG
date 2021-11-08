
module scene_mod

  use math_utils_mod
  ! Third-party modules
  use stringifor
  use mod_datetime

  implicit none

  !> A simple structure to keep the atmosphere data nice and tidy
  type atmosphere
     !> Number of levels in the model atmosphere
     integer :: num_levels
     !> Number of gases in the model atmosphere
     integer :: num_gases
     !> The name(s) of the gas(es), (gas number)
     type(string), allocatable :: gas_names(:)
     !> Gas mixing ratios (level, gas number)
     double precision, allocatable :: gas_vmr(:,:)
     !> To which spectroscopy data does this gas correspond to? (gas number)
     integer, allocatable :: gas_index(:)
     !> Surface pressure
     double precision :: psurf
     !> Model atmosphere temperature
     double precision, allocatable :: T(:)
     !> Model atmosphere pressure
     double precision, allocatable :: p(:)
     !> Model atmosphere mid-layer pressure
     double precision, allocatable :: p_layers(:)
     !> Model atmosphere specific humidity
     double precision, allocatable :: sh(:)
     !> Model atmosphere specific humidity at pressure layers
     double precision, allocatable :: sh_layers(:)
     !> Model atmosphere gravity at pressure levels
     double precision, allocatable :: grav(:)
     !> Model atmosphere gravity at pressure levels
     double precision, allocatable :: grav_layers(:)
     !> Model atmosphere altitude at pressure levels
     double precision, allocatable :: altitude_levels(:)
     !> Layer altitude
     double precision, allocatable :: altitude_layers(:)
     !> Model atmosphere dry air column at layers
     double precision, allocatable :: ndry(:)
     !> Level pressure weighting function for this scene
     double precision, allocatable :: pwgts(:)
     !> Layer pressure weighting function for this scene
     double precision, allocatable :: pwgts_layers(:)
     !> Layer pseudo spherical factors for solar angles
     double precision, allocatable :: ps_factors_solar(:)
     !> Layer pseudo spherical factors for viewing angles
     double precision, allocatable :: ps_factors_viewing(:)
     !> Layer Chapman factors (layer, layer)
     double precision, allocatable :: chapman_solar(:,:)
     double precision, allocatable :: chapman_viewing(:,:)
     !> Layer effective (slant) paths (1 / mu_eff)
     double precision, allocatable :: layer_slant_path_solar(:)
     double precision, allocatable :: layer_slant_path_viewing(:)

  end type atmosphere

  !> Holds optical gas/aerosol properties of the scene
  type optical_properties
     !> Wavelengths
     double precision, allocatable :: wl(:)
     !> Surface albedo
     double precision, allocatable :: albedo(:)

     !> Sphericity factor for the solar path
     double precision, allocatable :: sph_tau_factor_solar(:,:)
     !> Sphericity factor for the viewing / LOS path
     double precision, allocatable :: sph_tau_factor_viewing(:,:)
     !> Gas optical depth (spectral, layer, gas number)
     double precision, allocatable :: gas_tau(:,:,:)
     !> dtau / dtemp (spectral, layer, gas number)
     double precision, allocatable :: gas_tau_dtemp(:,:,:)
     !> dtau / dpsurf (spectral, layer, gas number)
     double precision, allocatable :: gas_tau_dpsurf(:,:,:)
     !> dtau / dvmr (spectral, layer, gas number, level-below-or-above)
     double precision, allocatable :: gas_tau_dvmr(:,:,:,:)
     !> Perturbed gas optical depth (spectral, layer, gas number)
     double precision, allocatable :: gas_tau_pert(:,:,:,:)

     !> Rayleigh extinction optical depth (spectral, layer)
     double precision, allocatable :: ray_tau(:,:)
     !> Rayleigh depolarization factor (spectral)
     double precision, allocatable :: ray_depolf(:)

     !> Total column optical depth (spectral)
     double precision, allocatable :: total_tau(:)
     ! Total layer optical properties go into the RT calculations
     !> Layer total optical depth (spectral, layer)
     double precision, allocatable :: layer_tau(:,:)
     !> Total single scatter albedo (spectral, layer)
     double precision, allocatable :: layer_omega(:,:)

     !> Map the aerosol "number" to the MCS entry (aerosol)
     ! aer_mcs_map(1) = 2 e.g. means that [aerosol-2] is the first
     ! aerosol we have in this configuration.
     integer, allocatable :: aer_mcs_map(:)
     !> Which left-wavelength index of the aerosol file are we using? (aerosol)
     integer, allocatable :: aer_wl_idx_l(:)
     !> Which right-wavelength index of the aerosol file are we using? (aerosol)
     integer, allocatable :: aer_wl_idx_r(:)
     !> Aerosol extinction cross section (spectral, aerosol)
     double precision, allocatable :: aer_ext_q(:,:)
     !> Aerosol scattering cross section (spectral, aerosol)
     double precision, allocatable :: aer_sca_q(:,:)
     !> Aerosol single scatter albedo (spectral, aerosol)
     double precision, allocatable :: aer_ssa(:,:)
     !> Aerosol interpolation factor (may be needed by PCA) (spectral)
     double precision, allocatable :: aer_frac(:)

     !> Aerosol reference AOD (aerosol)
     double precision, allocatable :: reference_aod(:)

     ! These here are calculated AFTER inferring the
     ! aerosol distribution (which layers?)

     !> Aerosol extinction optical depth (spectral, layer, aerosol)
     double precision, allocatable :: aer_ext_tau(:,:,:)
     !> Aerosol scattering optical depth (spectral, layer, aerosol)
     double precision, allocatable :: aer_sca_tau(:,:,:)

     ! These here are the same as the two above, but contain the
     ! extinction and scattering profiles at the wavelengths at
     ! which the miemom-type aerosol is actually calculated. This
     ! is needed for a quicker calculation of the phase function moments
     ! across the entire band

     !> Aerosol extinction optical depth (edge, layer, aerosol)
     double precision, allocatable :: aer_ext_tau_edge(:,:,:)
     !> Aerosol scattering optical depth (edge, layer, aerosol)
     double precision, allocatable :: aer_sca_tau_edge(:,:,:)

  end type optical_properties


  ! The "biggest" custom type that contains most per-scene
  ! information you'll ever need to run a physical retrieval.
  type, public :: scene

     !> Number of levels in the model atmosphere
     integer :: num_levels = 0
     !> Number of active levels currently used
     integer :: num_active_levels = 0
     !> Number of gases in the model atmosphere
     integer :: num_gases = 0
     !> Number of aerosols in the scene
     integer :: num_aerosols = 0
     !> Largest number of phase function moments, needed for allocation
     integer :: max_pfmom = 0
     !> Number of Stokes coefficients used in RT calculations
     integer :: num_stokes = 0

     !> The atmosphere of the scene
     type(atmosphere) :: atm
     !> Optical properties needed for RT
     type(optical_properties) :: op
     !> Date object
     type(datetime) :: date
     !> Epoch (Year, Month, Day, Hour, Minute, Second, Millisecond)
     integer :: epoch(7)
     !> Longitude
     double precision :: lon
     !> Latitude
     double precision :: lat
     !> Altitude at surface
     double precision :: alt
     !> Solar zenith angle
     double precision :: SZA
     !> Cosine of solar zenith angle
     double precision :: mu0
     !> Viewing zenith angle
     double precision :: VZA
     !> Cosine of viewing zenith angle
     double precision :: mu
     !> Solar azimuth angle
     double precision :: SAA
     !> Viewing azimuth angle
     double precision :: VAA
  end type scene


contains


  subroutine optical_depth_cleanup(scn)
    type(scene), intent(inout) :: scn

    ! Minimal value below which optical depths are clamped
    double precision :: MIN_VAL = 1d-10

    ! Set small values to some minimal value
    if (allocated(scn%op%gas_tau)) then
       where(scn%op%gas_tau < MIN_VAL) scn%op%gas_tau = MIN_VAL
    end if

    if (allocated(scn%op%ray_tau)) then
       where(scn%op%ray_tau < MIN_VAL) scn%op%ray_tau = MIN_VAL
    end if

    ! Total optical depth is calculated as sum of all gas ODs plus Rayleigh OD
    scn%op%layer_tau(:,:) = sum(scn%op%gas_tau, dim=3) + scn%op%ray_tau

    ! If there are aerosols in the scene, add them to the total OD

    if (allocated(scn%op%aer_ext_tau)) then
       scn%op%layer_tau(:,:) = scn%op%layer_tau(:,:) &
            + sum(scn%op%aer_ext_tau(:,:,:), dim=3)
    end if

    ! Total column optical depth is just the sum over all active layers
    ! (per wavelength still)
    scn%op%total_tau(:) = sum(scn%op%layer_tau(:, 1:scn%num_active_levels - 1), dim=2)

    ! The layer-resolved single scatter albedo is (Rayleigh + Aerosol) / (Total)
    ! extinctions.
    if (allocated(scn%op%aer_ext_tau)) then
       scn%op%layer_omega(:,:) = &
            (scn%op%ray_tau + sum(scn%op%aer_ext_tau(:,:,:), dim=3)) &
            / scn%op%layer_tau(:,:)
    else
       ! Or without aerosols: Rayleigh / (Rayleigh + Gas)
       scn%op%layer_omega(:,:) = scn%op%ray_tau(:,:) / scn%op%layer_tau(:,:)

    end if


  end subroutine optical_depth_cleanup


  !> @brief Calculates mid-layer pressures and humidities,
  !>        taking into account surface pressure
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
    scn%atm%sh_layers(:) = -1.0d0

    ! Loop over all layers and compute mid-layer pressure
    do l = 1, scn%num_active_levels - 1
       scn%atm%sh_layers(l) = 0.5d0 * (scn%atm%sh(l) + scn%atm%sh(l+1))
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
    allocate(scn%op%aer_frac(N_hires))

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
       allocate(scn%op%gas_tau_dvmr(2, N_hires, scn%num_levels-1, N_gases))
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
    deallocate(scn%op%aer_frac)

    if (allocated(scn%op%sph_tau_factor_solar)) deallocate(scn%op%sph_tau_factor_solar)
    if (allocated(scn%op%sph_tau_factor_viewing)) deallocate(scn%op%sph_tau_factor_viewing)
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
    double precision, intent(in) :: gdlat ! geodetic latitude [degrees]
    double precision, intent(in) :: altit ! geometric altitude [meters]
    double precision :: gravity ! gravitational acceleration [m/s^2]

    ! Local Variables
    double precision :: radius        ! radial distance (metres)
    double precision :: gclat         ! geocentric latitude [radians]
    double precision :: ff, hh, ge    ! scratch variables

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


  !> @brief Calculate Chapman factors per layer
  !> @details This code is inspired by MS3 (O'Brien, O'Dell et al.)
  !> @param scn Scene object
  subroutine chapman_factors(scn, mu, z0, chapman)
    type(scene), intent(inout) :: scn
    double precision, intent(in) :: mu ! either mu or mu0
    double precision, intent(in) :: z0 ! Earth radius at this scene
    double precision, intent(inout) :: chapman(:,:)

    integer :: i, j
    integer :: Nlay

    double precision :: a, b
    double precision :: r1, r2

    Nlay = scn%num_active_levels - 1
    chapman(:,:) = 0.0d0


    a = sin(acos(mu))
    a = a * a

    do i = Nlay, 1, -1

       b = z0 + scn%atm%altitude_levels(i+1)
       b = b * b * a

       do j = Nlay, i, -1

          r1 = z0 + scn%atm%altitude_levels(j+1)
          r2 = z0 + scn%atm%altitude_levels(j)

          !if (r1 == r2) then
          !   chapman(i,j) = chapman(i, j-1)
          !else
             ! This calculates 1 / mueff
          chapman(i,j) = (sqrt(r1*r1 - b) - sqrt(r2*r2 - b)) / (r1 - r2)
          !end if

       end do
    end do


  end subroutine chapman_factors


  !> @brief Calculate altitude levels, gravity levels and ndry
  !> @details This code is borrowed from MS3 (O'Brien, O'Dell et al.)
  !> @param scn Scene object
  subroutine scene_altitude(scn, N_sub)
    type(scene), intent(inout) :: scn
    integer, intent(in) :: N_sub

    ! Local stuff
    double precision :: SH_layer, g_layer, p_layer, T_layer, Tv, dP, dT, dQ
    double precision :: logratio, dz, constant, zlo, plo, tlo, tbar, qlo
    double precision :: alt_tmp
    double precision :: gclat, radius
    double precision, parameter  :: con = 0.006738d0


    ! Loop variable
    integer :: i, j

    scn%atm%altitude_layers(:) = 0.0d0
    ! Set altitudes to zero, and the lowest level to the scene altitude
    scn%atm%altitude_levels(:) = 0.0d0
    scn%atm%altitude_levels(scn%num_active_levels) = scn%alt
    scn%atm%grav(scn%num_active_levels) = jpl_gravity(scn%lat, scn%alt)

    scn%atm%ndry(:) = 0.0d0

    ! Loop through layers, starting with the bottom-most (surface) one
    do i = scn%num_active_levels, 2, -1

       zlo = scn%atm%altitude_levels(i)

       if (i == scn%num_active_levels) then
         ! Surface level
          dP = (scn%atm%psurf - scn%atm%p(i)) / N_sub
       else
          ! Any other level
          dP = (scn%atm%p(i-1) - scn%atm%p(i)) / N_sub
       end if

       dT = (scn%atm%T(i-1) - scn%atm%T(i)) / N_sub
       dQ = (scn%atm%sh(i-1) - scn%atm%sh(i)) / N_sub

       do j = 0, N_sub - 1

          if (i == scn%num_active_levels - 1) then
             plo = scn%atm%psurf + j * dP
          else
             plo = scn%atm%p(i) + j * dP
          end if

          logratio = log(plo / (plo + dP))

          tlo = scn%atm%T(i) + j * dT
          tbar = tlo + 0.5 * dT
          qlo = scn%atm%sh(i) + j * dQ

          g_layer = jpl_gravity(scn%lat, zlo)

          dz = (Rd * tbar / g_layer) * logratio
          zlo = zlo + dz

          scn%atm%ndry(i-1) = scn%atm%ndry(i-1) + dP / (DRY_AIR_MASS * g_layer) * (1.0d0 - qlo)

       end do

       !g_layer = jpl_gravity(scn%lat, scn%atm%altitude_levels(i+1))
       !SH_layer = (scn%atm%sh(i) + scn%atm%sh(i+1)) * 0.5d0
       !T_layer = (scn%atm%T(i) + scn%atm%T(i+1)) * 0.5d0
       !Tv = T_layer * (1.0d0 + SH_layer * (1.0d0 - EPSILON) / EPSILON)

       !dz = logratio * Tv * Rd / g_layer
       !g_layer = jpl_gravity(scn%lat, scn%atm%altitude_levels(i+1) + 0.5d0 * dz)
       !dz = logratio * Tv * Rd / g_layer
       !constant = dP / (DRY_AIR_MASS * g_layer)

       ! Write the important stuff back into the scene object
       !-----------------------------------------------------
       !scn%atm%ndry(i-1) = constant * (1.0d0 - SH_layer)
       !scn%atm%altitude_levels(i) = scn%atm%altitude_levels(i+1) + dz
       !scn%atm%altitude_levels(i) = scn%atm%altitude_levels(i) * 1.045

       scn%atm%altitude_levels(i-1) = zlo * 1.0d0
       scn%atm%grav(i-1) = jpl_gravity(scn%lat, scn%atm%altitude_levels(i-1))
       scn%atm%ndry(i-1) = (scn%atm%p(i) - scn%atm%p(i-1)) / (DRY_AIR_MASS * scn%atm%grav(i-1)) &
            * (1.0d0 - 0.5d0 * (scn%atm%sh(i) + scn%atm%sh(i-1)))
       !-----------------------------------------------------

    end do

    ! Some calculations want the layer altitude, so might as well compute them
    ! here and store them along with gravity on layers.
    do i = 1, scn%num_levels - 1
       scn%atm%altitude_layers(i) = 0.5d0 * (scn%atm%altitude_levels(i) + scn%atm%altitude_levels(i+1))
       scn%atm%grav_layers(i) = jpl_gravity(scn%lat, scn%atm%altitude_layers(i))
    end do

    ! Pre-calculate the factors needed to scale up optical depth
    ! to account for pseudo-spherical geometry. So any layer-tau
    ! should be multiplied by these factors to account for the in- and
    ! outcoming beams

    if (allocated(scn%atm%ps_factors_solar)) deallocate(scn%atm%ps_factors_solar)
    if (allocated(scn%atm%ps_factors_viewing)) deallocate(scn%atm%ps_factors_viewing)
    if (allocated(scn%atm%chapman_solar)) deallocate(scn%atm%chapman_solar)
    if (allocated(scn%atm%chapman_viewing)) deallocate(scn%atm%chapman_viewing)
    if (allocated(scn%atm%layer_slant_path_viewing)) deallocate(scn%atm%layer_slant_path_viewing)
    if (allocated(scn%atm%layer_slant_path_solar)) deallocate(scn%atm%layer_slant_path_solar)

    allocate(scn%atm%layer_slant_path_solar(scn%num_levels - 1))
    allocate(scn%atm%layer_slant_path_viewing(scn%num_levels - 1))
    allocate(scn%atm%ps_factors_solar(scn%num_levels - 1))
    allocate(scn%atm%ps_factors_viewing(scn%num_levels - 1))
    allocate(scn%atm%chapman_solar(scn%num_levels - 1, scn%num_levels - 1))
    allocate(scn%atm%chapman_viewing(scn%num_levels - 1, scn%num_levels - 1))



    ! Get a better estimate for radius of earth at this latitude
    gclat = atan(tan(DEG2RAD*scn%lat)/(1+con))
    radius = EARTH_EQUATORIAL_RADIUS/sqrt(1.0d0+con*sin(gclat)**2)
    !radius = EARTH_EQUATORIAL_RADIUS

    do i = 1, scn%num_levels - 1

       ! Maybe replace this by geoid radius..
       !alt_tmp = EARTH_EQUATORIAL_RADIUS / (EARTH_EQUATORIAL_RADIUS + scn%atm%altitude_levels(i))
       alt_tmp = radius / (radius + scn%atm%altitude_levels(i+1))
       scn%atm%ps_factors_solar(i) = cos(scn%SZA * DEG2RAD) / sqrt(1.0d0 - sin(scn%SZA * DEG2RAD)**2 * alt_tmp * alt_tmp)
       scn%atm%ps_factors_viewing(i) = cos(scn%VZA * DEG2RAD) / sqrt(1.0d0 - sin(scn%VZA * DEG2RAD)**2 * alt_tmp * alt_tmp)

    end do

    ! This provides effective 1 / mu0
    call chapman_factors(scn, scn%mu0, radius, scn%atm%chapman_solar)

    ! This provides effective 1 / mu
    call chapman_factors(scn, scn%mu, radius, scn%atm%chapman_viewing)

    do i = 1, scn%num_levels - 1
       !scn%atm%layer_slant_path_viewing(i) = scn%atm%chapman_viewing(i,i)
       !scn%atm%layer_slant_path_solar(i) = scn%atm%chapman_solar(i,i)

       scn%atm%layer_slant_path_viewing(i) = 1.0d0 / (cos(scn%VZA * DEG2RAD)) * scn%atm%ps_factors_viewing(i)
       scn%atm%layer_slant_path_solar(i) = 1.0d0 / (cos(scn%SZA * DEG2RAD)) * scn%atm%ps_factors_solar(i)
    end do

  end subroutine scene_altitude


end module scene_mod
