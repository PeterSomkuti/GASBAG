
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

  end type atmosphere

  !> Holds optical gas/aerosol properties of the scene
  type optical_properties
     !> Wavelengths
     double precision, allocatable :: wl(:)
     !> Surface albedo
     double precision, allocatable :: albedo(:)

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


  !> @brief Calculate altitude levels, gravity levels and ndry
  !> @details This code is borrowed from MS3 (O'Brien, O'Dell et al.)
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
    scn%atm%altitude_levels(scn%num_active_levels) = scn%alt
    scn%atm%grav(scn%num_active_levels) = jpl_gravity(scn%lat, scn%alt)
    scn%atm%ndry(:) = 0.0d0

    ! Loop through layers, starting with the bottom-most (surface) one
    do i = scn%num_active_levels - 1, 1, -1

       if (i == scn%num_active_levels - 1) then
          ! Surface level
          p_layer = (scn%atm%p(i) + scn%atm%psurf) * 0.5d0
          dP = scn%atm%psurf - scn%atm%p(i)
          logratio = log(scn%atm%psurf / scn%atm%p(i))
       else
          ! Any other level
          p_layer = (scn%atm%p(i) + scn%atm%p(i+1)) * 0.5d0
          dP = scn%atm%p(i+1) - scn%atm%p(i)
          logratio = log(scn%atm%p(i+1) / scn%atm%p(i))
       end if

       g_layer = jpl_gravity(scn%lat, scn%atm%altitude_levels(i+1))
       SH_layer = (scn%atm%sh(i) + scn%atm%sh(i+1)) * 0.5d0
       T_layer = (scn%atm%T(i) + scn%atm%T(i+1)) * 0.5d0
       Tv = T_layer * (1.0d0 + SH_layer * (1.0d0 - EPSILON) / EPSILON)

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
    ! here and store them along with gravity on layers.
    do i = 1, scn%num_levels - 1
       scn%atm%altitude_layers(i) = 0.5d0 * (scn%atm%altitude_levels(i) + scn%atm%altitude_levels(i+1))
       scn%atm%grav_layers(i) = jpl_gravity(scn%lat, scn%atm%altitude_layers(i))
    end do

  end subroutine scene_altitude


end module scene_mod
