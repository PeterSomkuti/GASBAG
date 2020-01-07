module scene_mod


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
     !> Model atmosphere specific humidity
     double precision, allocatable :: sh(:)
     !> Model atmosphere gravity at pressure levels
     double precision, allocatable :: grav(:)
     !> Model atmosphere altitude at pressure levels
     double precision, allocatable :: altitude_levels(:)
     !> Layer altitude
     double precision, allocatable :: altitude_layers(:)
     !> Model atmosphere dry air column at layers
     double precision, allocatable :: ndry(:)
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

end module scene_mod
