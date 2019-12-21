!> @brief Config file valid keyword list
!> @file Keywords.f90
!> @author Peter Somkuti
!>
!! This module contains the list of accepted keywords for the config file. They
!! are hard-coded in here, so the config_check subroutine can check against this
!! list. This is a bit on the annoying side, but we want to make sure that
!! no superfluous item is present in the config file.

module keywords_mod

  ! User modules
  use control_mod, only: MAX_WINDOWS, MAX_GASES, MAX_AEROSOLS

  ! Third party modules
  use stringifor, only: string

  implicit none

  !> Max number of sections allowed in the config file
  integer, parameter :: max_sections = 50
  !> Maximum number of options per section allowed in the config file
  integer, parameter :: max_options = 99
  !> String array to hold the names of valid sections
  type(string) :: valid_sections(max_sections)
  !> String array to hold the names of valid options for each section
  type(string) :: valid_options(max_sections, max_options)

  public :: initialize_valid_sections

contains

  !> @brief Populate the "valid sections/options" that the code accepts
  !> @detail If you want to extend the code by a section or option, the
  !> new names must be put into this list below. Upon starting, the code
  !> checks all entries in the config file against this list, and exits
  !> if there is a section or option that is not known. This should prevent
  !> (some) unwanted errors in the config file.
  !> Note that the order of the keywords is NOT meaningful, and mostly a
  !> record of when a function was implemented.
  subroutine initialize_valid_sections()

    implicit none
    !> Temporary char array to hold the windows/gases with variable names
    character(len=99) :: tmp_str
    !> At which index of valid_sections do we start the 'windows'?
    integer :: window_start
    !> Window loop counter
    integer :: window_nr
    !> At which index of valid_sectionsn do we start the 'gas'es?
    integer :: gas_start
    !> Gas loop counter
    integer :: gas_nr
    !> At which index of valid_sectionsn do we start the 'aerosol's?
    integer :: aerosol_start
    !> Aerosol loop counter
    integer :: aerosol_nr
    !> General loop counter
    integer :: this_idx

    ! Everything is lowercase here! When checking the various sections
    ! of the config file, the keys are turned into lowercase, which
    ! lets the user be a bit more liberal in terms of writing the files..

    ! Related to logging, messaging, verbosity
    valid_sections(1) = "logger"
    ! Where to write the logfile to?
    valid_options(1,1) = "logfile"
    ! What logging level are we using?
    valid_options(1,2) = "loglevel"


    ! Related to the algorithm general setup
    valid_sections(2) = "algorithm"
    ! Which algorithm(s) to use?
    valid_options(2,1) = "sif_algorithm"
    valid_options(2,2) = "n_basisfunctions"
    valid_options(2,3) = "observation_mode"
    valid_options(2,4) = "step_through"

    ! Related to the instrument in question
    valid_sections(3) = "instrument"
    ! Which one?
    valid_options(3,1) = "name"

    ! Related to input files
    valid_sections(4)  = "input"
    valid_options(4,1) = "l1b_file" ! L1b file location
    valid_options(4,2) = "met_file" ! MET file location

    ! Output file options
    valid_sections(5)  = "output"
    valid_options(5,1) = "output_file"
    valid_options(5,2) = "save_radiances"
    valid_options(5,3) = "overwrite_output"
    valid_options(5,4) = "pressure_weights"
    valid_options(5,5) = "gas_averaging_kernels"

    ! Solar model type and file path
    valid_sections(6)  = "solar"
    valid_options(6,1) = "solar_type"
    valid_options(6,2) = "solar_file"

    ! We have to define our windows manually here, and give them explicit
    ! numbers!
    window_start = 6
    do window_nr=1, MAX_WINDOWS
       write(tmp_str, '(A, G0.1)') "window-", window_nr
       this_idx = window_start + window_nr

       valid_sections(this_idx)   = trim(tmp_str)
       valid_options(this_idx,1)  = "name"
       valid_options(this_idx,2)  = "wl_min"
       valid_options(this_idx,3)  = "wl_max"
       valid_options(this_idx,4)  = "basisfunctions"
       valid_options(this_idx,5)  = "gases"
       valid_options(this_idx,6)  = "statevector"
       valid_options(this_idx,7)  = "albedo_apriori"
       valid_options(this_idx,8)  = "albedo_order"
       valid_options(this_idx,9)  = "dispersion_order"
       valid_options(this_idx,10) = "dispersion_perturbation"
       valid_options(this_idx,11) = "dispersion_covariance"
       valid_options(this_idx,12) = "atmosphere"
       valid_options(this_idx,13) = "algorithms"
       valid_options(this_idx,14) = "dsigma_scale"
       valid_options(this_idx,15) = "rt_strategy"
       valid_options(this_idx,16) = "wl_spacing"
       valid_options(this_idx,17) = "band"
       valid_options(this_idx,18) = "sublayers"
       valid_options(this_idx,19) = "max_iterations"
       valid_options(this_idx,20) = "solar_dispersion"
       valid_options(this_idx,21) = "lm_gamma"
       valid_options(this_idx,22) = "allow_divergences"
       valid_options(this_idx,23) = "frame_skip"
       valid_options(this_idx,24) = "footprint_skip"
       valid_options(this_idx,25) = "inverse_method"
       valid_options(this_idx,26) = "smart_scale_first_guess_wl_in"
       valid_options(this_idx,27) = "smart_scale_first_guess_wl_out"
       valid_options(this_idx,28) = "smart_scale_first_guess_delta_tau"
       valid_options(this_idx,29) = "ils_stretch_order"
       valid_options(this_idx,30) = "ils_stretch_perturbation"
       valid_options(this_idx,31) = "ils_stretch_covariance"
       valid_options(this_idx,32) = "rt_model"
       valid_options(this_idx,33) = "aerosols"
       valid_options(this_idx,34) = "xrtm_options"
       valid_options(this_idx,35) = "xrtm_solvers"
       valid_options(this_idx,36) = "polarization"
       valid_options(this_idx,37) = "gas_prior_type"
       valid_options(this_idx,38) = "rt_streams"

    end do

    ! Section for gases
    gas_start = this_idx
    do gas_nr=1, MAX_GASES
       write(tmp_str, '(A, G0.1)') "gas-", gas_nr
       this_idx = gas_start + gas_nr

       valid_sections(this_idx)  = trim(tmp_str)
       valid_options(this_idx,1) = "name"
       valid_options(this_idx,2) = "spectroscopy_type"
       valid_options(this_idx,3) = "spectroscopy_file"
       valid_options(this_idx,4) = "hitran_index"
    end do

    ! Section for aerosols
    aerosol_start = this_idx
    do aerosol_nr=1, MAX_AEROSOLS
       write(tmp_str, '(A, G0.1)') "aerosol-", aerosol_nr
       this_idx = aerosol_start + aerosol_nr

       valid_sections(this_idx)  = trim(tmp_str)
       valid_options(this_idx,1) = "aerosol_name"
       valid_options(this_idx,2) = "aerosol_type"
       valid_options(this_idx,3) = "mom_file"
       valid_options(this_idx,4) = "mie_file"
    end do


  end subroutine initialize_valid_sections


end module
