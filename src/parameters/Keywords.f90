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
  !> @details If you want to extend the code by a section or option, the
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
    ! Number of waveforms to fit
    valid_options(2,2) = "n_basisfunctions"
    ! What observation mode?
    valid_options(2,3) = "observation_mode"
    ! Step-through mode for debugging and exploring retrieval set-up
    valid_options(2,4) = "step_through"

    ! Related to the instrument in question
    valid_sections(3) = "instrument"
    ! Which one?
    valid_options(3,1) = "name"

    ! Related to input files
    valid_sections(4)  = "input"
    ! L1b file location
    valid_options(4,1) = "l1b_file"
    ! MET file location
    valid_options(4,2) = "met_file"
    ! Preload L1B spectra?
    valid_options(4,3) = "preload_spectra"

    ! Output file options
    valid_sections(5)  = "output"
    ! Where to save the output file
    valid_options(5,1) = "output_file"
    ! Do we want to save radiances, and noise (makes output file BIG)
    valid_options(5,2) = "save_radiances"
    ! Do we want to overwrite existing output files?
    valid_options(5,3) = "overwrite_output"
    ! Do we want to store pressure weights?
    valid_options(5,4) = "pressure_weights"
    ! Do we want to store gas averaging kernels? (Make the output a bit larger..)
    valid_options(5,5) = "gas_averaging_kernels"

    ! Solar model type and file path
    valid_sections(6)  = "solar"
    ! What solar model type are we using?
    valid_options(6,1) = "solar_type"
    ! Where is the solar model file?
    valid_options(6,2) = "solar_file"

    ! We have to define our windows manually here,
    ! and give them explicit numbers!
    window_start = 6
    do window_nr=1, MAX_WINDOWS

       write(tmp_str, '(A, G0.1)') "window-", window_nr
       this_idx = window_start + window_nr

       valid_sections(this_idx)   = trim(tmp_str)
       ! Name of the retrieval window/setup
       valid_options(this_idx,1)  = "name"
       ! Start of the window in terms of wavelength
       valid_options(this_idx,2)  = "wl_min"
       ! End of the window in terms of wavelength
       valid_options(this_idx,3)  = "wl_max"
       ! Where is the file that contains the waveforms for PCA fitting
       valid_options(this_idx,4)  = "basisfunctions"
       ! Which gases are in this atmosphere?
       valid_options(this_idx,5)  = "gases"
       ! What are the statevectors that we want to retrieve?
       valid_options(this_idx,6)  = "statevector"
       ! ??? UNUSED ???
       valid_options(this_idx,7)  = "albedo_apriori"
       ! What order polynomial to use for surface albedo (0 = constant, 1 = linear, etc.)
       valid_options(this_idx,8)  = "albedo_order"
       ! What order polynomial to use for dispersion retrieval
       valid_options(this_idx,9)  = "dispersion_order"
       ! What values to use for derivative calculation (finite differencing)
       valid_options(this_idx,10) = "dispersion_perturbation"
       ! What are the values for the dispersion covariance?
       valid_options(this_idx,11) = "dispersion_covariance"
       ! Where is the model atmosphere file?
       valid_options(this_idx,12) = "atmosphere"
       ! ??? UNUSED ???
       valid_options(this_idx,13) = "algorithms"
       ! Value to scale the "dsigma factor" for convergence determination
       valid_options(this_idx,14) = "dsigma_scale"
       ! Which RT strategy are we using, in combination with XRTM? (monochromatic, ..)
       valid_options(this_idx,15) = "rt_strategy"
       ! For the high-resolution spectrum, which interval are we using?
       valid_options(this_idx,16) = "wl_spacing"
       ! Which instrument band does this retrieval window lie in?
       valid_options(this_idx,17) = "band"
       ! How many sublayers are to be used for the gas property calculations?
       valid_options(this_idx,18) = "sublayers"
       ! How many iterations are allowed until the retrieval is considered non-converged
       valid_options(this_idx,19) = "max_iterations"
       ! ??? UNUSED ???
       valid_options(this_idx,20) = "solar_dispersion"
       ! What is the initial value of the Levenberg-Marquardt value?
       valid_options(this_idx,21) = "lm_gamma"
       ! Do we allow divergences in this retrieval?
       valid_options(this_idx,22) = "allow_divergences"
       ! GET RID OF THESE
       valid_options(this_idx,23) = "frame_skip"
       valid_options(this_idx,24) = "footprint_skip"
       ! Which inversion method to use?
       valid_options(this_idx,25) = "inverse_method"
       ! These three inform the retrieval about a smarter first guess for the
       ! scale factors. See documentation or ask Peter.
       valid_options(this_idx,26) = "smart_scale_first_guess_wl_in"
       valid_options(this_idx,27) = "smart_scale_first_guess_wl_out"
       valid_options(this_idx,28) = "smart_scale_first_guess_delta_tau"
       ! When retrieving ILS stretch, this is the polynomial order
       valid_options(this_idx,29) = "ils_stretch_order"
       ! Finite-difference jacobians for ILS stretch
       valid_options(this_idx,30) = "ils_stretch_perturbation"
       ! ILS stretch covariance values
       valid_options(this_idx,31) = "ils_stretch_covariance"
       ! Which RT Model to use (XRTM, Beer-Lambert)
       valid_options(this_idx,32) = "rt_model"
       ! What aerosols do we have in the atmosphere?
       valid_options(this_idx,33) = "aerosols"
       ! What options to supply to XRTM? (UNUSED)
       valid_options(this_idx,34) = "xrtm_options"
       ! Which XRTM solvers are to be used to calculate radiances (might be overruled in code)
       valid_options(this_idx,35) = "xrtm_solvers"
       ! Do we calculate more than just the first Stokes coefficient in the radiances?
       valid_options(this_idx,36) = "polarization"
       ! Which gas prior type do we use? (Maybe unused right now)
       valid_options(this_idx,37) = "gas_prior_type"
       ! For quadrature-type RT models, how many streams (both hemispheres)?
       valid_options(this_idx,38) = "rt_streams"
       ! When calculating phase function coefficients - do we keep them constant for the window?
       ! Or do we allow for wavelength dependence
       valid_options(this_idx,40) = "keep_scattering_constant"
       ! How do we distribute aerosols in the scene? Gaussian?
       valid_options(this_idx,41) = "aerosol_distribution_shape"
       ! Do we use another GASBAG result file to slip in as priors?
       valid_options(this_idx,42) = "gasbag_result_file_for_prior"
       ! Which data do we replace as priors?
       valid_options(this_idx,43) = "gasbag_priors"
       ! Solar irradiance scaling
       valid_options(this_idx,44) = "solar_irrad_scale_order"
       ! ILS stretch parameter priors
       valid_options(this_idx,45) = "ils_stretch_prior"

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
       valid_options(this_idx,5) = "default_height"
       valid_options(this_idx,6) = "default_aod"
       valid_options(this_idx,7) = "default_width"
    end do


  end subroutine initialize_valid_sections


end module
