!> @brief Main control_mod structure (MCS) of the program
!> @file Control.f90
!> @author Peter Somkuti
!>
!> This module is for easy access of important quantities throughout the
!> program, such as instrument name and retrieval settings, algorithm modes, ..
!> the whole shebang. A lot of other modules and subroutines access variables
!> from this module, which unfortunately makes those pieces of code very
!> dependent on the Control module - and thus less portable. On the other
!> hand, it makes some of the functions a bit cleaner, since you do not need
!> to drag every single parameter through a number of subroutines.
!> Note that this code does not do any heavy lifting, it really just passes
!> values from the FINER ini data into the various variables of the MCS structure.
!> E.g. the loading of ABSCO files etc. is all done somewhere else.

module control_mod

  ! User modules
  use file_utils_mod, only: check_config_files_exist, check_fini_error, &
       fini_extract, string_to_bool
  ! System modules
  use HDF5

  ! Third-party modules
  use stringifor
  use finer, only: file_ini
  use logger_mod, only: logger => master_logger


  implicit none

  ! These numbers sadly have to be hardcoded, as we cannot (easily) just
  ! extend the size of these arrays without a more complicated deallocation
  ! and reallocation procedure.

  !> Number of maximally allowed algorithms
  integer, parameter :: MAX_ALGORITHMS = 2
  !> Number of retrieval windows
  integer, parameter :: MAX_WINDOWS = 20
  !> Number of absorbers
  integer, parameter :: MAX_GASES = 10
  !> Number of aerosols
  integer, parameter :: MAX_AEROSOLS = 10

  type :: CS_general
     character(len=3) :: code_name = "gbg"
     !> Number of soundings to be processed
     integer :: N_soundings
     !> Number of frames and footprints
     integer :: N_fp, N_frame
     !> Number of bands
     integer :: N_bands
     !> Number of spectral points / channels per band
     integer, allocatable :: N_spec(:)
     !> Log level
     integer :: loglevel
  end type CS_general

  type :: CS_algorithm
     !> Name of the algorithm(s) used?
     type(string) :: name(MAX_ALGORITHMS)
     !> How many do we want to actually use?
     integer :: N_algorithms
     !> How mamy basis functions do we read in (and maybe use) for Guanter-type?
     integer :: N_basisfunctions
     !> Do we use the Guanter-type retrival?
     logical :: using_GK_SIF
     !> Do we use a physics-based retrieval?
     logical :: using_physical
     !> Path to the solar model file
     type(string) :: solar_file
     !> Which type of solar model?
     type(string) :: solar_type
     !> What is the observation mode? (downlooking, space-solar)
     type(string) :: observation_mode
     !> For solar observations, do we want to average all solar
     !> spectra within a L1B file (according to footprint)?
     logical :: solar_footprint_averaging
     !> Do we want to step through for debugging?
     logical :: step_through
  end type CS_algorithm

  type :: CS_window
     !> Is this CS_window structure used?
     logical :: used
     !> The name will be used in the output file
     type(string) :: name
     !> Window wavelength start (in microns)
     double precision :: wl_min
     !> Window wavelength end (in microns)
     double precision :: wl_max
     !> The high-resolution wavelength grid spacing
     double precision :: wl_spacing
     !> Which satellite instrument band are we using?
     integer :: band
     !> SV_string contains the space-separated state vector which will
     !> determine the state vector structure for the retrieval.
     type(string) :: SV_string
     !> Path to the basisfunction file
     type(string) :: basisfunction_file
     !> Order of albedo-polynomial to be retrieved (number + 1)
     integer :: albedo_order
     !> Number of dispersion coefficients to be retrieved
     integer :: dispersion_order
     !> Dispersion perturbation value for Jacobians
     double precision, allocatable :: dispersion_pert(:)
     !> Dispersion prior covariance value
     double precision, allocatable :: dispersion_cov(:)
     !> Number of ILS coefficients to be retrieved
     integer :: ils_stretch_order
     !> ILS perturbation value for Jacobians
     double precision, allocatable :: ils_stretch_pert(:)
     !> ILS prior covariance value
     double precision, allocatable :: ils_stretch_cov(:)
     !> Names of gases which are present in the window
     type(string), allocatable :: gases(:)
     !> What is the number of gases in this window?
     integer :: num_gases
     !> GAS prior string
     type(string) :: gas_prior_type_string
     !> Which type of priors are we using
     type(string), allocatable :: gas_prior_type(:)
     !> This gas_index variable holds the information about which gas-section (CS_gas)
     !> index corresponds to the gas that is stored in 'gases'
     integer, allocatable :: gas_index(:)
     !> Is this gas being retrieved?
     logical, allocatable :: gas_retrieved(:)
     !> Are we retrieving a gas profile?
     logical, allocatable :: gas_retrieve_profile(:)
     !> Are we retrieving scale factors?
     logical, allocatable :: gas_retrieve_scale(:)
     !> At which p/psurf does this partial column start?
     double precision, allocatable :: gas_retrieve_scale_start(:,:)
     !> At which p/psurf does this partial column stop?
     double precision, allocatable :: gas_retrieve_scale_stop(:,:)
     !> What is the covariance for this partial column?
     double precision, allocatable :: gas_retrieve_scale_cov(:,:)
     !> For a smart first guess of the gas scale factors, we store the
     !> wavelengths of line cores and pseudo-continuum, along with the
     !> expected delta-tau between those two.
     double precision, allocatable :: smart_scale_first_guess_wl_in(:)
     double precision, allocatable :: smart_scale_first_guess_wl_out(:)
     double precision, allocatable :: smart_scale_first_guess_delta_tau(:)
     !> Aerosols used in this window
     type(string), allocatable :: aerosol(:)
     !> Number of aerosols used in this window
     integer :: num_aerosols
     !> Aerosol distribution shape (gauss? other?) NOTE that we can only define the
     !> shape for the entire retrieval window (no mixing of different shapes)
     type(string) :: aerosol_distribution_shape
     !> This gas_index variable holds the information about which aerosol-section
     !> (CS_aerosol) index corresponds to the aerosol that is stored in 'aerosols'
     integer, allocatable :: aerosol_index(:)
     !> Do we retrieve AOD from this aerosol?
     logical, allocatable :: aerosol_retrieve_aod(:)
     !> Do we retrieve AOD in log-space?
     logical, allocatable :: aerosol_retrieve_aod_log(:)
     !> Prior AOD from SV string
     double precision, allocatable :: aerosol_prior_aod(:)
     !> Prior covariance from SV string
     double precision, allocatable :: aerosol_aod_cov(:)
     !> Do we retrieve height from this aerosol?
     logical, allocatable :: aerosol_retrieve_height(:)
     !> Do we retrieve  height in log-space?
     logical, allocatable :: aerosol_retrieve_height_log(:)
     !> Prior aerosol height from SV string
     double precision, allocatable :: aerosol_prior_height(:)
     !> Prior covariance from SV string
     double precision, allocatable :: aerosol_height_cov(:)
     !> Do we keep the scattering coefficients constant throughout the band? Speedup!
     logical :: constant_coef
     !> What is the number of sublayers to be used for gas OD calculations
     integer :: N_sublayers
     !> The dsigma_square factor to adjust convergence
     double precision :: dsigma_scale
     !> Number of iterations after which the rerieval is stopped
     integer :: max_iterations
     !> Location of the atmosphere file which must contain the gases mentioned
     !> in the 'gases' line
     type(string) :: atmosphere_file
     !> Initial value for the Levenberg-Marquart damping parameter
     double precision :: lm_gamma
     !> If this is true, then the physical retrieval will allow for divergent
     !> steps, where the LM-gamma parameter will be adjusted.
     logical :: allow_divergences
     !> Only retrieve every x'th frame
     integer :: frame_skip
     !> Only retrieve every x'th footprint
     integer :: footprint_skip
     !> Type of inverse method to use
     type(string) :: inverse_method
     !> Type of Radiative Transfer model to use
     type(string) :: RT_model
     !> What RT strategy to use
     type(string) :: RT_strategy
     !> How many quadrature streams to use
     integer, allocatable :: RT_streams(:)
     !> XRTM options to use
     type(string), allocatable :: XRTM_options(:)
     !> XRTM solvers to use
     type(string), allocatable :: XRTM_solvers(:)
     !> Do we account for polarization?
     logical :: do_polarization
  end type CS_window

  type :: CS_input
     !> Path to L1B file
     type(string) :: l1b_filename
     !> Path to MET file
     type(string) :: met_filename
     !> L1B HDF file ID
     integer(hid_t) :: l1b_file_id
     !> MET HDF file ID
     integer(hid_t) :: met_file_id
     !> Name of the instrument
     type(string) :: instrument_name
  end type CS_input

  type :: CS_output
     !> Where does the ouptut HDF file go?
     type(string) :: output_filename
     !> HDF File ID for the output file
     integer(hid_t) :: output_file_id
     !> Do we want to save radiances?
     logical :: save_radiances
     !> Do we want to override the output file?
     logical :: overwrite_output
     !> Do we want to store pressure weights?
     logical :: pressure_weights
     !> Do we want to store gas averaging kernels?
     logical :: gas_averaging_kernels
  end type CS_output

  type :: CS_gas
     !> Is this CS_gas used?
     logical :: used
     !> Name of the gas/asorber, this will be cross-referenced
     !> against the names in the CS_window%gases structure to find out which gases
     !> are to be used in the retrieval window.
     type(string) :: name
     !> Type of the cross section file (e.g. ABSCO)
     type(string) :: type
     !> Path to ABSCO file
     type(string) :: filename
     !> HITRAN index of the gas
     integer :: hitran_index

     ! At the time, we are using JPL ABSCO tables, which are 4-dimensional.
     ! Should we ever want to use different spectroscopy in the future, where the
     ! dimensionality is different, we sadly would need to either extend the array
     ! dimensions, or not use superfluous dimensions. During the calculation of
     ! optical properties, the type of spectroscopy will thus have to be referred to.
     ! The same holds true fo T,p,SH dimensions. E.g., T is 2-dimensional for ABSCO

     !> Does this spectroscopy have H2O broadening (ABSCO)?
     logical :: has_h2o
     !> Cross section data array
     double precision, allocatable :: cross_section(:,:,:,:)
     !> The wavelength dimension of the table
     double precision, allocatable :: wavelength(:)
     ! The temperature, pressure and water vapour dimensions of the cross sections
     !> The temperature dimension of the table
     double precision, allocatable :: T(:,:)
     !> The pressre dimension of the table
     double precision, allocatable :: p(:)
     !> The H2O dimension of the table
     double precision, allocatable :: H2O(:)
  end type CS_gas

  type :: CS_aerosol
     !> Is this aerosol used?
     logical :: used
     !> Name of the aerosol / identifier
     type(string) :: name
     !> Type of aerosol
     type(string) :: aer_type
     !> Name of the mom file
     type(string) :: mom_filename
     !> Name of the mie file
     type(string) :: mie_filename
     !> Default value for AOD
     double precision :: default_aod
     !> Default value for layer height
     double precision :: default_height
     !> Default value for layer width
     double precision :: default_width
     !> Coefficient array (coef, element, wavelength)
     double precision, allocatable :: coef(:,:,:)
     !> Max. no of coefs
     integer :: max_n_coef
     !> Wavelenghts
     double precision, allocatable :: wavelengths(:)
     !> Scattering efficiency (wavelength)
     double precision, allocatable :: qsca(:)
     !> Extinction efficiency (wavelength)
     double precision, allocatable :: qext(:)
     !> Single-scatter albedo (wavelength)
     double precision, allocatable :: ssa(:)
     !> Extinction cross section (wavelength)
     double precision, allocatable :: sigma_ext(:)
     !> Effective radius (wavelength)
     double precision, allocatable :: reff(:)

  end type CS_aerosol


  ! Main control_mod structure type
  type :: CS
     !> Algorithm/forwared model - related settings
     type(CS_algorithm) :: algorithm
     !> Retrieval windows
     type(CS_window) :: window(MAX_WINDOWS)
     !> Gas absorbers
     type(CS_gas) :: gas(MAX_GASES)
     !> Aerosols
     type(CS_aerosol) :: aerosol(MAX_AEROSOLS)
     !> Input files/handlers needed by the program
     type(CS_input) :: input
     !> Output settings
     type(CS_output) :: output
     !> General retrieval settings/info
     type(CS_general) :: general
  end type CS

  !> Main control structure
  type(CS), public :: MCS

  public populate_MCS

contains

  !> In here, the contents of the config file are being used to populate
  !> the main control structure of the program. It's mostly string/value
  !> parsing and making sure that the contents of the config file are
  !> in line with the expectation of the code. If something does not look
  !> right, the code will abort with error code 1, and a message stating
  !> what you did wrong.
  !
  !> @param fini FINER ini object
  subroutine populate_MCS(fini)

    implicit none

    type(file_ini), intent(in) :: fini

    ! Local variables

    ! Name of the function
    character(len=*), parameter :: fname = "populate_control_structure"
    ! Temporary character objection
    character(len=999) :: tmp_str
    character(len=999) :: win_str
    ! Number of algorithms (coming from the config)
    integer :: alg_count
    ! The string that will be turned into alg_count
    type(string), allocatable :: alg_strings(:)

    ! FINER stuff
    ! FINER error variable
    integer :: fini_error
    ! If we want to read a character array from FINER
    character(len=999) :: fini_char
    ! If we want to read a string from FINER
    type(string) :: fini_string
    ! If we want to read a double precision from FINER
    double precision :: fini_val
    ! If we want to read several double precisions from FINER
    double precision, allocatable :: fini_val_array(:)
    ! If we want to read several strings from FINER
    type(string), allocatable :: fini_string_array(:)
    ! If we want to read an integer from FINER
    integer :: fini_int

    ! Loop variables
    integer :: window_nr, gas_nr, aerosol_nr
    integer :: i
    ! Does a file exist?
    logical :: file_exists


    call logger%debug(fname, "Populating main program control structure..")

    ! ----------------------------------------------------------------------
    ! First, we set all those fields to -1 values, that are added/populated
    ! later in e.g. instrument-specific rountines.

    MCS%general%N_soundings = -1

    MCS%algorithm%using_GK_SIF = .false.
    MCS%algorithm%using_physical = .false.

    MCS%window(:)%name = ""
    MCS%window(:)%wl_min = 0.0d0
    MCS%window(:)%wl_max = 0.0d0

    call fini_extract(fini, 'logger', 'loglevel', .false., fini_int)
    if (fini_int > 0) then
       MCS%general%loglevel = fini_int
    else
       MCS%general%loglevel = 10
    end if
    
    ! ----------------------------------------------------------------------
    ! Check which algoirthms the user wants
    ! First make sure that the config file does not have more than the
    ! max. allowed number of algorithms

    alg_count = 0 ! Initialize with zero, otherwise we'll have garbage
    alg_count = fini%count_values(section_name="algorithm", &
         option_name="sif_algorithm")


    if (alg_count > MAX_ALGORITHMS) then
       write(tmp_str, '(A, I1.1, A, I3.3, A)') "We can only do ", MAX_ALGORITHMS, &
            " algorithms at most, but you want ", alg_count, '.'
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    else if (alg_count == 0) then
       call logger%warning(fname, "No SIF algorithms selected? " // &
            "Hope you know what you're doing!")
    else
       MCS%algorithm%N_algorithms = alg_count
       allocate(alg_strings(alg_count))
    end if

    ! Fortran and strings are annoying as always. First we have to have a
    ! character variable that is long enough to keep the contents of the
    ! option line, passed into it by fini%get. Then we need to cast that to
    ! a 'string' object fini_string, so we can perform the split operation
    ! where the results are going into a new string-array object.

    if (alg_count > 0) then
       call fini%get(section_name='algorithm', option_name='sif_algorithm', &
            val=fini_char, error=fini_error)
       if (fini_error /= 0) then
          call logger%fatal(fname, "Failure to get option value for " // &
               "algorithm/sif_algorithm")
          stop 1
       end if

       ! fini_string here is now hopefully space-delimited, i.e.
       ! ALG1 ALG2 ALG3
       fini_string = trim(fini_char)
       ! .. and is split and saved into alg_strings, such that
       ! alg_strings(1) = ALG1, alg_strings(2) = ALG2, etc.
       call fini_string%split(tokens=alg_strings, sep=' ', &
            max_tokens=alg_count)

       ! Stick names of algorithms into MCS
       do i=1, MCS%algorithm%N_algorithms
          MCS%algorithm%name(i) = alg_strings(i)
       end do

       ! And also check which one's we have to set the booleans correctly
       do i=1, MCS%algorithm%N_algorithms
          if (MCS%algorithm%name(i) == 'GK') then
             MCS%algorithm%using_GK_SIF = .true.
             call logger%trivia(fname, "Utilizing Guanter-type SIF retrieval!")
          else if(MCS%algorithm%name(i) == 'physical') then
             MCS%algorithm%using_physical = .true.
             call logger%trivia(fname, "Utilizing phyiscal retrieval!")
          end if
       end do
    end if

    tmp_str = "algorithm"
    if (fini%has_option(section_name=tmp_str, &
         option_name="n_basisfunctions")) then
       call fini%get(section_name='algorithm', option_name='n_basisfunctions', &
            val=fini_val, error=fini_error)

       MCS%algorithm%N_basisfunctions = int(fini_val)
    end if

    call fini_extract(fini, tmp_str, 'observation_mode', .false., fini_char)
    MCS%algorithm%observation_mode = trim(fini_char)
    ! Replace empty string (i.e. not supplied) by "downlooking" mode
    if (MCS%algorithm%observation_mode == "") then
       MCS%algorithm%observation_mode = "downlooking"
       call logger%trivia(fname, "No observation mode supplied - setting to: downlooking")
    else
       call logger%trivia(fname, "Observation mode - setting to: " &
            // MCS%algorithm%observation_mode%chars())
    end if

    call fini_extract(fini, tmp_str, 'solar_footprint_averaging', .false., fini_char)
    fini_string = fini_char
    if (fini_string == "") then
       ! If not supplied, default state is "no"
       MCS%algorithm%solar_footprint_averaging = .false.
    else
       MCS%algorithm%solar_footprint_averaging = string_to_bool(fini_string)
    end if

    call fini_extract(fini, tmp_str, 'step_through', .false., fini_char)
    fini_string = fini_char
    if (fini_string == "") then
       ! If not supplied, default state is "no"
       MCS%algorithm%step_through = .false.
    else
       MCS%algorithm%step_through = string_to_bool(fini_string)
       if (MCS%algorithm%step_through) then
          call logger%debug(fname, "Using step-through mode for debugging!")
       end if
    end if

    ! Algorithm section over------------------------------------------------

    ! Inputs section -------------------------------------------------------

    ! Check the L1b file input - this one is required
    call check_config_files_exist(fini, "input", "l1b_file", 1, file_exists)

    if(.not. file_exists) then
       call logger%fatal(fname, "L1b input check failed.")
       stop 1
    else
       ! All good? Stuff it into MCS
       call fini%get(section_name='input', option_name='l1b_file', &
            val=fini_char, error=fini_error)
       if (fini_error /= 0) then
          call logger%fatal(fname, "Error reading l1b_file string")
          stop 1
       end if
       MCS%input%l1b_filename = trim(fini_char)
    end if

    ! Do the same for the MET file
    ! If doing physical retrieval and downlooking mode, we MUST have the MET file
    if ((MCS%algorithm%using_physical .eqv. .true.) &
         .and. (MCS%algorithm%observation_mode == "downlooking")) then
       tmp_str = "input"
       if (.not. fini%has_option(section_name=tmp_str, &
            option_name="met_file")) then
          call logger%fatal(fname, "When using physical retrieval, you MUST supply a MET file.")
          stop 1
       end if

       call check_config_files_exist(fini, "input", "met_file", 1, file_exists)
       if (.not. file_exists) then
          call logger%fatal(fname, "MET file check failed.")
          stop 1
       end if

       call fini%get(section_name='input', option_name='met_file', &
            val=fini_char, error=fini_error)
       if (fini_error /= 0) then
          call logger%fatal(fname, "Error reading met_file string")
          stop 1
       end if

       MCS%input%met_filename = trim(fini_char)

    end if

    ! ----------------------------------------------------------------------

    ! Solar section ------------------------------------------------------
    ! If doing physical retrieval, we MUST have the solar section
    if (MCS%algorithm%using_physical .eqv. .true.) then
       if (.not. fini%has_section(section_name="solar")) then
          call logger%fatal(fname, "Need to have solar section when using physical retrieval.")
          stop 1
       else
          call fini%get(section_name='solar', option_name='solar_file', &
               val=fini_char, error=fini_error)
          if (fini_error /= 0) then
             call logger%fatal(fname, "Could not read solar model file name.")
             stop 1
          end if
          MCS%algorithm%solar_file = trim(fini_char)

          call fini%get(section_name='solar', option_name='solar_type', &
               val=fini_char, error=fini_error)
          if (fini_error /= 0) then
             call logger%fatal(fname, "Could not read solar model type.")
             stop 1
          end if
          MCS%algorithm%solar_type = trim(fini_char)

       end if
    end if


    ! Outputs section ------------------------------------------------------
    call fini_extract(fini, 'output', 'output_file', .true., fini_char)
    MCS%output%output_filename = trim(fini_char)

    call fini_extract(fini, 'output', 'save_radiances', .false., fini_char)
    fini_string = fini_char
    if (fini_string == "") then
       ! If not supplied, default state is "no"
       MCS%output%save_radiances = .false.
    else
       MCS%output%save_radiances = string_to_bool(fini_string)
    end if

    call fini_extract(fini, 'output', 'overwrite_output', .false., fini_char)
    fini_string = fini_char
    if (fini_string == "") then
       ! If not supplied, default state is "no"
       MCS%output%overwrite_output = .false.
    else
       MCS%output%overwrite_output = string_to_bool(fini_string)
    end if

    call fini_extract(fini, 'output', 'pressure_weights', .false., fini_char)
    fini_string = fini_char
    if (fini_string == "") then
       ! If not supplied, default state is "no"
       MCS%output%pressure_weights = .false.
    else
       MCS%output%pressure_weights = string_to_bool(fini_string)
    end if

    call fini_extract(fini, 'output', 'gas_averaging_kernels', .false., fini_char)
    fini_string = fini_char
    if (fini_string == "") then
       ! If not supplied, default state is "no"
       MCS%output%gas_averaging_kernels = .false.
    else
       MCS%output%gas_averaging_kernels = string_to_bool(fini_string)
    end if

    ! ----------------------------------------------------------------------

    ! Instrument section ---------------------------------------------------

    ! Get instrument name
    call fini_extract(fini, 'instrument', 'name', .true., fini_char)
    MCS%input%instrument_name = trim(fini_char)
    ! ----------------------------------------------------------------------


    ! Windows section ------------------------------------------------------
    ! The user might specify "window-2", and "window-5", so we need to check
    ! many possible windows here.

    do window_nr = 1, MAX_WINDOWS

       ! Is window "window_nr" in the config-file?
       write(win_str, '(A, G0.1)') "window-", window_nr
       win_str = trim(win_str)

       if (fini%has_section(section_name=win_str)) then

          ! Let's start with the required one's first!
          MCS%window(window_nr)%used = .true.

          ! Third argument in 'fini_extract' tells the function whether this
          ! is a required config setting or not. If a required setting is not found,
          ! the program will terminate with a useful error message.
          call fini_extract(fini, win_str, 'name', .true., fini_char)
          MCS%window(window_nr)%name = trim(fini_char)

          call fini_extract(fini, win_str, 'wl_min', .true., fini_val)
          MCS%window(window_nr)%wl_min = fini_val

          call fini_extract(fini, win_str, 'wl_max', .true., fini_val)
          MCS%window(window_nr)%wl_max = fini_val

          call fini_extract(fini, win_str, 'band', .true., fini_int)
          MCS%window(window_nr)%band = fini_int

          if (MCS%algorithm%using_physical) then
             ! The following are only required by the physical retrieval

             call fini_extract(fini, win_str, 'wl_spacing', .true., fini_val)
             MCS%window(window_nr)%wl_spacing = fini_val

             call fini_extract(fini, win_str, 'inverse_method', .true., fini_char)
             MCS%window(window_nr)%inverse_method = trim(fini_char)

             call fini_extract(fini, win_str, 'max_iterations', .true., fini_int)
             if (fini_int > 0) then
                MCS%window(window_nr)%max_iterations = fini_int
             else
                call logger%fatal(fname, "Max iterations has to be > 0")
                stop 1
             end if

             call fini_extract(fini, win_str, 'lm_gamma', .true., fini_val)
             if (fini_val >= 0) then
                MCS%window(window_nr)%lm_gamma = fini_val
             else
                call logger%fatal(fname, "LM-Gamma needs to be >= 0")
                stop 1
             end if

             call fini_extract(fini, win_str, 'statevector', .true., fini_char)
             MCS%window(window_nr)%SV_string = fini_char

             call fini_extract(fini, win_str, 'rt_strategy', .false., fini_char)
             MCS%window(window_nr)%RT_strategy = trim(fini_char)

             call fini_extract(fini, win_str, 'rt_model', .true., fini_char)
             MCS%window(window_nr)%RT_model = trim(fini_char)

             call fini_extract(fini, win_str, 'xrtm_options', .false., fini_string_array)
             if (allocated(fini_string_array)) then
                allocate(MCS%window(window_nr)%XRTM_options(size(fini_string_array)))
                MCS%window(window_nr)%XRTM_options(:) = fini_string_array(:)
                deallocate(fini_string_array)
             end if

             call fini_extract(fini, win_str, 'xrtm_solvers', .false., fini_string_array)
             if (allocated(fini_string_array)) then
                allocate(MCS%window(window_nr)%XRTM_solvers(size(fini_string_array)))
                MCS%window(window_nr)%XRTM_solvers(:) = fini_string_array(:)
                deallocate(fini_string_array)
             end if

             call fini_extract(fini, win_str, 'polarization', .false., fini_char)
             fini_string = fini_char
             if (fini_string == "") then
                ! If not supplied, default state is "no"
                MCS%window(window_nr)%do_polarization = .false.
             else
                MCS%window(window_nr)%do_polarization = string_to_bool(fini_string)
             end if

             call fini_extract(fini, win_str, 'rt_streams', .false., fini_string_array)
             ! If not supplied, RT_strems will not be allocated, and whether we can still
             ! carry on will be decided at a later point in the code.
             if (allocated(fini_string_array)) then

                allocate(MCS%window(window_nr)%RT_streams(size(fini_string_array)))

                call logger%info(fname, "RT stream numbers: ")
                ! Now we need to convert those strings to integers
                do i=1, size(fini_string_array)

                   fini_char = fini_string_array(i)%chars()
                   read(fini_char, '(I10)') MCS%window(window_nr)%RT_streams(i)
                   write(tmp_str, '(I2, A, I10)') i, ": ", MCS%window(window_nr)%RT_streams(i)
                   call logger%info(fname, trim(tmp_str))

                   if (mod(MCS%window(window_nr)%RT_streams(i), 2) /= 0) then
                      call logger%fatal(fname, "Sorry - Stream numbers must be > 2 and even.")
                      stop 1
                   end if

                end do

                deallocate(fini_string_array)
             end if


          end if

          ! The rest is potentially optional. Whether a certain option is
          ! required for a given retrieval setting, will be checked later
          ! on in the code, usually when it's needed the first time

          call fini_extract(fini, win_str, 'gas_prior_type', &
               .false., fini_char)
          fini_string = fini_char
          MCS%window(window_nr)%gas_prior_type_string = fini_string


          ! Arrays that are used for our super-duper smart first guess for the
          ! scalar gas retrieval
          call fini_extract(fini, win_str, 'smart_scale_first_guess_wl_in', &
               .false., fini_val_array)
          if (allocated(fini_val_array)) then
             allocate(MCS%window(window_nr)%smart_scale_first_guess_wl_in(size(fini_val_array)))
             do i=1, size(fini_val_array)
                MCS%window(window_nr)%smart_scale_first_guess_wl_in(i) = fini_val_array(i)
             end do
             deallocate(fini_val_array)
          end if

          call fini_extract(fini, win_str, 'smart_scale_first_guess_wl_out', &
               .false., fini_val_array)
          if (allocated(fini_val_array)) then
             allocate(MCS%window(window_nr)%smart_scale_first_guess_wl_out(size(fini_val_array)))
             do i=1, size(fini_val_array)
                MCS%window(window_nr)%smart_scale_first_guess_wl_out(i) = fini_val_array(i)
             end do
             deallocate(fini_val_array)
          end if

          call fini_extract(fini, win_str, 'smart_scale_first_guess_delta_tau', &
               .false., fini_val_array)
          if (allocated(fini_val_array)) then
             allocate(MCS%window(window_nr)%smart_scale_first_guess_delta_tau(size(fini_val_array)))
             do i=1, size(fini_val_array)
                MCS%window(window_nr)%smart_scale_first_guess_delta_tau(i) = fini_val_array(i)
             end do
             deallocate(fini_val_array)
          end if

          ! Now we need to check if they are the same size
          if ((size(MCS%window(window_nr)%smart_scale_first_guess_wl_in) /= &
               size(MCS%window(window_nr)%smart_scale_first_guess_wl_out)) .or. &
               (size(MCS%window(window_nr)%smart_scale_first_guess_wl_in) /= &
               size(MCS%window(window_nr)%smart_scale_first_guess_delta_tau))) then

             call logger%error(fname, "Error parsing smart gas scalar first guess!")
             call logger%error(fname, "I need to have the same number of values for each " &
                  // "of the three inputs!")
             stop 1
          end if

          call fini_extract(fini, win_str, 'allow_divergences', .false., fini_char)
          fini_string = fini_char
          if (fini_string == "") then
             ! If not supplied, default state is "no"
             MCS%window(window_nr)%allow_divergences = .false.
          else
             MCS%window(window_nr)%allow_divergences = string_to_bool(fini_string)
          end if

          call fini_extract(fini, win_str, 'frame_skip', .false., fini_int)
          if (fini_int < 1) then
             MCS%window(window_nr)%frame_skip = 1
          else
             MCS%window(window_nr)%frame_skip = fini_int
          end if

          call fini_extract(fini, win_str, 'footprint_skip', .false., fini_int)
          if (fini_int < 1) then
             MCS%window(window_nr)%footprint_skip = 1
          else
             MCS%window(window_nr)%footprint_skip = fini_int
          end if

          call fini_extract(fini, win_str, 'sublayers', .false., fini_int)
          ! We round the number of sublayers to the next odd value > 2
          if (fini_int < 2) then
             MCS%window(window_nr)%N_sublayers = 3
          else if (mod(fini_int, 2) == 0) then
             MCS%window(window_nr)%N_sublayers = fini_int + 1
          else
             MCS%window(window_nr)%N_sublayers = fini_int
          end if

          call fini_extract(fini, win_str, 'dsigma_scale', .false., fini_val)
          MCS%window(window_nr)%dsigma_scale = fini_val

          call fini_extract(fini, win_str, 'basisfunctions', .false., fini_char)
          MCS%window(window_nr)%basisfunction_file = trim(fini_char)

          call fini_extract(fini, win_str, 'albedo_order', .false., fini_int)
          MCS%window(window_nr)%albedo_order = fini_int

          call fini_extract(fini, win_str, 'ils_stretch_order', .false., fini_int)
          MCS%window(window_nr)%ils_stretch_order = fini_int

          call fini_extract(fini, win_str, 'ils_stretch_perturbation', .false., fini_val_array)
          if (allocated(fini_val_array)) then
             allocate(MCS%window(window_nr)%ils_stretch_pert(size(fini_val_array)))
             do i=1, size(fini_val_array)
                MCS%window(window_nr)%ils_stretch_pert(i) = fini_val_array(i)
             end do
             deallocate(fini_val_array)
          end if

          call fini_extract(fini, win_str, 'ils_stretch_covariance', .false., fini_val_array)
          if (allocated(fini_val_array)) then
             allocate(MCS%window(window_nr)%ils_stretch_cov(size(fini_val_array)))
             do i=1, size(fini_val_array)
                MCS%window(window_nr)%ils_stretch_cov(i) = fini_val_array(i)
             end do
             deallocate(fini_val_array)
          end if

          call fini_extract(fini, win_str, 'dispersion_order', .false., fini_int)
          MCS%window(window_nr)%dispersion_order = fini_int

          call fini_extract(fini, win_str, 'dispersion_perturbation', .false., fini_val_array)
          if (allocated(fini_val_array)) then
             allocate(MCS%window(window_nr)%dispersion_pert(size(fini_val_array)))
             do i=1, size(fini_val_array)
                MCS%window(window_nr)%dispersion_pert(i) = fini_val_array(i)
             end do
             deallocate(fini_val_array)
          end if

          call fini_extract(fini, win_str, 'dispersion_covariance', .false., fini_val_array)
          if (allocated(fini_val_array)) then
             allocate(MCS%window(window_nr)%dispersion_cov(size(fini_val_array)))
             do i=1, size(fini_val_array)
                MCS%window(window_nr)%dispersion_cov(i) = fini_val_array(i)
             end do
             deallocate(fini_val_array)
          end if

          MCS%window(window_nr)%num_gases = 0
          call fini_extract(fini, win_str, 'gases', .false., fini_string_array)

          if (allocated(fini_string_array)) then

             ! Depending on the number of gases supplied via the "gases" option, we
             ! allocate that number of fields for the various gas-related variables.
             allocate(MCS%window(window_nr)%gases(size(fini_string_array)))
             allocate(MCS%window(window_nr)%gas_prior_type(size(fini_string_array)))
             allocate(MCS%window(window_nr)%gas_retrieved(size(fini_string_array)))
             allocate(MCS%window(window_nr)%gas_retrieve_profile(size(fini_string_array)))
             allocate(MCS%window(window_nr)%gas_retrieve_scale(size(fini_string_array)))
             allocate(MCS%window(window_nr)%gas_retrieve_scale_start(size(fini_string_array), 99))
             allocate(MCS%window(window_nr)%gas_retrieve_scale_stop(size(fini_string_array), 99))
             allocate(MCS%window(window_nr)%gas_retrieve_scale_cov(size(fini_string_array), 99))
             allocate(MCS%window(window_nr)%gas_index(size(fini_string_array)))

             do i=1, size(fini_string_array)
                MCS%window(window_nr)%gas_prior_type(i) = ""
                MCS%window(window_nr)%gases(i) = fini_string_array(i)
                MCS%window(window_nr)%num_gases = MCS%window(window_nr)%num_gases + 1
             end do
             deallocate(fini_string_array)
          end if

          MCS%window(window_nr)%num_aerosols = 0
          call fini_extract(fini, win_str, 'aerosols', .false., fini_string_array)

          if (allocated(fini_string_array)) then

             ! If aerosols are considered in this window, one MUST specify an aerosol
             ! type.

             call fini_extract(fini, win_str, 'aerosol_distribution_shape', .true., fini_char)
             MCS%window(window_nr)%aerosol_distribution_shape = trim(fini_char)

             allocate(MCS%window(window_nr)%aerosol(size(fini_string_array)))
             allocate(MCS%window(window_nr)%aerosol_index(size(fini_string_array)))

             allocate(MCS%window(window_nr)%aerosol_retrieve_aod(size(fini_string_array)))
             allocate(MCS%window(window_nr)%aerosol_retrieve_aod_log(size(fini_string_array)))
             allocate(MCS%window(window_nr)%aerosol_prior_aod(size(fini_string_array)))
             allocate(MCS%window(window_nr)%aerosol_aod_cov(size(fini_string_array)))

             allocate(MCS%window(window_nr)%aerosol_retrieve_height(size(fini_string_array)))
             allocate(MCS%window(window_nr)%aerosol_retrieve_height_log(size(fini_string_array)))
             allocate(MCS%window(window_nr)%aerosol_prior_height(size(fini_string_array)))
             allocate(MCS%window(window_nr)%aerosol_height_cov(size(fini_string_array)))

             MCS%window(window_nr)%aerosol_retrieve_aod(:) = .false.
             MCS%window(window_nr)%aerosol_retrieve_aod_log(:) = .false.
             MCS%window(window_nr)%aerosol_retrieve_height(:) = .false.
             MCS%window(window_nr)%aerosol_retrieve_height_log(:) = .false.

             do i=1, size(fini_string_array)
                MCS%window(window_nr)%aerosol(i) = fini_string_array(i)
                MCS%window(window_nr)%num_aerosols = MCS%window(window_nr)%num_aerosols + 1
             end do

             deallocate(fini_string_array)
          end if
          

          call fini_extract(fini, win_str, 'atmosphere', .false., fini_char)
          MCS%window(window_nr)%atmosphere_file = trim(fini_char)

          call fini_extract(fini, win_str, 'keep_scattering_constant', .false., fini_char)
          fini_string = fini_char
          if (fini_string == "") then
             ! If not supplied, default state is "no"
             MCS%window(window_nr)%constant_coef = .false.
          else
             MCS%window(window_nr)%constant_coef = string_to_bool(fini_string)
          end if


       else
          MCS%window(window_nr)%used = .false.
       end if
    end do

    ! ----------------------------------------------------------------------

    ! Gases section --------------------------------------------------------
    ! This is done exactly the same as the windows section above.

    do gas_nr=1, MAX_GASES

       ! Is window "window_nr" in the config-file?
       write(tmp_str, '(A, G0.1)') "gas-", gas_nr
       tmp_str = trim(tmp_str)

       if (fini%has_section(section_name=tmp_str)) then

          MCS%gas(gas_nr)%used = .true.

          call fini_extract(fini, tmp_str, 'name', .true., fini_char)
          MCS%gas(gas_nr)%name = trim(fini_char)

          call fini_extract(fini, tmp_str, 'spectroscopy_type', .true., fini_char)
          MCS%gas(gas_nr)%type = trim(fini_char)

          call fini_extract(fini, tmp_str, 'spectroscopy_file', .true., fini_char)
          MCS%gas(gas_nr)%filename = trim(fini_char)

          call fini_extract(fini, tmp_str, 'hitran_index', .false., fini_int)
          if (fini_int == -9999) then
             call logger%debug(fname, "No HITRAN index supplied.")
             fini_int = -1
          end if
          MCS%gas(gas_nr)%hitran_index = fini_int

       else
          MCS%gas(gas_nr)%used = .false.
       end if
    end do

    ! Aerosols section ---------------------------------------------------
    ! This is done exactly the same as the windows section above.
    do aerosol_nr=1, MAX_AEROSOLS

       ! Is window "window_nr" in the config-file?
       write(tmp_str, '(A, G0.1)') "aerosol-", aerosol_nr
       tmp_str = trim(tmp_str)

       if (fini%has_section(section_name=tmp_str)) then

          MCS%aerosol(aerosol_nr)%used = .true.

          call fini_extract(fini, tmp_str, 'aerosol_name', .true., fini_char)
          MCS%aerosol(aerosol_nr)%name = trim(fini_char)

          call fini_extract(fini, tmp_str, 'aerosol_type', .true., fini_char)
          MCS%aerosol(aerosol_nr)%aer_type = trim(fini_char)

          call fini_extract(fini, tmp_str, 'mom_file', .true., fini_char)
          MCS%aerosol(aerosol_nr)%mom_filename = trim(fini_char)

          call fini_extract(fini, tmp_str, 'mie_file', .true., fini_char)
          MCS%aerosol(aerosol_nr)%mie_filename = trim(fini_char)

          call fini_extract(fini, tmp_str, 'default_aod', .true., fini_val)
          MCS%aerosol(aerosol_nr)%default_aod = fini_val

          call fini_extract(fini, tmp_str, 'default_width', .true., fini_val)
          MCS%aerosol(aerosol_nr)%default_width = fini_val

          call fini_extract(fini, tmp_str, 'default_height', .true., fini_val)
          MCS%aerosol(aerosol_nr)%default_height = fini_val

       else
          MCS%aerosol(aerosol_nr)%used = .false.
       end if

    end do


  end subroutine populate_MCS

  !> Subroutine to match gases from the window object to the
  !> gases / spectroscopy sections.
  !>
  !> @param window CS_window object array
  !> @param gas CS_gas object array
  !> @param i_win CS_window index
  subroutine MCS_find_gases(window, gas, i_win)
    implicit none

    type(CS_window), intent(inout) :: window(:)
    type(CS_gas), intent(inout) :: gas(:)
    integer, intent(in) :: i_win

    ! Local variables

    ! Loop variables
    integer :: i, j
    ! Set if gas is found
    logical :: gas_found
    ! Function name
    character(len=*), parameter :: fname = "MCS_find_gases"
    ! Temp character array for conversions
    character(len=999) :: tmp_str

    ! If we have no gases in this window, there's nothing to do
    ! so might as well just return.
    if (window(i_win)%num_gases == 0) return

    do i=1, size(window(i_win)%gases)
       ! Loop over all gases specified in the retrieval window

       ! Skip unused retrieval windows
       if (.not. window(i_win)%used) cycle

       gas_found = .false.
       do j=1, MAX_GASES
          if (window(i_win)%gases(i) == gas(j)%name) then
             gas_found = .true.
             write(tmp_str, '(A, A, A, G0.1, A)')  "Gas found: ",  &
                  window(i_win)%gases(i)%chars(), " (gas-", j, ")"
             call logger%info(fname, trim(tmp_str))
             ! And also store which gas section corresponds to this particular gas
             ! in the window gas definition.
             window(i_win)%gas_index(i) = j
             ! Gas was found, step out of loop
             exit
          end if
       end do

       ! If this specific gas was not found, kill the program immediately. There's no use-case
       ! for a gas being specified in a retrieval window, and that gas then not being defined
       ! in a 'gas'-section.
       if (.not. gas_found) then
          write(tmp_str, '(A, A, A, G0.1)') "Sorry - gas '", window(i_win)%gases(i)%chars(), &
               "' was not found in window-", dble(i_win)
          call logger%fatal(fname, trim(tmp_str))
          stop 1
       end if

    end do ! Finish first loop to find/match gases with window gases

  end subroutine MCS_find_gases

  !> Subroutine to match aerosols from the window object to the
  !> aerosol sections.
  !>
  !> @param window CS_window object array
  !> @param aerosol CS_aerosol object array
  !> @param i_win CS_window index
  subroutine MCS_find_aerosols(window, aerosol, i_win)
    implicit none

    type(CS_window), intent(inout) :: window(:)
    type(CS_aerosol), intent(inout) :: aerosol(:)
    integer, intent(in) :: i_win

    ! Local variables

    ! Loop variables
    integer :: i, j
    ! Set if gas is found
    logical :: aerosol_found
    ! Function name
    character(len=*), parameter :: fname = "MCS_find_aerosols"
    ! Temp character array for conversions
    character(len=999) :: tmp_str

    ! If we have no aerosols in this window, there's nothing to do
    ! so might as well just return.
    if (window(i_win)%num_aerosols == 0) return

    do i=1, size(window(i_win)%aerosol)
       ! Loop over all gases specified in the retrieval window

       ! Skip unused retrieval windows
       if (.not. window(i_win)%used) cycle

       aerosol_found = .false.
       do j=1, MAX_AEROSOLS
          if (window(i_win)%aerosol(i) == aerosol(j)%name) then
             aerosol_found = .true.
             write(tmp_str, '(A, A, A, G0.1, A)')  "Aerosol found: ",  &
                  window(i_win)%aerosol(i)%chars(), " (aerosol-", j, ")"
             call logger%info(fname, trim(tmp_str))
             ! And also store which gas section corresponds to this particular
             ! aerosol in the window aerosol definition.
             window(i_win)%aerosol_index(i) = j
             ! Gas was found, step out of loop
             exit
          end if
       end do

       ! If this specific gas was not found, kill the program immediately. There's no use-case
       ! for a gas being specified in a retrieval window, and that gas then not being defined
       ! in a 'gas'-section.
       if (.not. aerosol_found) then
          write(tmp_str, '(A, A, A, G0.1)') "Sorry - aerosol '", window(i_win)%aerosol(i)%chars(), &
               "' was not found in window-", dble(i_win)
          call logger%fatal(fname, trim(tmp_str))
          stop 1
       end if

    end do ! Finish first loop to find/match gases with window gases

  end subroutine MCS_find_aerosols


  subroutine MCS_find_gas_priors(window, gas, i_win)

    implicit none

    type(CS_window), intent(inout) :: window(:)
    type(CS_gas), intent(in) :: gas(:)
    integer, intent(in) :: i_win

    ! Function name
    character(len=*), parameter :: fname = "MCS_find_gas_priors"
    character(len=999) :: tmp_str
    type(string), allocatable :: split_string(:), prior_strings(:)
    integer :: i, j
    integer :: gas_idx, MCS_gas_pos
    logical :: found_gas

    call window(i_win)%gas_prior_type_string%split(&
         tokens=split_string, sep=' ', &
         max_tokens=MAX_GASES)

    ! Split the gas_prior_type string into substrings, separated by spaces
    do i=1, size(split_string)

       ! And now separate them at the colon ':'
       call split_string(i)%split(tokens=prior_strings, sep=':', &
            max_tokens=2)

       write(tmp_str, '(A,A)') "Matching up gas prior for requested gas ", prior_strings(1)%chars()
       call logger%debug(fname, trim(tmp_str))

       ! Going through the gases, we have to find which one matches the gas
       ! specified in this prior type
       found_gas = .false.
       do j=1, window(i_win)%num_gases
          ! Skip unused
          if (.not. gas(j)%used) cycle

          ! Found the gas!
          if (gas(j)%name == prior_strings(1)) then
             write(tmp_str, '(A, G0.1, A, A)') "Gas matching gas prior type found at index ", &
                  j, " with name ", gas(j)%name%chars()
             call logger%debug(fname, trim(tmp_str))
             found_gas = .true.
             gas_idx = j
             exit
          end if
       end do

       ! If gas wasn't found - end it here. We don't want the default behavior to
       ! be a fallback solution. The user did something wrong, so flag it up here!!
       if (.not. found_gas) then
          write(tmp_str, '(A,A)') "No gas matching gas prior type found for requested gas: ", &
               prior_strings(1)%chars()
          call logger%fatal(fname, trim(tmp_str))
          stop 1
       end if

       ! Write the gas prior type string into the corresponding position in
       ! MCS%window
       MCS_gas_pos = MCS%window(i_win)%gas_index(gas_idx)
       MCS%window(i_win)%gas_prior_type(MCS_gas_pos) = prior_strings(2)

    end do

  end subroutine MCS_find_gas_priors


end module control_mod
