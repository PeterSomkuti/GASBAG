!! Main control_mod structure (MCS) of the program, for easy access of important
!! quantities throughout the program, such as instrument name and retrieval
!! settings, algorithm modes, .. the whole shebang.
!! This is designed to be instrument-independent, so apart from the name of the
!! instrument, no instrument-specifics should be stored here.

module control_mod

    use stringifor
    use finer, only: file_ini
    use logger_mod, only: logger => master_logger
    use file_utils_mod, only: check_config_files_exist, check_fini_error, fini_extract
    use HDF5

    implicit none

    ! These numbers sadly have to be hardcoded, as we cannot (easily) just
    ! extend the size of these arrays without a more complicated deallocation
    ! and reallocation procedure.

    ! At the moment, we only plan the Guanter and Frankenberg methods
    integer, parameter :: MAX_ALGORITHMS = 2
    ! Number of retrieval windows
    integer, parameter :: MAX_WINDOWS = 10
    ! Number of absorbers
    integer, parameter :: MAX_GASES = 10

    type, private :: CS_general
        integer(8) :: N_soundings ! Number of soundings to be processed
    end type

    type, private :: CS_algorithm
       type(string) :: name(MAX_ALGORITHMS) ! Name of the algorithm(s) used?
       integer :: N_algorithms ! How many do we want to actually use?
       integer :: N_basisfunctions ! How mamy basis functions do we read in (and maybe use)?
       logical :: using_GK_SIF ! Do we use the Guanter-type retrival?
       logical :: using_physical ! Do we use a physics-based retrieval?
       type(string) :: solar_file ! Path to the solar model file_exists
       type(string) :: solar_type ! Which type of solar model?
    end type CS_algorithm

    type, private :: CS_window
       logical :: used ! Is this CS_window structure used?
       type(string) :: name ! The name will be used in the output file
       double precision :: wl_min ! Window wavelength start and end
       double precision :: wl_max
       double precision :: wl_spacing ! The high-resolution wavelength grid spacing
       integer :: band ! Which satellite instrument band are we using? 
       ! SV_string contains the space-separated state vector which will
       ! determine the state vector structure for the retrieval.
       type(string) :: SV_string
       type(string) :: basisfunction_file
       ! How many parameters of various state vector elements do we want to
       ! retrieve (physical retrieval only)
       integer :: albedo_order, dispersion_order
       double precision, allocatable :: dispersion_pert(:), dispersion_cov(:)
       type(string), allocatable :: gases(:)
       ! This gas_index variable holds the information about which gas-section (CS_gas)
       ! index corresponds to the gas that is stored in 'gases'
       integer, allocatable :: gas_index(:)
       ! Is this gas being retrieved?
       logical, allocatable :: gas_retrieved(:)
       ! What type of gas retrieval is used here?
       logical, allocatable :: gas_retrieve_profile(:), gas_retrieve_scale(:)
       double precision, allocatable :: gas_retrieve_scale_start(:,:), &
            gas_retrieve_scale_stop(:,:), gas_retrieve_scale_cov(:,:)
       integer :: num_gases
       integer :: N_sublayers
       double precision :: dsigma_scale ! Convergence scaling factor
       integer :: max_iterations
       ! Do we use the less-accurate, but faster FFT convolution with an
       ! averaged ILS kernel?
       logical :: fft_convolution 
       ! Location of the atmosphere file which must contain the gases mentioned
       ! in the 'gases' line
       type(string) :: atmosphere_file
       ! Initial value for the Levenberg-Marquart damping parameter
       double precision :: lm_gamma
    end type CS_window

    type, private :: CS_input
       type(string) :: l1b_filename , met_filename! Path to the L1B/MET file
       integer(hid_t) :: l1b_file_id, met_file_id ! HDF file handler of L1B/MET file
       type(string) :: instrument_name ! Name of the instrument
    end type CS_input

    type, private :: CS_output
       type(string) :: output_filename ! Where does the ouptut HDF file go?
       integer(hid_t) :: output_file_id
    end type CS_output

    type :: CS_gas
       logical :: used ! Is this CS_gas used?
       type(string) :: name ! Name of the gas/absorber, this will be cross-referenced
       ! against the names in the CS_window%gases structure to find out which gases
       ! are to be used in the retrieval window.
       type(string) :: type
       type(string) :: filename

       ! At the time, we are using JPL ABSCO tables, which are 4-dimensional.
       ! Should we ever want to use different spectroscopy in the future, where the
       ! dimensionality is different, we sadly would need to either extend the array
       ! dimensions, or not use superfluous dimensions. During the calculation of
       ! optical properties, the type of spectroscopy will thus have to be referred to.
       ! The same holds true fo T,p,SH dimensions. E.g., T is 2-dimensional for ABSCO
       logical :: has_h2o ! Does this spectroscopy have H2O broadening?
       double precision, allocatable :: cross_section(:,:,:,:)
       ! The wavelength dimension
       double precision, allocatable :: wavelength(:)
       ! The temperature, pressure and water vapour dimensions of the cross sections
       double precision, allocatable :: T(:,:), p(:), H2O(:)
    end type CS_gas


    ! Main control_mod structure type
    type, private :: CS
        type(CS_algorithm) :: algorithm ! Algorithm/forwared model - related settings
        type(CS_window) :: window(MAX_WINDOWS) ! Retrieval windows
        type(CS_gas) :: gas(MAX_GASES) ! Absorbers
        type(CS_input) :: input ! Input files needed by the program
        type(CS_output) :: output ! Output file path(s)
        type(CS_general) :: general
    end type

    ! Define it here. Rest of the code should be allowed to change data?
    type(CS), public :: MCS

    public populate_MCS

contains

    subroutine populate_MCS(fini)
        !! In here, the contents of the config file are being used to populate
        !! the main control structure of the program. It's mostly string/value
        !! parsing and making sure that the contents of the config file are
        !! in line with the expectation of the code. If something does not look
        !! right, the code will abort with error code 1, and a message stating
        !! what you did wrong.

        implicit none

        type(file_ini), intent(in) :: fini
        character(len=*), parameter :: fname = "populate_control_structure"
        character(len=999) :: tmp_str
        integer :: alg_count
        type(string), allocatable :: alg_strings(:)

        ! FINER stuff
        integer :: fini_error
        character(len=999) :: fini_char
        type(string) :: fini_string
        double precision :: fini_val
        double precision, allocatable :: fini_val_array(:)
        type(string), allocatable :: fini_string_array(:)
        integer :: fini_int

        integer :: window_nr, gas_nr
        integer :: i
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
        ! If doing physical retrieval, we MUST have the MET file
        if (MCS%algorithm%using_physical .eqv. .true.) then
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
            write(tmp_str, '(A, G0.1)') "window-", window_nr
            tmp_str = trim(tmp_str)

            if (fini%has_section(section_name=tmp_str)) then

                ! Let's start with the required one's first!
                MCS%window(window_nr)%used = .true.

                ! Third argument in 'fini_extract' tells the function whether this
                ! is a required config setting or not. If a required setting is not found,
                ! the program will terminate with a useful error message.
                call fini_extract(fini, tmp_str, 'name', .true., fini_char)
                MCS%window(window_nr)%name = trim(fini_char)

                call fini_extract(fini, tmp_str, 'wl_min', .true., fini_val)
                MCS%window(window_nr)%wl_min = fini_val

                call fini_extract(fini, tmp_str, 'wl_max', .true., fini_val)
                MCS%window(window_nr)%wl_max = fini_val

                call fini_extract(fini, tmp_str, 'wl_spacing', .true., fini_val)
                MCS%window(window_nr)%wl_spacing = fini_val

                call fini_extract(fini, tmp_str, 'band', .true., fini_int)
                MCS%window(window_nr)%band = fini_int

                call fini_extract(fini, tmp_str, 'max_iterations', .true., fini_int)
                if (fini_int > 0) then
                   MCS%window(window_nr)%max_iterations = fini_int
                else
                   call logger%fatal(fname, "Max iterations has to be > 0")
                   stop 1
                end if

                call fini_extract(fini, tmp_str, 'lm_gamma', .true., fini_val)
                if (fini_val >= 0) then
                   MCS%window(window_nr)%lm_gamma = fini_val
                else
                   call logger%fatal(fname, "LM-Gamma needs to be >= 0")
                   stop 1
                end if

                call fini_extract(fini, tmp_str, 'statevector', .true., fini_char)
                MCS%window(window_nr)%SV_string = fini_char

                ! The rest is potentially optional. Whether a certain option is
                ! required for a given retrieval setting, will be checked later
                ! on in the code, usually when it's needed the first time

                call fini_extract(fini, tmp_str, 'fft_convolution', .false., fini_char)
                fini_string = fini_char
                if (fini_string == "") then
                   ! If not supplied, default state is "no"
                   MCS%window(window_nr)%fft_convolution = .false.
                else
                   if (fini_string%lower() == "true") then
                      MCS%window(window_nr)%fft_convolution = .true.
                   else if (fini_string%lower() == "false") then
                      MCS%window(window_nr)%fft_convolution = .false.
                   else
                      call logger%fatal(fname, "Sorry, -fft_convolution- option accepts " &
                           // "only T/true or F/false.")
                      stop 1
                   end if
                end if


                call fini_extract(fini, tmp_str, 'sublayers', .false., fini_int)
                ! We round the number of sublayers to the next odd value > 2
                if (fini_int < 2) then
                   MCS%window(window_nr)%N_sublayers = 3
                else if (mod(fini_int, 2) == 0) then
                   MCS%window(window_nr)%N_sublayers = fini_int + 1
                else
                   MCS%window(window_nr)%N_sublayers = fini_int
                end if

                call fini_extract(fini, tmp_str, 'dsigma_scale', .false., fini_val)
                MCS%window(window_nr)%dsigma_scale = fini_val

                call fini_extract(fini, tmp_str, 'basisfunctions', .false., fini_char)
                MCS%window(window_nr)%basisfunction_file = trim(fini_char)

                call fini_extract(fini, tmp_str, 'albedo_order', .false., fini_int)
                MCS%window(window_nr)%albedo_order = fini_int

                call fini_extract(fini, tmp_str, 'dispersion_order', .false., fini_int)
                MCS%window(window_nr)%dispersion_order = fini_int

                call fini_extract(fini, tmp_str, 'dispersion_perturbation', .false., fini_val_array)
                if (allocated(fini_val_array)) then
                    allocate(MCS%window(window_nr)%dispersion_pert(size(fini_val_array)))
                    do i=1, size(fini_val_array)
                        MCS%window(window_nr)%dispersion_pert(i) = fini_val_array(i)
                    end do
                    deallocate(fini_val_array)
                end if

                call fini_extract(fini, tmp_str, 'dispersion_covariance', .false., fini_val_array)
                if (allocated(fini_val_array)) then
                    allocate(MCS%window(window_nr)%dispersion_cov(size(fini_val_array)))
                    do i=1, size(fini_val_array)
                       MCS%window(window_nr)%dispersion_cov(i) = fini_val_array(i)
                     end do
                    deallocate(fini_val_array)
                 end if

                 MCS%window(window_nr)%num_gases = 0
                 call fini_extract(fini, tmp_str, 'gases', .false., fini_string_array)
                 if (allocated(fini_string_array)) then

                    allocate(MCS%window(window_nr)%gases(size(fini_string_array)))
                    allocate(MCS%window(window_nr)%gas_retrieved(size(fini_string_array)))
                    allocate(MCS%window(window_nr)%gas_retrieve_profile(size(fini_string_array)))
                    allocate(MCS%window(window_nr)%gas_retrieve_scale(size(fini_string_array)))
                    allocate(MCS%window(window_nr)%gas_retrieve_scale_start(size(fini_string_array), 99))
                    allocate(MCS%window(window_nr)%gas_retrieve_scale_stop(size(fini_string_array), 99))
                    allocate(MCS%window(window_nr)%gas_retrieve_scale_cov(size(fini_string_array), 99))
                    allocate(MCS%window(window_nr)%gas_index(size(fini_string_array)))


                    do i=1, size(fini_string_array)
                       MCS%window(window_nr)%gases(i) = fini_string_array(i)
                       MCS%window(window_nr)%num_gases = MCS%window(window_nr)%num_gases + 1
                    end do
                    deallocate(fini_string_array)
                 end if

                 call fini_extract(fini, tmp_str, 'atmosphere', .false., fini_char)
                 MCS%window(window_nr)%atmosphere_file = trim(fini_char)

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

           else
              MCS%gas(gas_nr)%used = .false.
           end if
        end do

    end subroutine

    subroutine MCS_find_gases(window, gas, i_win)
      implicit none
      type(CS_window), intent(inout) :: window(:)
      type(CS_gas), intent(inout) :: gas(:)
      integer, intent(in) :: i_win

      integer :: i, j
      logical :: gas_found
      character(len=*), parameter :: fname = "MCS_find_gases"
      character(len=999) :: tmp_str

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
               call logger%trivia(fname, trim(tmp_str))
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



end module
