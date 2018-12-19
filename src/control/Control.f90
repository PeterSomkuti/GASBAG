!! Main control structure (MCS) of the program, for easy access of important
!! quantities throughout the program, such as instrument name and retrieval
!! settings, algorithm modes, .. the whole shebang.
!! This is designed to be instrument-independent, so apart from the name of the
!! instrument, no instrument-specifics should be stored here.

module control

    use stringifor
    use finer, only: file_ini
    use logger_mod, only: logger => master_logger
    use file_utils, only: check_config_files_exist
    use HDF5

    implicit none

    ! At the moment, we only plan the Guanter and Frankenberg methods
    integer, parameter :: MAX_ALGORITHMS = 2
    ! Why would you even need more than 2?
    integer, parameter :: MAX_WINDOWS = 5

    type, private :: CS_general
        integer(8) :: N_soundings ! Number of soundings to be processed
    end type

    type, private :: CS_algorithm
        type(string) :: name(MAX_ALGORITHMS) ! Name of the algorithm(s) used?
        integer :: N_algorithms ! How many do we want to actually use?
        integer :: N_basisfunctions ! How mamy basis functions do we read in (and maybe use)?
        logical :: using_GK_SIF
        logical :: using_physical
    end type

    type, private :: CS_window
        logical :: used ! Is this CS_window structure used?
        type(string) :: name
        double precision :: wl_min
        double precision :: wl_max
        type(string) :: basisfunction_file
    end type

    type, private :: CS_input
        type(string) :: L1B_filename ! Path to the L1B file
        integer(hid_t) :: l1b_file_id ! HDF file handler of L1B file
        type(string) :: MET_filename ! Path to MET file
        type(string) :: instrument_name ! Name of the instrument
    end type

    type, private :: CS_output
        type(string) :: output_filename ! Where does the ouptut HDF file go?
        integer(hid_t) :: output_file_id
    end type

    ! Main control structure type
    type, private :: CS
        type(CS_algorithm) :: algorithm ! Algorithm/forwared model - related settings
        type(CS_window), dimension(MAX_WINDOWS) :: window
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

        integer :: window_nr
        integer :: i
        logical :: file_exists


        call logger%debug(fname, "Populating main program control structure..")


        ! ----------------------------------------------------------------------
        ! First, we set all those fields to -1 values, that are added/populated
        ! later in e.g. instrument-specific rountines.

        MCS%general%N_soundings = -1

        MCS%algorithm%using_GK_SIF = .false.
        MCS%algorithm%using_physical = .false.

        MCS%window%name = ""
        MCS%window%wl_min = 0.0d0
        MCS%window%wl_max = 0.0d0

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
            MCS%input%L1B_filename = trim(fini_char)
        end if

        ! ----------------------------------------------------------------------

        ! Outputs section -------------------------------------------------------
        call fini%get(section_name='output', option_name='output_file', &
                      val=fini_char, error=fini_error)
        if (fini_error /= 0) then
            call logger%fatal(fname, "Error reading output file name")
            stop 1
        end if

        MCS%output%output_filename = trim(fini_char)


        ! ----------------------------------------------------------------------

        ! Instrument section ---------------------------------------------------

        ! Get instrument name
        call fini%get(section_name='instrument', option_name='name', &
                      val=fini_char, error=fini_error)
        if (fini_error /= 0) then
            call logger%fatal(fname, "Error reading instrument name")
            stop 1
        end if

        MCS%input%instrument_name = trim(fini_char)
        ! ----------------------------------------------------------------------


        ! Windows section ------------------------------------------------------
        ! The user might specify "window-2", and "window-5", so we need to check
        ! many possible windows here.

        do window_nr = 1, MAX_WINDOWS

            ! Is window "window_nr" in the config-file?
            write(tmp_str, '(A, G0.1)') "window-", window_nr
            if (fini%has_section(section_name=tmp_str)) then

                MCS%window(window_nr)%used = .true.

                call fini%get(section_name=tmp_str, option_name='name', &
                              val=fini_char, error=fini_error)
                if (fini_error /= 0) then
                    call logger%fatal(fname, "Could not read name in window section!")
                    stop 1
                else
                    MCS%window(window_nr)%name = trim(fini_char)
                end if

                call fini%get(section_name=tmp_str, option_name='wl_min', &
                              val=fini_val, error=fini_error)
                if (fini_error /= 0) then
                    call logger%fatal(fname, "Could not read wl_min in window section!")
                    stop 1
                else
                    MCS%window(window_nr)%wl_min = fini_val
                end if

                call fini%get(section_name=tmp_str, option_name='wl_max', &
                              val=fini_val, error=fini_error)
                if (fini_error /= 0) then
                    call logger%fatal(fname, "Could not read wl_ax in window section!")
                    stop 1
                else
                    MCS%window(window_nr)%wl_max = fini_val
                end if

                call fini%get(section_name=tmp_str, option_name='basisfunctions', &
                              val=fini_char, error=fini_error)
                if (fini_error /= 0) then
                    call logger%fatal(fname, "Could not read wl_ax in window section!")
                    stop 1
                else
                    MCS%window(window_nr)%basisfunction_file = trim(fini_char)
                end if
            else
                MCS%window(window_nr)%used = .false.
            end if
        end do

        ! ----------------------------------------------------------------------

    end subroutine




end module
