!! Main control structure (CS) of the program, for easy access of important
!! quantities throughout the program, such as instrument and retrieval
!! settings, algorithm modes, .. the whole shebang.

module control

    use stringifor
    use finer, only: file_ini
    use logger_mod, only: logger => master_logger
    use file_utils, only: check_config_files_exist
    implicit none

    ! At the moment, we only plan the Guanter and Frankenberg methods
    integer, parameter :: MAX_ALGORITHMS = 2

    type, private :: CS_algorithm
        type(string) :: name(MAX_ALGORITHMS) ! Name of the algorithm(s) used?
        integer :: N_algorithms ! Now many do we want to actually use?
    end type

    type, private :: CS_input
        type(string) :: L1B_filename ! Path to the L1B file
        type(string) :: MET_filename ! Path to MET file
    end type

    type, private :: CS_output
        type(string) :: output_filename ! Where does the ouptut HDF file go?
    end type

    ! Main control structure type
    type, private :: CS
        type(CS_algorithm) :: algorithm ! Algorithm-related settings
        type(CS_input) :: input ! Input files needed by the program
        type(CS_output) :: output ! Output file path(s)
    end type

    ! Define it here. Rest of the code should be allowed to change data?
    type(CS), public :: MCS

    public populate_control_structure

contains

    subroutine populate_control_structure(fini)
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

        integer :: i
        logical :: file_exists

        call logger%debug(fname, "Populating main program control structure..")

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
            MCS%input%L1B_filename = trim(fini_char)
        end if

        ! ----------------------------------------------------------------------



    end subroutine




end module