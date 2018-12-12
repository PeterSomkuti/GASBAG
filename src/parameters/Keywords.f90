!! This module contains the list of accepted keywords for the config file. They
!! are hard-coded in here, so the config_check subroutine can check against this
!! list. This is a bit on the annoying side, but we want to make sure that
!! no superfluous item is present in the config file.

module keywords

    use stringifor, only: string

    implicit none

    ! Just parameters/dimensions to generate the string arrays
    integer, parameter :: max_sections = 10
    integer, parameter :: max_options = 99

    type(string) :: valid_sections(max_sections)
    type(string) :: valid_options(max_sections, max_options)

    public :: initialize_valid_sections

contains

    subroutine initialize_valid_sections()
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        ! Everything is lowercase here!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

        ! Related to logging, messaging, verbosity
        valid_sections(1) = "logger"
            ! Where to write the logfile to?
            valid_options(1,1) = "logfile"
            ! What logging level are we using?
            valid_options(1,2) = "loglevel"


        ! Related to the algorithm general setup
        valid_sections(2) = "algorithm"
            ! Which SIF algorithm(s) to use?
            valid_options(2,1) = "sif_algorithm"
            valid_options(2,2) = "n_basisfunctions"

        ! Related to the instrument in question
        valid_sections(3) = "instrument"
            ! Which one?
            valid_options(3,1) = "name"

        ! Related to input files
        valid_sections(4) = "input"
            valid_options(4,1) = "l1b_file" ! L1b file location
            valid_options(4,2) = "met_file" ! MET file location

        valid_sections(5) = "window"
            valid_options(5,1) = "name"
            valid_options(5,2) = "wl_min"
            valid_options(5,3) = "wl_max"
            valid_options(5,4) = "basisfunctions"


    end subroutine


end module
