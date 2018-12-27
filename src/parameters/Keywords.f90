!! This module contains the list of accepted keywords for the config file. They
!! are hard-coded in here, so the config_check subroutine can check against this
!! list. This is a bit on the annoying side, but we want to make sure that
!! no superfluous item is present in the config file.

module keywords

    use stringifor, only: string
    use control, only: MAX_WINDOWS

    implicit none

    ! Just parameters/dimensions to generate the string arrays
    integer, parameter :: max_sections = 50
    integer, parameter :: max_options = 99

    type(string) :: valid_sections(max_sections)
    type(string) :: valid_options(max_sections, max_options)

    public :: initialize_valid_sections

contains

    subroutine initialize_valid_sections()

        implicit none
        character(len=99) :: tmp_str
        integer :: window_start, window_nr, this_idx

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
            valid_options(2,3) = "solar_dopper_shift"

        ! Related to the instrument in question
        valid_sections(3) = "instrument"
            ! Which one?
            valid_options(3,1) = "name"

        ! Related to input files
        valid_sections(4) = "input"
            valid_options(4,1) = "l1b_file" ! L1b file location
            valid_options(4,2) = "met_file" ! MET file location

        ! Output file options
        valid_sections(5) = "output"
            valid_options(5,1) = "output_file"
            valid_options(5,2) = "save_radiances"

        ! Solar model type and file path
        valid_sections(6) = "solar"
            valid_options(6,1) = "solar_type"
            valid_options(6,2) = "solar_file"

        ! We have to define our windows manually here, and give them explicit
        ! numbers!
        window_start = 6
        do window_nr=1, MAX_WINDOWS
            write(tmp_str, '(A, G0.1)') "window-", window_nr
            this_idx = window_start + window_nr
                valid_sections(this_idx) = trim(tmp_str)
                valid_options(this_idx,1) = "name"
                valid_options(this_idx,2) = "wl_min"
                valid_options(this_idx,3) = "wl_max"
                valid_options(this_idx,4) = "basisfunctions"
                valid_options(this_idx,5) = "gases"
        end do


    end subroutine


end module
