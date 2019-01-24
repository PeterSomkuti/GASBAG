!! This module contains the list of accepted keywords for the config file. They
!! are hard-coded in here, so the config_check subroutine can check against this
!! list. This is a bit on the annoying side, but we want to make sure that
!! no superfluous item is present in the config file.

module keywords_mod

    use stringifor, only: string
    use control_mod, only: MAX_WINDOWS, MAX_GASES

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
        integer :: window_start, window_nr, gas_start, gas_nr, this_idx

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
                valid_options(this_idx,6) = "statevector"
                valid_options(this_idx,7) = "albedo_apriori"
                valid_options(this_idx,8) = "albedo_order"
                valid_options(this_idx,9) = "dispersion_order"
                valid_options(this_idx,10) = "dispersion_perturbation"
                valid_options(this_idx,11) = "dispersion_covariance"
                valid_options(this_idx,12) = "atmosphere"
                valid_options(this_idx,13) = "algorithms"
                valid_options(this_idx,14) = "dsigma_scale"
                valid_options(this_idx,15) = "fft_convolution"
                valid_options(this_idx,16) = "wl_spacing"
                valid_options(this_idx,17) = "band"
        end do

        ! Section for gases
        gas_start = this_idx
        do gas_nr=1, MAX_GASES
           write(tmp_str, '(A, G0.1)') "gas-", gas_nr
           this_idx = gas_start + gas_nr

           valid_sections(this_idx) = trim(tmp_str)
               valid_options(this_idx,1) = "name"
               valid_options(this_idx,2) = "spectroscopy_type"
               valid_options(this_idx,3) = "spectroscopy_file"
        end do

    end subroutine


end module
