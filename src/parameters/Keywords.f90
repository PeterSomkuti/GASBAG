!! This module contains the list of accepted keywords for the config file. They
!! are hard-coded in here, so the config_check subroutine can check against this
!! list. This is a bit on the annoying side, but we want to make sure that
!! no superfluous item is present in the config file.

module keywords

    use stringifor

    implicit none

    type(string) :: valid_sections(10)

    public :: initialize_valid_sections

contains

    subroutine initialize_valid_sections()
        ! Everything is lowercase

        valid_sections(1) = "logger"
        valid_sections(2) = "algorithm"
    end subroutine


end module
