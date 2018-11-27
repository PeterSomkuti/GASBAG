module Keywords

    use stringifor

    implicit none

    type(string) :: valid_sections(10)

    public :: initialize_valid_sections

contains

    subroutine initialize_valid_sections()
        ! Everything is lowercase

        valid_sections(1) = "logger"
    end subroutine


end module
