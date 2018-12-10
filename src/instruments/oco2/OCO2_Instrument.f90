module oco2

    use stringifor
    use instruments, only: generic_instrument

    use HDF5
    use H5LT

    implicit none

    type, extends(generic_instrument), public :: oco2_instrument

        ! dispersion coefficients, noise model, ILS array, ...
        ! all of that should be read in immediately at the start of the
        ! program, and kept in memory for quick access

    contains
        procedure, nopass :: check_l1b_validity
    end type


contains

    subroutine check_l1b_validity(l1b_file)
        !! Checks whether this is indeed a valid OCO-2 l1b file that has all
        !! the required fields and variables..

        type(string) :: l1b_file

        integer(hid_t) :: file_id
        integer :: hdferr

        ! Open the HDF file
        call h5fopen_f(l1b_file%chars(), H5F_ACC_RDONLY_F, file_id, hdferr)

    end subroutine
end module
