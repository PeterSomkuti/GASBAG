!! Contains various helpers to deal with file-related matters, such as
!! checking whether files exist

module file_utils

    use HDF5
    use logger_mod, only: logger => master_logger

    implicit none
    private

    interface get_HDF5_dset_dims
        module procedure get_HDF5_dset_dims_2d
    end interface


    public check_config_files_exist, get_HDF5_dset_dims

contains


    subroutine get_HDF5_dset_dims_2d(file_id, dset_name, dset_dims)

        ! Takes a HDF5 file_id, looks for the dataset with dset_name, and
        ! reads the shape of the dataset array and puts it into dset_dims

        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: dset_name
        integer(kind=8), intent(out) :: dset_dims(2)

        character(len=*), parameter :: fname = "get_HDF5_dset_dims_2d"

        integer(hid_t) :: dset_id, dspace_id
        integer(hsize_t) :: dim_dset_2d(2), max_dim_dset_2d(2)
        integer :: hdferr

        ! First establish how many measurements we have in total by interr-
        ! ogating the sounding_id list
        call h5dopen_f(file_id, dset_name, dset_id, hdferr)
        if (hdferr /= 0) then
          call logger%fatal(fname, "Error reading field: " // trim(dset_name))
          stop 1
        end if

        call h5dget_space_f(dset_id, dspace_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error getting dataspace field for " // trim(dset_name))
            stop 1
        end if

        call h5sget_simple_extent_dims_f(dspace_id, dim_dset_2d, max_dim_dset_2d, hdferr)
        if (hdferr == -1) then
            ! Note, this function returns the non-negative dataspace rank as the
            ! error code on success. -1 on a fail.
            call logger%fatal(fname, "Error getting dataspace extents for " // trim(dset_name))
            stop 1
        end if

        dset_dims(:) = dim_dset_2d(:)

    end subroutine




    subroutine check_config_files_exist(fini, section_name, option_name, &
                                        num_max_files, result)
    !! With this function, you can very easily check if the files
    !! supplied by a config file option actually exist or not.

        use finer, only: file_ini
        use stringifor
        use logger_mod, only: logger => master_logger

        type(file_ini), intent(in) :: fini
        character(len=*), intent(in) :: section_name, option_name
        integer, intent(in) :: num_max_files
        logical, intent(out) :: result

        character(len=*), parameter :: fname = "check_config_files_exist"

        ! FINER stuff
        integer :: fini_error, fini_count
        character(len=999) :: fini_char, tmp_str
        type(string) :: fini_string
        type(string), allocatable :: file_strings(:)
        logical :: file_exists

        integer :: i

        result = .true.

        ! First, let's see if the number of entries in the config option
        ! are NOT more than the maximally allowed number
        write(tmp_str, '(I3.1)') num_max_files
        fini_count = fini%count_values(section_name=section_name, &
                                       option_name=option_name)
        if (fini_count > num_max_files) then
            call logger%error(fname, "Number of values in [" // &
                              section_name // "/" // &
                              option_name // "] is greater than " // trim(tmp_str))
            result = .false.
        end if

        ! Check if the file(s) is/are actually exist(s)
        ! Grab filename. But ONLY if we actually have values. There
        ! might be cases, where these options are allowed to be left
        ! blank (empty).
        if (fini_count > 0) then
            call fini%get(section_name=section_name, &
                          option_name=option_name, &
                          val=fini_char, error=fini_error)
            if (fini_error /= 0) then
                call logger%fatal(fname, "Failure to get option value for " // &
                                         "input/l1b_file")
                stop 1
            end if
            ! We can now split the option value (space delimited) and go
            ! through every entry, and see if they are proper files.
            allocate(file_strings(fini_count))
            fini_string = trim(fini_char)
            call fini_string%split(tokens=file_strings, sep=' ', &
                                   max_tokens=fini_count)
            do i=1, fini_count
                ! Loop through all supplied file paths..
                inquire(file=trim(file_strings(i)%chars()), exist=file_exists)
                ! .. and check if they exist
                if (.not. file_exists) then
                    ! No? Throw an error message to the screen and
                    ! set result accordinly.
                    result = .false.
                    call logger%error(fname, trim(file_strings(i)%chars()) &
                                      // " does not exist!")
                else
                    call logger%trivia(fname, trim(file_strings(i)%chars()) &
                                       // " was found.")
                end if
            end do
        end if

    end subroutine



end module
