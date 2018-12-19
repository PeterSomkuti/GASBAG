!! Contains various helpers to deal with file-related matters, such as
!! checking whether files exist

module file_utils

    use HDF5
    use logger_mod, only: logger => master_logger

    implicit none


    ! Interface for writing double precision arrays into HDF files, subroutines
    ! exist for several dimensionalities. The output is still always single
    ! precision, and the conversion is done within the subroutines.
    interface write_DP_hdf_dataset
        module procedure write_2D_DP_hdf_dataset
        module procedure write_3D_DP_hdf_dataset
    end interface write_DP_hdf_dataset

    public check_hdf_error, check_config_files_exist, get_HDF5_dset_dims, &
           write_DP_hdf_dataset

contains


    subroutine write_2D_DP_hdf_dataset(file_id, dset_name, array, dims, fill_value)

        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: dset_name
        double precision, dimension(:,:), intent(in) :: array
        integer(hsize_t), dimension(:), intent(in) :: dims
        double precision, optional, intent(in) :: fill_value

        real, dimension(:,:), allocatable :: conv_array
        character(len=*), parameter :: fname = "write_2D_DP_hdf_dataset"
        integer :: hdferr
        integer(hid_t) :: dspace_id, dset_id, dcpl

        allocate(conv_array(size(array, 1), size(array, 2)))
        conv_array = real(array)

        include "HDF5_write_DP_array.inc"

    end subroutine

    subroutine write_3D_DP_hdf_dataset(file_id, dset_name, array, dims, fill_value)

        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: dset_name
        double precision, dimension(:,:,:), intent(in) :: array
        integer(hsize_t), dimension(:), intent(in) :: dims
        double precision, optional, intent(in) :: fill_value

        real, dimension(:,:,:), allocatable :: conv_array
        character(len=*), parameter :: fname = "write_3D_DP_hdf_dataset"
        integer :: hdferr
        integer(hid_t) :: dspace_id, dset_id, dcpl

        allocate(conv_array(size(array, 1), size(array, 2), size(array, 3)))
        conv_array = real(array)

        include "HDF5_write_DP_array.inc"

    end subroutine



    subroutine check_hdf_error(hdferr, fname, msg)
        integer, intent(in) :: hdferr
        character(len=*), intent(in) :: fname
        character(len=*), intent(in) :: msg

        if (hdferr < 0) then
            call logger%fatal(fname, msg)
            stop 1
        end if
    end subroutine


    subroutine get_HDF5_dset_dims(file_id, dset_name, dset_dims)

        ! Takes a HDF5 file_id, looks for the dataset with dset_name, and
        ! reads the shape of the dataset array and puts it into the allocatable
        ! dset_dims variable

        integer(hid_t), intent(in) :: file_id
        character(len=*), intent(in) :: dset_name
        integer(kind=8), allocatable, intent(out) :: dset_dims(:)

        character(len=*), parameter :: fname = "get_HDF5_dset_dims"

        integer(hid_t) :: dset_id, dspace_id
        integer :: ndims
        integer(hsize_t), allocatable :: dim_dset(:), max_dim_dset(:)
        integer :: hdferr


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

        call h5sget_simple_extent_ndims_f(dspace_id, ndims, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error getting ndims for " // trim(dset_name))
            stop 1
        end if

        allocate(dim_dset(ndims))
        allocate(max_dim_dset(ndims))
        allocate(dset_dims(ndims))

        call h5sget_simple_extent_dims_f(dspace_id, dim_dset, max_dim_dset, hdferr)
        if (hdferr == -1) then
            ! Note, this function returns the non-negative dataspace rank as the
            ! error code on success. -1 on a fail.
            call logger%fatal(fname, "Error getting dataspace extents for " // trim(dset_name))
            stop 1
        end if

        dset_dims(:) = dim_dset(:)

        call h5dclose_f(dset_id, hdferr)
        if (hdferr /= 0) then
            call logger%fatal(fname, "Error closing dataset " // trim(dset_name))
        end if

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
