!! Contains various helpers to deal with file-related matters, such as
!! checking whether files exist

module file_utils_mod

  use HDF5
  use logger_mod, only: logger => master_logger
  use finer, only: file_ini
  use stringifor, only: string
  use physical_model_mod, only: atmosphere

  use iso_fortran_env

  implicit none


  ! Interface for writing double precision arrays into HDF files, subroutines
  ! exist for several dimensionalities. The output is still always single
  ! precision, and the conversion is done within the subroutines.
  interface write_DP_hdf_dataset
     module procedure write_2D_DP_hdf_dataset
     module procedure write_3D_DP_hdf_dataset
  end interface write_DP_hdf_dataset

  interface read_DP_hdf_dataset
     module procedure read_1D_DP_hdf_dataset
     module procedure read_2D_DP_hdf_dataset
     module procedure read_3D_DP_hdf_dataset
     module procedure read_4D_DP_hdf_dataset
  end interface read_DP_hdf_dataset

  interface write_INT_hdf_dataset
     module procedure write_2D_INT_hdf_dataset
  end interface write_INT_hdf_dataset

  interface read_INT_hdf_dataset
     module procedure read_3D_INT_hdf_dataset
  end interface read_INT_hdf_dataset


  interface fini_extract
     module procedure fini_extract_DP
     module procedure fini_extract_DP_array
     module procedure fini_extract_INT
     module procedure fini_extract_CHAR
     module procedure fini_extract_STRING_array
  end interface fini_extract

  public check_hdf_error, check_config_files_exist, get_HDF5_dset_dims, &
       check_fini_error, fini_extract, &
       write_DP_hdf_dataset, read_DP_hdf_dataset, &
       write_INT_hdf_dataset, read_INT_hdf_dataset

contains

  subroutine write_2D_INT_hdf_dataset(file_id, dset_name, array, dims, fill_value)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    integer, dimension(:,:), intent(in) :: array
    integer(hsize_t), dimension(:), intent(in) :: dims
    integer, optional, intent(in) :: fill_value

    integer, dimension(:,:), allocatable :: conv_array
    integer :: conv_fill_value
    character(len=*), parameter :: fname = "write_2D_DP_hdf_dataset"
    integer :: hdferr
    integer(hid_t) :: dspace_id, dset_id, dcpl

    allocate(conv_array(size(array, 1), size(array, 2)))
    conv_array = int(array)
    if (present(fill_value)) then
       conv_fill_value = int(fill_value)
    end if

    include "HDF5_write_INT_array.inc"

  end subroutine write_2D_INT_hdf_dataset


  subroutine read_3D_INT_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    integer, dimension(:,:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_3D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3)))

    call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_3D_INT_hdf_dataset

  subroutine read_4D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, allocatable, intent(inout) :: array(:,:,:,:)
    integer(hsize_t), allocatable, intent(inout) :: dset_dims(:)

    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_4D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))

    call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_4D_DP_hdf_dataset

  subroutine read_3D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:,:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_3D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3)))

    call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_3D_DP_hdf_dataset

  subroutine read_2D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_2D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    allocate(array(dset_dims(1), dset_dims(2)))

    call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_2D_DP_hdf_dataset

  subroutine read_1D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_1D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    allocate(array(dset_dims(1)))

    call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_1D_DP_hdf_dataset


  subroutine write_2D_DP_hdf_dataset(file_id, dset_name, array, dims, fill_value)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:,:), intent(in) :: array
    integer(hsize_t), dimension(:), intent(in) :: dims
    double precision, optional, intent(in) :: fill_value

    real, dimension(:,:), allocatable :: conv_array
    real :: conv_fill_value
    character(len=*), parameter :: fname = "write_2D_DP_hdf_dataset"
    integer :: hdferr
    integer(hid_t) :: dspace_id, dset_id, dcpl

    allocate(conv_array(size(array, 1), size(array, 2)))
    conv_array = real(array)
    if (present(fill_value)) then
       conv_fill_value = real(fill_value)
    end if

    include "HDF5_write_DP_array.inc"

  end subroutine write_2D_DP_hdf_dataset

  subroutine write_3D_DP_hdf_dataset(file_id, dset_name, array, dims, fill_value)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:,:,:), intent(in) :: array
    integer(hsize_t), dimension(:), intent(in) :: dims
    double precision, optional, intent(in) :: fill_value

    real, dimension(:,:,:), allocatable :: conv_array
    real :: conv_fill_value
    character(len=*), parameter :: fname = "write_3D_DP_hdf_dataset"
    integer :: hdferr
    integer(hid_t) :: dspace_id, dset_id, dcpl

    allocate(conv_array(size(array, 1), size(array, 2), size(array, 3)))
    conv_array = real(array)
    if (present(fill_value)) then
       conv_fill_value = real(fill_value)
    end if

    include "HDF5_write_DP_array.inc"

  end subroutine write_3D_DP_hdf_dataset



  subroutine check_hdf_error(hdferr, fname, msg)
    integer, intent(in) :: hdferr
    character(len=*), intent(in) :: fname, msg

    if (hdferr < 0) then
       call logger%fatal(fname, msg)
       stop 1
    end if
  end subroutine check_hdf_error

  subroutine check_fini_error(fini_error, fname, msg)
    integer, intent(in) :: fini_error
    character(len=*), intent(in) :: fname, msg

    if (fini_error /= 0) then
       call logger%fatal(fname, msg)
       stop 1
    end if
  end subroutine check_fini_error


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

  end subroutine get_HDF5_dset_dims


  subroutine fini_extract_INT(fini, section_name, option_name, required, value)
    integer, intent(inout) :: value
!!!!
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: section_name
    character(len=*), intent(in) :: option_name

    logical, intent(in) :: required
!!!
    character(len=len(section_name)) :: section_name_fini
    character(len=len(option_name)) :: option_name_fini

    logical :: found_option
    character(len=*), parameter :: fname = "fini_extract"
    integer :: fini_error, section_index

    value = -9999
    include "Fini_extract.inc"

  end subroutine fini_extract_INT

  subroutine fini_extract_DP(fini, section_name, option_name, required, value)
    double precision, intent(inout) :: value
!!!!
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: section_name
    character(len=*), intent(in) :: option_name

    logical, intent(in) :: required
!!!
    character(len=len(section_name)) :: section_name_fini
    character(len=len(option_name)) :: option_name_fini

    logical :: found_option
    character(len=*), parameter :: fname = "fini_extract"
    integer :: fini_error, section_index

    value = -9999.99d0
    include "Fini_extract.inc"

  end subroutine fini_extract_DP

  subroutine fini_extract_DP_array(fini, section_name, option_name, required, &
       value_array)
    double precision, allocatable, intent(inout) :: value_array(:)
!!!!
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: section_name
    character(len=*), intent(in) :: option_name
    character(len=999) :: value

    logical, intent(in) :: required
!!!
    character(len=len(section_name)) :: section_name_fini
    character(len=len(option_name)) :: option_name_fini

    logical :: found_option
    character(len=*), parameter :: fname = "fini_extract"
    character(len=999) :: tmp_str
    integer :: fini_error, section_index
    type(string) :: fini_string
    type(string), allocatable :: split_strings(:)
    integer :: i

    include "Fini_extract.inc"
    ! This gives us a character variable with all the numbers in it

    ! Convert to a string object
    fini_string = trim(value)
    ! And use the split function to get the single items out
    call fini_string%split(tokens=split_strings, sep=' ')
    ! Now allocate a double precision object according to the number of
    ! splitted strings
    allocate(value_array(size(split_strings)))
    ! And cast the objects into double precision values
    do i=1, size(split_strings)
       tmp_str = trim(split_strings(i)%chars())
       read(tmp_str, *) value_array(i)
    end do

  end subroutine fini_extract_DP_array

  subroutine fini_extract_CHAR(fini, section_name, option_name, required, value)
    character(len=*), intent(inout) :: value
!!!!
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: section_name
    character(len=*), intent(in) :: option_name

    logical, intent(in) :: required
!!!
    character(len=len(section_name)) :: section_name_fini
    character(len=len(option_name)) :: option_name_fini

    logical :: found_option
    character(len=*), parameter :: fname = "fini_extract"
    integer :: fini_error

    value = ""
    include "Fini_extract.inc"

  end subroutine fini_extract_CHAR

  subroutine fini_extract_STRING_array(fini, section_name, option_name, required, &
       value_array)
    type(string), allocatable, intent(inout) :: value_array(:)
!!!!
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: section_name
    character(len=*), intent(in) :: option_name
    character(len=999) :: value

    logical, intent(in) :: required
!!!
    character(len=len(section_name)) :: section_name_fini
    character(len=len(option_name)) :: option_name_fini

    logical :: found_option
    character(len=*), parameter :: fname = "fini_extract"
    character(len=999) :: tmp_str
    integer :: fini_error, section_index
    type(string) :: fini_string
    type(string), allocatable :: split_strings(:)
    integer :: i

    include "Fini_extract.inc"
    ! This gives us a character variable with all the substrings in it

    ! Convert to a string object
    fini_string = trim(value)
    ! And use the split function to get the single items out
    call fini_string%split(tokens=split_strings, sep=' ')
    ! Now allocate a double precision object according to the number of
    ! splitted strings
    allocate(value_array(size(split_strings)))
    ! And cast the objects into double precision values
    do i=1, size(split_strings)
       value_array(i) = trim(split_strings(i)%chars())
    end do

  end subroutine fini_extract_STRING_array


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

  end subroutine check_config_files_exist

  subroutine read_atmosphere_file(filename, gas_indices, gas_strings, atm)

    implicit none
    character(len=*), intent(in) :: filename
    integer, intent(in) :: gas_indices(:)
    type(string), intent(in) :: gas_strings(:)
    type(atmosphere), intent(inout) :: atm

    character(len=*), parameter :: fname = "read_atmosphere_file"
    integer :: funit, iostat
    logical :: file_exist
    integer :: line_count, nonempty_line_count, file_start, level_count
    integer :: idx_p, idx_t
    character(len=999) :: dummy, tmp_str
    double precision :: dummy_dp
    type(string) :: dummy_string
    type(string), allocatable :: split_string(:)

    integer :: i,j,cnt
    integer, allocatable :: this_gas_index(:)

    integer :: num_gases

    inquire(file=filename, exist=file_exist)
    if (.not. file_exist) then
       call logger%fatal(fname, "Atmosphere file does not exist: " // filename)
       stop 1
    else
       call logger%debug(fname, "File does exist.")
    end if

    ! First pass: we scan the file to see how many levels our atmosphere has
    open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')
    rewind(unit=funit, iostat=iostat)

    line_count = 0
    level_count = 0
    nonempty_line_count = 0
    file_start = -1

    ! Loop through the file until we have reached EOF
    do
       read(funit, '(A)', iostat=iostat) dummy

       if (iostat == iostat_end) then
          ! End of file?
          exit
       end if

       line_count = line_count + 1

       if (scan(dummy, "!#;") > 0) then
          ! Skip this line, as it's commented
          cycle
       else if (trim(dummy) == "") then
          ! Skip empty lines
          cycle
       end if

       nonempty_line_count = nonempty_line_count + 1

       ! First non-empty line should be saved here
       if (file_start == -1) then
          file_start = line_count
       else
          if (nonempty_line_count > 1) then
             level_count = level_count + 1
          end if
       end if

    end do

    !close(funit)
    !open(newunit=funit, file=filename, iostat=iostat, action='read', status='old')

    ! Go back to the top of the file.
    rewind(unit=funit, iostat=iostat)

    line_count = 0
    do
       read(funit, '(A)', iostat=iostat) dummy

       line_count = line_count + 1

       if (line_count < file_start) cycle

       if (line_count == file_start) then
          ! This is the proper atmosphere header that should contain
          ! the information about the gases. So first, let's check if
          ! the numbers match
          dummy_string = dummy
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now that we know both the number of levels and gases, we can allocate the
          ! arrays in the atmosphere structure.

          num_gases = 0
          idx_p = -1
          idx_t = -1
          do j=1, size(split_string)
             ! Skip temp or pressure
             if (split_string(j)%lower() == "p") then
                idx_p = j
                cycle
             end if
             if (split_string(j)%lower() == "t") then
                idx_t = j
                cycle
             end if
             num_gases = num_gases + 1
          end do

          write(tmp_str, '(A,G0.1,A,A)') "There seem to be ", num_gases, " gases in ", filename
          call logger%info(fname, trim(tmp_str))

          if (allocated(atm%p)) deallocate(atm%p)
          if (allocated(atm%T)) deallocate(atm%T)
          if (allocated(atm%sh)) deallocate(atm%sh)
          if (allocated(atm%gas_names)) deallocate(atm%gas_names)
          if (allocated(atm%gas_index)) deallocate(atm%gas_index)
          if (allocated(atm%gas_vmr)) deallocate(atm%gas_vmr)

          allocate(atm%T(level_count))
          allocate(atm%p(level_count))
          allocate(atm%sh(level_count))
          allocate(atm%gas_names(num_gases))
          allocate(atm%gas_vmr(level_count, num_gases))
          allocate(atm%gas_index(num_gases))



          allocate(this_gas_index(size(gas_strings)))

          ! But we also want to know what gas index to use for storage
          do i=1, size(gas_strings)
             this_gas_index(i) = -1
             cnt = 1
             do j=1, size(split_string)
                ! Skip temp or pressure
                if (split_string(j)%lower() == "p") cycle
                if (split_string(j)%lower() == "t") cycle
                ! Check if gas description matches gases we know
                if (split_string(j) == gas_strings(i)) then
                   this_gas_index(i) = j
                   write(tmp_str, '(A,A,A,G0.1)') "Index for atmosphere gas ", &
                        split_string(j)%chars(), ": ", j
                   call logger%debug(fname, trim(tmp_str))
                   exit
                end if
                cnt = cnt + 1
             end do
          end do

          ! And last check - do all required gases have a VMR column in the
          ! atmosphere file.

          do i=1, size(gas_strings)
             if (this_gas_index(i) == -1) then
                ! Uh-oh, gas that was speficied in the "gases" option of the window
                ! could not be found in the atmosphere! Hard exit.
                write(tmp_str, "(A,A)") "The following gas was not found in the atmosphere file: " &
                     // gas_strings(i)%chars()
                call logger%fatal(fname, trim(tmp_str))
                stop 1
             end if
          end do

       end if


       if (line_count > file_start) then
          ! Right after the header, we should have the data in rows. We
          ! use the string split option to split the row string into substrings,
          ! and then convert each into double precision values and feed them into
          ! the atmosphere structure - at the right position!

          dummy_string = dummy
          ! Need to deallocate the split_string object first
          if (allocated(split_string)) deallocate(split_string)
          call dummy_string%split(tokens=split_string, sep=' ')

          ! Now here we need to check again whether a certain line has more than
          ! num_gases+1 columns.
          if (size(split_string) /= (num_gases + 2)) then
             write(tmp_str, '(A, G0.1)') "Too many values in line ", line_count
             call logger%fatal(fname, trim(tmp_str))
             stop 1
          end if

          ! Get the pressure value
          tmp_str = split_string(idx_p)%chars()
          read(tmp_str, *) dummy_dp
          atm%p(line_count - file_start) = dummy_dp

          ! Get the temperature value
          tmp_str = split_string(idx_t)%chars()
          read(tmp_str, *) dummy_dp
          atm%T(line_count - file_start) = dummy_dp

          ! And the gas value(s) from the other column(s)
          do i=1, size(gas_strings)
             tmp_str = split_string(this_gas_index(i))%chars()
             read(tmp_str, *) dummy_dp
             atm%gas_vmr(line_count - file_start, i) = dummy_dp
          end do

       end if

       if (line_count == (file_start + level_count)) exit

    end do

    close(funit)

  end subroutine read_atmosphere_file


end module file_utils_mod
