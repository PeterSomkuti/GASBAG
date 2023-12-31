!! Contains various helpers to deal with file-related matters, such as
!! checking whether files exist

module file_utils_mod

  use HDF5
  use logger_mod, only: logger => master_logger
  use finer, only: file_ini
  use stringifor, only: string

  use iso_c_binding
  use iso_fortran_env
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

  implicit none

  integer :: MAX_HDF_ATTEMPTS = 10

  !> Interface for writing double precision arrays into HDF files
  interface write_DP_hdf_dataset
     module procedure write_2D_DP_hdf_dataset
     module procedure write_3D_DP_hdf_dataset
  end interface write_DP_hdf_dataset

  !> Interface for reading double precision arrays from an HDF file
  interface read_DP_hdf_dataset
     module procedure read_1D_DP_hdf_dataset
     module procedure read_2D_DP_hdf_dataset
     module procedure read_3D_DP_hdf_dataset
     module procedure read_4D_DP_hdf_dataset
     module procedure read_5D_DP_hdf_dataset
  end interface read_DP_hdf_dataset

  !> Interface for reading single precision (real) arrays from an HDF file
  interface read_SP_hdf_dataset
     module procedure read_2D_SP_hdf_dataset
     module procedure read_4D_SP_hdf_dataset
  end interface read_SP_hdf_dataset

  !> Interface for writing integer arrays into an HDF file
  interface write_INT_hdf_dataset
     module procedure write_2D_INT_hdf_dataset
  end interface write_INT_hdf_dataset

  !> Interface for reading integer arrays from an HDF file
  interface read_INT_hdf_dataset
     module procedure read_2D_INT_hdf_dataset
     module procedure read_3D_INT_hdf_dataset
  end interface read_INT_hdf_dataset

  !> Interface to extract value(s) out of the config file
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
       write_INT_hdf_dataset, read_INT_hdf_dataset, &
       read_one_arbitrary_value_dp, write_string_hdf_dataset

contains


  subroutine write_config_into_groups(file_id, fini, root_group)
    integer(hid_t), intent(in) :: file_id
    type(file_ini), intent(in) :: fini
    character(len=*), intent(in) :: root_group

    character(len=*), parameter :: fname = "write_config_into_groups"
    integer :: i
    logical :: group_exists
    integer(hid_t) :: gid
    integer :: hdferr
    character(len=:), allocatable :: item(:)
    character(len=:), allocatable :: sections(:)
    character(len=:), allocatable :: this_group
    character(len=:), allocatable :: this_dataset
    type(string) :: tmp_str


    ! Get the sections from the config option
    call fini%get_sections_list(sections)

    ! Check if the root group exists in the file
    call h5lexists_f(file_id, root_group, group_exists, hdferr)

    ! Create new group if it doesn't exist (it really shouldn't)
    if (.not. group_exists) then
       call h5gcreate_f(file_id, root_group, gid, hdferr)
       call check_hdf_error(hdferr, fname, "Error. Could not create group: " // trim(root_group))
    end if

    call logger%info(fname, "Writing configuration keys into file.")

    ! Loop through each section
    do i = 1, size(sections)

       ! Current group will be this:
       this_group = root_group // "/" // sections(i)

       ! Again, check for the existence for the group
       call h5lexists_f(file_id, this_group, group_exists, hdferr)
       if (group_exists) then
          ! This really should not happen, so flag it as warning
          call logger%warning(fname, "Group " // this_group // " already exists. ")
       else
          call h5gcreate_f(file_id, this_group, gid, hdferr)
          call check_hdf_error(hdferr, fname, "Error. Could not create group: " // trim(root_group))
       end if

       ! Now loop through each option in this section and save every key/value
       ! pair as a string dataset
       do while (fini%loop(section_name=sections(i), option_pairs=item))
          ! Construct group name and value
          this_dataset = trim(this_group) // "/" // item(1)
          tmp_str = item(2)
          ! Write it ito string dataset
          call write_string_hdf_dataset(file_id, this_dataset, tmp_str)
       end do
    end do

  end subroutine write_config_into_groups


  subroutine write_string_hdf_dataset(file_id, dset_name, str)
    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    type(string), intent(in) :: str

    character(len=:), allocatable, target :: c_arr
    integer(hid_t) :: dims(1) = 1
    integer(hid_t) :: c_len
    integer :: hdferr

    type(c_ptr) :: f_ptr
    integer(hid_t) :: filetype, dspace_id, dset_id

    c_arr = str%chars()
    c_len = len(c_arr)

    call h5tcopy_f(H5T_FORTRAN_S1, filetype, hdferr)
    call h5tset_size_f(filetype, c_len, hdferr)
    call h5screate_simple_f(1, dims, dspace_id, hdferr)
    call h5dcreate_f(file_id, dset_name, filetype, dspace_id, dset_id, hdferr)

    f_ptr = c_loc(c_arr)

    call h5dwrite_f(dset_id, filetype, f_ptr, hdferr)

    call h5dclose_f(dset_id, hdferr)
    call h5sclose_f(dspace_id, hdferr)
    call h5tclose_f(filetype, hdferr)

  end subroutine write_string_hdf_dataset

  subroutine write_2D_INT_hdf_dataset(file_id, dset_name, array, dims, fill_value)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    integer, dimension(:,:), intent(in) :: array
    integer(hsize_t), dimension(:), intent(in) :: dims
    integer, optional, intent(in) :: fill_value

    integer, dimension(:,:), allocatable :: conv_array
    integer :: conv_fill_value
    character(len=*), parameter :: fname = "write_2D_INT_hdf_dataset"
    integer :: hdferr
    integer(hid_t) :: dspace_id, dset_id, dcpl

    allocate(conv_array(size(array, 1), size(array, 2)))
    conv_array = int(array)
    if (present(fill_value)) then
       conv_fill_value = int(fill_value)
    end if

    include "HDF5_write_INT_array.inc"

  end subroutine write_2D_INT_hdf_dataset


  subroutine read_2D_INT_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    integer, dimension(:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_2D_INT_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    allocate(array(dset_dims(1), dset_dims(2)))

    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do

    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_2D_INT_hdf_dataset

  subroutine read_3D_INT_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    integer, dimension(:,:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_3D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3)))

    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do

    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_INTEGER, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_3D_INT_hdf_dataset

  subroutine read_5D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, allocatable, intent(inout) :: array(:,:,:,:,:)
    integer(hsize_t), allocatable, intent(inout) :: dset_dims(:)

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "read_5D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)

    if (size(dset_dims) /= 5) then
       write(tmp_str, '(A, G0.1)') "Error reading " // trim(dset_name) &
            // ". Expected 5 dimensions, but got: ", size(dset_dims)
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if

    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4), dset_dims(5)))

    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do

    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_5D_DP_hdf_dataset


  subroutine read_4D_SP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    real, allocatable, intent(inout) :: array(:,:,:,:)
    integer(hsize_t), allocatable, intent(inout) :: dset_dims(:)

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "read_4D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)

    if (size(dset_dims) /= 4) then
       write(tmp_str, '(A, G0.1)') "Error reading " // trim(dset_name) &
            // ". Expected 4 dimensions, but got: ", size(dset_dims)
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if

    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))

    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_4D_SP_hdf_dataset

  subroutine read_4D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, allocatable, intent(inout) :: array(:,:,:,:)
    integer(hsize_t), allocatable, intent(inout) :: dset_dims(:)

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "read_4D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)

    if (size(dset_dims) /= 4) then
       write(tmp_str, '(A, G0.1)') "Error reading " // trim(dset_name) &
            // ". Expected 4 dimensions, but got: ", size(dset_dims)
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if

    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3), dset_dims(4)))

    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_4D_DP_hdf_dataset

  subroutine read_3D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:,:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "read_3D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)

    if (size(dset_dims) /= 3) then
       write(tmp_str, '(A, G0.1)') "Error reading " // trim(dset_name) &
            // ". Expected 3 dimensions, but got: ", size(dset_dims)
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if

    allocate(array(dset_dims(1), dset_dims(2), dset_dims(3)))
    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_3D_DP_hdf_dataset


  subroutine read_2D_SP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    real, dimension(:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_2D_SP_dataset"
    character(len=999) :: tmp_str

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    if (size(dset_dims) /= 2) then
       write(tmp_str, '(A, G0.1)') "Error reading " // trim(dset_name) &
            // ". Expected 2 dimensions, but got: ", size(dset_dims)
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if

    allocate(array(dset_dims(1), dset_dims(2)))

    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_REAL, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_2D_SP_hdf_dataset

  subroutine read_2D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:,:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=*), parameter :: fname = "read_2D_DP_dataset"
    character(len=999) :: tmp_str

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)
    if (size(dset_dims) /= 2) then
       write(tmp_str, '(A, G0.1)') "Error reading " // trim(dset_name) &
            // ". Expected 2 dimensions, but got: ", size(dset_dims)
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if

    allocate(array(dset_dims(1), dset_dims(2)))

    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do
    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, array, dset_dims, hdferr)
    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

  end subroutine read_2D_DP_hdf_dataset

  subroutine read_1D_DP_hdf_dataset(file_id, dset_name, array, dset_dims)

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    double precision, dimension(:), allocatable, intent(inout) :: array
    integer(hsize_t), dimension(:), allocatable, intent(inout) :: dset_dims

    integer :: i
    integer :: hdferr
    integer(hid_t) :: dset_id
    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "read_1D_DP_dataset"

    if (allocated(array)) deallocate(array)

    call get_HDF5_dset_dims(file_id, trim(dset_name), dset_dims)

    if (size(dset_dims) /= 1) then
       write(tmp_str, '(A, G0.1)') "Error reading " // trim(dset_name) &
            // ". Expected 1 dimension, but got: ", size(dset_dims)
       call logger%fatal(fname, trim(tmp_str))
       stop 1
    end if

    allocate(array(dset_dims(1)))
    hdferr = -1
    i = 1
    do while ((hdferr < 0) .and. (i <= MAX_HDF_ATTEMPTS))
       call h5dopen_f(file_id, trim(dset_name), dset_id, hdferr)
       i = i + 1
       if (hdferr < 0) then
          call sleep(3)
       else
          exit
       end if
    end do
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
    integer :: fini_error

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
    integer :: fini_error

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
    integer :: fini_error
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
    integer :: fini_error
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


  function string_to_bool(input_string) result(value)
    implicit none
    type(string) :: input_string
    logical :: value

    if (input_string%lower() == 'true') then
       value = .true.
    else if (input_string%lower() == 'false') then
       value = .false.
    else
       call logger%fatal("string_to_bool",&
            "Sorry - we only accept T/true and F/false.")
       stop 1
    end if

  end function string_to_bool


  !> @brief Reads an old-school aerosol moment file
  !> @param filename Momfile location
  !> @param wavelengths Number of wavelengths
  !> @param coefs The phase matrix Legendre coefficients
  subroutine read_mom_file(filename, wavelengths, coefs, max_coef, success)

    implicit none
    character(len=*) :: filename
    double precision, allocatable, intent(inout) :: wavelengths(:)
    ! Coefs: coef, matrix element, wavelength
    double precision, allocatable, intent(inout) :: coefs(:,:,:)
    ! Largest number of coeffs in file, needed for array allocation
    integer, intent(inout) :: max_coef
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "read_mom_file"
    ! Does this file exist?
    logical :: file_exist
    ! Total number of wavelengths contained in the file
    integer :: wl_count
    ! Temporary value
    integer :: tmp_value

    ! Loop variables
    integer :: i,j
    ! Other loop variables
    integer :: cnt_coef, cnt_wl, this_n_coef
    ! Dummy var needed to convert string to char
    character(len=999) :: dummy
    ! Strings for splitting
    type(string), allocatable :: split_string(:)
    ! String object to hold the entire file
    type(string) :: file_string
    ! String object to hold the file per-line
    type(string), allocatable :: line_strings(:)

    success = .false.

    ! Try and open the file - see if it exists:
    inquire(file=filename, exist=file_exist)
    if (.not. file_exist) then
       call logger%fatal(fname, "Mom file does not exist: " // trim(filename))
       stop 1
    end if

    ! Read the file if it exists
    call file_string%read_file(file=filename)
    call file_string%split(tokens=line_strings, sep=new_line('a'))

    ! Find out how many wavelengths this mom file contains, and how many
    ! coefficients each wavelength has. We can also already infer what the
    ! size of the full array needs to be by keeping track of the largest
    ! number of coefficients.
    wl_count = 0
    max_coef = 0
    do i=1, size(line_strings, dim=1)

       ! Skip comments
       if (scan(line_strings(i)%chars(), "!#;") > 0) then
          cycle
       end if

       if (allocated(split_string)) deallocate(split_string)
       call line_strings(i)%split(tokens=split_string, sep=' ')

       ! Some files have lines like:
       ! "wl mom" as their first line (no comment sign)
       ! so we need to take that into account too
       if (split_string(1) == "wl") cycle

       if (size(split_string) == 2) then
          dummy = split_string(2)%chars()
          read(dummy, *) tmp_value
          if (tmp_value > max_coef) max_coef = tmp_value
          wl_count = wl_count + 1
       end if
    end do

    max_coef = max_coef

    write(dummy, '(A, G0.1)') "Mom file has this many wavelengths: ", wl_count
    call logger%debug(fname, trim(dummy))
    write(dummy, '(A, G0.1)') "Maximum number of coefficients: ", max_coef
    call logger%debug(fname, trim(dummy))


    ! Now we know how many elements our coef array needs to have
    ! "6" is hardcoded here, but I guess we don't expect physics
    ! to change all that much..
    allocate(coefs(max_coef, 6, wl_count))
    allocate(wavelengths(wl_count))


    ! NOTE
    ! MOM files are organized in a strange way: for a band, the
    ! number of expansion coefficients at one band edge might
    ! different from the coefficients at the other band edge for
    ! the same band! Here, we fill them up with some small value,
    ! hoping that it will not impact the calculations. Otherwise,
    ! the RT code is prone to produce garbage if you leave these
    ! values as NaNs.

    coefs(:,:,:) = 1.0d-10 !IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    wavelengths(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

    ! Do a second sweep through the file, this time
    ! stick the values into the coef array.

    cnt_coef = 0
    cnt_wl = 0
    do i=1, size(line_strings, dim=1)

       ! Skip comments
       if (scan(line_strings(i)%chars(), "!#;") > 0) then
          cycle
       end if

       if (allocated(split_string)) deallocate(split_string)
       call line_strings(i)%split(tokens=split_string, sep=' ')

       ! Some files have lines like:
       ! "wl mom" as their first line (no comment sign)
       ! so we need to take that into account too
       if (split_string(1) == "wl") then
          cycle
       end if

       ! This line indicates a beginning of a new section, so we need to figure out
       ! the positions in the array etc.
       if (size(split_string) == 2) then
          cnt_wl = cnt_wl + 1

          ! Read wavelength
          dummy = split_string(1)%chars()
          read(dummy, *) wavelengths(cnt_wl)
          ! Read number of coefs for this section
          dummy = split_string(2)%chars()
          read(dummy, *) this_n_coef

          cnt_coef = 1
       else if (size(split_string) == 6) then
          ! This is a line with 6 coefficients
          do j=1, 6
             dummy = split_string(j)%chars()
             read(dummy, *) coefs(cnt_coef, j, cnt_wl)
          end do

          cnt_coef = cnt_coef + 1
       else
          call logger%debug(fname, "Unexpected line in moment file: " // filename)
          return
       end if

    end do

    success = .true.

  end subroutine read_mom_file


  !> @brief Reads an old-school aerosol Mie file
  !> @param filename Miefile location
  !> @param wavelengths Wavelengths
  !> @param qext Extinction efficiency per wavelength
  !> @param qsca Scattering efficiency per wavelength
  !> @param ssa Single-scatter albedo per wavelength
  !> @param sigma_ext Extinction cross section per wavelength [um^2]
  !> @param reff Effective radius per wavelength [um]
  subroutine read_mie_file(filename, &
       wavelengths, qext, qsca, &
       ssa, sigma_ext, reff, &
       success)

    implicit none
    character(len=*), intent(in) :: filename
    double precision, allocatable, intent(inout) :: wavelengths(:)
    double precision, allocatable, intent(inout) :: qext(:)
    double precision, allocatable, intent(inout) :: qsca(:)
    double precision, allocatable, intent(inout) :: ssa(:)
    double precision, allocatable, intent(inout) :: sigma_ext(:)
    double precision, allocatable, intent(inout) :: reff(:)

    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "read_mie_file"
    ! Does the file exist?
    logical :: file_exist
    ! Dummy variable
    character(len=999) :: dummy
    ! Strings for splitting
    type(string), allocatable :: split_string(:)
    ! String object to hold the entire file
    type(string) :: file_string
    ! String object to hold the file per-line
    type(string), allocatable :: line_strings(:)
    integer :: cnt_wl
    integer :: i

    success = .false.

    ! Try and open the file - see if it exists:
    inquire(file=filename, exist=file_exist)
    if (.not. file_exist) then
       call logger%fatal(fname, "Mie file does not exist: " // trim(filename))
       stop 1
    end if

    ! Read the file if it exists
    call file_string%read_file(file=filename)
    call file_string%split(tokens=line_strings, sep=new_line('a'))

    cnt_wl = 0
    do i=1, size(line_strings, dim=1)

       ! Skip comments
       if (scan(line_strings(i)%chars(), "!#;") > 0) then
          cycle
       end if

       if (allocated(split_string)) deallocate(split_string)
       call line_strings(i)%split(tokens=split_string, sep=' ')

       ! Some files have lines like:
       ! "wl mom" as their first line (no comment sign)
       ! so we need to take that into account too
       if (split_string(1) == "wl") cycle

       cnt_wl = cnt_wl + 1
    end do

    write(dummy, '(A, G0.1)') "Mie file has this many wavelengths: ", cnt_wl
    call logger%debug(fname, trim(dummy))

    ! Allocate containers since we know now how many..
    allocate(wavelengths(cnt_wl))
    wavelengths(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    allocate(qext(cnt_wl))
    qext(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    allocate(qsca(cnt_wl))
    qsca(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    allocate(ssa(cnt_wl))
    ssa(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    allocate(sigma_ext(cnt_wl))
    sigma_ext(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
    allocate(reff(cnt_wl))
    reff(:) = IEEE_VALUE(1D0, IEEE_QUIET_NAN)


    cnt_wl = 0
    do i=1, size(line_strings, dim=1)

       ! Skip comments
       if (scan(line_strings(i)%chars(), "!#;") > 0) then
          cycle
       end if

       if (allocated(split_string)) deallocate(split_string)
       call line_strings(i)%split(tokens=split_string, sep=' ')

       ! Some files have lines like:
       ! "wl mom" as their first line (no comment sign)
       ! so we need to take that into account too
       if (split_string(1) == "wl") cycle

       cnt_wl = cnt_wl + 1
       ! Read wavelength
       dummy = split_string(1)%chars()
       read(dummy, *) wavelengths(cnt_wl)
       ! Read Qext
       dummy = split_string(2)%chars()
       read(dummy, *) qext(cnt_wl)
       ! Read Qsca
       dummy = split_string(3)%chars()
       read(dummy, *) qsca(cnt_wl)
       ! Read SSA
       dummy = split_string(4)%chars()
       read(dummy, *) ssa(cnt_wl)

       ! Apparently mie files come in at least two
       ! variants - some have 4 columns, other have 6
       ! and include sigma_ext and effective radius
       if (size(split_string) == 6) then
          ! Read Sigma_ext
          dummy = split_string(5)%chars()
          read(dummy, *) sigma_ext(cnt_wl)
          ! Read Reff
          dummy = split_string(6)%chars()
          read(dummy, *) reff(cnt_wl)
       end if
    end do

    success = .true.

  end subroutine read_mie_file

  !> @brief Takes a single double precision value from an HDF file/variable
  !> @param file_id HDF5 file ID
  !> @param dset_name Dataset name
  !> @param i_fp Footprint index
  !> @param i_fr Frame index
  !> @param value Return value
  function read_one_arbitrary_value_dp(file_id, dset_name, i_fp, i_fr) result(value)

    implicit none

    integer(hid_t), intent(in) :: file_id
    character(len=*), intent(in) :: dset_name
    integer, intent(in) :: i_fp
    integer, intent(in) :: i_fr

    double precision :: tmp_value(1)
    double precision :: value

    character(len=*), parameter :: fname = "read_one_arbitrary_value"
    integer(hid_t) :: dspace_id
    integer(hid_t) :: dset_id
    integer(hid_t) :: memspace_id
    integer(hsize_t) :: hs_offset(2), hs_count(2)
    integer(hsize_t) :: dim_mem(1)

    integer :: hdferr

    hs_offset(1) = i_fp - 1
    hs_offset(2) = i_fr - 1

    hs_count(1) = 1
    hs_count(2) = 1

    dim_mem(1) = 1

!$OMP CRITICAL
    call h5dopen_f(file_id, dset_name, dset_id, hdferr)
!$OMP END CRITICAL
    call check_hdf_error(hdferr, fname, "Error opening dataset at: " // trim(dset_name))
!$OMP CRITICAL
    call h5dget_space_f(dset_id, dspace_id, hdferr)
!$OMP END CRITICAL
    call check_hdf_error(hdferr, fname, "Error getting dataspace id for " // trim(dset_name))

!$OMP CRITICAL
    call h5sselect_hyperslab_f(dspace_id, H5S_SELECT_SET_F, hs_offset, hs_count, hdferr)
!$OMP END CRITICAL
    call check_hdf_error(hdferr, fname, "Error performing hyperslab selection.")

!$OMP CRITICAL
    call h5screate_simple_f(1, dim_mem, memspace_id, hdferr)
!$OMP END CRITICAL
    call check_hdf_error(hdferr, fname, "Error creating simple memory space.")
!$OMP CRITICAL
    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp_value, dim_mem, &
         hdferr, memspace_id, dspace_id)
!$OMP END CRITICAL
    call check_hdf_error(hdferr, fname, "Error reading data from " // trim(dset_name))
!$OMP CRITICAL
    call h5dclose_f(dset_id, hdferr)
!$OMP END CRITICAL
    call check_hdf_error(hdferr, fname, "Error closing dataset id for " // trim(dset_name))

    value = tmp_value(1)

  end function read_one_arbitrary_value_dp

end module file_utils_mod
