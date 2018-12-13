subroutine guanter_retrieval(l1b_file_id, my_instrument)

    use control, only: MCS
    use instruments, only: generic_instrument
    use oco2
    use logger_mod, only: logger => master_logger
    use file_utils, only: get_HDF5_dset_dims
    use HDF5

    implicit none

    integer(hid_t), intent(in) :: l1b_file_id
    class(generic_instrument), intent(in) :: my_instrument

    integer(hid_t) :: basisfunction_file_id

    character(len=*), parameter :: fname = "guanter_retrieval"

    double precision, allocatable :: tmp_data(:)
    integer(hsize_t), allocatable :: num_pixels(:)

    double precision, allocatable :: basisfunctions(:,:,:)
    double precision, allocatable :: dispersion_coefs(:,:,:)
    double precision, allocatable :: snr_coefs(:,:,:,:)

    integer :: i_fr, i_fp ! Indices for frames and footprints
    integer :: num_frames(1)

    character(len=999) :: dset_name
    integer(hid_t) :: dset_id
    integer :: hdferr


    integer :: i,j


    ! The strategy of this retrieval is:
    ! 1) Read necessary quantities from HDF files (basisfunctions, dispersion, noise)



    ! a) Read-in of the radiance basis functions
    call h5fopen_f(MCS%window%basisfunction_file%chars(), H5F_ACC_RDONLY_F, &
                   basisfunction_file_id, hdferr)
    if (hdferr /= 0) then
        call logger%fatal(fname, "Error opening HDF file: " // trim(MCS%window%basisfunction_file%chars()))
        stop 1
    end if

    ! And then start reading values into our own arrays
    ! Get the array dimensions by inquiring the first one..
    call get_HDF5_dset_dims(basisfunction_file_id, "/BasisFunction_SV1_FP1", num_pixels)

    ! Allocate the basisfunction array
    allocate(basisfunctions(num_pixels(1), my_instrument%num_fp, MCS%algorithm%n_basisfunctions))
    allocate(tmp_data(num_pixels(1)))

    ! And read them in, one by one - the number of footprints is conviently
    ! stored in my_instrument, so we don't need to invoke any type checking YET
    do i=1, my_instrument%num_fp
        do j=1, MCS%algorithm%n_basisfunctions

            ! For now, I've decided the naming pattern for the basisfunctions
            ! are BasisFunction_SVX_FPY, for basisfunction number X and footprint Y
            write(dset_name, '(A,G0.1,A,G0.1)') "/BasisFunction_SV", j, "_FP", i

            call logger%debug(fname, "Looking for " // trim(dset_name))
            call h5dopen_f(basisfunction_file_id, dset_name, dset_id, hdferr)
            if (hdferr /= 0) then
                call logger%fatal(fname, "Error. Could not open " // trim(dset_name))
                stop 1
            end if

            call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp_data, num_pixels, hdferr)
            if (hdferr /= 0) then
                call logger%fatal(fname, "Error. Could not read " // trim(dset_name))
                stop 1
            end if
            call logger%debug(fname, "Read in " // trim(dset_name))

            ! Copy over data to basisfunctions array. Note!! This might have a
            ! good number of NaNs in there. Does not mean the file is bad!
            basisfunctions(:,i,j) = tmp_data(:)

        end do
    end do

    ! b) Read-in of dispersion and noise coefficients (THIS IS INSTRUMENT SPECIFIC!)
    select type(my_instrument)
        type is (oco2_instrument)
            call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
            call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
            call my_instrument%read_num_frames(l1b_file_id, num_frames)
    end select


    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    do i_fp=1, my_instrument%num_fp
        do i_fr=1, num_frames(1)




        end do
    end do
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!


end subroutine


subroutine guanter_FM()

end subroutine
