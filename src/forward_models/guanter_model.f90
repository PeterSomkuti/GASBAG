module guanter_model

    use control, only: MCS
    use instruments, only: generic_instrument
    use oco2
    use logger_mod, only: logger => master_logger
    use file_utils, only: get_HDF5_dset_dims

    implicit none

    ! In the Guanter-scheme, dispersion is not really touched, hence we can just
    ! keep it as a module-wide fixed set of numbers (pixel, footprint)
    double precision, allocatable :: dispersion(:,:)

    double precision, allocatable :: basisfunctions(:,:,:)

    public guanter_retrieval


contains

    subroutine guanter_retrieval(l1b_file_id, my_instrument)

        implicit none

        integer(hid_t), intent(in) :: l1b_file_id
        class(generic_instrument), intent(in) :: my_instrument

        integer(hid_t) :: basisfunction_file_id

        character(len=*), parameter :: fname = "guanter_retrieval"

        double precision, allocatable :: tmp_data(:)
        integer(hsize_t), allocatable :: num_pixels(:)


        double precision, allocatable :: dispersion_coefs(:,:,:)
        double precision, allocatable :: dispersion(:,:,:)
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
                ! Read dispersion coefficients and create dispersion array
                call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
                call my_instrument%calculate_dispersion(dispersion_coefs, dispersion, num_pixels(1))

                ! Grab the SNR coefficients for noise calculations
                call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
                ! How many frames do we have in this file again?
                call my_instrument%read_num_frames(l1b_file_id, num_frames)

        end select

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i_fp=1, 1!my_instrument%num_fp
            do i_fr=1, 1!num_frames(1)

                call guanter_FM(my_instrument, l1b_file_id, i_fr, i_fp, &
                                snr_coefs)

            end do
        end do
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    end subroutine


    subroutine guanter_FM(my_instrument, l1b_file_id, i_fr, i_fp, &
                          snr_coefs)

        implicit none
        class(generic_instrument), intent(in) :: my_instrument
        integer(hid_t) :: l1b_file_id
        integer :: i_fr, i_fp
        double precision :: snr_coefs(:,:,:,:)

        !! This is where most of the juice is. The retrieval concept is fairly
        !! simple, but there are still some tasks to do to get it done. i_fr and
        !! i_fp tell us the location of the sounding we will process here.

        ! The measured radiance coming straight from the L1B file, and then to
        ! be further processed (slope-corrected)
        double precision, allocatable :: radiance_work(:)
        double precision, allocatable :: radiance_l1b(:)


        select type(my_instrument)
            type is (oco2_instrument)

                call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, 1, radiance_l1b)

        end select





    end subroutine


    subroutine slope_correction(radiance, perc)

        implicit none
        double precision, intent(in) :: radiance(:)
        integer, intent(in) :: perc

        double precision, allocatable :: tmp_rad(:)
        ! Slope normalization works like this: first, we grab the top X percent
        ! all radiance points (to make sure we mostly have continuum). X ~ 60
        ! After that, we fit a linear function through these points to get
        ! slope and intersect - by which linear function we then divide the
        ! radiance.



    end subroutine



end module
