module guanter_model

    use control, only: MCS
    use instruments, only: generic_instrument
    use oco2
    use logger_mod, only: logger => master_logger
    use file_utils, only: get_HDF5_dset_dims

    use math_utils

    implicit none

    public guanter_retrieval

    private

    ! In the Guanter-scheme, dispersion is not really touched, hence we can just
    ! keep it as a module-wide fixed set of numbers (pixel, footprint)
    double precision, save, allocatable :: dispersion(:,:,:)
    ! Basisfunctions are stored as (spectral pixel, footprint, #SV)
    double precision, save, allocatable :: basisfunctions(:,:,:)


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
                call my_instrument%calculate_dispersion(dispersion_coefs, dispersion)

                ! Grab the SNR coefficients for noise calculations
                call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
                ! How many frames do we have in this file again?
                call my_instrument%read_num_frames(l1b_file_id, num_frames)

        end select

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        do i_fp=3, 3 !my_instrument%num_fp
            do i_fr=2345, 2345 !num_frames(1)

                write(*,*) i_fp, i_fr

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
        integer :: l1b_wl_idx_min, l1b_wl_idx_max
        double precision, dimension(:,:,:,:) :: snr_coefs
        integer :: i, j, funit

        !! This is where most of the juice is. The retrieval concept is fairly
        !! simple, but there are still some tasks to do to get it done. i_fr and
        !! i_fp tell us the location of the sounding we will process here.

        ! The measured radiance coming straight from the L1B file, and then to
        ! be further processed (slope-corrected)
        double precision, dimension(:), allocatable :: radiance_work, noise_work, tmp_rad
        double precision, dimension(:), allocatable :: radiance_l1b


        !! Here are all the retrieval scheme matrices, quantities, etc., and
        !! temporary matrices

        integer :: N_sv ! Number of statevector elements

        double precision, dimension(:), allocatable :: xhat
        double precision, dimension(:,:), allocatable :: K, Se_inv, Shat, Shat_inv
        double precision, dimension(:,:), allocatable :: m_tmp1, m_tmp2, m_tmp3


        select type(my_instrument)
            type is (oco2_instrument)
                call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, 1, radiance_l1b)
        end select

        !! We now have the full L1b spectrum of a given index, we need to grab
        !! the relevant portion of the spectrum.
        l1b_wl_idx_min = 0
        l1b_wl_idx_max = 0

        do i=1, size(dispersion, 1)
            if (dispersion(i, i_fp, 1) < MCS%window%wl_min) then
                l1b_wl_idx_min = i
            end if
            if (dispersion(i, i_fp, 1) < MCS%window%wl_max) then
                l1b_wl_idx_max = i
            end if
        end do

        if (l1b_wl_idx_max < size(dispersion, 1)) then
            l1b_wl_idx_max = l1b_wl_idx_max + 1
        end if

        if (l1b_wl_idx_max > size(dispersion, 1)) then
            l1b_wl_idx_max = size(dispersion, 1)
        end if

        !write(*,*) "Dispersion min value: ", dispersion(l1b_wl_idx_min, i_fp, 1), l1b_wl_idx_min
        !write(*,*) "Dispersion max value: ", dispersion(l1b_wl_idx_max, i_fp, 1), l1b_wl_idx_max

        ! Allocate some space for the new work array which we can do the
        ! retrieval with, and copy over the corrsponding part of the spectrum
        allocate(radiance_work(l1b_wl_idx_max - l1b_wl_idx_min + 1))
        allocate(noise_work(l1b_wl_idx_max - l1b_wl_idx_min + 1))

        ! Copy the relevant part of the spectrum to radiance_work
        radiance_work(:) = radiance_l1b(l1b_wl_idx_min:l1b_wl_idx_max)

        ! New calculate the noise-equivalent radiances
        select type(my_instrument)
            type is (oco2_instrument)
                call my_instrument%calculate_noise(snr_coefs, radiance_work, &
                noise_work, i_fp, 1, l1b_wl_idx_min, l1b_wl_idx_max)
        end select

        ! Do the slope correction for the selected spectrum
        call slope_correction(radiance_work, 40.0d0)
        ! And we have our "y" mesurement vector ready!!

        ! NOTE: HERE, RADIANCE_L1B IS IS IN PHYSICAL UNITS, AS IS NOISE_WORK.
        ! RADIANCE_WORK, HOWEVER, HAS BEEN SLOPE-NORMALIZED AND IS OF UNIT 1.




        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! And here is the whole retrieval magic. Simple stuff, linear retrieval
        !! see Rodgers (2000) for the usual details.

        N_sv = MCS%algorithm%n_basisfunctions + 1
        allocate(xhat(N_sv))
        allocate(Shat(N_sv, N_sv))
        allocate(Shat_inv(N_sv, N_sv))

        !! Stick basisfunctions into Jacobian matrix, which has following
        !! structure: first M entries are the M basisfunctions, and then there
        !! is the fluorescence, so we have M+1

        allocate(K(size(radiance_work), N_sv))
        do i=1, MCS%algorithm%n_basisfunctions
            K(:,i) = basisfunctions(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i)
        end do
        K(:, N_sv) = 1.0d0 ! fluorescence jacobian is flat!

        !! Calculate the inverse(!) noise matrix Se_inv
        allocate(Se_inv(size(radiance_work), size(radiance_work)))

        !noise_work(:) = noise_work(:) / radiance_l1b(l1b_wl_idx_min)
        Se_inv(:,:) = 0.0d0
        do i=1, size(radiance_work)
            Se_inv(i,i) = 1 / (noise_work(i) ** 2)
        end do



        !! Start by calculating posterior covariance for linear case
        !! Shat = (K^T S_e^-1 K)^-1 (Rodgers 2.27 when S_a is large)
        allocate(m_tmp1, mold=Shat)
        m_tmp1 = matmul(matmul(transpose(K), Se_inv), K) ! this is K^T Se_inv K
        Shat_inv(:,:) = m_tmp1

        allocate(m_tmp2, mold=m_tmp1)
        call invert_matrix(m_tmp1, m_tmp2)
        ! m_tmp2 is know (K^T Se_inv K)^-1
        !allocate(Shat, mold=m_tmp2)
        !allocate(Shat_inv, mold=m_tmp2)

        Shat(:,:) = m_tmp2(:,:)
        !write(*,*) shape(Shat), shape(K), shape(transpose(K))
        !write(*,*) matmul(Shat, transpose(K))

        xhat = matmul(matmul(matmul(Shat, transpose(K)), Se_inv), radiance_work)

        allocate(tmp_rad, mold=radiance_work)

        do i=1, N_sv
            write(*,*) i, xhat(i)
            tmp_rad(:) = tmp_rad(:) + K(:,i) * xhat(i)
        end do

        open(newunit=funit, file='K.dat')
        do i=1, size(K, 1)
            write(funit, *) (K(i,j), j=1,size(K,2))
        end do


        open(newunit=funit, file="test.dat")
        do i=1, size(radiance_work)
            write(funit, *) radiance_work(i), &
                            tmp_rad(i), &
                            radiance_l1b(l1b_wl_idx_min+i-1), &
                            noise_work(i), &
                            basisfunctions(i-1+l1b_wl_idx_min,i_fp,1), &
                            basisfunctions(i-1+l1b_wl_idx_min,i_fp,2), &
                            basisfunctions(i-1+l1b_wl_idx_min,i_fp,3)
        enddo
        close(funit)

    end subroutine






    subroutine slope_correction(radiance, perc)

        implicit none
        double precision, intent(inout) :: radiance(:)
        double precision, intent(in) :: perc

        double precision :: perc_value
        integer :: num_used
        double precision, allocatable :: adata(:,:), bdata(:,:)
        double precision, allocatable :: work(:)
        integer :: sgels_info
        character(len=*), parameter :: fname = "slope_correction"

        integer :: i, cnt
        ! Slope normalization works like this: first, we grab the top X percent
        ! all radiance points (to make sure we mostly have continuum). X ~ 60
        ! After that, we fit a linear function through these points to get
        ! slope and intersect - by which linear function we then divide the
        ! radiance.

        ! This is the percentile value from the radiance array
        perc_value = percentile(radiance, perc)
        num_used = count(radiance > perc_value)

        ! We need a new container for values ABOVE the percentile one
        allocate(adata(num_used, 2))
        allocate(bdata(num_used, 1))
        allocate(work(2 * num_used))

        ! Take out all the radiance values (with pixel positions) that are
        ! larger than the percentile value, and stick them into tmp_rad and tmp_coord
        adata(:, 1) = 1.0 ! For calculating intercepts
        cnt = 1
        do i=1, size(radiance)
            if (radiance(i) > perc_value) then
                bdata(cnt, 1) = radiance(i)
                adata(cnt, 2)= i
                cnt = cnt + 1
            end if
        end do

        !! Now we need to perform a linear regression y = ax + b.
        !! Use the handy LAPACK routine DGELS here.

        call DGELS('N', &          ! TRANS
                   num_used, &     ! M
                   2, &            ! N
                   1, &            ! NRHS
                   adata, &        ! A
                   num_used, &     ! LDA
                   bdata, &        ! B
                   num_used, &     ! LDB
                   work, &         ! WORK
                   2*num_used, &   ! LWORK
                   sgels_info)

        if (sgels_info /= 0) then
            call logger%fatal(fname, "Error from DGELS.")
            stop 1
        end if

        ! We now have the intercept in bdata(1,1), and the slope in bdata(2,1)
        ! and can thus slope-correct the spectrum.
        do i=1, size(radiance)
            radiance(i) = radiance(i) / (bdata(1,1) + (i * bdata(2,1)))
        end do


    end subroutine



end module
