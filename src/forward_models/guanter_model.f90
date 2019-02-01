module guanter_model_mod

    use control_mod, only: MCS, MAX_WINDOWS
    use instruments, only: generic_instrument
    use oco2_mod
    use logger_mod, only: logger => master_logger
    use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset
    use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

    use math_utils_mod

    implicit none

    public guanter_retrieval

    private

    ! In the Guanter-scheme, dispersion is not really touched, hence we can just
    ! keep it as a module-wide fixed set of numbers (pixel, footprint)
    double precision, dimension(:,:,:), allocatable :: dispersion
    double precision, allocatable :: dispersion_coefs(:,:,:)
    double precision, allocatable :: snr_coefs(:,:,:,:)

    ! Basisfunctions are stored as (file, spectral pixel, footprint, #SV)
    double precision, dimension(:,:,:,:), allocatable :: basisfunctions
    ! Sounding_ids
    integer(8), dimension(:,:), allocatable :: sounding_ids

    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!! Retrieval results !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
    !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

    ! Final statevector and post-uncert for each measurement
    double precision, dimension(:,:,:), allocatable :: final_SV, final_SV_uncert
    ! Final SIF value in physical radiance units, chi2
    double precision, dimension(:,:), allocatable :: retrieved_SIF_abs, &
                                                     retrieved_SIF_rel, chi2

    double precision, dimension(:,:,:), allocatable :: final_radiance, &
                                                       measured_radiance, &
                                                       noise_radiance

contains

    subroutine guanter_retrieval(my_instrument)

        implicit none
        class(generic_instrument), intent(in) :: my_instrument


        integer(hid_t) :: l1b_file_id, basisfunction_file_id, output_file_id

        character(len=*), parameter :: fname = "guanter_retrieval"
        character(len=999) :: tmp_str

        double precision, allocatable :: tmp_data(:)
        integer(hsize_t), dimension(2) :: out_dims2d
        integer(hsize_t), dimension(3) :: out_dims3d
        integer(hsize_t), dimension(:), allocatable :: num_pixels


        integer :: i_fr, i_fp, i_win ! Indices for frames, footprints, microwindows
        integer :: num_frames, num_fp, num_total_soundings, cnt, band, fp

        character(len=999) :: dset_name
        integer(hid_t) :: dset_id, sif_result_gid
        integer :: hdferr

        integer :: i,j


        l1b_file_id = MCS%input%l1b_file_id
        output_file_id = MCS%output%output_file_id

        ! And also create the result group in the output file
        call h5gcreate_f(output_file_id, "linear_fluorescence_results", &
                         sif_result_gid, hdferr)
        call check_hdf_error(hdferr, fname, "Error. Could not create group: linear_fluorescence_results")

        ! The strategy of this retrieval is:
        ! 1) Read necessary quantities from HDF files (basisfunctions, dispersion, noise)
        ! a) Read-in of the radiance basis functions

        do i_win=1, MAX_WINDOWS

            if (MCS%window(i_win)%used .eqv. .false.) then
                cycle ! This window is not used!
            end if

            call h5fopen_f(MCS%window(i_win)%basisfunction_file%chars(), &
                           H5F_ACC_RDONLY_F, basisfunction_file_id, hdferr)
            call check_hdf_error(hdferr, fname, "Error opening HDF file: " &
                                 // trim(MCS%window(i_win)%basisfunction_file%chars()))

            if (.not. allocated(basisfunctions)) then
                ! And then start reading values into our own arrays
                ! Get the array dimensions by inquiring the first one..
                call get_HDF5_dset_dims(basisfunction_file_id, "/BasisFunction_SV1_FP1", num_pixels)
                ! Allocate the basisfunction array
                allocate(basisfunctions(MAX_WINDOWS, num_pixels(1), &
                                        my_instrument%num_fp, MCS%algorithm%n_basisfunctions))
                allocate(tmp_data(num_pixels(1)))
            end if

            ! And read them in, one by one - the number of footprints is conviently
            ! stored in my_instrument, so we don't need to invoke any type checking YET
            do i=1, my_instrument%num_fp
                do j=1, MCS%algorithm%n_basisfunctions

                    ! For now, I've decided the naming pattern for the basisfunctions
                    ! are BasisFunction_SVX_FPY, for basisfunction number X and footprint Y
                    write(dset_name, '(A,G0.1,A,G0.1)') "/BasisFunction_SV", j, "_FP", i

                    call logger%debug(fname, "Looking for " // trim(dset_name))
                    call h5dopen_f(basisfunction_file_id, dset_name, dset_id, hdferr)
                    call check_hdf_error(hdferr, fname, "Error. Could not open " // trim(dset_name))

                    call h5dread_f(dset_id, H5T_NATIVE_DOUBLE, tmp_data, num_pixels, hdferr)
                    call check_hdf_error(hdferr, fname, "Error. Could not read " // trim(dset_name))

                    call logger%debug(fname, "Read in " // trim(dset_name))

                    ! Copy over data to basisfunctions array. Note!! This might have a
                    ! good number of NaNs in there. Does not mean the file is bad!
                    basisfunctions(i_win, :, i, j) = tmp_data(:)

                end do
            end do

            call h5fclose_f(basisfunction_file_id, hdferr)
            call check_hdf_error(hdferr, fname, "Could not close basisfunction file " &
                                 // trim(MCS%window(i_win)%basisfunction_file%chars()))

        end do

        ! b) Read-in of dispersion and noise coefficients (THIS IS INSTRUMENT SPECIFIC!)
        select type(my_instrument)
            type is (oco2_instrument)
                ! Read dispersion coefficients and create dispersion array
                call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)

                ! Grab the SNR coefficients for noise calculations
                call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
                ! How many frames do we have in this file again?
                call my_instrument%read_num_frames_and_fp(l1b_file_id, num_frames, num_fp)
                ! Read in the sounding id's
                call my_instrument%read_sounding_ids(l1b_file_id, sounding_ids)

                allocate(dispersion(1016,8,3))
                do band=1,3
                    do fp=1, num_fp
                        call my_instrument%calculate_dispersion(dispersion_coefs(:,fp,band), dispersion(:,fp,band), band, fp)
                    end do
                end do



        end select

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Allocate the (module-wide) result arrays here
        num_total_soundings = my_instrument%num_fp * num_frames

        allocate(retrieved_SIF_abs(my_instrument%num_fp, num_frames))
        allocate(retrieved_SIF_rel(my_instrument%num_fp, num_frames))
        allocate(chi2(my_instrument%num_fp, num_frames))
        allocate(final_SV(my_instrument%num_fp, num_frames, MCS%algorithm%n_basisfunctions + 1))
        allocate(final_SV_uncert(my_instrument%num_fp, num_frames, MCS%algorithm%n_basisfunctions + 1))

        allocate(final_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
        allocate(measured_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
        allocate(noise_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! Loop through all frames and footprints, and perform the retrieval
        !! The outermose loop simply does it for the various microwindows

        do i_win=1, MAX_WINDOWS

            if (MCS%window(i_win)%used .eqv. .false.) then
                cycle ! This window is not used!
            end if

            call logger%info(fname, "Processing window: " // trim(MCS%window(i_win)%name))

            retrieved_SIF_abs(:,:) = -9999.99d0
            retrieved_SIF_rel(:,:) = -9999.99d0
            final_SV(:,:,:) = -9999.99d0
            final_SV_uncert(:,:,:) = -9999.99d0
            chi2(:,:) = -9999.99d0

            final_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
            measured_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)
            noise_radiance = IEEE_VALUE(1D0, IEEE_QUIET_NAN)

            cnt = 0
            do i_fr=1, num_frames
                do i_fp=1, my_instrument%num_fp

                    write(tmp_str,'(A,I7.1,A,I1.1)') "Frame: ", i_fr, ", FP: ", i_fp
                    call logger%trivia(fname, trim(tmp_str))

                    ! Retrieval time!
                    call guanter_FM(my_instrument, i_fr, i_fp, i_win)

                    ! Print out progress..
                    if (MODULO(cnt, floor(0.05 * num_total_soundings)) == 0) then
                        write(tmp_str, "(A, F6.2, A)") "Progress: ", (1.0 * cnt) / num_total_soundings * 100.0,  "%"
                        call logger%info(fname, trim(tmp_str))
                    end if
                    cnt = cnt + 1
                end do
            end do
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

            !! Write the results into the output HDF5 file

            call logger%info(fname, "Writing out results.")
            ! And then just write out the datasets / arrays
            out_dims3d = shape(final_SV)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/final_statevector_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      final_SV, out_dims3d, -9999.99d0)

            out_dims3d = shape(final_SV_uncert)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/final_statevector_uncert_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      final_SV_uncert, out_dims3d, -9999.99d0)

            out_dims2d = shape(retrieved_SIF_rel)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/retrieved_sif_rel_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      retrieved_SIF_rel, out_dims2d, -9999.99d0)

            out_dims2d = shape(retrieved_SIF_abs)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/retrieved_sif_abs_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      retrieved_SIF_abs, out_dims2d, -9999.99d0)

            out_dims2d = shape(chi2)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/retrieved_chi2_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      chi2, out_dims2d, -9999.99d0)

            out_dims3d = shape(final_radiance)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/modelled_radiance_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      final_radiance, out_dims3d)

            out_dims3d = shape(measured_radiance)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/measured_radiance_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      measured_radiance, out_dims3d)

            out_dims3d = shape(noise_radiance)
            write(tmp_str, '(A,A)') "/linear_fluorescence_results/noise_radiance_" // MCS%window(i_win)%name
            call write_DP_hdf_dataset(output_file_id, &
                                      trim(tmp_str), &
                                      noise_radiance, out_dims3d)

            call logger%info(fname, "Finished writing out results.")
        end do ! microwindow loop
    end subroutine


    subroutine guanter_FM(my_instrument, i_fr, i_fp, i_win)

        implicit none

        class(generic_instrument), intent(in) :: my_instrument
        integer, intent(in) :: i_fr, i_fp, i_win

        integer :: l1b_wl_idx_min, l1b_wl_idx_max
        integer :: i, j, i_all
        integer(hid_t) :: l1b_file_id

        character(len=*), parameter :: fname = "guanter_FM"
        character(len=999) :: tmp_str

        !! This is where most of the juice is. The retrieval concept is fairly
        !! simple, but there are still some tasks to do to get it done. i_fr and
        !! i_fp tell us the location of the sounding we will process here.

        ! The measured radiance coming straight from the L1B file, and then to
        ! be further processed (slope-corrected)
        double precision, dimension(:), allocatable :: radiance_work, noise_work, rad_conv, residual
        double precision, dimension(:), allocatable :: radiance_l1b

        double precision :: cont_rad, fit_intercept
        ! Are the radiances that we process sensible?
        logical :: radiance_OK

        logical :: success_inv_mat

        double precision :: tmp_chi2, BIC
        logical :: BIC_converged

        !! Here are all the retrieval scheme matrices, quantities, etc., and
        !! temporary matrices

        integer :: N_sv ! Number of statevector elements
        integer :: N_spec

        double precision, dimension(:), allocatable :: xhat
        double precision, dimension(:,:), allocatable :: K, Se_inv, Shat, Shat_inv
        double precision, dimension(:,:), allocatable :: m_tmp1, m_tmp2


        l1b_file_id = MCS%input%l1b_file_id

        ! Turn frame/FP index into a flat index for saving results
        i_all = (i_fr - 1) * 8 + i_fp

        select type(my_instrument)
            type is (oco2_instrument)
                call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, 1, radiance_l1b)
        end select

        cont_rad = percentile(radiance_l1b, 99.9d0)

        !! We now have the full L1b spectrum of a given index, we need to grab
        !! the relevant portion of the spectrum.
        l1b_wl_idx_min = 0
        l1b_wl_idx_max = 0

        do i=1, size(dispersion, 1)
            if (dispersion(i, i_fp, 1) < MCS%window(i_win)%wl_min) then
                l1b_wl_idx_min = i
            end if
            if (dispersion(i, i_fp, 1) < MCS%window(i_win)%wl_max) then
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
        N_spec = size(radiance_work)
        allocate(residual, mold=radiance_work)


        ! Here we check the radiance
        select type(my_instrument)
            type is (oco2_instrument)
                call my_instrument%check_radiance_valid(l1b_file_id, radiance_work, &
                                                        l1b_wl_idx_min, l1b_wl_idx_max, &
                                                        radiance_OK)
        end select

        if (radiance_OK .eqv. .false.) then
            write(tmp_str, '(A,G0.1,A,G0.1,A)') "Sounding (", i_fr, ",", i_fp, &
                                                ") has invalid radiances. Skipping."
            call logger%warning(fname, trim(tmp_str))
            return
        end if


        ! Now calculate the noise-equivalent radiances
        select type(my_instrument)
            type is (oco2_instrument)
                call my_instrument%calculate_noise(snr_coefs, radiance_work, &
                                                   noise_work, i_fp, 1, &
                                                   l1b_wl_idx_min, l1b_wl_idx_max)
        end select

        ! Do the slope correction for the selected spectrum
        call slope_correction(radiance_work, 40.0d0, fit_intercept)
        ! And we have our "y" mesurement vector ready!!

        ! NOTE: HERE, RADIANCE_L1B IS IS IN PHYSICAL UNITS, AS IS NOISE_WORK.
        ! RADIANCE_WORK, HOWEVER, HAS BEEN SLOPE-NORMALIZED AND IS OF UNIT 1.
        ! WE THUS NORMALIZE THE NOISE DATA AS WELL (see lines below)

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !! And here is the whole retrieval magic. Simple stuff, linear retrieval
        !! see Rodgers (2000) for the usual details.

        ! Initial statevector number: amount of basisfunctions plus 1 for SIF
        N_sv = MCS%algorithm%n_basisfunctions + 1

        !! Calculate the inverse(!) noise matrix Se_inv
        allocate(Se_inv(N_spec, N_spec))
        ! Since we work in slope-normalized space, we also need to scale our
        ! noise levels.
        noise_work(:) = noise_work(:) / fit_intercept

        Se_inv(:,:) = 0.0d0
        do i=1, size(radiance_work)
            Se_inv(i,i) = 1 / (noise_work(i) ** 2)
        end do

        allocate(xhat(N_sv))
        allocate(Shat(N_sv, N_sv))
        allocate(Shat_inv(N_sv, N_sv))

        !! Stick basisfunctions into Jacobian matrix, which has following
        !! structure: first M entries are the M basisfunctions, and then there
        !! is the fluorescence as the last entry
        allocate(K(N_spec, N_sv))
        do i=1, MCS%algorithm%n_basisfunctions
            K(:,i) = basisfunctions(i_win, l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i)
        end do
        K(:, N_sv) = 1.0d0 ! fluorescence jacobian is flat!

        !! Easy check if we did everything right - are there NaNs in K?
        do i=1, size(K, 1)
            do j=1, size(K, 2)
                if (ieee_is_nan(K(i,j))) then
                    call logger%fatal(fname, "There are NaN's in the Jacobian.")
                    stop 1
                end if
            end do
        end do


        !! Start by calculating posterior covariance for linear case
        !! Shat = (K^T S_e^-1 K)^-1 (Rodgers 2.27 when S_a is large)
        allocate(m_tmp1, mold=Shat)
        m_tmp1 = matmul(matmul(transpose(K), Se_inv), K) ! this is K^T Se_inv K
        Shat_inv(:,:) = m_tmp1

        allocate(m_tmp2, mold=m_tmp1)
        call invert_matrix(m_tmp1, m_tmp2, success_inv_mat)
        if (.not. success_inv_mat) then
           call logger%error(fname, "Failed to invert K^T Se^-1 K")
           return
        end if

        ! m_tmp2 is know (K^T Se_inv K)^-1 = Shat
        Shat(:,:) = m_tmp2(:,:)

        ! Report the state vector posteriori uncertainty
        do i=1, N_sv
            final_SV_uncert(i_fp, i_fr, i) = sqrt(Shat(i,i))
        end do

        !! Calculate xhat!
        xhat = matmul(matmul(matmul(Shat, transpose(K)), Se_inv), radiance_work)
        final_SV(i_fp, i_fr, :) = xhat(:)

        allocate(rad_conv, mold=radiance_work)
        rad_conv(:) = 0.0d0
        !! Calculate the modeled radiance via xhat
        do i=1, N_sv
            !write(*,*) i, xhat(i), sqrt(Shat(i,i))
            rad_conv(:) = rad_conv(:) + K(:,i) * xhat(i)
        end do


        noise_work(:) = noise_work(:) * fit_intercept
        retrieved_SIF_abs(i_fp, i_fr) = xhat(N_sv) * fit_intercept !cont_rad !radiance_l1b(l1b_wl_idx_min)
        write(tmp_str, '(A, E12.5)') "SIF abs: ", xhat(N_sv) * fit_intercept !cont_rad ! radiance_l1b(l1b_wl_idx_min)
        call logger%trivia(fname, trim(tmp_str))

        retrieved_SIF_rel(i_fp, i_fr) = xhat(N_sv)
        write(tmp_str, '(A, F12.3,A)') "SIF rel: ", xhat(N_sv) * 100, "%"
        call logger%trivia(fname, trim(tmp_str))

        !! Now we still are in a unitless world, so let's get back to the
        !! realm of physical units by back-scaling our result here with the first
        !! pixel of the microwindow.
        radiance_work(:) = radiance_work(:) * fit_intercept !cont_rad !radiance_l1b(l1b_wl_idx_min)
        rad_conv(:) = rad_conv(:) * fit_intercept !cont_rad !radiance_l1b(l1b_wl_idx_min)
        ! Get the spectral residual
        residual(:) = rad_conv(:) - radiance_work(:)

        final_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = rad_conv(:)
        measured_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = radiance_work(:)
        noise_radiance(l1b_wl_idx_min:l1b_wl_idx_max, i_fp, i_fr) = noise_work(:)

        ! reduced chi2
        tmp_chi2 = SUM((residual ** 2) / (noise_work ** 2))
        tmp_chi2 = tmp_chi2 / (N_spec - N_sv)
        chi2(i_fp, i_fr) = tmp_chi2

        BIC = N_spec * log(SUM(residual ** 2) / N_spec) + N_sv * log(1.0 * N_spec)

        write(tmp_str, '(A, I3.1)') "N_sv: ", N_sv
        call logger%trivia(fname, trim(tmp_str))
        write(tmp_str, '(A,F7.3)') "Chi2: ", tmp_chi2
        call logger%trivia(fname, trim(tmp_str))

    end subroutine






    subroutine slope_correction(radiance, perc, intercept)

        implicit none
        double precision, intent(inout) :: radiance(:)
        double precision, intent(in) :: perc
        double precision, intent(inout) :: intercept

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
                adata(cnt, 2)= i-1
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
            radiance(i) = radiance(i) / (bdata(1,1) + ((i-1) * bdata(2,1)))
        end do

        intercept = bdata(1,1)

    end subroutine



end module
