module physical_model_mod

    use control_mod, only: MCS, MAX_WINDOWS
    use instruments, only: generic_instrument
    use oco2_mod
    use logger_mod, only: logger => master_logger
    use file_utils_mod, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset, &
                          read_DP_hdf_dataset
    use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

    use solar_model_mod
    use math_utils_mod
    use statevector_mod
    use mod_datetime

    implicit none

    public physical_retrieval

    private

    ! Dispersion/wavelength array
    double precision, dimension(:,:,:), allocatable :: dispersion
    ! Sounding_ids
    integer(8), dimension(:,:), allocatable :: sounding_ids
    ! Sounding time strings
    character(len=25), dimension(:,:), allocatable :: sounding_time_strings
    ! Sounding scene geometry (solar and viewing zenith and azimuth)
    double precision, dimension(:,:), allocatable :: SZA, SAA, VZA, VAA
    ! Sounding scene location (lon, lat, altitude)
    double precision, dimension(:,:), allocatable :: lon, lat, altitude, &
                                                     relative_velocity, &
                                                     relative_solar_velocity

    ! We grab the full MET profiles from the MET data file, for quick access
    double precision, allocatable, dimension(:,:,:) :: met_T_profiles, &
                                                       met_P_levels, &
                                                       met_SH_profiles
    ! MET surface pressure
    double precision, allocatable, dimension(:,:) :: met_psurf

    ! The solar spectrum (wavelength, transmission)
<<<<<<< HEAD
    double precision, allocatable, dimension(:,:) :: solar_spectrum
    ! And the broadband spectrum
    double precision, dimension(:), allocatable :: solar_irrad
=======
    double precision, dimension(:,:), allocatable :: solar_spectrum
    double precision, dimension(:), allocatable :: solar_irrad

>>>>>>> ba4753e55eb1ec7c4c974954f0612d7cb6c7cbd6
    ! Radiances
    double precision, dimension(:,:,:), allocatable :: final_radiance, &
                                                       measured_radiance, &
                                                       noise_radiance

    type(statevector) :: SV


contains

    subroutine physical_retrieval(my_instrument)

        implicit none
        class(generic_instrument), intent(in) :: my_instrument

        integer(hid_t) :: l1b_file_id, met_file_id, output_file_id, dset_id
        integer(hsize_t), dimension(:), allocatable :: dset_dims
        integer :: hdferr
        character(len=999) :: dset_name, tmp_str
        character(len=*), parameter :: fname = "physical_retrieval"

        double precision, allocatable :: dispersion_coefs(:,:,:)
        double precision, allocatable :: snr_coefs(:,:,:,:)
        integer :: num_frames

        integer :: i_fp, i_fr, band, fp
        integer :: i, retr_count

        double precision :: cpu_time_start, cpu_time_stop, mean_duration

        !! Open up the MET file
        call h5fopen_f(MCS%input%met_filename%chars(), &
                       H5F_ACC_RDONLY_F, MCS%input%met_file_id, hdferr)
        call check_hdf_error(hdferr, fname, "Error opening MET file: " &
                             // trim(MCS%input%met_filename%chars()))

        !! Store HDF file handler for more convenient access
        met_file_id = MCS%input%met_file_id
        l1b_file_id = MCS%input%l1b_file_id

        !! Get the necesary MET data profiles. This probably needs expanding for
        !! other instrument types / file structures. (will Geocarb be different?)
        dset_name = "/Meteorology/vector_pressure_levels_met"
        call read_DP_hdf_dataset(met_file_id, dset_name, met_P_levels, dset_dims)
        call logger%trivia(fname, "Finished reading in pressure levels.")

        dset_name = "/Meteorology/temperature_profile_met"
        call read_DP_hdf_dataset(met_file_id, dset_name, met_T_profiles, dset_dims)
        call logger%trivia(fname, "Finished reading in temperature profiles.")

        dset_name = "/Meteorology/specific_humidity_profile_met"
        call read_DP_hdf_dataset(met_file_id, dset_name, met_SH_profiles, dset_dims)
        call logger%trivia(fname, "Finished reading in specific humidity profiles.")

        dset_name = "/Meteorology/surface_pressure_met"
        call read_DP_hdf_dataset(met_file_id, dset_name, met_psurf, dset_dims)
        call logger%trivia(fname, "Finished reading in surface pressure.")

        !do i=1, size(met_P_levels, 1)
        !    write(*,*) i, met_P_levels(i,1,1), met_T_profiles(i,1,1), met_SH_profiles(i,1,1)
        !end do
        !write(*,*) "Psurf: ", met_psurf(1,1)


        ! Read-in of dispersion and noise coefficients (THIS IS INSTRUMENT SPECIFIC!)
        select type(my_instrument)
            type is (oco2_instrument)
                ! Read dispersion coefficients and create dispersion array
                call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
                allocate(dispersion(1016,8,3))

                do band=1,3
                    do fp=1,8
                        call my_instrument%calculate_dispersion(dispersion_coefs, dispersion(:,fp,band), band, fp)
                    end do
                end do

                ! Grab the SNR coefficients for noise calculations
                call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
                ! How many frames do we have in this file again?
                call my_instrument%read_num_frames(l1b_file_id, num_frames)
                ! Read in the sounding id's
                call my_instrument%read_sounding_ids(l1b_file_id, sounding_ids)
                ! Read the time strings
                call my_instrument%read_time_strings(l1b_file_id, sounding_time_strings)
                ! Read in the measurement geometries
                call my_instrument%read_sounding_geometry(l1b_file_id, 1, SZA, SAA, VZA, VAA)
                ! Read in the measurement location
                call my_instrument%read_sounding_location(l1b_file_id, 1, lon, lat, &
                                                          altitude, relative_velocity, &
                                                          relative_solar_velocity)

        end select

        allocate(final_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
        allocate(measured_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
        allocate(noise_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))


        ! Read in the solar model
        if (MCS%algorithm%solar_type == "toon") then
            call read_toon_spectrum(MCS%algorithm%solar_file%chars(), &
                                    solar_spectrum)
<<<<<<< HEAD
            allocate(solar_irrad(size(solar_spectrum, 1)))
            call calculate_solar_planck_function(5800.d0, solar_spectrum(:,1), solar_irrad)
=======

            allocate(solar_irrad(size(solar_spectrum, 1)))
            call calculate_solar_planck_function(5800.d0, solar_spectrum(:,1), solar_irrad)

>>>>>>> ba4753e55eb1ec7c4c974954f0612d7cb6c7cbd6
        else
            call logger%fatal(fname, "Sorry, solar model type " &
                                     // MCS%algorithm%solar_type%chars() &
                                     // " is not known.")
        end if

        !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
        !!!!!! Set up state vector structure here

        ! We can do this once, and simply clear the contents of the state vectors
        ! inside the loop later on. Might not save an awful lot, but every little
        ! helps, I guess..

        ! At the moment, we can do following state vectors:
        !
        ! Lambertian surface albedo, arbitrary (?) polynomial order
        ! Zero-level offset / SIF (constant)
        ! Instrument disperison

        ! TODO: read in state vector stuff from config file
        call initialize_statevector(SV, 2, 1, 0)


        ! Main loop over all soundings to perform actual retrieval
        ! footprint, frame, microwindow #, band

        retr_count = 0
        mean_duration = 0.0d0
        do i_fp=1, my_instrument%num_fp
            do i_fr=1, num_frames

                call cpu_time(cpu_time_start)
                call physical_FM(my_instrument, i_fp, i_fr, 1, 1)
                call cpu_time(cpu_time_stop)

                retr_count = retr_count + 1
                mean_duration = ((mean_duration * retr_count) + (cpu_time_stop - cpu_time_start)) / (retr_count + 1)

                if (mod(retr_count, 100) == 0) then
                    write(tmp_str, '(A, F7.3)') "Duration: ", mean_duration
                    call logger%debug(fname, trim(tmp_str))
                    !read(*,*)
                end if
            end do
        end do


    end subroutine


    subroutine physical_FM(my_instrument, i_fp, i_fr, i_win, band)

        implicit none

        class(generic_instrument), intent(in) :: my_instrument
        integer, intent(in) :: i_fr, i_fp, i_win, band

        double precision, parameter :: SPEED_OF_LIGHT = 299792458d0

        !!
        integer(hid_t) :: l1b_file_id

        !!
        double precision, dimension(:), allocatable :: radiance_l1b, radiance_work
        integer :: l1b_wl_idx_min, l1b_wl_idx_max, solar_idx_min, solar_idx_max
        integer :: funit

        !! Sounding location /time stuff
        type(datetime) :: date ! Datetime object for sounding date/time
        double precision :: doy_dp ! Day of year as double precision

        !! Instrument stuff
        double precision :: instrument_doppler
        double precision, dimension(:), allocatable :: this_dispersion

        !! Solar stuff
        double precision :: solar_dist, solar_rv, earth_rv, solar_doppler
<<<<<<< HEAD
        double precision, dimension(:,:), allocatable :: this_solar

        !! Albedo
        double precision :: albedo_apriori


=======
>>>>>>> ba4753e55eb1ec7c4c974954f0612d7cb6c7cbd6

        !! Misc
        integer :: i

        l1b_file_id = MCS%input%l1b_file_id


        ! Create a datetime object for the current measurement
        select type(my_instrument)
            type is (oco2_instrument)
                call my_instrument%convert_time_string_to_date(sounding_time_strings(i_fp, i_fr), date)
        end select
        doy_dp = dble(date%yearday()) + dble(date%getHour()) / 24.0


        !write(*,*) "Isoformat: ", date%isoformat()
        !write(*,*) "Day of the year: ", doy_dp

        !write(*,*) solar_rv * 1000.0d0 + earth_rv, solar_doppler
        !write(*,*) "SZA: ", SZA(i_fp, i_fr)
        !write(*,*) "SAA: ", SAA(i_fp, i_fr)
        !write(*,*) "VZA: ", VZA(i_fp, i_fr)
        !write(*,*) "VAA: ", VAA(i_fp, i_fr)
        !write(*,*) "Lon: ", lon(i_fp, i_fr), " Lat: ", lat(i_fp, i_fr), "Altitude: ", altitude(i_fp, i_fr)

        ! Grab the radiance for this particular sounding
        select type(my_instrument)
            type is (oco2_instrument)
                call my_instrument%read_one_spectrum(l1b_file_id, i_fr, i_fp, 1, radiance_l1b)
        end select

        ! Here we grab the index limits for the radiance and solar spectra for
        ! the choice of our microwindow.
        l1b_wl_idx_min = 0
        l1b_wl_idx_max = 0
        solar_idx_min = 0
        solar_idx_max = 0

        do i=1, size(dispersion, 1)
            if (dispersion(i, i_fp, band) < MCS%window(i_win)%wl_min) then
                l1b_wl_idx_min = i
            end if
            if (dispersion(i, i_fp, band) < MCS%window(i_win)%wl_max) then
                l1b_wl_idx_max = i
            end if
            if (dispersion(i, i_fp, band) > MCS%window(i_win)%wl_max) then
                exit
            end if
        end do

        if (l1b_wl_idx_max < size(dispersion, 1)) then
            l1b_wl_idx_max = l1b_wl_idx_max + 1
        end if

        if (l1b_wl_idx_max > size(dispersion, 1)) then
            l1b_wl_idx_max = size(dispersion, 1)
        end if

        do i=1, size(solar_spectrum, 1)
            if (solar_spectrum(i,1) < MCS%window(i_win)%wl_min) then
                solar_idx_min = i
            end if
            if (solar_spectrum(i,1) < MCS%window(i_win)%wl_max) then
                solar_idx_max = i
            end if
            if (solar_spectrum(i,1) > MCS%window(i_win)%wl_max) then
                exit
            end if
        end do

<<<<<<< HEAD
        allocate(this_solar(solar_idx_max - solar_idx_min + 1, 2))
        allocate(radiance_work(l1b_wl_idx_max - l1b_wl_idx_min + 1))

        ! Grab a copy of the L1b radiances
        radiance_work(:) = radiance_l1b(l1b_wl_idx_min:l1b_wl_idx_max)

        ! If solar doppler-shift is needed, calculate the distance and relative
        ! velocity between point of measurement and the (fixed) sun

        !call calculate_solar_distance_and_rv(doy_dp, solar_dist, solar_rv)
        !call calculate_rel_velocity_earth_sun(lat(i_fp, i_fr), SZA(i_fp, i_fr), SAA(i_fp, i_fr), altitude(i_fp, i_fr), earth_rv)
        !solar_doppler =  (solar_rv * 1000.0d0 + earth_rv) / SPEED_OF_LIGHT

        solar_doppler = relative_solar_velocity(i_fp, i_fr) / SPEED_OF_LIGHT
        ! Take a copy of the solar spectrum and re-adjust the solar spectrum wavelength grid
        this_solar(:,1) = solar_spectrum(solar_idx_min:solar_idx_max, 1) / (1.0d0 + solar_doppler)
        this_solar(:,2) = solar_spectrum(solar_idx_min:solar_idx_max, 2) * solar_irrad(solar_idx_min:solar_idx_max)
=======

>>>>>>> ba4753e55eb1ec7c4c974954f0612d7cb6c7cbd6

        ! Correct for instrument doppler shift
        instrument_doppler = relative_velocity(i_fp, i_fr) / SPEED_OF_LIGHT
        ! Grab a copy of the dispersion relation that's shifted according to
        ! spacecraft velocity.
        allocate(this_dispersion(size(dispersion, 1)))
        this_dispersion = dispersion(:, i_fp, band) / (1.0d0 + instrument_doppler)


        select type(my_instrument)
            type is (oco2_instrument)
        !        ! OCO-2 has Stokes coefficient 0.5 for intensity, so we need to
        !        ! take that into account for the incoming solar irradiance
                albedo_apriori = PI * maxval(radiance_work(:)) / &
                (0.5 * maxval(this_solar(:,2)) * cos(DEG2RAD * SZA(i_fp, i_fr)))
        end select

        ! TODO: What do we do in case of invalid albedo (<0, >1)?

        ! Populate state vector priors (this should also be )
        SV%svap(SV%idx_albedo(1)) = albedo_apriori
        do i=2, SV%num_albedo
            SV%svap(SV%idx_albedo(i)) = 0.0d0
        end do
        SV%svap(SV%idx_sif(1)) = 0.0d0

        !write(*,*) i_fp, i_fr, "Albedo: ", albedo_apriori, "Altitude: ", altitude(i_fp, i_fr), "Latitude: ", lat(i_fp, i_fr)
        if (albedo_apriori > 1) then
            write(*,*) "albedo too large"
        !    read(*,*)
        end if




        !open(file="l1b_spec.dat", newunit=funit)
        !do i=l1b_wl_idx_min, l1b_wl_idx_max
        !    write(funit,*) dispersion(i, i_fp, band), radiance_l1b(i)
        !end do
        !close(funit)

        !open(file="solar_spec.dat", newunit=funit)
        !do i=1, size(this_solar, 1)
        !    write(funit,*) this_solar(i,1), this_solar(i,2)
        !end do
        !close(funit)


    end subroutine physical_FM


end module
