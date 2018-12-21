module physical_model

    use control, only: MCS, MAX_WINDOWS
    use instruments, only: generic_instrument
    use oco2
    use logger_mod, only: logger => master_logger
    use file_utils, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset, &
                          read_DP_hdf_dataset
    use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

    use solar_model
    use math_utils

    implicit none

    public physical_retrieval

    private

    ! Dispersion/wavelength array
    double precision, dimension(:,:,:), allocatable :: dispersion
    ! Sounding_ids
    integer(8), dimension(:,:), allocatable :: sounding_ids
    ! Radiances
    double precision, dimension(:,:,:), allocatable :: final_radiance, &
                                                       measured_radiance, &
                                                       noise_radiance

    ! We grab the full MET profiles from the MET data file, for quick access
    double precision, allocatable, dimension(:,:,:) :: met_T_profiles, &
                                                       met_P_levels, &
                                                       met_SH_profiles
    double precision, allocatable, dimension(:,:) :: met_psurf

    ! The solar spectrum
    double precision, allocatable, dimension(:,:) :: solar_spectrum

contains

    subroutine physical_retrieval(my_instrument)

        implicit none
        class(generic_instrument), intent(in) :: my_instrument

        integer(hid_t) :: l1b_file_id, met_file_id, output_file_id, dset_id
        integer(hsize_t), dimension(:), allocatable :: dset_dims
        integer :: hdferr
        character(len=999) :: dset_name
        character(len=*), parameter :: fname = "physical_retrieval"

        double precision, allocatable :: dispersion_coefs(:,:,:)
        double precision, allocatable :: snr_coefs(:,:,:,:)
        integer :: num_frames

        integer :: i

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

        do i=1, size(met_P_levels, 1)
            write(*,*) i, met_P_levels(i,1,1), met_T_profiles(i,1,1), met_SH_profiles(i,1,1)
        end do
        write(*,*) "Psurf: ", met_psurf(1,1)


        ! Read-in of dispersion and noise coefficients (THIS IS INSTRUMENT SPECIFIC!)
        select type(my_instrument)
            type is (oco2_instrument)
                ! Read dispersion coefficients and create dispersion array
                call my_instrument%read_l1b_dispersion(l1b_file_id, dispersion_coefs)
                call my_instrument%calculate_dispersion(dispersion_coefs, dispersion)

                ! Grab the SNR coefficients for noise calculations
                call my_instrument%read_l1b_snr_coef(l1b_file_id, snr_coefs)
                ! How many frames do we have in this file again?
                call my_instrument%read_num_frames(l1b_file_id, num_frames)
                ! Read in the sounding id's
                call my_instrument%read_sounding_ids(l1b_file_id, sounding_ids)

        end select

        allocate(final_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
        allocate(measured_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))
        allocate(noise_radiance(size(dispersion, 1), my_instrument%num_fp, num_frames))


        ! Read in the solar model
        if (MCS%algorithm%solar_type == "toon") then
            call read_toon_spectrum(MCS%algorithm%solar_file%chars(), &
                                    solar_spectrum)
        else
            call logger%fatal(fname, "Sorry, solar model type " &
                                     // MCS%algorithm%solar_type%chars() &
                                     // " is not known.")
        end if



    end subroutine


end module
