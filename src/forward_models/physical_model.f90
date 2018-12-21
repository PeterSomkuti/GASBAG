module physical_model

    use control, only: MCS, MAX_WINDOWS
    use instruments, only: generic_instrument
    use oco2
    use logger_mod, only: logger => master_logger
    use file_utils, only: get_HDF5_dset_dims, check_hdf_error, write_DP_hdf_dataset, &
                          read_DP_hdf_dataset
    use, intrinsic:: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan

    use math_utils

    implicit none

    public physical_retrieval

    private

    ! We grab the full MET profiles from the MET data file, for quick access
    double precision, allocatable, dimension(:,:,:) :: met_T_profiles, &
                                                       met_P_levels, &
                                                       met_SH_profiles
    double precision, allocatable, dimension(:,:) :: met_psurf

contains

    subroutine physical_retrieval(my_instrument)

        implicit none
        class(generic_instrument), intent(in) :: my_instrument

        integer(hid_t) :: l1b_file_id, met_file_id, output_file_id, dset_id
        integer(hsize_t), dimension(:), allocatable :: dset_dims
        integer :: hdferr
        character(len=999) :: dset_name
        character(len=*), parameter :: fname = "physical_retrieval"

        integer :: i

        !! Open up the MET file
        call h5fopen_f(MCS%input%met_filename%chars(), &
                       H5F_ACC_RDONLY_F, MCS%input%met_file_id, hdferr)
        call check_hdf_error(hdferr, fname, "Error opening MET file: " &
                             // trim(MCS%input%met_filename%chars()))

        !! Store HDF file handler for more convenient access
        met_file_id = MCS%input%met_file_id
        l1b_file_id = MCS%input%l1b_file_id

        !! Get the necesary MET data profiles
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


    end subroutine


end module
