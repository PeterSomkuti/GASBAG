!> @brief Main GeoCARB SIF Retrieval Program
!> @file GeoCARBSIF_main.f90
!> @author Peter Somkuti
!>
!! This is where it all starts. This main program calls functions to read in
!! the command-line input via FLAP as well as setting up the configuration through
!! through the text file via FINER.

program GeoCARBSIF

  use logger_mod, only: logger => master_logger
  use startup_mod, only: initialize_config
  use version_mod, only: git_branch, git_commit_hash, git_rev_no
  use control_mod, only: MCS, populate_MCS
  use finer, only: file_ini
  use file_utils_mod, only: check_hdf_error
  use instruments
  use oco2_mod


  use iso_fortran_env
  use HDF5

  implicit none


  type(file_ini) :: fini ! The config file structure
  class(generic_instrument), allocatable :: my_instrument ! The used instrument type
  integer :: hdferr ! HDF error variable


  ! Initilize the HDF5 library program-wide
  call h5open_f(hdferr)
  if (hdferr /= 0) then
     write(*, '(A)') "Error initializing HDF5 library."
     stop 1
  end if

  ! Greet the user and display information about the build itself.

  write(*,'(A)') "------------------------------"
  write(*,'(A)') "╔═╗┌─┐┌─┐╔═╗╔═╗╦═╗╔╗   ╔═╗╦╔═╗"
  write(*,'(A)') "║ ╦├┤ │ │║  ╠═╣╠╦╝╠╩╗  ╚═╗║╠╣ "
  write(*,'(A)') "╚═╝└─┘└─┘╚═╝╩ ╩╩╚═╚═╝  ╚═╝╩╚  "
  write(*,'(A)') "------------------------------"

  write(*,'(A)') "Version [" // git_branch // " " // git_commit_hash // &
       " #" // git_rev_no // "]"

  ! Initialize the whole thing by reading the configuration file
  call initialize_config(fini)

  ! Initialize the program control_mod structure (MCS) with the settings
  ! from the config file. MCS is designed to be read program-wide, so you
  ! can simply read any user-settings from any module.
  call populate_MCS(fini)

  ! This is where the my_insturment type is properly allocated using one of
  ! the derived types, depending on the instrument specified in the config.
  ! From here on, all derived-type bound subroutines MUST be called through
  ! a (sadly maybe cumbersome) SELECT TYPE statement - but this is just how
  ! Fortran works with run-time polymorphism.

  if (MCS%input%instrument_name == 'oco2') then
     allocate(oco2_instrument :: my_instrument)
     call logger%info("Main", "Using instrument: OCO-2")
  else
     call logger%fatal("Main", "Unknown instrument " // MCS%input%instrument_name)
     stop 1
  end if

  select type(my_instrument)
  type is (oco2_instrument)
     ! Scan the L1b file - we need some info from there, mostly the
     ! number of frames, footprints, bands and spectral points
     call my_instrument%scan_l1b_file(MCS%input%l1b_filename)
  end select


  ! Open up the output HDF5 file for writing and save the HDF5 file handler in
  ! MCS to be used all over the program.
  call h5fcreate_f(MCS%output%output_filename%chars(), H5F_ACC_TRUNC_F, &
       MCS%output%output_file_id, hdferr)
  call check_hdf_error(hdferr, "Main", "Error creating output HDF5 file at: " &
       // trim(MCS%output%output_filename%chars()))

  ! Go and perform the retrieval process. At this stage, all information from
  ! the config file should have been passed onto the MCS, hence why the main
  ! retrieval function needs no arguments apart from the choice of instrumentm,
  ! and also does not return anything back really.

  call perform_retrievals(my_instrument)


  ! Finishing touches

  ! Close the HDF5 files
  call h5fclose_f(MCS%input%l1b_file_id, hdferr)
  call check_hdf_error(hdferr, "Main", "Error closing input HDF5 file")

  call h5fclose_f(MCS%output%output_file_id, hdferr)
  call check_hdf_error(hdferr, "Main", "Error closing output HDF5 file")

  ! Close the HDF5 library
  call h5close_f(hdferr)
  call check_hdf_error(hdferr, "Main", "Error closing HDF5 library")

  call logger%info("Main", "That's all, folks!")

end program GeoCARBSIF
