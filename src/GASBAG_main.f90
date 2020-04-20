!> @brief Main GASBAG Retrieval Program
!> @file GASBAG_main.f90
!> @author Peter Somkuti
!>
!> @details
!> This is where it all starts. This main program calls functions to read in
!> the command-line input via FLAP as well as setting up the configuration through
!> the text file via FINER. The instrument is set and initialized right here, and
!> the program goes straight into the perform_retrieval subroutine.
!> After the retrievals are done, the remaining open HDF files are closed, and the
!> program terminates with a zero exit code.

program GASBAG

  ! User modules
  use startup_mod, only: initialize_config
  use version_mod, only: git_branch, git_commit_hash, git_rev_no
  use control_mod, only: CS_t, populate_MCS
  use instruments_mod, only: generic_instrument
  use file_utils_mod, only: check_hdf_error
  use oco2_mod

  ! Third party modules
  use logger_mod, only: logger => master_logger
  use finer, only: file_ini
  use datetime_module

  ! System modules
  use iso_fortran_env
  use HDF5

  implicit none

  ! Local variables
  type(CS_t) :: CS
  type(file_ini) :: fini ! The config file structure
  class(generic_instrument), allocatable :: my_instrument ! The used instrument type
  integer :: hdferr ! HDF error variable


  type(datetime) :: gasbag_start
  type(datetime) :: gasbag_stop
  type(timedelta) :: duration


  ! Start measuring the total program execution time in real
  ! machine clock.
  gasbag_start = gasbag_start%now()
  write(*, '(A, A)') "Started GASBAG at " // gasbag_start%isoformat()

  ! Initilize the HDF5 library program-wide
  call h5open_f(hdferr)
  if (hdferr /= 0) then
     write(*, '(A)') "Error initializing HDF5 library."
     stop 1
  end if

  ! Greet the user and display information about the build itself.
  write(*,'(A)') "========================================="
  write(*,'(A)') " / ___|   / \  / ___|| __ )  / \  / ___| "
  write(*,'(A)') " | |  _  / _ \ \___ \|  _ \ / _ \| |  _  "
  write(*,'(A)') " | |_| |/ ___ \ ___) | |_) / ___ | |_| | "
  write(*,'(A)') " \_____/_/   \_|____/|____/_/   \_\____| "
  write(*,'(A)') "========================================="

  write(*,'(A)') "Version [" // git_branch // " " // git_commit_hash // &
       " #" // git_rev_no // "]"

  ! Initialize the whole thing by reading the configuration file
  call initialize_config(fini)

  ! Initialize the program control_mod structure (MCS) with the settings
  ! from the config file.
  call populate_MCS(fini, CS)

  ! This is where the my_insturment type is properly allocated using one of
  ! the derived types, depending on the instrument specified in the config.
  ! From here on, all derived-type bound subroutines MUST be called through
  ! a (sadly maybe cumbersome) SELECT TYPE statement - but this is just how
  ! Fortran works with run-time polymorphism.

  if (CS%input%instrument_name%lower() == 'oco2') then
     allocate(oco2_instrument :: my_instrument)
     call logger%info("Main", "Using instrument: OCO-2")
  else
     call logger%fatal("Main", "Unknown instrument " // CS%input%instrument_name)
     stop 1
  end if

  ! If the user chooses to overwrite the output, use the H5F_ACC_TRUNC_F key
  ! which does exactly that. If not, use H5F_ACC_EXCL_F, which will return
  ! with an error that we catch afterwards (and terminate the program).
  if (CS%output%overwrite_output) then
     call h5fcreate_f(CS%output%output_filename%chars(), H5F_ACC_TRUNC_F, &
          CS%output%output_file_id, hdferr)
  else
     call h5fcreate_f(CS%output%output_filename%chars(), H5F_ACC_EXCL_F, &
          CS%output%output_file_id, hdferr)

     call check_hdf_error(hdferr, "Main", "Error creating output HDF5 file at: " &
          // trim(CS%output%output_filename%chars()))
  end if

  ! For traceability of results, we store the full config file contents
  ! in the output file itself. (THIS MIGHT NEED REWORKING WHEN MOVING TO NC4)
  call h5gcreate_f(CS%output%output_file_id, "/Metadata", &
       CS%output%metadata_gid, hdferr)
  call check_hdf_error(hdferr, "Main", "Error creating /Metadata group.")
  

  call fini%print(6, retain_comments=.true.)



  select type(my_instrument)
  type is (oco2_instrument)
     ! Scan the L1b file - we need some info from there, mostly the
     ! number of frames, footprints, bands and spectral points
     call my_instrument%scan_l1b_file(CS%input%l1b_filename, CS%general)
  end select


  ! Go and perform the retrieval process. At this stage, all information from
  ! the config file should have been passed onto the MCS, hence why the main
  ! retrieval function needs no arguments apart from the choice of instrumentm,
  ! and also does not return anything back really.
  call perform_retrievals(my_instrument, CS)


  ! Finishing touches

  ! Close the output HDF5 file
  call h5fclose_f(CS%output%output_file_id, hdferr)
  call check_hdf_error(hdferr, "Main", "Error closing output HDF5 file")

  ! Close the HDF5 library
  call h5close_f(hdferr)
  call check_hdf_error(hdferr, "Main", "Error closing HDF5 library")

  ! Say goodbye
  call logger%info("Main", "That's all, folks!")

  ! Let the user know how long it took
  gasbag_stop = gasbag_stop%now()
  write(*, '(A, A)') "Finished GASBAG at " // gasbag_stop%isoformat()

  duration = gasbag_stop - gasbag_start
  write(*,'(A, ES15.5)') "Total duration in seconds: ", duration%total_seconds()
  write(*, '(G0.1, A, G0.1, A, G0.1, A, G0.1, A)') &
       duration%getDays(), "d ", duration%getHours(), "h ", &
       duration%getMinutes(), "m ", duration%getSeconds(), "s"

end program GASBAG
