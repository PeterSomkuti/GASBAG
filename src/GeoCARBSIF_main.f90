!
! Main GeoCARB SIF Retrieval program. This is where it all starts. This main
! program handles the command-line input via FLAP as well as setting up the
! configuration through the text file via FINER.
!
program GeoCARBSIF

    ! Pick up the logging module for status messages
    use logger_mod, only: logger_init, logger => master_logger
    use startup, only: initialize_config
    use version

    use iso_fortran_env

    implicit none

    ! Before we can read the contents of the logfile, cast it all into /dev/null
    ! logger_init REQUIRES a logfile as an argument
    call logger_init('/dev/null')

    ! Greet the user and display information about the build itself.
    call logger%info("", "Welcome to GeoCARBSIF!")
    call logger%info("", "GIT Branch: " // git_branch)
    call logger%info("", "GIT commit hash: " // git_commit_hash)

    ! Initialize the whole thing by reading the configuration file
    call initialize_config()

end program
