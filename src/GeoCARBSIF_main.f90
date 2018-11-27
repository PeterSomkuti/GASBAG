!
! Main GeoCARB SIF Retrieval program. This is where it all starts. This main
! program handles the command-line input via FLAP as well as setting up the
! configuration through the text file via FINER.
!
program GeoCARBSIF

    use startup, only: initialize_config
    use version, only: git_branch, git_commit_hash

    use iso_fortran_env

    implicit none

    ! Greet the user and display information about the build itself.
    write(*,'(A)') "Welcome to GeoCARBSIF!"
    write(*, '(A)') "Version [" // git_branch // " " // git_commit_hash // "]"

    ! Initialize the whole thing by reading the configuration file
    call initialize_config()

end program
