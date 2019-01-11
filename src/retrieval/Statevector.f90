module statevector_mod

    use logger_mod, only: logger => master_logger

    type statevector
        ! Number of state vector elements per type
        integer :: num_albedo, num_sif, num_dispersion, num_psurf
        integer, dimension(:), allocatable :: idx_albedo, idx_sif, idx_dispersion, idx_psurf
        ! State vector (current), state vector a priori
        double precision, dimension(:), allocatable :: svsv, svap, sver
        double precision, dimension(:,:), allocatable :: sv_ap_cov, sv_post_cov
    end type

    public initialize_statevector

contains

  subroutine initialize_statevector(sv, count_albedo, count_sif, &
       count_dispersion, count_psurf)

    implicit none
    type(statevector), intent(inout) :: sv
    integer, intent(in) :: count_albedo, count_sif, count_dispersion, count_psurf

    character(len=*), parameter :: fname = "initialize_statevector"
    character(len=999) :: tmp_str
    integer :: sv_count
    integer :: i

    ! First, check if we've got silly values
    if (count_albedo < 0) then
       call logger%fatal(fname, "Albedo SV count < 0.")
       stop 1
    end if

    if (count_sif < 0) then
       call logger%fatal(fname, "SIF SV count < 0.")
       stop 1
    end if

    if (count_dispersion < 0) then
       call logger%fatal(fname, "Albedo SV count < 0.")
       stop 1
    end if

    if (count_dispersion < 0) then
       call logger%fatal(fname, "Albedo SV count < 0.")
       stop 1
    end if

    if ((count_psurf < 0) .or. (count_psurf > 1)) then
       call logger%fatal(fname, "Surface pressure count < 0 or > 1")
       stop 1
    end if

    ! Set The Number of Statevector parameters
    sv%num_albedo = count_albedo
    sv%num_sif = count_sif
    sv%num_dispersion = count_dispersion
    sv%num_psurf = count_psurf

    sv_count = 0
    ! And determine the position (indices) of the state vecto elements within
    ! the state vector

    ! Albedo: arbitrary number of parameters allowed
    if (sv%num_albedo > 0) then

       write(tmp_str, '(A, G0.1)') "Number of albedo SV elements: ", sv%num_albedo
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_albedo(sv%num_albedo))
       do i=1, sv%num_albedo
          sv_count = sv_count + 1
          sv%idx_albedo(i) = sv_count
       end do
    else
       allocate(sv%idx_albedo(1))
       sv%idx_albedo(1) = -1
    end if

    ! SIF: we can do only two things here, SIF magnitude and slope
    if (sv%num_sif > 2) then
       call logger%fatal(fname, "Sorry! Only up to 2 SIF SV elements supported!")
       stop 1
    end if

    if (sv%num_sif > 0) then

       write(tmp_str, '(A, G0.1)') "Number of SIF SV elements: ", sv%num_sif
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_sif(sv%num_sif))
       do i=1, sv%num_sif
          sv_count = sv_count + 1
          sv%idx_sif(i) = sv_count
       end do
    else
       allocate(sv%idx_sif(1))
       sv%idx_sif(1) = -1
    end if

    ! Dispersion: we do arbitrary coefficients here, and we'll have to check
    ! in the retrieval subroutines whether the number passed into this function
    ! actually makes sense.

    if (sv%num_dispersion > 0) then

       write(tmp_str, '(A, G0.1)') "Number of dispersion SV elements: ", sv%num_dispersion
       call logger%info(fname, trim(tmp_str))

       allocate(sv%idx_dispersion(sv%num_dispersion))
       do i=1, sv%num_dispersion
          sv_count = sv_count + 1
          sv%idx_dispersion(i) = sv_count
       end do
    else
       allocate(sv%idx_dispersion(1))
       sv%idx_dispersion(1) = -1
    end if


    ! Surface pressure
    if (sv%num_psurf == 1) then
       write(tmp_str, '(A)') "Number of surface pressure SV elements: 1"
       call logger%info(fname, trim(tmp_str))
       allocate(sv%idx_psurf(1))
       sv_count = sv_count + 1
       sv%idx_psurf = sv_count
    end if

    write(tmp_str, '(A, G0.1, A)') "We have ", sv_count, " SV elements."
    call logger%debug(fname, trim(tmp_str))

    allocate(sv%svap(sv_count))
    allocate(sv%svsv(sv_count))
    allocate(sv%sver(sv_count))

    sv%svap(:) = -9999.99
    sv%svsv(:) = -9999.99
    sv%sver(:) = -9999.99

    allocate(sv%sv_ap_cov(sv_count, sv_count))
    allocate(sv%sv_post_cov(sv_count, sv_count))

  end subroutine initialize_statevector

end module
