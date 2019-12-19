  !> @brief Aerosol module
  !> @file Aerosols.f90
  !> @author Peter Somkuti


module aerosols_mod

  ! User modules
  use file_utils_mod, only: read_mom_file, read_mie_file
  use physical_model_addon_mod
  use control_mod, only: MCS, CS_aerosol
  use math_utils_mod
  
  ! Third-party modules
  use stringifor
  use logger_mod, only: logger => master_logger

  public :: ingest_aerosol_files

contains


  !> @brief Feed in the aerosol data, given a mom/mie file combination
  !> @param name Identifier
  !> @param aerosol Aerosol object
  !> @param mom_filename Location of mom file
  !> @param mie_filename Location of mie file
  subroutine ingest_aerosol_files(aerosol)

    ! For a CS_aerosol object, which is part of MCS, we have already the
    ! filenames for the mie/mom files that belong to this particular aerosol,
    ! so all this subroutine has to do, is to load the files and stick the
    ! coefficients etc. into the various arrays.

    implicit none
    type(CS_aerosol), intent(inout) :: aerosol
    logical :: success

    ! Function name
    character(len=*), parameter :: fname = "create_aerosol"
    ! Temporary holder for the wavelengths, we need that to check whether the
    ! mie and mom files 'match'
    double precision, allocatable :: wavelengths_tmp(:)
    integer :: i ! Loop variable
    character(len=999) :: tmp_str

    call logger%debug(fname, "Reading in data for aerosol: " // aerosol%name%chars())

    ! Attempt at read-in of mom file
    call logger%debug(fname, "Reading in mom file: " // aerosol%mom_filename%chars())
    call read_mom_file(aerosol%mom_filename%chars(), &
         aerosol%wavelengths, aerosol%coef, success)

    ! Attempt at read-in of mie file
    call logger%debug(fname, "Reading in mie file: " // aerosol%mie_filename%chars())
    call read_mie_file(aerosol%mie_filename%chars(), &
         wavelengths_tmp, aerosol%qext, &
         aerosol%qsca, aerosol%ssa, &
         aerosol%sigma_ext, aerosol%reff, success)

    ! Do some checks to confirm that the mie and mom file
    ! 'belong' together by checking if the wavelengths are the same.
    if (size(wavelengths_tmp) /= size(aerosol%wavelengths)) then
       call logger%fatal(fname, "Mie and Mom files have different numbers of wavelengths.")
       stop 1
    end if

    ! Check if the wavelength values in both mom and mie file correspond. There are obviously
    ! rounding errors, and probably depend on the code that wrote the files, so rather than
    ! terminating the program, just spit out a warning. The user should know whether this
    ! is acceptable or not.
    do i=1, size(aerosol%wavelengths)
       if (abs(aerosol%wavelengths(i) - wavelengths_tmp(i)) > 1.0d-3) then

          call logger%warning(fname, "Mie and Mom files have different wavelength values (>1nm).")
          write(tmp_str, '(A, ES10.6)') "Mie file: ", wavelengths_tmp(i)
          call logger%warning(fname, trim(tmp_str))
          write(tmp_str, '(A, ES10.6)') "Mom file: ", aerosol%wavelengths(i)
          call logger%warning(fname, trim(tmp_str))
          
       end if
    end do

  end subroutine ingest_aerosol_files


  !> @brief Initialize the required aerosol data
  !> @param scn Scene object
  subroutine aerosol_init(scn)

    type(scene), intent(inout) :: scn

    ! Local
    character(len=*), parameter :: fname = "aerosol_init"
    character(len=999) :: tmp_str
    integer :: i, j, l
    integer :: count_aer
    integer :: idx_l, idx_r

    !> Exctinction Angstrom coefficient
    double precision :: alpha_ext
    !> Scattering Angstrom coefficient
    double precision :: alpha_sca

    ! ----------------------------------
    ! Find out how many aerosols we have
    ! ----------------------------------

    count_aer = 0
    do i = 1, size(MCS%aerosol)
       if (MCS%aerosol(i)%used) then
          count_aer = count_aer + 1
       end if
    end do

    scn%num_aerosols = count_aer
    scn%max_pfmom = 0

    allocate(scn%op%aer_ext_q(size(scn%op%wl), scn%num_aerosols))
    allocate(scn%op%aer_sca_q(size(scn%op%wl), scn%num_aerosols))
    allocate(scn%op%aer_ssa(size(scn%op%wl), scn%num_aerosols))

    allocate(scn%op%aer_ext_tau(size(scn%op%wl), &
         scn%num_levels-1, scn%num_aerosols))
    allocate(scn%op%aer_sca_tau(size(scn%op%wl), &
         scn%num_levels-1, scn%num_aerosols))

    ! -------------------------------
    ! Loop through all used aerosols
    ! -------------------------------

    j = 0
    do i = 1, size(MCS%aerosol)
       if (MCS%aerosol(i)%used) then
          j = j + 1

          ! Find out which wavelength regions of the
          ! aerosol files are needed for this band
          idx_l = minloc(abs(MCS%aerosol(i)%wavelengths(:) - scn%op%wl(1)), 1)
          idx_r = idx_l + 1

          ! Now check the values against the available data
          if (idx_r > size(MCS%aerosol(i)%wavelengths)) then
             call logger%fatal(fname, "Problem initializing aerosols!")
             write(tmp_str, '(A, F15.5, A, F15.5)') "First wavelength of band at ", &
                  scn%op%wl(1), " is closest to the highest-wavelength " &
                  // "value of the aerosol files at "
             call logger%fatal(fname, trim(tmp_str))
             call logger%fatal(fname, "Maybe review the aerosol files?")
             stop 1
          end if

          if (scn%max_pfmom < MCS%aerosol(i)%max_n_coef ) then
             scn%max_pfmom = MCS%aerosol(i)%max_n_coef
          end if

          ! -------------------------------------------------
          !
          ! Calculate the aerosol extinction and scattering
          ! cross sections for every spectral point via a
          ! simple Angstrom exponent Ansatz. Also store the
          ! aerosol SSA's for convenient access.
          !
          ! -------------------------------------------------

          alpha_ext = -log(MCS%aerosol(i)%qext(idx_l) / MCS%aerosol(i)%qext(idx_r)) &
               / log(MCS%aerosol(i)%wavelengths(idx_l) / MCS%aerosol(i)%wavelengths(idx_r))
          
          alpha_sca = -log(MCS%aerosol(i)%qsca(idx_l) / MCS%aerosol(i)%qsca(idx_r)) &
               / log(MCS%aerosol(i)%wavelengths(idx_l) / MCS%aerosol(i)%wavelengths(idx_r))

          do l=1, size(scn%op%wl)

             scn%op%aer_ext_q(l, j) = MCS%aerosol(i)%qext(idx_l) &
                  * (scn%op%wl(l) / MCS%aerosol(i)%wavelengths(idx_l)) ** (-alpha_ext)

             scn%op%aer_sca_q(l, j) = MCS%aerosol(i)%qsca(idx_l) &
                  * (scn%op%wl(l) / MCS%aerosol(i)%wavelengths(idx_l)) ** (-alpha_ext)

             scn%op%aer_ssa(l, j) = scn%op%aer_sca_q(l, j) / scn%op%aer_ext_q(l, j)

          end do

       end if
    end do

  end subroutine aerosol_init


  subroutine aerosol_gauss_shape(scn)

    type(scene), intent(inout) :: scn

    double precision, parameter :: aero_height = 4000.0
    double precision, parameter :: aero_aod = 1.0
    double precision, parameter :: aero_width = 1000.0


    double precision, allocatable :: layer_height(:)
    double precision :: aod_norm
    integer :: aer, wl, lay

    allocate(layer_height(scn%num_levels - 1))

    ! Get layer heights in advance (mean height between two levels)
    do lay = 1, scn%num_levels - 1
       layer_height(lay) = 0.5d0 * (scn%atm%altitude_levels(lay) &
            + scn%atm%altitude_levels(lay+1))
    end do

    do aer = 1, scn%num_aerosols
       do wl = 1, size(scn%op%wl)
          do lay = 1, scn%num_levels - 1

             scn%op%aer_ext_tau(wl,lay,aer) = exp(-((layer_height(lay) - aero_height)**2) &
                  / (2 * aero_width * aero_width))

          end do

          if (wl == 1) then
             aod_norm = sum(scn%op%aer_ext_tau(wl,:,aer)) / aero_aod
          end if

          scn%op%aer_ext_tau(wl,:,aer) = scn%op%aer_ext_tau(wl,:,aer) / aod_norm
          scn%op%aer_sca_tau(wl,:,aer) = scn%op%aer_ext_tau(wl,:,aer) * scn%op%aer_ssa(wl, aer)

          ! Bump up tiny values to some lower threshold
          ! For these very small values, it doesn't matter if the aerosol SSA
          ! ends up being 1.0
          where(scn%op%aer_ext_tau(wl,:,aer) < 1.0d-10) scn%op%aer_ext_tau(wl,:,aer) = 1d-10
          where(scn%op%aer_sca_tau(wl,:,aer) < 1.0d-10) scn%op%aer_sca_tau(wl,:,aer) = 1d-10

          if (wl == 1) then
             do lay = 1, scn%num_levels - 1
                write(*,*) lay, layer_height(lay), scn%op%aer_ext_tau(wl,lay,aer), scn%op%aer_sca_tau(wl,lay,aer)
             end do
          end if

       end do

    end do
    


  end subroutine aerosol_gauss_shape


  subroutine precompute_all_coef(scn, n_stokes, coef, lcoef)

    type(scene), intent(in) :: scn
    integer, intent(in) :: n_stokes
    double precision, allocatable, intent(inout) :: coef(:,:,:,:)
    double precision, allocatable, intent(inout) :: lcoef(:,:,:,:,:)

    integer :: nmom
    integer :: n_layers

    n_layers = scn%num_active_levels - 1

    ! Depending on whether we use polarization or not,
    ! we only need to do a certain number of phase matrix
    ! elements.

    if (n_stokes == 1) then
       nmom = 1
    else if (n_stokes == 3) then
       nmom = 6
    end if

    allocate(coef(size(scn%op%wl), scn%max_pfmom, nmom, n_layers))
    allocate(lcoef(size(scn%op%wl), scn%max_pfmom, nmom, 1, n_layers))

    coef(:,:,:,:) = 0.0d0
    lcoef(:,:,:,:,:) = 0.0d0


  end subroutine precompute_all_coef


end module aerosols_mod
