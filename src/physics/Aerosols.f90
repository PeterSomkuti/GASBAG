  !> @brief Aerosol module
  !> @file Aerosols.f90
  !> @author Peter Somkuti


module aerosols_mod

  ! User modules
  use scene_mod
  use file_utils_mod, only: read_mom_file, read_mie_file
  use control_mod, only: CS_aerosol_t
  use math_utils_mod
  use statevector_mod

  ! Third-party modules
  use stringifor
  use logger_mod, only: logger => master_logger

  !> Generic aerosol type
  type, abstract :: generic_aerosol
     ! Aerosol name (corresponding to MCS%aerosol(..)%name)
     type(string) :: name
  end type generic_aerosol

  !> Gaussian aerosol type
  type, extends(generic_aerosol) :: gauss_aerosol
     !> Total column aerosol optical depth, can be either in
     !> log or linear space.
     double precision :: AOD
     !> Center height of the aerosol Gaussian (in Pascal)
     double precision :: height
     !> Width of the aerosol Gaussian (in Pascal)
     double precision :: width
   contains
     procedure, pass(self) :: print_info => print_my_gauss_info
  end type gauss_aerosol

contains

  !> @brief Print a status message about which type of aerosol this is.
  !> @param aer gauss_aerosol type object
  subroutine print_my_gauss_info(self)
    class(gauss_aerosol), intent(in) :: self

    character(len=*), parameter :: fname = "print_my_gauss_info"
    character(len=999) :: tmp_str

    call logger%debug(fname, "___________________________________")
    call logger%debug(fname, "| Gauss-type aerosol: " // self%name)
    write(tmp_str, '(A, F10.4)') "| Total column OD: ", self%AOD
    call logger%debug(fname, trim(tmp_str))
    write(tmp_str, '(A, F10.1, A)') "| Layer height: ", self%height, " Pa"
    call logger%debug(fname, trim(tmp_str))
    write(tmp_str, '(A, F10.1, A)') "| Layer width : ", self%width, " Pa"
    call logger%debug(fname, trim(tmp_str))
    call logger%debug(fname, "___________________________________") 

  end subroutine print_my_gauss_info


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
    type(CS_aerosol_t), intent(inout) :: aerosol
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
         aerosol%wavelengths, aerosol%coef, aerosol%max_n_coef, success)

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
  !> @param i_win Retrieval window index
  subroutine aerosol_init(scn, i_win, CS_win, CS_aerosol)

    type(scene), intent(inout) :: scn
    integer, intent(in) :: i_win
    type(CS_window_t), intent(in) :: CS_win
    type(CS_aerosol_t), intent(inout) :: CS_aerosol(:)
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


    allocate(scn%op%aer_wl_idx_l(scn%num_aerosols))
    allocate(scn%op%aer_wl_idx_r(scn%num_aerosols))
    allocate(scn%op%aer_mcs_map(scn%num_aerosols))
    allocate(scn%op%aer_ext_q(size(scn%op%wl), scn%num_aerosols))
    allocate(scn%op%aer_sca_q(size(scn%op%wl), scn%num_aerosols))
    allocate(scn%op%aer_ssa(size(scn%op%wl), scn%num_aerosols))

    allocate(scn%op%aer_ext_tau(size(scn%op%wl), &
         scn%num_levels-1, scn%num_aerosols))
    allocate(scn%op%aer_sca_tau(size(scn%op%wl), &
         scn%num_levels-1, scn%num_aerosols))

    allocate(scn%op%aer_ext_tau_edge(2, &
         scn%num_levels-1, scn%num_aerosols))
    allocate(scn%op%aer_sca_tau_edge(2, &
         scn%num_levels-1, scn%num_aerosols))

    allocate(scn%op%reference_aod(scn%num_aerosols))

    ! ------------------------------------------------------
    ! Loop through all aerosols specified in the config file
    ! ------------------------------------------------------

    do j = 1, scn%num_aerosols

       do i = 1, size(CS_aerosol)

          if (CS_win%aerosol(j) /= CS_aerosol(i)%name) cycle


          ! Find out which wavelength regions of the
          ! aerosol files are needed for this band
          idx_l = minloc(abs(CS_aerosol(i)%wavelengths(:) - scn%op%wl(1)), 1)
          idx_r = idx_l + 1

          ! Store these wavelength indices for later use
          scn%op%aer_wl_idx_l(j) = idx_l
          scn%op%aer_wl_idx_r(j) = idx_r

          ! Now check the values against the available data
          if (idx_r > size(CS_aerosol(i)%wavelengths)) then
             call logger%fatal(fname, "Problem initializing aerosols!")
             write(tmp_str, '(A, F15.5, A, F15.5)') "First wavelength of band at ", &
                  scn%op%wl(1), " is closest to the highest-wavelength " &
                  // "value of the aerosol files at "
             call logger%fatal(fname, trim(tmp_str))
             call logger%fatal(fname, "Maybe review the aerosol files?")
             stop 1
          end if

          if (scn%max_pfmom < CS_aerosol(i)%max_n_coef) then
             scn%max_pfmom = CS_aerosol(i)%max_n_coef
          end if

          ! Store the relevant index that allows us to access the aerosol
          scn%op%aer_mcs_map(j) = i

          ! -------------------------------------------------
          !
          ! Calculate the aerosol extinction and scattering
          ! cross sections for every spectral point via a
          ! simple Angstrom exponent Ansatz. Also store the
          ! aerosol SSA's for convenient access.
          !
          ! -------------------------------------------------

          alpha_ext = -log(CS_aerosol(i)%qext(idx_l) / CS_aerosol(i)%qext(idx_r)) &
               / log(CS_aerosol(i)%wavelengths(idx_l) / CS_aerosol(i)%wavelengths(idx_r))

          alpha_sca = -log(CS_aerosol(i)%qsca(idx_l) / CS_aerosol(i)%qsca(idx_r)) &
               / log(CS_aerosol(i)%wavelengths(idx_l) / CS_aerosol(i)%wavelengths(idx_r))

          do l=1, size(scn%op%wl)

             scn%op%aer_ext_q(l, j) = CS_aerosol(i)%qext(idx_l) &
                  * (scn%op%wl(l) / CS_aerosol(i)%wavelengths(idx_l)) ** (-alpha_ext)

             scn%op%aer_sca_q(l, j) = CS_aerosol(i)%qsca(idx_l) &
                  * (scn%op%wl(l) / CS_aerosol(i)%wavelengths(idx_l)) ** (-alpha_ext)

             scn%op%aer_ssa(l, j) = scn%op%aer_sca_q(l, j) / scn%op%aer_ext_q(l, j)

          end do

       end do

    end do

  end subroutine aerosol_init

  subroutine destroy_aerosol(scn)

    type(scene), intent(inout) :: scn

    deallocate(scn%op%aer_wl_idx_l)
    deallocate(scn%op%aer_wl_idx_r)
    deallocate(scn%op%aer_mcs_map)
    deallocate(scn%op%aer_ext_q)
    deallocate(scn%op%aer_sca_q)
    deallocate(scn%op%aer_ssa)

    deallocate(scn%op%aer_ext_tau)
    deallocate(scn%op%aer_sca_tau)

    deallocate(scn%op%aer_ext_tau_edge)
    deallocate(scn%op%aer_sca_tau_edge)

    deallocate(scn%op%reference_aod)

  end subroutine destroy_aerosol


  subroutine insert_aerosols_in_scene(SV, CS_win, CS_aerosol, scn, scn_aer)

    type(statevector), intent(in) :: SV
    type(CS_window_t), intent(in) :: CS_win
    type(CS_aerosol_t), intent(inout) :: CS_aerosol(:)
    type(scene), intent(inout) :: scn
    class(generic_aerosol), intent(inout) :: scn_aer(:)

    integer :: i_win

    character(len=*), parameter :: fname = "insert_aerosols_in_scene"
    integer :: i

    select type(scn_aer)

    !------------------------------------------
    !
    !   GAUSSIAN AEROSOL DISTRIBUTION
    !
    !------------------------------------------
    type is (gauss_aerosol)

       do i = 1, scn%num_aerosols

          scn_aer(i)%name = CS_win%aerosol(i)

          ! Feed in the default values into the aerosol distribution objects
          ! Note that these values are supposed to be in real-space (no log!)
          scn_aer(i)%AOD = CS_aerosol(scn%op%aer_mcs_map(i))%default_aod
          scn_aer(i)%height = CS_aerosol(scn%op%aer_mcs_map(i))%default_height * scn%atm%psurf
          scn_aer(i)%width = CS_aerosol(scn%op%aer_mcs_map(i))%default_width

       end do

       ! If any of the aerosol parameters are retrieved, replace those values here
       ! with those coming from the state vector.

       do i = 1, SV%num_aerosol_aod
          if (CS_win%aerosol_retrieve_aod_log(SV%aerosol_aod_idx_lookup(i))) then
             ! AOD supplied in log-space
             scn_aer(SV%aerosol_aod_idx_lookup(i))%AOD = exp(SV%svsv(SV%idx_aerosol_aod(i)))
          else
             ! AOD supplied in linear space
             scn_aer(SV%aerosol_aod_idx_lookup(i))%AOD = SV%svsv(SV%idx_aerosol_aod(i))
          end if
       end do

       do i = 1, SV%num_aerosol_height
          if (CS_win%aerosol_retrieve_height_log(SV%aerosol_height_idx_lookup(i))) then
             ! Aerosol height supplied in log-space
             scn_aer(SV%aerosol_height_idx_lookup(i))%height = exp(SV%svsv(SV%idx_aerosol_height(i)))
          else
             ! Aerosol height supplied in linear space
             scn_aer(SV%aerosol_height_idx_lookup(i))%height = SV%svsv(SV%idx_aerosol_height(i))
          end if
          ! Aerosol height is given as fraction of psurf
          scn_aer(SV%aerosol_height_idx_lookup(i))%height = &
               scn_aer(SV%aerosol_height_idx_lookup(i))%height * scn%atm%psurf
       end do

       ! This function distributes the aerosols according to the data given in
       ! the scn_aer object.
       call aerosol_gauss_shape(CS_aerosol, scn, scn_aer)

       ! Print some debug information on each aerosol used in the scene.
       do i = 1, scn%num_aerosols
          call scn_aer(i)%print_info()
       end do

    !------------------------------------------
    !
    !   DEFAULT CASE WHICH TERMINATES THE PROGRAM
    !
    !------------------------------------------
    class default
       call logger%fatal(fname, "NOT SUPPORTED AEROSOL DISTRIBUTION TYPE. MAJOR CODE PROBLEM.")
       stop 1
    end select


  end subroutine insert_aerosols_in_scene


  subroutine aerosol_gauss_shape(CS_aerosol, scn, scn_aer)

    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    type(scene), intent(inout) :: scn
    type(gauss_aerosol), intent(inout) :: scn_aer(:)
    
    double precision :: aod_norm
    integer :: aer
    integer :: wl
    integer :: lay, j

    character(len=*), parameter :: fname = "aerosol_gauss_shape"
    character(len=999) :: tmp_str


    do aer = 1, scn%num_aerosols

       ! First calculate ext/sca for the wavelenghts for which the aerosols are
       ! actually specified / calculated. This is usually at two wavelengths near
       ! the band edges.

       ! We also store the reference AOD for this scene (iteration) for this aerosol type.
       ! Comes in handy when we need to calculate jacobians.
       scn%op%reference_aod(aer) = scn_aer(aer)%AOD

       aod_norm = sum(exp(-((scn%atm%p_layers(:) - scn_aer(aer)%height)**2) &
                  / (2 * scn_aer(aer)%width * scn_aer(aer)%width)))

       do wl = 1, 2

          do lay = 1, scn%num_active_levels - 1
             scn%op%aer_ext_tau_edge(wl,lay,aer) = exp( &
                  -((scn%atm%p_layers(lay) - scn_aer(aer)%height)**2) &
                  / (2 * scn_aer(aer)%width * scn_aer(aer)%width)) &
                  * scn_aer(aer)%AOD / aod_norm
          end do

          ! We are grabbing the left and right hand side SSAs from the file to compute
          ! the scattering contribution to the extinction.
          if (wl == 1) then
             scn%op%aer_sca_tau_edge(wl,:,aer) = scn%op%aer_ext_tau_edge(wl,:,aer) &
                  * CS_aerosol(scn%op%aer_mcs_map(aer))%ssa(scn%op%aer_wl_idx_l(aer))
          else if (wl == 2) then
             scn%op%aer_sca_tau_edge(wl,:,aer) = scn%op%aer_ext_tau_edge(wl,:,aer) &
                  * CS_aerosol(scn%op%aer_mcs_map(aer))%ssa(scn%op%aer_wl_idx_r(aer))
          end if

          where(scn%op%aer_ext_tau_edge(wl,:,aer) < 1.0d-10) scn%op%aer_ext_tau_edge(wl,:,aer) = 1d-10
          where(scn%op%aer_sca_tau_edge(wl,:,aer) < 1.0d-10) scn%op%aer_sca_tau_edge(wl,:,aer) = 1d-10

          ! -----------------
          ! Debug information
          ! -----------------

          if (wl == 1) then

             write(tmp_str, '(A,A)') "Aerosol extinction, scattering and SSA for left edge and aerosol: ", &
                  CS_aerosol(scn%op%aer_mcs_map(aer))%name%chars()
             call logger%debug(fname, trim(tmp_str))
             do lay = 1, scn%num_levels - 1

                write(tmp_str,'(A, I3, ES15.5, ES15.5, ES15.5)') "Layer: ", lay, &
                     scn%op%aer_ext_tau_edge(1, lay, aer), scn%op%aer_sca_tau_edge(1, lay, aer), &
                     scn%op%aer_sca_tau_edge(1, lay, aer) / scn%op%aer_ext_tau_edge(1, lay, aer)
                call logger%debug(fname, trim(tmp_str))
             end do

          end if

       end do

       ! Now do it for all wavelengths in our hires grid
       do wl = 1, size(scn%op%wl)

          do lay = 1, scn%num_active_levels - 1
             scn%op%aer_ext_tau(wl,lay,aer) = exp( &
                  -((scn%atm%p_layers(lay) - scn_aer(aer)%height)**2) &
                  / (2 * scn_aer(aer)%width * scn_aer(aer)%width)) &
                  * scn_aer(aer)%AOD / aod_norm
          end do

          scn%op%aer_sca_tau(wl,:,aer) = scn%op%aer_ext_tau(wl,:,aer) * scn%op%aer_ssa(wl, aer)

          ! Bump up tiny values to some lower threshold
          ! For these very small values, it doesn't
          ! matter if the aerosol SSA ends up being 1.0
          where(scn%op%aer_ext_tau(wl,:,aer) < 1.0d-10) scn%op%aer_ext_tau(wl,:,aer) = 1d-10
          where(scn%op%aer_sca_tau(wl,:,aer) < 1.0d-10) scn%op%aer_sca_tau(wl,:,aer) = 1d-10

       end do

    end do

  end subroutine aerosol_gauss_shape


  subroutine calculate_aero_height_factors(layer_p, aero_height, aero_width, &
       factor)

    double precision, intent(in) :: layer_p(:)
    double precision, intent(in) :: aero_height
    double precision, intent(in) :: aero_width
    double precision, intent(inout) :: factor(:)

    double precision, allocatable :: aero_shape(:)
    integer :: i, j

    allocate(aero_shape(size(layer_p)))

    aero_shape(:) = 0.0d0
    factor(:) = 0.0d0

    aero_shape(:) = exp(-(layer_p(:) - aero_height)**2 / (2 * aero_width**2))
    factor(:) = (layer_p(:) - aero_height) / (aero_width**2)

    do i = 1, size(layer_p)
       factor(i) = factor(i) - sum(aero_shape(:) * (layer_p(:) - aero_height) / (aero_width**2)) / sum(aero_shape)
    end do

  end subroutine calculate_aero_height_factors


end module aerosols_mod
