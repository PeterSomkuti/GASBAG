!> @brief Module that contains functions for the PCA-based RT acceleration
!> @file PCART.f90
!> @author Peter Somkuti

module PCART_mod

  ! User modules
  use control_mod, only: CS_aerosol_t
  use scene_mod
  use math_utils_mod

  ! System modules
  use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_is_nan, ieee_is_finite

  implicit none

  !> Contains all the relevant data for a bin
  type PCA_bin_t
     !> Number of spectral points in this bin
     integer :: N_spect
     !> Number of used EOFs in this bin
     integer :: N_EOF
     !> Matrix of optical parameters (N_spect, N_opt)
     double precision, allocatable :: F(:,:)
     !> F_centered is the mean-removed F (N_spect, N_opt)
     double precision, allocatable :: F_centered(:,:)
     !> Correlation matrix
     double precision, allocatable :: C(:,:)
     !> Eigenvectors, Eigenvalues, Principal Components and scaled EOFs
     double precision, allocatable :: EigVec(:,:)
     !> Scaled EOFs
     double precision, allocatable :: sEOF(:,:)
     !> Eigenvalues
     double precision, allocatable :: EigVal(:)
     !> Principal components
     double precision, allocatable :: PC(:,:)
     !> Mean optical properties (N_opt)
     double precision, allocatable :: mean_opt(:)
     !> Perturbed optical states
     double precision, allocatable :: pert_opt_states(:,:)
  end type PCA_bin_t

  !> User-type that contains data for PCA-based fast RT
  type PCA_handler_t
     !> Assignment of spectral index to PCA bin
     integer, allocatable :: assigned_bin(:)
     !> Assignment of PCA bin index to assigned bin index
     integer, allocatable :: PCA_bin_idx_map(:)
     !> Number of non-empty bins
     integer :: N_bin
     !> Number of optical properties used
     integer :: N_opt
     !> Start and end indices for optical properties
     integer :: s_tau, e_tau, &
          s_ray, e_ray, &
          s_aer, e_aer, &
          s_fra, e_fra, &
          s_sur, e_sur

     !> PCA bin objects
     type(PCA_bin_t), allocatable :: PCA_bin(:)
  end type PCA_handler_t

contains


  !> @brief Initilizes the full PCA routine
  !> @param scn Scene object
  !> @param PCA_handler PCA_handler object
  subroutine initialize_PCA(scn, PCA_handler)
    type(scene), intent(in) :: scn
    type(PCA_handler_t), intent(inout) :: PCA_handler

    ! Initialize with the total number of spectral points
    allocate(PCA_handler%assigned_bin(size(scn%op%wl)))
    PCA_handler%assigned_bin(:) = -1
    ! Calculate the number of optical properties that enter
    ! the calculations.
    ! At the moment:
    ! 2 * Number of layers
    ! + number of aerosols
    ! + number of surface parameters
    ! + 1 for wavelength interpolation coefficient

    PCA_handler%N_opt = 2 * (scn%num_active_levels - 1)

    ! Add the number of aerosols
    ! (regardless of number of retrieved aerosol parameters)
    PCA_handler%N_opt = PCA_handler%N_opt + scn%num_aerosols

    ! Add one for wavelength interpolation coefficient
    PCA_handler%N_opt = PCA_handler%N_opt + 1

    ! Add the number surface parameters
    ! (regardless of number of retrieved surface parameters)
    ! (for Lambertian-type surface, this is 1)
    PCA_handler%N_opt = PCA_handler%N_opt + 1

    ! Determine the start and end indices for the various
    ! portions of the optical state matrix F

    ! Total optical depth
    PCA_handler%s_tau = 1
    PCA_handler%e_tau = scn%num_active_levels - 1 ! nlay

    ! Rayleigh optical depth
    PCA_handler%s_ray = PCA_handler%e_tau + 1
    PCA_handler%e_ray = PCA_handler%s_ray + scn%num_active_levels - 2! + nlay - 1

    ! Surface parameteres (so far only Lambertian surface here)
    PCA_handler%s_sur = PCA_handler%e_ray + 1
    PCA_handler%e_sur = PCA_handler%s_sur

    ! Aerosol/wavelength interpolation coefficient (only one needed)
    PCA_handler%s_fra = PCA_handler%e_sur + 1
    PCA_handler%e_fra = PCA_handler%s_fra

    if (allocated(scn%op%reference_aod)) then
       ! Aerosol scattering ratio/coefficient
       PCA_handler%s_aer = PCA_handler%e_fra + 1
       PCA_handler%e_aer =  PCA_handler%s_aer + scn%num_aerosols - 1
    else
       ! Set these indices to negative values, so the code
       ! will throw an error if elements of F are accessed
       ! that shouldn't be accessed.
       PCA_handler%s_aer = -1
       PCA_handler%e_aer = -1
    end if

  end subroutine initialize_PCA

  !> @brief Binning method according to Kopparla et al.
  !> @param scn Scene object
  !> @param PCA_handler PCA handler object
  subroutine assign_pushkar_bins(scn, PCA_handler)
    type(scene), intent(in) :: scn
    type(PCA_handler_t), intent(inout) :: PCA_handler

    ! Manually set bin boundaries in total tau_gas space
    ! huge(0.0d0) is largest double precision number
    !double precision, parameter :: bounds(12) = &
    !     [0d0, 0.01d0, 0.025d0, 0.05d0, &
    !      0.1d0, 0.25d0, 0.5d0, 0.625d0, &
    !      0.75d0, 1.0d0, 5.0d0, huge(0.0d0)]

    double precision, parameter :: bounds(12) = &
         [0d0, 0.005d0, 0.01d0, 0.015d0, 0.025d0, 0.05d0, &
         0.1d0, 0.25d0, 0.5d0, &
         0.75d0, 1.0d0, huge(0.0d0)]

    integer :: nspect, i, j

    integer, allocatable :: single_assigned_bins(:)

    double precision, allocatable :: gas_tot(:)
    double precision, allocatable :: ssa_tot(:)
    double precision :: ssa_median
    double precision :: total_aod


    character(len=*), parameter :: fname = "assign_pushkar_bins"
    character(len=99) :: tmp_str

    nspect = size(PCA_handler%assigned_bin)

    allocate(ssa_tot(nspect))
    allocate(gas_tot(nspect))

    allocate(single_assigned_bins(nspect))
    single_assigned_bins(:) = 0

    ! Sum up the gas opds to get the total taug
    ! taug_tot(:) = sum(taug, 1)
    ssa_tot(:) = sum(scn%op%layer_omega(:, 1:scn%num_active_levels - 1), 2)
    gas_tot(:) = sum(sum(scn%op%gas_tau(:,:,:), 2), 2)

    if (scn%num_aerosols > 0) then
       total_aod = sum(scn%op%reference_aod(:))
    else
       total_aod = 0.0d0
    end if

    ! Assign each spectral point to a bin according to taug
    ! is there a neat way to do an equivalent to np.searchsorted
    do i = 1, nspect
       do j = 1, size(bounds) - 1
          if ( &
               (max(gas_tot(i), 0.0d0) .ge. bounds(j)) .and. &
               (max(gas_tot(i), 0.0d0) .lt. bounds(j+1))) then
             single_assigned_bins(i) = j
          endif
       enddo
    enddo

    ! Now split each bin along median SSA
    ! Note that the total possible bins now is 2*size(bounds) - 2

    j = 1
    do i = 1, size(bounds) - 1

       if (count(single_assigned_bins == i) == 0) then
          cycle
       end if

       ! calculate the median of total SSA within bin "i"
       ssa_median = percentile(pack(ssa_tot, single_assigned_bins == i), 50.0d0)

       ! Find all scenes with a total SSA <= than median and assign it to bin index "j"
       where ((single_assigned_bins == i) .and. (ssa_tot <= ssa_median))
          PCA_handler%assigned_bin = j
       endwhere
       j = j + 1

       ! Find all scenes with a total SSA > than median and assign it to bin index "j+1"
       where ((single_assigned_bins == i) .and. (ssa_tot > ssa_median))
          PCA_handler%assigned_bin = j
       endwhere
       j = j + 1

    enddo

    ! Now, if there is only a single wavelength in a bin (which can happen),
    ! we need to move that spectral point over to some other bin

    do i = 1, maxval(PCA_handler%assigned_bin)

       if (count(PCA_handler%assigned_bin == i) == 1) then
          ! Loop through next higher bin and see if that one is non-empty
          do j = i + 1, maxval(PCA_handler%assigned_bin)
             if (count(PCA_handler%assigned_bin == j) > 0) then
                where(PCA_handler%assigned_bin == i) PCA_handler%assigned_bin = j
                exit
             end if
          end do

       end if

    end do



    PCA_handler%N_bin = 0
    ! Debug: print the number of points in each bin
    i = 1
    do j = 1, 2 * size(bounds) - 2

       write(tmp_str, "(A, G3.1, A, E10.5, A, E10.5, A, G7.1, A)") &
            "Bin ", j, " [", bounds(i), " - ", bounds(i+1), "]: ", &
            count(PCA_handler%assigned_bin(:) == j), " points"
       call logger%debug(fname, trim(tmp_str))

       if (count(PCA_handler%assigned_bin(:) == j) > 0) then
          PCA_handler%N_bin = PCA_handler%N_bin + 1
       end if

       if ((modulo(j, 2) == 0) .and. j > 1) i = i + 1
    enddo

    ! Also let the user know if there are any unassigned spectral points
    if (count(PCA_handler%assigned_bin(:) == -1) > 0) then
       write(tmp_str, "(A, G0.1)") "Unassigned points: ", count(PCA_handler%assigned_bin(:) == -1)
       call logger%error(fname, trim(tmp_str))
    end if


    ! How many non-empty bins do we have?
    write(tmp_str, "(A, G0.1)") "Number of non-empty bins: ", PCA_handler%N_bin
    call logger%debug(fname, trim(tmp_str))

    ! Map the 1-based PCA bin idx with the spectral bin assigment
    allocate(PCA_handler%PCA_bin_idx_map(PCA_handler%N_bin))
    i = 1
    do j = 1, 2 * size(bounds) - 2
       if (count(PCA_handler%assigned_bin(:) == j) > 0) then
          PCA_handler%PCA_bin_idx_map(i) = j
          i = i + 1
       end if
    end do

  end subroutine assign_pushkar_bins


  
  !> @brief Allocates bin-related quantities (needs to be done
  !> AFTER assignment into bins using some scheme)
  !> @param PCA_handler PCA handler object
  subroutine allocate_first_PCA(PCA_handler)
    type(PCA_handler_t), intent(inout) :: PCA_handler

    character(len=99) :: tmp_str
    character(len=*), parameter :: fname = "allocate_PCA"

    integer :: i
    integer :: nspect, nopt

    nopt = PCA_handler%N_opt

    ! Allocate the main object for nonempty bins
    call logger%debug(fname, "Allocating PCA bin objects and associated arrays.")
    allocate(PCA_handler%PCA_bin(PCA_handler%N_bin))

    ! Allocate some of the arrays within, and set some variables
    ! NOTE: the number of required EOFs will be calculated in the
    !       next step, so matrices related to that can only be
    !       allocated in a second pass
    do i = 1, PCA_handler%N_bin

       nspect = count(PCA_handler%assigned_bin == PCA_handler%PCA_bin_idx_map(i))

       allocate(PCA_handler%PCA_bin(i)%F(nspect, nopt))
       allocate(PCA_handler%PCA_bin(i)%F_centered(nspect, nopt))
       allocate(PCA_handler%PCA_bin(i)%mean_opt(nopt))

       PCA_handler%PCA_bin(i)%N_spect = nspect
       PCA_handler%PCA_bin(i)%N_EOF = nopt

    end do
  end subroutine allocate_first_PCA

  !> @brief Fills F-matrix (optical states/properties) with values
  !> @param scn Scene object
  !> @param PCA_handler PCA handler object
  subroutine create_F_matrix(scn, PCA_handler)
    type(scene), intent(in) :: scn
    type(PCA_handler_t), intent(inout) :: PCA_handler

    integer :: i, j, k, l

    integer :: n_lay

    integer :: s_tau, e_tau
    integer :: s_ray, e_ray
    integer :: s_sur, e_sur
    integer :: s_aer, e_aer
    integer :: s_fra, e_fra

    n_lay = scn%num_active_levels - 1

    s_tau = PCA_handler%s_tau
    e_tau = PCA_handler%e_tau

    s_ray = PCA_handler%s_ray
    e_ray = PCA_handler%e_ray

    s_sur = PCA_handler%s_sur
    e_sur = PCA_handler%e_sur

    s_aer = PCA_handler%s_aer
    e_aer = PCA_handler%e_aer

    s_fra = PCA_handler%s_fra
    e_fra = PCA_handler%e_fra

    do j = 1, PCA_handler%N_bin
       k = 1
       do i = 1, size(scn%op%wl)

          if (PCA_handler%assigned_bin(i) == PCA_handler%PCA_bin_idx_map(j)) then
             ! spectral point "i" belongs to bin "PCA_bin_idx_map(j)"

             ! Copies layer-resolved total optical depth
             PCA_handler%PCA_bin(j)%F(k, s_tau:e_tau) = scn%op%layer_tau(i, 1:n_lay)
             ! Copies layer-resolved Rayleigh optical depth
             PCA_handler%PCA_bin(j)%F(k, s_ray:e_ray) = scn%op%ray_tau(i, 1:n_lay)
             ! Copies surface parameter
             PCA_handler%PCA_bin(j)%F(k, s_sur:e_sur) = scn%op%albedo(i)
             ! Copies aerosol interpolation fraction (needed for wavelength)
             PCA_handler%PCA_bin(j)%F(k, s_fra:e_fra) = scn%op%aer_frac(i)

             ! Aerosols present?
             if (scn%num_aerosols > 0) then
                PCA_handler%PCA_bin(j)%F(k, s_aer:e_aer) = scn%op%aer_sca_q(i, :)
             end if

             ! Increment wavelength counter
             k = k + 1

          end if

       end do
    end do

  end subroutine create_F_matrix


  !> @brief Peform forward transform of optical property matrix
  !> @param PCA_handler PCA handler object
  subroutine forward_transform_F_matrix(PCA_handler)
    type(PCA_handler_t), intent(inout) :: PCA_handler
    integer :: wl, bin, opt

    do bin = 1, PCA_handler%N_bin
          do opt = 1, PCA_handler%N_opt

             if (opt == PCA_handler%s_fra) then
                PCA_handler%PCA_bin(bin)%F(:, opt) = &
                     PCA_handler%PCA_bin(bin)%F(:, opt)
             else
                PCA_handler%PCA_bin(bin)%F(:, opt) = &
                     log(PCA_handler%PCA_bin(bin)%F(:, opt))
             end if

          end do
    end do

  end subroutine forward_transform_F_matrix


  !> @brief Performs mean-removal of optical property matrix
  !> @param PCA_handler PCA handler object
  subroutine F_matrix_mean_removal(PCA_handler)
    type(PCA_handler_t) :: PCA_handler
    integer :: bin, opt

    do bin=1, PCA_handler%N_bin
       do opt=1, PCA_handler%N_opt

          ! Compute the mean optical property and subsetted into bins
          PCA_Handler%PCA_bin(bin)%mean_opt(opt) = sum(PCA_handler%PCA_bin(bin)%F(:, opt)) / &
               PCA_handler%PCA_bin(bin)%N_spect
          PCA_Handler%PCA_bin(bin)%F_centered(:, opt) = PCA_Handler%PCA_bin(bin)%F(:, opt) - &
               PCA_Handler%PCA_bin(bin)%mean_opt(opt)

       end do
    end do

  end subroutine F_matrix_mean_removal


  !> @brief Peforms the eigenproblem calculations for each bin
  !> @param PCA_handler PCA handler object
  !> @param success Was the calculation successful?
  subroutine perform_PCA(PCA_handler, success)
    type(PCA_handler_t), intent(inout) :: PCA_handler
    logical, intent(inout) :: success

    character(len=*), parameter :: fname = "perform_PCA"
    character(len=99) :: tmp_str

    external DSYEVD
    external DSYEV

    integer :: bin, opt, eof, wl

    ! covariance matrix (need a temporary container)
    double precision, allocatable :: C(:,:)
    ! array to hold eigenvalues temporarily
    double precision, allocatable :: tmp_eigval(:)

    double precision, allocatable :: work(:), iwork(:)
    integer :: lwork, liwork, info, info2

    integer :: N_opt
    integer :: N_spect
    integer :: N_EOF

    success = .false.

    N_opt = PCA_handler%N_opt

    LWORK = 10 + 6*PCA_handler%N_opt + 2*(PCA_handler%N_opt**2)
    LIWORK = 10 + 5*PCA_handler%N_opt

    ! Required workspaces for DSYEVD and DSYEV
    allocate(work(lwork))
    allocate(iwork(liwork))
    allocate(tmp_eigval(N_opt))
    allocate(C(N_opt, N_opt))

    do bin = 1, PCA_handler%N_bin

       N_spect = PCA_handler%PCA_bin(bin)%N_spect
       !N_EOF = PCA_handler%PCA_bin(bin)%N_EOF
       info = -999
       info2 = -999

       allocate(PCA_handler%PCA_bin(bin)%C(N_opt, N_opt))
       allocate(PCA_handler%PCA_bin(bin)%EigVal(N_opt))
       allocate(PCA_handler%PCA_bin(bin)%EigVec(N_opt, N_opt))
       allocate(PCA_handler%PCA_bin(bin)%sEOF(N_opt, N_opt))
       allocate(PCA_handler%PCA_bin(bin)%PC(N_spect, N_opt))

       PCA_handler%PCA_bin(bin)%C = 0.0d0
       PCA_handler%PCA_bin(bin)%EigVal = 0.0d0
       PCA_handler%PCA_bin(bin)%EigVec = 0.0d0
       PCA_handler%PCA_bin(bin)%sEOF = 0.0d0
       PCA_handler%PCA_bin(bin)%PC = 0.0d0


       ! Calculate correlation matrix F^T dot T / (N - 1)
       C = matmul( &
            transpose(PCA_handler%PCA_bin(bin)%F_centered), &
            PCA_handler%PCA_bin(bin)%F_centered) / (N_spect - 1)

       PCA_handler%PCA_bin(bin)%C(:,:) = C(:,:)

       ! Perform eigenvalue decomposition on C
       call DSYEVD('V', 'U', N_opt, C, N_opt, &
            tmp_eigval, work, lwork, &
            iwork, liwork, info)

       if (info /= 0) then
          ! Occasionally, DSYEVD will fail (numerical instabilities?)
          ! However, we can still try to perform the eigenvalue decomposition
          ! using DSYEV

          ! Copy back C in case it was modified by DSYEVD
          C(:,:) = PCA_handler%PCA_bin(bin)%C(:,:)
          ! Try again
          call DSYEV('V', 'U', N_opt, C, N_opt, &
               tmp_eigval, work, lwork, info2)

          if (info2 /= 0) then
             call logger%error(fname, "Eigenvalue decomposition failed with both DSYEVD and DSYEV!")
             return
          else
             call logger%debug(fname, "Eigenvalue decomposition with DSYEV was successful!")
          end if

       else
          call logger%debug(fname, "Eigenvalue decomposition with DSYEVD was successful!")
       end if

       ! Store Eigenvalues, but need to reverse the order (descending)!
       PCA_handler%PCA_bin(bin)%EigVal(:) = tmp_eigval(N_opt:1:-1)

       ! Store Eigenvectors, and reverse order as well
       PCA_handler%PCA_bin(bin)%EigVec(:,:) = C(:, N_opt:1:-1)

       ! If some eigenvalue drops below a tiny value, just set them to zero
       ! along with the corresponding eigenvectors and EOFs
       do opt = 1, N_opt
          if (PCA_handler%PCA_bin(bin)%EigVal(opt) < 1E-15) then

             PCA_handler%PCA_bin(bin)%EigVal(opt) = 0.0
             PCA_handler%PCA_bin(bin)%EigVec(:,opt) = 0.0
             PCA_handler%PCA_bin(bin)%sEOF(:,opt) = 0.0

          else
             ! Scaled EOF is eigenvector times sqrt of eigenvalue
             PCA_handler%PCA_bin(bin)%sEOF(:, opt) = PCA_handler%PCA_bin(bin)%EigVec(:, opt) * &
                  sqrt(PCA_handler%PCA_bin(bin)%EigVal(opt))

          end if
       end do

       ! Determine the number of EOFs used via the explained fraction
       do opt = 1, N_opt
          !write(tmp_str, '(A, G0.1, A, G0.1, A, F6.2, A)') &
          !     "Explained variance for bin ", bin, ", #", opt, ": ", &
          !     100.0 * sum(PCA_handler%PCA_bin(bin)%EigVal(1:opt)) / sum(PCA_handler%PCA_bin(bin)%EigVal(:)), &
          !     "%"
          !call logger%debug(fname, trim(tmp_str))

          if (100.0 * sum(PCA_handler%PCA_bin(bin)%EigVal(1:opt)) / &
               sum(PCA_handler%PCA_bin(bin)%EigVal(:)) > 96.0d0) then

             PCA_handler%PCA_bin(bin)%N_EOF = opt

             write(tmp_str, '(A, G0.1, A, G0.1, A), ') "Chosen ", opt, " PCAs for bin # ", bin, "."
             call logger%debug(fname, trim(tmp_str))
             exit

          end if

       end do

       ! Since we now know how many EOFs we need per bin, allocate
       ! the structure to hold the perturbed states
       N_EOF = PCA_handler%PCA_bin(bin)%N_EOF
       allocate(PCA_handler%PCA_bin(bin)%pert_opt_states(-N_EOF:N_EOF, N_opt))
       PCA_handler%PCA_bin(bin)%pert_opt_states = 0.0d0

       ! Calculate and store the optical state vectors for every perturbation
       do eof = -N_EOF, N_EOF
             do opt = 1, N_opt

                if (eof == 0) then
                   PCA_handler%PCA_bin(bin)%pert_opt_states(eof, opt) = &
                        PCA_handler%PCA_bin(bin)%mean_opt(opt)
                else
                   PCA_handler%PCA_bin(bin)%pert_opt_states(eof, opt) = &
                        PCA_handler%PCA_bin(bin)%mean_opt(opt) + &
                        sign(1, eof) * PCA_handler%PCA_bin(bin)%sEOF(opt, abs(eof))
                end if

             end do
       end do

       ! Calculate principal components (PCs)
       do opt = 1, N_opt

          if (PCA_handler%PCA_bin(bin)%EigVal(opt) > 1E-15) then

             PCA_handler%PCA_bin(bin)%PC(:, opt) = matmul(&
                  PCA_handler%PCA_bin(bin)%F_centered, &
                  PCA_handler%PCA_bin(bin)%sEOF(:, opt)) / &
                  PCA_handler%PCA_bin(bin)%EigVal(opt)

          end if

       end do

    end do

    deallocate(C, work, iwork, tmp_eigval)

    success = .true.

  end subroutine perform_PCA

  !> @brief Creates new scene objects for each PCA bin
  !> @param PCA_handler PCA handler object
  !> @param CS_aerosol Control-aerosol object (array)
  !> @param scn Scene object (main scene)
  !> @param PCA_scn Array of new scenes for every bin (bin, eof)
  !>
  !> @detail
  !> This step here is required so that we end up with a scene object
  !> for every bin and eof combination. As the PCA process produces
  !> new optical property profiles for each bin/eof combination, we need
  !> to drive RT calculations with these new profiles. Creating new scene
  !> objects makes this process somewhat managable.
  subroutine create_scenes_from_PCA(PCA_handler, CS_aerosol, scn, PCA_scn)
    type(PCA_handler_t), intent(in) :: PCA_handler
    type(CS_aerosol_t), intent(in) :: CS_aerosol(:)
    type(scene), intent(in) :: scn
    type(scene), allocatable, intent(inout) :: PCA_scn(:, :)

    type(optical_properties) :: op
    integer :: N_EOF, max_EOF
    integer :: i, j, aer
    integer :: bin, opt, wl, eof
    integer :: n_lay
    integer :: idx

    ! Count the total number of binned calculations
    ! N_BIN * (2*N_EOF + 1) holds only when all bins have the same
    ! number of EOFs, but we might change that at some point..

    max_EOF = maxval(PCA_handler%PCA_bin(:)%N_EOF)
    allocate(PCA_scn(PCA_handler%N_bin, -max_EOF:max_EOF))

    ! Go through all scene objects and populate them with
    ! the constants that go unchanged from the original
    ! scene. (pretty much everything apart from the optical props)
    ! Each scene will only consist of a single wavelength point.

    do bin = 1, PCA_handler%N_bin
       N_EOF = PCA_handler%PCA_bin(bin)%N_EOF

       do eof = -N_EOF, N_EOF
          PCA_scn(bin, eof)%num_levels = scn%num_levels
          PCA_scn(bin, eof)%num_active_levels = scn%num_active_levels
          PCA_scn(bin, eof)%num_gases = 1 !scn%num_gases
          PCA_scn(bin, eof)%num_aerosols = scn%num_aerosols
          PCA_scn(bin, eof)%max_pfmom = scn%max_pfmom
          PCA_scn(bin, eof)%num_stokes = scn%num_stokes

          PCA_scn(bin, eof)%atm = scn%atm

          PCA_scn(bin, eof)%date = scn%date
          PCA_scn(bin, eof)%epoch = scn%epoch
          PCA_scn(bin, eof)%lon = scn%lon
          PCA_scn(bin, eof)%lat = scn%lat
          PCA_scn(bin, eof)%SZA = scn%SZA
          PCA_scn(bin, eof)%mu0 = scn%mu0
          PCA_scn(bin, eof)%VZA = scn%VZA
          PCA_scn(bin, eof)%mu = scn%mu
          PCA_scn(bin, eof)%SAA = scn%SAA
          PCA_scn(bin, eof)%VAA = scn%VAA
       end do
    end do

    n_lay = scn%num_active_levels - 1

    ! Now construct the optical properties for each bin, and for each perturbed
    ! state within each bin.
    do bin = 1, PCA_handler%N_bin
       N_EOF = PCA_handler%PCA_bin(bin)%N_EOF
       do eof = -N_EOF, N_EOF

          ! Allocate the optical property arrays inside the scene.
          ! Each scene has only one wavelength
          call allocate_optical_properties(PCA_scn(bin, eof), 1, 1)

          ! Aerosol related objects require a separate allocation
          ! NOTE:
          ! We don't want to make use of the "aerosol_init" function,
          ! as it ends up doing more than we need right now. So just allocate
          ! arrays manually and copy over some values from the "mother" scn object
          allocate(PCA_scn(bin, eof)%op%aer_wl_idx_l(scn%num_aerosols))
          allocate(PCA_scn(bin, eof)%op%aer_wl_idx_r(scn%num_aerosols))
          allocate(PCA_scn(bin, eof)%op%aer_mcs_map(scn%num_aerosols))
          allocate(PCA_scn(bin, eof)%op%aer_ext_q(1, scn%num_aerosols))
          allocate(PCA_scn(bin, eof)%op%aer_sca_q(1, scn%num_aerosols))
          allocate(PCA_scn(bin, eof)%op%aer_ssa(1, scn%num_aerosols))
          allocate(PCA_scn(bin, eof)%op%aer_ext_tau(1, &
               scn%num_levels-1, scn%num_aerosols))
          allocate(PCA_scn(bin, eof)%op%aer_sca_tau(1, &
               scn%num_levels-1, scn%num_aerosols))

          ! Copy over the needed bookkeeping variables

          if (scn%num_aerosols > 0) then
             PCA_scn(bin, eof)%op%reference_aod = scn%op%reference_aod
             PCA_scn(bin, eof)%op%aer_mcs_map = scn%op%aer_mcs_map
             PCA_scn(bin, eof)%op%aer_wl_idx_l = scn%op%aer_wl_idx_l
             PCA_scn(bin, eof)%op%aer_wl_idx_r = scn%op%aer_wl_idx_r
          end if

          ! Reconstruct total optical depth profiles
          PCA_scn(bin, eof)%op%layer_tau(1, 1:n_lay) = &
               exp(PCA_handler%PCA_bin(bin)%pert_opt_states(eof, &
               PCA_handler%s_tau:PCA_handler%e_tau))
          where(PCA_scn(bin, eof)%op%gas_tau < 1d-10) PCA_scn(bin, eof)%op%gas_tau = 1d-10


          ! Reconstruct Rayleigh extinction optical depth
          PCA_scn(bin, eof)%op%ray_tau(1, 1:n_lay) = &
               exp(PCA_handler%PCA_bin(bin)%pert_opt_states(eof, &
               PCA_handler%s_ray:PCA_handler%e_ray))
          where(PCA_scn(bin, eof)%op%ray_tau < 1d-10) PCA_scn(bin, eof)%op%ray_tau = 1d-10

          ! Reconstruct surface albedo
          PCA_scn(bin, eof)%op%albedo(1) = &
               exp(PCA_handler%PCA_bin(bin)%pert_opt_states(eof, &
               PCA_handler%s_sur))

          if (PCA_scn(bin, eof)%op%albedo(1) > 0.999999) then
             PCA_scn(bin, eof)%op%albedo(1) = 0.999999
          end if

          if (PCA_scn(bin, eof)%op%albedo(1) < 1d-6) then
             PCA_scn(bin, eof)%op%albedo(1) = 1d-6
          end if

          PCA_scn(bin, eof)%op%aer_frac(1) = &
               PCA_handler%PCA_bin(bin)%pert_opt_states(eof, &
               PCA_handler%s_fra)

          ! Reconstruct aerosol scattering efficiency, one per aerosol
          do j = 1, PCA_scn(bin, eof)%num_aerosols
             PCA_scn(bin, eof)%op%aer_sca_q(1, j) = &
                  exp(PCA_handler%PCA_bin(bin)%pert_opt_states(eof, &
                  PCA_handler%s_aer + j - 1))
          end do

          ! Rayleigh scattering matrix needs a notion of wavelength
          PCA_scn(bin, eof)%op%wl(1) = (1.0d0 - PCA_scn(bin, eof)%op%aer_frac(1)) * &
               scn%op%wl(1) + PCA_scn(bin, eof)%op%aer_frac(1) * scn%op%wl(size(scn%op%wl))

          ! --------------------------------------------------------
          ! From here on we can derive the other required quantities
          ! --------------------------------------------------------

          ! First, we must reconstruct the aerosol scattering optical depth,
          ! which is done via the aerosol scattering efficiency, that
          ! relates the aerosol scattering optical depth at the band edges
          ! to whatever spectral point we have.

          ! Step 1) Reconstruct the single scattering albedo as well as
          !         the aerosol scattering optical depth.

          ! Single scattering albedo:

          ! Rayleigh component
          PCA_scn(bin, eof)%op%layer_omega(1,1:n_lay) = &
               PCA_scn(bin, eof)%op%ray_tau(1,1:n_lay) / PCA_scn(bin, eof)%op%layer_tau(1,1:n_lay)

          ! Aerosol component(s)
          do aer = 1, scn%num_aerosols
             PCA_scn(bin, eof)%op%aer_sca_tau(1,1:n_lay, aer) = &
                  scn%op%aer_ext_tau_edge(1,1:n_lay, aer) * PCA_scn(bin, eof)%op%aer_sca_q(1,aer) &
                  / CS_aerosol(scn%op%aer_mcs_map(aer))%qext(scn%op%aer_wl_idx_l(aer))

             WHERE(PCA_scn(bin, eof)%op%aer_sca_tau(1,1:n_lay,aer) < 1d-10) &
                  PCA_scn(bin, eof)%op%aer_sca_tau(1,1:n_lay,aer) = 1d-10

             PCA_scn(bin, eof)%op%layer_omega(1,1:n_lay) = &
                  PCA_scn(bin, eof)%op%layer_omega(1,1:n_lay) + &
                  PCA_scn(bin, eof)%op%aer_sca_tau(1,1:n_lay,aer) / &
                  PCA_scn(bin, eof)%op%layer_tau(1,1:n_lay)

          end do

          ! Trim the layer single-scattering albedos back to 0 > ssa > 0.999999
          WHERE(PCA_scn(bin, eof)%op%layer_omega > 0.999999) PCA_scn(bin, eof)%op%layer_omega = 0.999999
          WHERE(PCA_scn(bin, eof)%op%layer_omega < 1d-10) PCA_scn(bin, eof)%op%layer_omega = 1d-10

       end do
    end do

  end subroutine create_scenes_from_PCA

  !> @brief Calculates linearized inputs for the RT models
  !> @param PCA_handler PCA handler object
  !> @param bin Current bin
  !> @param ltau Linearized total optical depth for main scene
  !> @param lomega Linearized total single scattering albedo for main scene
  !> @param lsurf Linearized surface parameter(s) for main scene
  !> @param ltau_pca Linearized total optical depth for PCA
  !> @param lomega_pca Linearized total single scattering albedo for PCA
  !> @param lsurf_pca Linearized surface parameter(s) for PCA
  subroutine extract_linputs_for_PCA(PCA_handler, bin, &
       ltau, lomega, lsurf, &
       ltau_pca, lomega_pca, lsurf_pca)

    type(PCA_handler_t), intent(in) :: PCA_handler
    integer, intent(in) :: bin
    double precision, intent(in) :: ltau(:,:,:)
    double precision, intent(in) :: lomega(:,:,:)
    double precision, intent(in) :: lsurf(:,:)
    double precision, intent(inout) :: ltau_pca(:,:,:)
    double precision, intent(inout) :: lomega_pca(:,:,:)
    double precision, intent(inout) :: lsurf_pca(:,:)

    integer :: i
    integer :: j
    integer :: wl

    ! NOTE
    ! Keep in mind, these *_pca arrays have to conform to the general
    ! structure of linearized inputs, where the first array dimension
    ! corresponds to wavelength. In our PCA-related scene objects, we
    ! do *not* have wavelengths explicitly, hence there will always be
    ! one element to that dimension.

    ltau_pca(1,:,:) = 0.0d0
    lomega_pca(1,:,:) = 0.0d0
    lsurf_pca(1,:) = 0.0d0

    do wl = 1, size(ltau, 1)

       if (PCA_handler%PCA_bin_idx_map(bin) == PCA_handler%assigned_bin(wl)) then
          ltau_pca(1,:,:) = ltau_pca(1,:,:) + ltau(wl,:,:)
          lomega_pca(1,:,:) = lomega_pca(1,:,:) + lomega(wl,:,:)
          lsurf_pca(1,:) = lsurf_pca(1,:) + lsurf(wl,:)
       end if

    end do

    ltau_pca(1,:,:) = ltau_pca(1,:,:) / PCA_handler%PCA_bin(bin)%N_spect
    lomega_pca(1,:,:) = lomega_pca(1,:,:) / PCA_handler%PCA_bin(bin)%N_spect
    lsurf_pca(1,:) = lsurf_pca(1,:) / PCA_handler%PCA_bin(bin)%N_spect

  end subroutine extract_linputs_for_PCA


  !> @brief Maps radiances obtained through PCA back to the monochromatic grid
  !> @param PCA_handler PCA handler object
  !> @param binned_lo Results from the binned low-accuracy calculations
  !> @param binned_hi Results from the binned high-accuracy calculations
  !> @param monochromatic Monochromatic array with radiance results
  subroutine map_PCA_radiances(PCA_handler, binned_lo, binned_hi, monochromatic)

    type(PCA_handler_t), intent(in) :: PCA_handler
    double precision, intent(inout) :: binned_lo(:,-maxval(PCA_handler%PCA_bin(:)%N_EOF):,:)
    double precision, intent(inout) :: binned_hi(:,-maxval(PCA_handler%PCA_bin(:)%N_EOF):,:)
    double precision, intent(inout) :: monochromatic(:,:)

    integer :: bin, eof, wl, q, i
    integer :: n_stokes

    double precision, allocatable :: delta_minus(:), delta_plus(:), delta_zero(:)
    double precision, allocatable :: monochromatic_tmp(:,:)
    double precision :: zero_ratio, plus_ratio, minus_ratio


    n_stokes = size(monochromatic, 2)

    allocate(delta_minus(n_stokes))
    allocate(delta_plus(n_stokes))
    allocate(delta_zero(n_stokes))


    do bin = 1, PCA_handler%N_bin

       allocate(monochromatic_tmp(PCA_handler%PCA_bin(bin)%N_spect, n_stokes))
       monochromatic_tmp(:,:) = 0.0d0

       zero_ratio = binned_hi(bin, 0, 1) / binned_lo(bin, 0, 1)
       delta_zero(1) = log(zero_ratio)

       if (ieee_is_nan(delta_zero(1))) delta_zero(1) = 0.0d0
       monochromatic_tmp(:,1) = delta_zero(1)

       do q = 2, n_stokes
          delta_zero(q) = binned_hi(bin, 0, q) - binned_lo(bin, 0, q)
       end do

       do eof = 1, PCA_handler%PCA_bin(bin)%N_EOF

          plus_ratio = binned_hi(bin, eof, 1) / binned_lo(bin, eof, 1)
          delta_plus(1) = log(plus_ratio)

          minus_ratio = binned_hi(bin, -eof, 1) / binned_lo(bin, -eof, 1)
          delta_minus(1) = log(minus_ratio)

          do q = 2, n_stokes
             delta_plus(q) = binned_hi(bin, eof, q) - binned_lo(bin, eof, q)
             delta_minus(q) = binned_hi(bin, -eof, q) - binned_lo(bin, -eof, q)
          end do

          do q = 1, n_stokes

             if (ieee_is_nan(delta_zero(q)) .or. ieee_is_nan(delta_plus(q)) .or. ieee_is_nan(delta_minus(q))) then
                delta_plus(q) = 0.0
                delta_minus(q) = 0.0
                delta_zero(q) = 0.0
             end if

             do wl = 1, PCA_handler%PCA_bin(bin)%N_spect

                monochromatic_tmp(wl, q) = monochromatic_tmp(wl, q) + &
                     0.5d0 * (delta_plus(q) - delta_minus(q)) * &
                     PCA_handler%PCA_bin(bin)%PC(wl, eof)

                monochromatic_tmp(wl, q) = monochromatic_tmp(wl, q) + &
                     0.5d0 * (delta_plus(q) - 2.0d0 * delta_zero(q) + delta_minus(q)) * &
                     (PCA_handler%PCA_bin(bin)%PC(wl, eof) ** 2)

             end do
          end do

       end do

       ! Exponentiate radiance (intensity only) correction factors
       monochromatic_tmp(:, 1) = exp(monochromatic_tmp(:, 1))


       ! Back-assign and add to low-accuracy result
       i = 1
       do wl = 1, size(monochromatic, 1)

          if (PCA_handler%PCA_bin_idx_map(bin) == PCA_handler%assigned_bin(wl)) then

             monochromatic(wl, 1) = monochromatic(wl, 1) * monochromatic_tmp(i, 1)

             do q = 2, n_stokes
                monochromatic(wl, q) = monochromatic(wl, q) + monochromatic_tmp(i, q)
             end do

             ! Wavelength counter for binned spectral points
             i = i + 1
          end if
       end do



       deallocate(monochromatic_tmp)

    end do


  end subroutine map_PCA_radiances


  !> @brief Maps Jacobians obtained through PCA back to the monochromatic grid
  !> @param PCA_handler PCA handler object
  !> @param wfunctions_lo Results from the binned low-accuracy calculations
  !> @param wfunctions_hi Results from the binned high-accuracy calculations
  !> @param wfunctions_mono Monochromatic array with Jacobian results
  subroutine map_PCA_jacobians(PCA_handler, wfunctions_lo, wfunctions_hi, wfunctions_mono)
    type(PCA_handler_t), intent(in) :: PCA_handler
    double precision, intent(in) :: wfunctions_lo(:,-maxval(PCA_handler%PCA_bin(:)%N_EOF):,:,:)
    double precision, intent(in) :: wfunctions_hi(:,-maxval(PCA_handler%PCA_bin(:)%N_EOF):,:,:)
    double precision, intent(inout) :: wfunctions_mono(:,:,:)

    integer :: n_stokes
    integer :: n_deriv
    integer :: n_wl
    integer :: q, l, i, wl, eof, bin, deriv

    double precision :: delta_minus, delta_plus, delta_zero
    double precision :: delta_o1, delta_o2
    double precision, allocatable :: corr_wfunctions(:,:,:)


    n_wl = size(wfunctions_mono, 1)
    n_stokes = size(wfunctions_mono, 2)
    n_deriv = size(wfunctions_mono, 3)

    allocate(corr_wfunctions, mold=wfunctions_mono)
    corr_wfunctions = 0.0d0

    ! Jacobians will be corrected just like non-radiance Stokes components,
    ! i.e. via differences and not through a product.

    do bin = 1, PCA_handler%N_bin
       do deriv = 1, n_deriv
          do q = 1, n_stokes

             delta_zero = wfunctions_hi(bin, 0, q, deriv) - wfunctions_lo(bin, 0, q, deriv)

             do eof = 0, 0 !PCA_handler%PCA_bin(bin)%N_EOF

                delta_plus = wfunctions_hi(bin, eof, q, deriv) - wfunctions_lo(bin, eof, q, deriv)
                delta_minus = wfunctions_hi(bin, -eof, q, deriv) - wfunctions_lo(bin, -eof, q, deriv)

                delta_o1 = 0.5d0 * (delta_plus - delta_minus)
                delta_o2 = 0.5d0 * (delta_plus + delta_minus - 2 * delta_zero)

                i = 1
                do wl = 1, n_wl
                   if (PCA_handler%PCA_bin_idx_map(bin) == PCA_handler%assigned_bin(wl)) then

                      if (eof == 0) then
                         corr_wfunctions(wl,q,deriv) = delta_zero
                      end if

                      if (eof > 0) then
                         corr_wfunctions(wl,q,deriv) = corr_wfunctions(wl,q,deriv) + &
                              delta_o1 * PCA_handler%PCA_bin(bin)%PC(i, eof) + &
                              delta_o2 * PCA_handler%PCA_bin(bin)%PC(i, eof)**2
                      end if

                      i = i + 1
                   end if
                end do ! wavelength loop

             end do ! EOF loop


          end do ! Stokes component loop
       end do ! derivative loop
    end do ! bin loop

    wfunctions_mono(:,:,:) = wfunctions_mono(:,:,:) + corr_wfunctions(:,:,:)

  end subroutine map_PCA_jacobians



end module PCART_mod
