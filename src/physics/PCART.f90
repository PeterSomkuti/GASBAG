!> @brief Module that contains functions for the PCA-based RT acceleration
!> @file PCART.f90
!> @author Peter Somkuti

module PCART_mod

  ! User modules
  use scene_mod
  use math_utils_mod

  implicit none

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

  !> @brief User-type that contains data for PCA-based fast RT
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
    ! + 1 for aerosol interpolation coefficient

    PCA_handler%N_opt = 2 * (scn%num_active_levels - 1)

    ! Add the number of aerosols plus one for the interpolation coefficient
    ! (regardless of number of retrieved aerosol parameters)
    if (allocated(scn%op%reference_aod)) then
       PCA_handler%N_opt = PCA_handler%N_opt + size(scn%op%reference_aod)
       PCA_handler%N_opt = PCA_handler%N_opt + 1
    end if

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
    PCA_handler%e_ray = PCA_handler%s_ray + scn%num_active_levels ! + nlay - 1

    ! Surface parameteres (so far only Lambertian surface here)
    PCA_handler%s_sur = PCA_handler%e_ray + 1
    PCA_handler%e_sur = PCA_handler%s_sur

    ! Aerosol scattering ratio/coefficient
    PCA_handler%s_aer = PCA_handler%e_sur + 1
    PCA_handler%e_aer =  PCA_handler%s_aer + size(scn%op%reference_aod) - 1

    ! Aerosol interpolation coefficient (one one needed)
    PCA_handler%s_fra = PCA_handler%e_aer + 1
    PCA_handler%e_fra = PCA_handler%s_fra

  end subroutine initialize_PCA

  !> @brief Binning method according to Kopparla et al.
  subroutine assign_pushkar_bins(scn, PCA_handler)
    type(scene), intent(in) :: scn
    type(PCA_handler_t), intent(inout) :: PCA_handler

    ! Manually set bin boundaries in total tau_gas space
    ! huge(0.0d0) is largest double precision number
    double precision, parameter :: bounds(12) = &
         [0d0, 0.01d0, 0.025d0, 0.05d0, &
          0.1d0, 0.25d0, 0.5d0, 0.625d0, &
          0.75d0, 1.0d0, 5.0d0, huge(0.0d0)]

    integer :: nspect, i, j

    integer, allocatable :: single_assigned_bins(:)

    double precision, allocatable :: ssa_tot(:)
    double precision :: ssa_median

    character(len=*), parameter :: fname = "assign_pushkar_bins"
    character(len=99) :: tmp_str

    nspect = size(PCA_handler%assigned_bin)

    allocate(ssa_tot(nspect))

    allocate(single_assigned_bins(nspect))
    single_assigned_bins(:) = 0

    ! Sum up the gas opds to get the total taug
    ! taug_tot(:) = sum(taug, 1)
    ssa_tot(:) = sum(scn%op%layer_omega, 2)

    ! Assign each spectral point to a bin according to taug
    ! is there a neat way to do an equivalent to np.searchsorted
    do i = 1, nspect

       do j = 1, size(bounds) - 1
          if ( &
               (scn%op%total_tau(i) .ge. bounds(j)) .and. &
               (scn%op%total_tau(i) .lt. bounds(j+1))) then
             single_assigned_bins(i) = j
          endif
       enddo
    enddo

    ! Now split each bin along median SSA
    ! Note that the total possible bins now is 2*size(bounds) - 2

    j = 1
    do i = 1, size(bounds) - 1

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
    write(tmp_str, "(A, G0.1)") "Unassigned points: ", count(PCA_handler%assigned_bin(:) == -1)
    call logger%debug(fname, trim(tmp_str))

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

    end do
  end subroutine allocate_first_PCA

  !> @brief Fills F-matrix with values
  subroutine create_F_matrix(scn, PCA_handler)
    type(scene), intent(in) :: scn
    type(PCA_handler_t), intent(inout) :: PCA_handler

    integer :: i, j, k

    integer :: s_tau
    integer :: e_tau
    integer :: s_ray
    integer :: e_ray
    integer :: s_sur
    integer :: e_sur

    s_tau = PCA_handler%s_tau
    e_tau = PCA_handler%e_tau
    s_ray = PCA_handler%s_ray
    e_ray = PCA_handler%e_ray
    s_sur = PCA_handler%s_sur
    e_sur = PCA_handler%e_sur

    do j = 1, PCA_handler%N_bin
       k = 1
       do i = 1, size(scn%op%wl)

          if (PCA_handler%assigned_bin(i) == PCA_handler%PCA_bin_idx_map(j)) then
             ! spectral point "i" belongs to bin "PCA_bin_idx_map(j)"

             ! Copies layer-resolved total optical depth
             PCA_handler%PCA_bin(j)%F(k, s_tau:e_tau) = scn%op%layer_tau(i, :)
             ! Copies layer-resolved Rayleigh optical depth
             PCA_handler%PCA_bin(j)%F(k, s_ray:e_ray) = scn%op%ray_tau(i, :)
             ! Copies surface parameter
             PCA_handler%PCA_bin(j)%F(k, s_sur:e_sur) = scn%op%albedo(i)

             ! Aerosols present?
             if (allocated(scn%op%reference_aod)) then

             end if

          end if

       end do
    end do




  end subroutine create_F_matrix




end module PCART_mod
