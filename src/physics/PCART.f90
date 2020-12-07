!> @brief Module that contains functions for the PCA-based RT acceleration
!> @file PCART.f90
!> @author Peter Somkuti

module PCART_mod

  ! User modules
  use scene_mod
  use math_utils_mod

  implicit none

  type PCA_bin_t
     !> Is this bin used?
     logical :: used
     !> Number of spectral points in this bin
     integer :: nspect
     !> Number of used EOFs in this bin
     integer :: neof
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
     !> Mean optical properties
     double precision, allocatable :: mean_opt(:)
     !> Perturbed optical states
     double precision, allocatable :: pert_opt_states(:,:)
  end type PCA_bin_t

  !> @brief User-type that contains data for PCA-based fast RT
  type PCA_handler_t
     !> Assignment of spectral index to PCA bin
     integer, allocatable :: assigned_bin(:)
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

    double precision, allocatable :: ssa_tot(:), ssa_tot_sorted(:)
    double precision :: ssa_median

    character(len=*), parameter :: fname = "assign_pushkar_bins"
    character(len=99) :: tmp_str

    nspect = size(PCA_handler%assigned_bin)

    allocate(ssa_tot(nspect))
    allocate(ssa_tot_sorted(nspect))

    allocate(single_assigned_bins(nspect))
    single_assigned_bins(:) = 0

    ! Sum up the gas opds to get the total taug
    ! taug_tot(:) = sum(taug, 1)
    ssa_tot(:) = sum(scn%op%layer_omega, 1)
    ssa_tot_sorted(:) = ssa_tot(:)

    ! Sort ssa_tot
    call combsort(ssa_tot_sorted)

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

    j = 1
    do i = 1, size(bounds) - 1

       ! calculate the median of total SSA
       ssa_median = percentile(pack(ssa_tot, single_assigned_bins == i), 50.0d0)

       ! Find all scenes with a total SSA <= than median and assign it to bin "j"
       where ((single_assigned_bins == i) .and. (ssa_tot <= ssa_median))
          PCA_handler%assigned_bin = j
       endwhere
       j = j + 1

       ! Find all scenes with a total SSA > than median and assign it to bin "j+1"
       where ((single_assigned_bins == i) .and. (ssa_tot > ssa_median))
          PCA_handler%assigned_bin = j
       endwhere
       j = j + 1

    enddo

    ! Debug: print the number of points in each bin
    do j = 1, size(bounds) - 1
       write(tmp_str, "(A, G3.1, A, E10.5, A, E10.5, A, G7.1, A)") &
            "Bin ", j, " [", bounds(j), " - ", bounds(j+1), "]: ", &
            count(PCA_handler%assigned_bin(:) == j), " points"
       call logger%debug(fname, trim(tmp_str))
    enddo
    ! Also let the user know if there are any unassigned spectral points
    write(tmp_str, "(A, G0.1)") "Unassigned points: ", count(PCA_handler%assigned_bin(:) == -1)
    call logger%debug(fname, trim(tmp_str))

  end subroutine assign_pushkar_bins

  subroutine allocate_PCA(PCA_handler)
    type(PCA_handler_t), intent(inout) :: PCA_handler


  end subroutine allocate_PCA




end module PCART_mod
