!> @brief Gas optical depth calculations
!> @file gas_tau.f90
!> @author Peter Somkuti
!>
!> @details This module provides functions to calculate the gas optical depths for a given
!> model atmosphere for one specific gas. There are essentially two functions:
!> one to extract the cross section value from a given (wavelength,H2O,p,T)-coordinate
!> via linear interpolation (get_CS_value_at). The other function (calculate_gas_tau)
!> does the physics to calculate the layer and sub-layer ODs, and calls the other
!> to grab the cross section values needed for the calculations. Some notation in this
!> code is more explicit than in other parts of the program, in hopes that the compiler
!> can optimize this performance-critical section of the code a bit better.

module gas_tau_mod

  ! User modules
  use control_mod, only: CS_gas_t
  use math_utils_mod

  public :: calculate_gas_tau, get_CS_value_at

contains


  !> @brief Calculate the per-layer and gas optical depths
  !> @param pre_gridded Do we use (wavelength) pre-gridded spectroscopy?
  !> @param is_H2O Is this gas water vapour?
  !> @param wl Wavelength array
  !> @param gas_vmr Gas VMR array
  !> @param psurf Surface pressure
  !> @param p Pressure levels array
  !> @param T Temperature levels array
  !> @param sh Specific humidity levels array
  !> @param gas MCS%gas object
  !> @param N_sub Number of sublayers for the sublayer integration
  !> @param need_psurf_jac Do we want to return surface pressure jacobians?
  !> @param gas_tau Per-layer gas optical depths
  !> @param gas_tau_dpsurf Jacobian: dtau/dpsurf
  !> @param gas_tau_dvmr Jacobain: dtau/dvmr (level)
  !> @param ndry Number of dry air molecules
  !> @param success Was the calculation successful?
  !>
  !> @details
  !> "calculate_gas_tau" uses Gauss-Kronrod integration to calculate per-layer
  !> gas optical depths for a given gas. This is probably quite similar to the
  !> OCO full-physics way of doing it, but somewhat optimized for speed, since
  !> all variables are loaded into memory. Using pre-gridded spectroscopy makes
  !> this a bit faster than having to interpolate into the ABSCO wavelength grid.
  subroutine calculate_gas_tau(pre_gridded, &
       is_H2O, &
       N_levels, &
       N_active_levels, &
       N_wl, &
       wl, &
       gas_vmr, &
       psurf, &
       p, &
       T, &
       sh, &
       grav, &
       gas, &
       N_sub, &
       need_psurf_jac, &
       gas_tau, &
       gas_tau_dpsurf, &
       gas_tau_dvmr, &
       success)

    implicit none
    logical, intent(in) :: pre_gridded
    logical, intent(in) :: is_H2O
    integer, intent(in) :: N_levels
    integer, intent(in) :: N_active_levels
    integer, intent(in) :: N_wl
    double precision, intent(in) :: wl(N_wl)
    double precision, intent(in) :: gas_vmr(N_levels)
    double precision, intent(in) :: psurf
    double precision, intent(in) :: p(N_levels)
    double precision, intent(in) :: T(N_levels)
    double precision, intent(in) :: sh(N_levels)
    double precision, intent(in) :: grav(N_levels)
    type(CS_gas_t), intent(in) :: gas
    integer, intent(in) :: N_sub
    logical, intent(in) :: need_psurf_jac

    double precision, intent(inout) :: gas_tau(N_wl, N_levels-1)
    double precision, intent(inout) :: gas_tau_dpsurf(N_wl, N_levels-1)
    double precision, intent(inout) :: gas_tau_dvmr(2, N_wl, N_levels-1)

    logical, intent(inout) :: success


    ! character(len=999) :: tmp_str ! Temporary string
    character(len=*), parameter :: fname = "calculate_gas_tau" ! Function name
    logical :: log_scaling ! Are we interpolating in log-p space? (rather than linear p space)
    ! Placeholders to store the lower and higher (in altitude, i.e. layer boundaries)
    ! values of p,T,sh, etc.
    double precision :: p_lower, p_higher
    double precision :: T_lower, T_higher
    double precision :: sh_lower, sh_higher
    double precision :: VMR_lower, VMR_higher
    double precision :: grav_lower, grav_higher
    ! For H2O, we need a different factor to calculate the ODs (no SH-adjustment)
    double precision :: H2O_corr
    ! Values for this particular layer, and perturbations for psurf Jacobian
    double precision :: this_H2O, this_H2O_pert
    double precision :: this_p, p_fac, this_p_fac
    double precision :: this_T, this_sh, this_VMR, this_M
    double precision :: this_p_pert, p_fac_pert, this_p_fac_pert, p_lower_pert
    double precision :: this_grav, this_grav_pert
    double precision :: T_lower_pert, this_T_pert
    double precision :: VMR_lower_pert, this_VMR_pert
    double precision :: sh_lower_pert, this_sh_pert
    double precision :: this_M_pert
    double precision :: C_tmp
    double precision :: ndry
    double precision :: gas_tmp(N_wl)
    double precision :: this_CS_value(N_wl)
    integer :: wl_left_indices(N_wl)

    ! Variables related to the Gauss-Kronrod weights
    double precision, allocatable :: GK_abscissae(:), GK_weights(:), G_weights(:)
    double precision, allocatable :: GK_abscissae_f(:), GK_weights_f(:), G_weights_f(:)
    double precision :: GK_abscissae_f_pert(N_sub), GK_weights_f_pert(N_sub), G_weights_f_pert(N_sub)

    integer :: j,k,l


    ! Notes to Gauss-Kronrod integration:
    ! According to GK-rules, the knots ("x-coords") and weights are symmetric
    ! around 0. So the GK algorithm used here (third_party/kronrod/kronrod.f90)
    ! will calculate N+1 knots and weights, where N is the number you supply to
    ! the function "kronrod". To perform the integration, we then just mirror those
    ! points to eventually end up with 2N+1 knots and weights. The user here supplies
    ! "N_sub", the number of sublayers used for the integration. We don't want the
    ! user to have to think too much, we convert these numbers accordingly.
    ! If the user says: Nsub=5, then we supply N=(N_sub+1)/2=3 to the function,
    ! which gives us 3 GK knots (including 0), which is then later on mirrored
    ! to obtain back our requested 5.

    ! These three are the placeholders for GK knots (abscissae) and weights
    ! that come out of the kronrod function.
    allocate(GK_abscissae((N_sub+1)/2))
    allocate(GK_weights((N_sub+1)/2))
    allocate(G_weights((N_sub+1)/2))

    ! These are the fully-back-mirrored GK knots and weights
    allocate(GK_abscissae_f(N_sub))
    allocate(GK_weights_f(N_sub))
    allocate(G_weights_f(N_sub))

    ! Maybe some time later we will let the user decide if log-scaling should be
    ! used or not, but I don't see a reason why you would ever NOT want to use this.
    log_scaling = .true.

    ! Just to make sure there's nothing in there already..
    gas_tau(:,:) = 0.0d0
    gas_tau_dpsurf(:,:) = 0.0d0
    gas_tau_dvmr(:,:,:) = 0.0d0

    ! Pre-compute the indicies of the ABSCO wavelength dimension at which
    ! we can find the wavelengths supplied to this function and at which we
    ! want to calculate the optical depths.
    ! If the spectroscopy comes pre-gridded already, this step is not
    ! needed.

    if (.not. pre_gridded) then
       ! A binary search seems to be similarly quick as the odd thing we're doing above..
       do j=1, size(wl)
          wl_left_indices(j) = searchsorted_dp(gas%wavelength(:), wl(j), .true.)
       end do
    else
       if (size(gas%wavelength) /= size(wl)) then
          ! This is a sanity check - the pre-gridded spectroscopy needs to have
          ! the same size as the high-resolution wavelength grid, obviously.
          call logger%fatal(fname, "Spectroscopy wavelength size not equal to hires grid!")
          stop 1
       end if
    end if

    ! Given the number of sublayers, we precompute the Gauss-Kronrod weights for
    ! the integral approximation.
    call kronrod((N_sub-1)/2, 1d-6, GK_abscissae, GK_weights, G_weights)

    ! Traverse the atmosphere layers(!), starting from the bottom to the top
    do l=N_active_levels, 2, -1

       ! First, grab the lower (altitude-wise) and higher values for
       ! p,T,sh from the atmosphere profiles for whatever layer l we're in.

       if (l == N_active_levels) then
          ! BOA layer - psurf should be between lowermost level and
          ! the level above (ideally). But this will extrapolate linearly
          ! anyway, which might cause some trouble.
          p_lower = psurf

          if (log_scaling) then
             p_fac = (log(p_lower) - log(p(l-1))) / (log(p(l)) - log(p(l-1)))
          else
             p_fac = (p_lower - p(l-1)) / (p(l) - p(l-1))
          end if

          T_lower = (1.0d0 - p_fac) * T(l-1) + p_fac * T(l)
          sh_lower = (1.0d0 - p_fac) * sh(l-1) + p_fac * sh(l)
          VMR_lower = (1.0d0 - p_fac) * gas_vmr(l-1) + p_fac * gas_vmr(l)
          grav_lower = (1.0d0 - p_fac) * grav(l-1) + p_fac * grav(l)

          if (need_psurf_jac) then
             ! Calculate perturbed properties when psurf jacobian is needed.
             ! We perturb psurf towards lower values, making it essentially a
             ! backwards-differentiation, so we don't accidentally step outside
             ! the lower atmsophere bound. It's then 'reversed' via a minus
             ! sign when calculating the proper forward Jacobians.
             p_lower_pert = psurf - PSURF_PERTURB

             if (log_scaling) then
                p_fac_pert = (log(p_lower_pert) - log(p(l-1))) / (log(p(l)) - log(p(l-1)))
             else
                p_fac_pert = (p_lower_pert - p(l-1)) / (p(l) - p(l-1))
             end if

             T_lower_pert = (1.0d0 - p_fac_pert) * T(l-1) + p_fac_pert * T(l)
             sh_lower_pert = (1.0d0 - p_fac_pert) * sh(l-1) + p_fac_pert * sh(l)
             VMR_lower_pert = (1.0d0 - p_fac_pert) * gas_vmr(l-1) + p_fac_pert * gas_vmr(l)
          end if
       else
          p_lower = p(l)
          T_lower = T(l)
          sh_lower = sh(l)
          VMR_lower = gas_vmr(l)
          grav_lower = grav(l)
       end if

       p_higher = p(l-1)
       T_higher = T(l-1)
       sh_higher = sh(l-1)
       VMR_higher = gas_vmr(l-1)
       grav_higher = grav(l-1)

       ! Should the perturbed surface pressure actually fall onto the next-higher level,
       ! we need to make a small adjustment, otherwise, we end up dividing by zero later on.
       if (need_psurf_jac) then
          if (abs(p_lower_pert - p_higher) < 1d-3) then
             p_lower_pert = p_lower_pert + PSURF_PERTURB/10.0d0
          end if
       end if


       ! -------------------------
       !
       ! Gauss - Kronrod portion
       !
       ! -------------------------

       ! Map the GK abscissae and weights from [0,1] to symmetric [-1,1] by mirroring. This section
       ! seems a bit excessive - maybe a shorter way of doing it?

       ! For odd number of sublayers
       GK_abscissae_f((N_sub+1)/2) = GK_abscissae(size(GK_abscissae))
       GK_weights_f((N_sub+1)/2) = GK_weights(size(GK_weights))
       G_weights_f((N_sub+1)/2) = G_weights(size(G_weights))

       do k=1, (N_sub-1)/2
          GK_abscissae_f(k) = -GK_abscissae(k)
          GK_weights_f(k) = GK_weights(k)
          G_weights_f(k) = G_weights(k)

          GK_abscissae_f(k+(N_sub+1)/2) = GK_abscissae(size(GK_abscissae)-k)
          GK_weights_f(k+(N_sub+1)/2) = GK_weights(size(GK_weights)-k)
          G_weights_f(k+(N_sub+1)/2) = G_weights(size(GK_weights)-k)
       end do

       ! And adjust the GK abscissae and weights to our pressure interval
       ! between p_lower and p_higher. This way, we don't have to re-scale the
       ! result after integration and can obtain it simply by multiplying with
       ! the GK weights. NOTE: kronrod_adjust requires the order or the rule as
       ! the third argument, but this is N_sub-1 (rather than just N_sub)!

       if (need_psurf_jac) then
          ! Since we altered the integration limit, we also need to re-scale the
          ! GK weights accordingly for the Jacobian
          GK_abscissae_f_pert = GK_abscissae_f
          GK_weights_f_pert = GK_weights_f
          G_weights_f_pert = G_weights_f

          if (log_scaling) then
             call kronrod_adjust(log(p_higher), log(p_lower_pert), N_sub-1, &
                  GK_abscissae_f_pert, GK_weights_f_pert, G_weights_f_pert)
          else
             call kronrod_adjust(p_higher, p_lower_pert, N_sub-1, &
                  GK_abscissae_f_pert, GK_weights_f_pert, G_weights_f_pert)
          end if
       end if

       if (log_scaling) then
          call kronrod_adjust(log(p_higher), log(p_lower), N_sub-1, &
               GK_abscissae_f, GK_weights_f, G_weights_f)
       else
          call kronrod_adjust(p_higher, p_lower, N_sub-1, &
               GK_abscissae_f, GK_weights_f, G_weights_f)
       end if


       ! -------------------------
       !
       ! Sub-layer loop
       !
       ! -------------------------
       !do k = 1, N_sub
       do k = 0, N_sub - 1

          ! Get the value of pressure according to the Gauss-Kronrod abscissa
          ! value, at which we evaluate the cross section.
          !this_p = GK_abscissae_f(k)
          if (log_scaling) then
             this_p = log(p_lower) + k * (log(p_higher) - log(p_lower)) / N_sub
          else
             this_p = p_lower + k * (p_higher - p_lower) / N_sub
          end if

          ! And get corresponding values for T,sh and the VMR at this (scaled)
          ! pressure value, which is obtained via simple linear interpolation.
          if (log_scaling) then
             this_p_fac = (this_p - log(p_higher)) / (log(p_lower) - log(p_higher))
             this_p = exp(this_p)
          else
             this_p_fac = (this_p - p_higher) / (p_lower - p_higher)
          end if

          this_T = (1.0d0 - this_p_fac) * T_lower + this_p_fac * T_higher
          this_sh = (1.0d0 - this_p_fac) * sh_lower + this_p_fac * sh_higher
          this_VMR = (1.0d0 - this_p_fac) * VMR_lower + this_p_fac * VMR_higher
          this_grav = (1.0d0 - this_p_fac) * grav_lower + this_p_fac * grav_higher
          this_M = 1.0d3 * DRY_AIR_MASS

          ! Gas CS routine works in H2O VMR rather than SH
          ! this_H2O = this_sh / (1.0d0 - this_sh) * SH_H2O_CONV
          this_H2O = this_sh * DRY_AIR_MASS / (H2Om + this_sh * (DRY_AIR_MASS - H2Om))

          if (is_H2O) then
             H2O_corr = 1.0d0
          else
             H2O_corr = (1.0d0 - this_sh)
          end if

          this_CS_value = get_CS_value_at( &
               pre_gridded, &
               gas, &
               N_wl, &
               wl(1:N_wl), &
               this_p, &
               this_T, &
               this_H2O, &
               wl_left_indices(1:N_wl) &
               )

          ! Constant value ~ almost the same as Ndry but without GK weights
          ! and potential factor from log-scale (d(exp(p)) = p * dp)
          C_tmp = 1.0d0 / this_grav * NA * 0.1d0 * H2O_corr / this_M

          ! Tau for this sublayer
          if (log_scaling) then
             !gas_tmp(:) = GK_weights_f(k) * this_CS_value(:) * this_VMR * C_tmp * this_p
             !ndry = GK_weights_f(k) * C_tmp * this_p
             ndry = C_tmp * ( &
                  exp(log(p_lower) + k * (log(p_higher) - log(p_lower)) / N_sub) &
                  - exp(log(p_lower) + (k + 1) * (log(p_higher) - log(p_lower)) / N_sub) &
                  )
             ! The expression in the brackets is delta pressure in log-spaced intervals

          else
             !gas_tmp(:) = GK_weights_f(k) * this_CS_value(:) * this_VMR * C_tmp
             !ndry = GK_weights_f(k) * C_tmp
             ndry = C_tmp * abs((p_higher - p_lower) / N_sub)
          end if

          ! ----------------------------------------------------------------------------------------------
          ! Add sublayer contribution to full layer tau
          ! ----------------------------------------------------------------------------------------------
          gas_tau(1:N_wl,l-1) = gas_tau(1:N_wl,l-1) + (ndry * this_CS_value(1:N_wl) * this_VMR)
          ! Jacobians of tau w.r.t. changes in level VMR. Layer gas tau can change in two
          ! ways: either change the VMR at the level above (index 1, higher), or below (index 2, lower)

          ! "Higher" (i.e. closer to TOA)
          gas_tau_dvmr(1,1:N_wl,l-1) = gas_tau_dvmr(1,1:N_wl,l-1) + (ndry * this_CS_value(1:N_wl) * (this_p_fac))
          ! "Lower" (i.e. closer to surface)
          gas_tau_dvmr(2,1:N_wl,l-1) = gas_tau_dvmr(2,1:N_wl,l-1) + (ndry * this_CS_value(1:N_wl) * (1.0d0 - this_p_fac))
          ! ----------------------------------------------------------------------------------------------

          ! The same is required in the case of surface pressure jacobians,
          ! but we obviously only do this for the BOA layer

          if (need_psurf_jac .and. (l == N_active_levels)) then

             this_p_pert = GK_abscissae_f_pert(k)

             if (log_scaling) then
                this_p_fac_pert = (this_p_pert - log(p_higher)) / (log(p_lower_pert) - log(p_higher))
             else
                this_p_fac_pert = ((this_p_pert - p_higher) / (p_lower_pert - p_higher))
             end if

             this_T_pert = (1.0d0 - this_p_fac_pert) * T_lower_pert + this_p_fac_pert * T_higher
             this_sh_pert = (1.0d0 - this_p_fac_pert) * sh_lower_pert + this_p_fac_pert * sh_higher
             this_VMR_pert = (1.0d0 - this_p_fac_pert) * VMR_lower_pert + this_p_fac_pert * VMR_higher
             this_grav_pert = (1.0d0 - this_p_fac_pert) * grav_lower + this_p_fac_pert * grav_higher
             this_M_pert = 1.0d3 * DRY_AIR_MASS

             !this_H2O_pert = this_sh_pert / (1.0d0 - this_sh_pert) * SH_H2O_CONV
             this_H2O_pert = this_sh_pert * DRY_AIR_MASS / (H2Om + this_sh_pert * (DRY_AIR_MASS - H2Om))

             if (log_scaling) this_p_pert = exp(this_p_pert)

             this_CS_value = get_CS_value_at( &
                  pre_gridded, &
                  gas, &
                  N_wl, &
                  wl(1:N_wl), &
                  this_p_pert, &
                  this_T_pert, &
                  this_H2O_pert, &
                  wl_left_indices(1:N_wl) &
                  )

             if (is_H2O) then
                H2O_corr = 1.0d0
             else
                H2O_corr = (1.0d0 - this_sh_pert)
             end if

             ! This calculates the gas OD, as a result of a perturbed surface pressure

             C_tmp = 1.0d0 / this_grav_pert * NA * 0.1d0 * H2O_corr / this_M_pert

             if (log_scaling) then
                ! this_p_pert is already exponentiated here!!
                !gas_tmp(1:N_wl) = GK_weights_f_pert(k) * C_tmp * this_p_pert &
                !     * this_CS_value(1:N_wl) * this_VMR_pert
                gas_tmp(1:N_wl) = C_tmp * this_CS_value(1:N_wl) * this_VMR_pert &
                     * ( &
                     exp(log(p_lower) + k * (log(p_higher) - log(p_lower)) / N_sub) &
                     - exp(log(p_lower) + (k + 1) * (log(p_higher) - log(p_lower)) / N_sub) &
                     )
             else
                !gas_tmp(1:N_wl) = GK_weights_f_pert(k) * C_tmp &
                !     * this_CS_value(1:N_wl) * this_VMR_pert
                gas_tmp(1:N_wl) = C_tmp * this_CS_value(1:N_wl) * this_VMR_pert &
                     * abs((p_higher - p_lower) / N_sub)

             end if

             gas_tau_dpsurf(1:N_wl,l-1) = gas_tau_dpsurf(1:N_wl,l-1) + gas_tmp(1:N_wl)

          end if
       end do

       if (need_psurf_jac .and. (l == N_active_levels)) then
          ! Get the difference: tau(psurf - psurb_perturb) - tau(psurf)
          gas_tau_dpsurf(1:N_wl,l-1) = (gas_tau_dpsurf(1:N_wl,l-1) - gas_tau(1:N_wl,l-1)) / PSURF_PERTURB
       end if

       if (need_psurf_jac .and. (l < N_active_levels)) then
          ! Copy the non-BOA layer ODs to the gas_tau_dpsurf array, as the
          ! surface pressure Jacobian merely affects the BOA layer ODs
          gas_tau_dpsurf(1:N_wl,l-1) = 0.0d0
       end if
    end do

    ! When all completed nicely, mark it as successful
    success = .true.

  end subroutine calculate_gas_tau


  !> @brief Grab cross section value at given wl, p, T, H2O
  !> @param pre_gridded Is the spectroscopy pre-gridded already?
  !> @param gas CS_gas object
  !> @param wl Wavelength array [um]
  !> @param p Pressure [Pa]
  !> @param T Temperature [K]
  !> @param H2O H2O VMR
  !> @param wl_left_idx Positions of wl in gas%wl array
  !> @param CS_value Cross section values for wavelengths wl
  pure function get_CS_value_at(pre_gridded, gas, N_wl, wl, p, T, H2O, wl_left_idx) result(CS_value)

    implicit none
    logical, intent(in) :: pre_gridded
    type(CS_gas_t), intent(in) :: gas
    integer, intent(in) :: N_wl
    double precision, intent(in) :: wl(N_wl)
    double precision, intent(in) :: p, T, H2O
    integer, intent(in) :: wl_left_idx(N_wl)

    double precision :: CS_value(N_wl)

    ! character(len=*), parameter :: fname = "get_CS_value_at" ! Function name

    ! These here are the indices between which the CS array will
    ! be linearly interpolated: e.g. idx_(left, right)_(pressure)
    integer :: idx_l_p, idx_r_p ! Pressure
    integer :: idx_lr_T, idx_ll_T, idx_rl_T, idx_rr_T ! Temperature has two dimensions
    integer :: idx_l_H2O, idx_r_H2O
    integer :: idx_l_wl, idx_r_wl
    double precision :: C3(0:1, 0:1, 0:1) ! 0 is 'left', 1 is 'right'
    double precision :: C2(0:1, 0:1), C1(0:1)
    double precision :: wl_d, p_d, T_d_l, T_d_r, H2O_d
    integer :: i


    ! This next section here "just" finds the left and right indices of p,T,H2O,wl values within the
    ! grids of the cross section tables. Unfortunately, sourcing this bit out into a PURE FUNCTION still
    ! comes with such a large overhead that I had to inline them explicitly. Maybe there's a way around
    ! it, but for now the performance increase is totally worth it!

    ! Get the pressure indices
    if (p <= gas%p(1)) then
       idx_l_p = 1
    else if (p >= gas%p(size(gas%p) - 1)) then
       idx_l_p = size(gas%p) - 1
    else
       do i=1, size(gas%p)-1
          if ((p >= gas%p(i)) .and. (p < gas%p(i+1))) then
             idx_l_p = i
             exit
          end if
       end do
    end if
!!$   idx_l_p = searchsorted_dp(gas%p(:), p)
    idx_r_p = idx_l_p + 1

    if ((idx_l_p < 1) .or. (idx_l_p > size(gas%p) - 1)) then
    !   call logger%error(fname, "idx_l_p out of range.")
       return
    end if

    ! Get the temperature indices, which depend on P - so we have two
    ! indices for the temperature dimension. One set of T indices (idx_l(l,r)_T)
    ! corresponds to the left index of pressure (idx_l_p), the other set
    ! (idx_r(l,r)_T) corresponds to the right index of pressure (idx_r_p).

    if (T <= gas%T(1, idx_l_p)) then
       idx_ll_T = 1
    else if (T >= gas%T(size(gas%T, 1) - 1, idx_l_p)) then
       idx_ll_T = size(gas%T, 1) - 1
    else
       do i=1, size(gas%T, 1) - 1
          if ((T > gas%T(i, idx_l_p)) .and. (T <= gas%T(i+1, idx_l_p))) then
             idx_ll_T = i
             exit
          end if
       end do
    end if
!!$   idx_ll_T = searchsorted_dp(gas%T(:, idx_l_p), T)
    idx_lr_T = idx_ll_T + 1

    if ((idx_ll_T < 1) .or. (idx_ll_T > size(gas%T, 1) - 1)) then
    !   call logger%error(fname, "idx_ll_T out of range.")
       return
    end if

    if (T <= gas%T(1, idx_r_p)) then
       idx_rl_T = 1
    else if (T >= gas%T(size(gas%T, 1) - 1, idx_r_p)) then
       idx_rl_T = size(gas%T, 1) - 1
    else
       do i=1, size(gas%T, 1) - 1
          if ((T >= gas%T(i, idx_r_p)) .and. (T < gas%T(i+1, idx_r_p))) then
             idx_rl_T = i
             exit
          end if
       end do
    end if
!!$   idx_rl_T = searchsorted_dp(gas%T(:, idx_r_p), T)
    idx_rr_T = idx_rl_T + 1

    if ((idx_rl_T < 1) .or. (idx_rl_T > size(gas%T, 1) - 1)) then
    !   call logger%error(fname, "idx_rl_T out of range.")
       return
    end if

    ! Get the water vapor indices
    if (gas%has_H2O) then
       if (H2O <= gas%H2O(1)) then
          idx_l_H2O = 1
       else if (H2O >= gas%H2O(size(gas%H2O, 1) - 1)) then
          idx_l_H2O = size(gas%H2O, 1) - 1
       else
          do i=1, size(gas%H2O) - 1
             if ((H2O > gas%H2O(i)) .and. (H2O <= gas%H2O(i+1))) then
                idx_l_H2O = i
                exit
             end if
          end do
       end if
!!$      idx_l_H2O = searchsorted_dp(gas%H2O, H2O)
       idx_r_H2O = idx_l_H2O + 1
    else
       ! If the cross sections have no water vapour dependence,
       ! we just set both left and right indices to 1. I guess there
       ! is some overhead from this - compared to writing a separate
       ! function that skips the WV completely.
       idx_l_H2O = 1
       idx_r_H2O = 1
    end if

    if (gas%has_H2O) then
       ! Check if H2O is out of range, but this check only makes sense
       ! if we have spectroscopy with an H2O dimension.
       if ((idx_l_H2O < 1) .or. (idx_l_H2O > size(gas%H2O) - 1)) then
       !   call logger%error(fname, "idx_l_H2O out of range.")
          return
       end if
    end if


    ! And perform the linear interpolation in now 3 or 4 dimensions!
    ! Remember, we store the CS grid in the following way:
    ! wavelength, h2o, temperature, pressure
    ! .. so we order the vertices of our polytope in the same way.
    p_d = (p - gas%p(idx_l_p)) / (gas%p(idx_r_p) - gas%p(idx_l_p))
    if (gas%has_H2O) then
       H2O_d = (H2O - gas%H2O(idx_l_H2O)) / (gas%H2O(idx_r_H2O) - gas%H2O(idx_l_H2O))
    else
       H2O_d = 0.0d0
    end if

    ! Normalised temperature when pressure is at lower, and higher index
    T_d_l = (T - gas%T(idx_ll_T, idx_l_p)) / &
         (gas%T(idx_lr_T, idx_l_p) - gas%T(idx_ll_T, idx_l_p))
    T_d_r = (T - gas%T(idx_rl_T, idx_r_p)) / &
         (gas%T(idx_rr_T, idx_r_p) - gas%T(idx_rl_T, idx_r_p))


    do i=1, N_wl

       if (.not. pre_gridded) then

          ! Non-pre-gridded spectroscopy has to be sampled at the requested
          ! wavelengths first.

          if (wl_left_idx(i) == -1) then
             CS_value(i) = 0.0d0
             cycle
          end if

          idx_l_wl = wl_left_idx(i)
          idx_r_wl = idx_l_wl + 1


          wl_d = (wl(i) - gas%wavelength(idx_l_wl)) / &
               (gas%wavelength(idx_r_wl) - gas%wavelength(idx_l_wl))

          ! Start interpolating along wavelength, meaning that C3 = C3(h2o, temperature, pressure)
          C3(0,0,0) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_ll_T, idx_l_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_l_H2O, idx_ll_T, idx_l_p) * wl_d

          C3(1,0,0) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_ll_T, idx_l_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_r_H2O, idx_ll_T, idx_l_p) * wl_d

          C3(0,1,0) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_lr_T, idx_l_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_l_H2O, idx_lr_T, idx_l_p) * wl_d

          C3(0,0,1) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_rl_T, idx_r_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_l_H2O, idx_rl_T, idx_r_p) * wl_d

          C3(1,1,0) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_lr_T, idx_l_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_r_H2O, idx_lr_T, idx_l_p) * wl_d

          C3(0,1,1) = gas%cross_section(idx_l_wl, idx_l_H2O, idx_rr_T, idx_r_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_l_H2O, idx_rr_T, idx_r_p) * wl_d

          C3(1,0,1) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_rl_T, idx_r_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_r_H2O, idx_rl_T, idx_r_p) * wl_d

          C3(1,1,1) = gas%cross_section(idx_l_wl, idx_r_H2O, idx_rr_T, idx_r_p) * (1.0d0 - wl_d) &
               + gas%cross_section(idx_r_wl, idx_r_H2O, idx_rr_T, idx_r_p) * wl_d

       else

          ! For pre-gridded cross sections, we do not need to interpolate along
          ! the wavelength dimension, and the requested wavelength wl(i) corresponds
          ! to the cross section wavelength at that index i.
          C3(0,0,0) = gas%cross_section(i, idx_l_H2O, idx_ll_T, idx_l_p)
          C3(1,0,0) = gas%cross_section(i, idx_r_H2O, idx_ll_T, idx_l_p)
          C3(0,1,0) = gas%cross_section(i, idx_l_H2O, idx_lr_T, idx_l_p)
          C3(0,0,1) = gas%cross_section(i, idx_l_H2O, idx_rl_T, idx_r_p)
          C3(1,1,0) = gas%cross_section(i, idx_r_H2O, idx_lr_T, idx_l_p)
          C3(0,1,1) = gas%cross_section(i, idx_l_H2O, idx_rr_T, idx_r_p)
          C3(1,0,1) = gas%cross_section(i, idx_r_H2O, idx_rl_T, idx_r_p)
          C3(1,1,1) = gas%cross_section(i, idx_r_H2O, idx_rr_T, idx_r_p)

       end if

       ! Next, go across the H2O dimension, such that C2 = C2(temperature, pressure).
       C2(0,0) = C3(0,0,0) * (1.0d0 - H2O_d) + C3(1,0,0) * H2O_d
       C2(0,1) = C3(0,0,1) * (1.0d0 - H2O_d) + C3(1,0,1) * H2O_d
       C2(1,0) = C3(0,1,0) * (1.0d0 - H2O_d) + C3(1,1,0) * H2O_d
       C2(1,1) = C3(0,1,1) * (1.0d0 - H2O_d) + C3(1,1,1) * H2O_d

       ! Next, go across the T dimension - remember that we have two normalised
       ! temperatures, one is for the lower, one is for the higher pressure index!
       ! C1 = C1(pressure)
       C1(0) = C2(0,0) * (1.0d0 - T_d_l) + C2(1,0) * T_d_l
       C1(1) = C2(0,1) * (1.0d0 - T_d_r) + C2(1,1) * T_d_r

       ! Finally, we interpolate along pressure, which is the final result that
       ! we pass back. If below zero, set to zero.
       CS_value(i) = max(C1(0) * (1.0d0 - p_d) + C1(1) * p_d, 0.0d0)

    end do

  end function get_CS_value_at


end module gas_tau_mod
