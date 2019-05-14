module gas_tau_mod

  use control_mod, only: CS_gas
  use math_utils_mod

  public calculate_gas_tau, get_CS_value_at

contains


  subroutine calculate_gas_tau(pre_gridded, is_H2O, &
       wl, gas_vmr, psurf, p, T, sh, &
       gas, N_sub, need_psurf_jac, &
       gas_tau, gas_tau_dpsurf, &
       ndry, success)

    implicit none
    logical, intent(in) :: pre_gridded ! Is the spectroscopy pre-gridded?
    logical, intent(in) :: is_H2O ! Is this gas water vapor?
    double precision, intent(in) :: wl(:) ! Wavelength array
    double precision, intent(in) :: gas_vmr(:) ! Gas volume mixing ratio
    double precision, intent(in) :: psurf ! Surface pressure
    ! Pressure, temperature and specific humidity profiles
    double precision, intent(in) :: p(:), T(:), sh(:)
    type(CS_gas), intent(in) :: gas ! Gas structure which contains the cross sections
    integer, intent(in) :: N_sub ! How many sublayers for more accurate calculation?
    logical, intent(in) :: need_psurf_jac
    ! Gas optical depth - wavelength, layer
    double precision, intent(inout) :: gas_tau(:,:)
    ! Gas optical depth after psurf pert. - wavelength, layer
    double precision, intent(inout) :: gas_tau_dpsurf(:,:)
    ! dry air mass per layer
    double precision, intent(inout) :: ndry(:)
    ! Success?
    logical, intent(inout) :: success


    character(len=999) :: tmp_str
    character(len=*), parameter :: fname = "calculate_gas_tau"
    logical :: log_scaling
    double precision :: p_lower, p_higher, T_lower, T_higher, sh_lower, sh_higher
    double precision :: VMR_lower, VMR_higher
    double precision :: H2O_corr
    double precision :: this_H2O, this_H2O_pert
    double precision :: this_p, p_fac, this_p_fac, this_T, this_sh, this_VMR, this_M
    double precision :: this_p_pert, p_fac_pert, this_p_fac_pert, p_lower_pert, T_lower_pert, &
         sh_lower_pert, VMR_lower_pert, VMR_sigma_pert, this_VMR_pert, this_M_pert, &
         this_T_pert, this_sh_pert
    double precision :: C_tmp, p_lower_gas_pert, this_p_fac_lower
    double precision, allocatable :: gas_tmp(:), this_CS_value(:)
    integer, allocatable :: wl_left_indices(:)

    double precision, allocatable :: GK_abscissae(:), GK_weights(:), G_weights(:)
    double precision, allocatable :: GK_abscissae_f(:), GK_weights_f(:), G_weights_f(:)
    double precision :: GK_abscissae_f_pert(N_sub), GK_weights_f_pert(N_sub), G_weights_f_pert(N_sub)

    integer :: N_lay, N_lev, N_wl, num_active_levels
    integer :: i,j,k,l
    integer :: funit

    double precision :: gas_wl_step_avg, gas_wl_start
    integer :: gas_idx_fg

    double precision :: cpu_start, cpu_end

    allocate(this_CS_value(size(wl)))
    allocate(wl_left_indices(size(wl)))

    allocate(GK_abscissae((N_sub+1)/2))
    allocate(GK_weights((N_sub+1)/2))
    allocate(G_weights((N_sub+1)/2))

    allocate(GK_abscissae_f(N_sub))
    allocate(GK_weights_f(N_sub))
    allocate(G_weights_f(N_sub))

    N_lev = size(gas_vmr)
    N_lay = N_lev - 1
    N_wl = size(wl)

    do j=1, N_lev
       if (psurf >= p(j)) then
          num_active_levels = j+1
       end if
    end do

    if (num_active_levels > N_lev) then
       write(tmp_str, '(A,F12.2,A,F12.2)') "Psurf: ", psurf, " larger than BOA level: ", p(N_lev)
       call logger%error(fname, trim(tmp_str))
       success = .false.
       return
    end if

    log_scaling = .true.

    allocate(gas_tmp(N_wl))

    ! Just to make sure there's nothing in there already..
    gas_tau(:,:) = 0.0d0
    gas_tau_dpsurf(:,:) = 0.0d0
    ndry(:) = 0.0d0

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
    do l=num_active_levels, 2, -1

       ! First, grab the lower (altitude-wise) and higher values for
       ! p,T,sh from the atmosphere profiles for whatever layer l we're in.

       if (l == num_active_levels) then
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
       end if

       p_higher = p(l-1)
       T_higher = T(l-1)
       sh_higher = sh(l-1)
       VMR_higher = gas_vmr(l-1)

       ! Should the perturbed surface pressure actually fall onto the next-higher level,
       ! we need to make a small adjustment, otherwise, we end up dividing by zero later on.
       if (need_psurf_jac) then
          if (abs(p_lower_pert - p_higher) < 1d-3) then
             p_lower_pert = p_lower_pert + PSURF_PERTURB/10.0d0
          end if
       end if


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

       do k=1, N_sub

          ! Get the value of pressure according to the Gauss-Kronrod abscissa
          ! value, at which we evaluate the cross section.
          this_p = GK_abscissae_f(k)

          ! And get corresponding values for T,sh and the VMR at this (scaled)
          ! pressure value, which is obtained via simple linear interpolation.
          if (log_scaling) then
             this_p_fac = (this_p - log(p_higher)) / (log(p_lower) - log(p_higher))
          else
             this_p_fac = (this_p - p_higher) / (p_lower - p_higher)
          end if

          this_T = (1.0d0 - this_p_fac) * T_lower + this_p_fac * T_higher
          this_sh = (1.0d0 - this_p_fac) * sh_lower + this_p_fac * sh_higher
          this_VMR = (1.0d0 - this_p_fac) * VMR_lower + this_p_fac * VMR_higher
          this_M = 1.0d3 * DRY_AIR_MASS

          ! Gas CS routine works in H2O VMR rather than SH
          this_H2O = this_sh / (1.0d0 - this_sh) * SH_H2O_CONV

          if (log_scaling) this_p = exp(this_p)

          this_CS_value = get_CS_value_at(pre_gridded, gas, wl(:), this_p, &
               this_T, this_H2O, wl_left_indices(:))

          if (is_H2O) then
             H2O_corr = 1.0d0
          else
             H2O_corr = (1.0d0 - this_sh)
          end if

          C_tmp = 1.0d0 / 9.80665d0 * NA * 0.1d0 * H2O_corr / this_M

          ! Tau for this sublayer
          if (log_scaling) then
             !gas_tmp(:) = GK_weights_f(k) * this_CS_value(:) * this_VMR * C_tmp * this_p
             ndry(l-1) = GK_weights_f(k) * C_tmp * this_p
          else
             !gas_tmp(:) = GK_weights_f(k) * this_CS_value(:) * this_VMR * C_tmp
             ndry(l-1) = GK_weights_f(k) * C_tmp
          end if

          ! Add sublayer contribution to full layer tau
          gas_tau(:,l-1) = gas_tau(:,l-1) + (ndry(l-1) * this_CS_value(:) * this_VMR)

          ! The same is required in the case of surface pressure jacobians,
          ! but we obviously only do this for the BOA layer
          if (need_psurf_jac .and. (l == num_active_levels)) then

             this_p_pert = GK_abscissae_f_pert(k)

             if (log_scaling) then
                this_p_fac_pert = (this_p_pert - log(p_higher)) / (log(p_lower_pert) - log(p_higher))
             else
                this_p_fac_pert = ((this_p_pert - p_higher) / (p_lower_pert - p_higher))
             end if

             this_T_pert = (1.0d0 - this_p_fac_pert) * T_lower_pert + this_p_fac_pert * T_higher
             this_sh_pert = (1.0d0 - this_p_fac_pert) * sh_lower_pert + this_p_fac_pert * sh_higher
             this_VMR_pert = (1.0d0 - this_p_fac_pert) * VMR_lower_pert + this_p_fac_pert * VMR_higher
             this_M_pert = 1.0d3 * DRY_AIR_MASS

             this_H2O_pert = this_sh_pert / (1.0d0 - this_sh_pert) * SH_H2O_CONV

             if (log_scaling) this_p_pert = exp(this_p_pert)

             this_CS_value = get_CS_value_at(pre_gridded, gas, wl(:), this_p_pert, &
                  this_T_pert, this_H2O_pert, wl_left_indices(:))

             if (is_H2O) then
                H2O_corr = 1.0d0
             else
                H2O_corr = (1.0d0 - this_sh_pert)
             end if
             ! This calculates the gas OD, as a result of a perturbed surface pressure
             C_tmp = 1.0d0 / (9.80665d0 * this_M) * NA * 0.1d0 * H2O_corr

             if (log_scaling) then
                gas_tmp(:) = GK_weights_f_pert(k) * this_CS_value(:) &
                * this_VMR_pert * C_tmp * exp(this_p_pert)
             else
                gas_tmp(:) = GK_weights_f_pert(k) * this_CS_value(:) &
                     * this_VMR_pert * C_tmp
             end if

             gas_tau_dpsurf(:,l-1) = gas_tau_dpsurf(:,l-1) + gas_tmp(:)

          end if
       end do

       if (need_psurf_jac .and. (l == num_active_levels)) then
          ! Get the difference: tau(psurf - psurb_perturb) - tau(psurf)
          gas_tau_dpsurf(:,l-1) = (gas_tau_dpsurf(:,l-1) - gas_tau(:,l-1)) / PSURF_PERTURB
       end if

       if (need_psurf_jac .and. (l < num_active_levels)) then
          ! Copy the non-BOA layer ODs to the gas_tau_dpsurf array, as the
          ! surface pressure Jacobian merely affects the BOA layer ODs
          gas_tau_dpsurf(:,l-1) = 0.0d0
       end if
    end do

    success = .true.

  end subroutine calculate_gas_tau



 pure function get_CS_value_at(pre_gridded, gas, wl, p, T, H2O, wl_left_idx) result(CS_value)

   implicit none
   logical, intent(in) :: pre_gridded
   type(CS_gas), intent(in) :: gas
   double precision, intent(in) :: wl(:), p, T, H2O
   integer, intent(in) :: wl_left_idx(:)

   double precision :: CS_value(size(wl))

   ! These here are the indices between which the CS array will
   ! be linearly interpolated: e.g. idx_(left, right)_(pressure)
   integer :: idx_closest
   integer :: idx_l_p, idx_r_p ! Pressure
   integer :: idx_lr_T, idx_ll_T, idx_rl_T, idx_rr_T ! Temperature has two dimensions
   integer :: idx_l_H2O, idx_r_H2O
   integer :: idx_l_wl, idx_r_wl
   double precision :: C3(0:1, 0:1, 0:1) ! 0 is 'left', 1 is 'right'
   double precision :: C2(0:1, 0:1), C1(0:1)
   double precision :: wl_d, p_d, T_d_l, T_d_r, H2O_d
   integer :: i,j,k


   ! This next section here "just" finds the left and right indices of p,T,H2O,wl values within the
   ! grids of the cross section tables. Unfortunately, sourcing this bit out into a PURE FUNCTION still
   ! comes with such a large overhead that I had to inline them explicitly. Maybe there's a way around
   ! it, but for now the performance increase is totally worth it!

   ! Get the pressure indices
   if (p <= gas%p(1)) then
      idx_l_p = 1
   else if (p >= gas%p(size(gas%p))) then
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

   ! Get the temperature indices, which depend on P - so we have two
   ! indices for the temperature dimension. One set of T indices (idx_l(l,r)_T)
   ! corresponds to the left index of pressure (idx_l_p), the other set
   ! (idx_r(l,r)_T) corresponds to the right index of pressure (idx_r_p).

   if (T <= gas%T(1, idx_l_p)) then
      idx_ll_T = 1
   else if (T >= gas%T(size(gas%T, 1), idx_l_p)) then
      idx_ll_T = size(gas%T, 1) - 1
   else
      do i=1, size(gas%T, 1)-1
         if ((T > gas%T(i, idx_l_p)) .and. (T <= gas%T(i+1, idx_l_p))) then
            idx_ll_T = i
            exit
         end if
      end do
   end if
!!$   idx_ll_T = searchsorted_dp(gas%T(:, idx_l_p), T)
   idx_lr_T = idx_ll_T + 1

   if (T <= gas%T(1, idx_r_p)) then
      idx_rl_T = 1
   else if (T >= gas%T(size(gas%T, 1), idx_r_p)) then
      idx_rl_T = size(gas%T, 1) - 1
   else
      do i=1, size(gas%T, 1)-1
         if ((T >= gas%T(i, idx_r_p)) .and. (T < gas%T(i+1, idx_r_p))) then
            idx_rl_T = i
            exit
         end if
      end do
   end if
!!$   idx_rl_T = searchsorted_dp(gas%T(:, idx_r_p), T)
   idx_rr_T = idx_rl_T + 1

   ! Get the water vapor indices
   if (gas%has_H2O) then
      if (H2O <= gas%H2O(1)) then
         idx_l_H2O = 1
      else if (H2O >= gas%H2O(size(gas%H2O, 1))) then
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

   do i=1, size(wl)

      if (.not. pre_gridded) then
         wl_d = (wl(i) - gas%wavelength(idx_l_wl)) / &
              (gas%wavelength(idx_r_wl) - gas%wavelength(idx_l_wl))

         if (wl_left_idx(i) == -1) then
            CS_value(i) = 0.0d0
            cycle
         end if

         idx_l_wl = wl_left_idx(i)
         idx_r_wl = idx_l_wl + 1

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
      !if (CS_value(i) < 0.0d0) CS_value(i) = 0.0d0

   end do

 end function get_CS_value_at


end module gas_tau_mod
