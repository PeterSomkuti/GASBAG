  !> @brief Common Jacobian calculations (not RT dependent)
  !> @file Jacobians.f90
  !> @author Peter Somkuti

module jacobians_mod


  use math_utils_mod, only: pwl_value_1d_v2, oco_type_convolution
  use control_mod, only: MCS
  use instruments_mod, only: generic_instrument
  use oco2_mod
  use statevector_mod, only: statevector

contains


  !> @brief Calculate solar shift/stretch jacobians
  !> @param solar_wl Solar spectrum wavelength array
  !> @param solar_trans Solar transmittance spectrum array (same size as solar_wl)
  !> @param solar_irrad Solar irradiance (same size as solar_wl)
  !> @param hires_grid Hi-res wavelength grid
  !> @param new_solar_trans Solar transmittance at hires grid
  !> @param dsolar_dlambda dSolar/dLambda weighting function
  subroutine calculate_solar_jacobian(solar_wl, solar_trans, solar_irrad, &
       hires_grid, new_solar_trans, dsolar_dlambda)

    implicit none

    double precision, intent(in) :: solar_wl(:)
    double precision, intent(in) :: solar_trans(:)
    double precision, intent(in), allocatable :: solar_irrad(:)
    double precision, intent(in) :: hires_grid(:)
    double precision, intent(inout) :: new_solar_trans(:)
    double precision, intent(inout) :: dsolar_dlambda(:)

    ! Size of hires wavelength array
    integer :: N_hires

    double precision, allocatable :: solar_tmp(:)
    integer :: i

    N_hires = size(solar_wl)
    allocate(solar_tmp(N_hires))

    ! If we retrieve a solar stretch or shift,
    ! "solar_shift" and/or "solar_stretch" are /= 1.0.
    ! Thus, the solar spectrum needs to be re-sampled at the
    ! hires-grid given the new shifted solar wavelength array.

    call pwl_value_1d_v2( &
         N_hires, &
         solar_wl, &
         solar_trans, &
         N_hires, &
         hires_grid, solar_tmp(:))

    new_solar_trans(:) = solar_tmp(:)
    deallocate(solar_tmp)

    ! Calculate dSolar / dlambda using central differences
    ! apart from the first and last point obviously. This is required for
    ! the solar parameter Jacobians.

    dsolar_dlambda(1) = (solar_trans(2) - solar_trans(1)) / &
         (solar_wl(2) - solar_wl(1))
    dsolar_dlambda(N_hires) = (solar_trans(N_hires-1) - solar_trans(N_hires)) / &
         (solar_wl(N_hires-1) -  solar_wl(N_hires))

    do i=2, N_hires-1
       dsolar_dlambda(i) = (solar_trans(i+1) - solar_trans(i-1)) / &
            (solar_wl(i+1) - solar_wl(i-1))
    end do


  end subroutine calculate_solar_jacobian



  !> @brief Dispersion Jacobian calculation
  !> @param my_instrument Instrument instance
  !> @param disp_idx Dispersion order index
  !> @param dispersion_coefs Dispersion polynomial coefficients for this band/FP
  !> @param band This band
  !> @param i_win This retrieval window
  !> @param i_fp This footprint
  !> @param N_spec Number of spectral points in this retrieval window
  !> @param ILS_delta_lambda ILS wavelength array
  !> @param ILS_relative_response ILS values
  !> @param l1b_wl_idx_min Lower WL index in total detector pixels
  !> @param l1b_wl_idx_max Higher WL index in total detector pixels
  !> @param instrument_doppler Instrument doppler factor
  !> @param hires_grid Hires WL grid
  !> @param radiance_calc_work Modelled detector-resolution radiances
  !> @param radiance_calc_work_hi Modelled high-resolution radiances
  !> @param dispersion_jacobian Output dIntensity/dDispersionCoefficient
  subroutine calculate_dispersion_jacobian(my_instrument, disp_idx, &
       dispersion_coefs, band, i_win, i_fp, N_spec, &
       ILS_delta_lambda, ILS_relative_response, &
       l1b_wl_idx_min, l1b_wl_idx_max, &
       instrument_doppler, hires_grid, radiance_calc_work, &
       radiance_calc_work_hi, dispersion_jacobian)

    implicit none
    class(generic_instrument), intent(in) :: my_instrument
    double precision, intent(in) :: dispersion_coefs(:)
    integer, intent(in) :: disp_idx
    integer, intent(in) :: band
    integer, intent(in) :: i_win
    integer, intent(in) :: i_fp
    integer, intent(in) :: N_spec
    double precision, intent(in) :: ILS_delta_lambda(:,:)
    double precision, intent(in) :: ILS_relative_response(:,:,:,:)
    integer, intent(in) :: l1b_wl_idx_min
    integer, intent(in) :: l1b_wl_idx_max
    double precision, intent(in) :: instrument_doppler
    double precision, intent(in) :: hires_grid(:)
    double precision, intent(in) :: radiance_calc_work(:)
    double precision, intent(in) :: radiance_calc_work_hi(:)
    double precision, intent(inout) :: dispersion_jacobian(:)

    ! Function name for logging
    character(len=99), parameter :: fname = "calculate_dispersion_jacobian"
    ! Perturbed dispersion coefficients
    double precision, allocatable :: dispersion_pert(:)
    ! New / perturbed dispersion array
    double precision, allocatable :: this_dispersion(:)
    ! Radiances after perturbation of dispersion coeffs
    double precision, allocatable :: radiance_tmp_work(:)
    ! Convolution with ILS successful?
    logical :: ILS_success

    allocate(dispersion_pert(size(dispersion_coefs, 1)))
    allocate(this_dispersion(MCS%general%N_spec(band)))
    allocate(radiance_tmp_work(N_spec))

    ! Perturb the dispersion coefficient (given the index disp_idx)
    dispersion_pert(:) = dispersion_coefs(:)
    dispersion_pert(disp_idx) = dispersion_pert(disp_idx) &
         + MCS%window(i_win)%dispersion_pert(disp_idx)

    select type(my_instrument)
    type is (oco2_instrument)

       ! Calculate new dispersion grid using perturbed coefficients
       call my_instrument%calculate_dispersion(dispersion_pert, &
            this_dispersion(:), band, i_fp)

       ! Apply instrument Doppler shift
       this_dispersion = this_dispersion / (1.0d0 - instrument_doppler)

       ! Convolve the perturbed TOA radiance
       call oco_type_convolution(hires_grid, radiance_calc_work_hi, &
            ILS_delta_lambda(:,:), &
            ILS_relative_response(:,l1b_wl_idx_min:l1b_wl_idx_max,i_fp,band), &
            this_dispersion(l1b_wl_idx_min:l1b_wl_idx_max), radiance_tmp_work, &
            ILS_success)

       if (.not. ILS_success) then
          call logger%error(fname, "ILS convolution error.")
          return
       end if

       ! Calculate the dispersion jacobian via finite differencing
       dispersion_jacobian(:) = (radiance_tmp_work - radiance_calc_work) &
            / MCS%window(i_win)%dispersion_pert(disp_idx)

    end select

  end subroutine calculate_dispersion_jacobian


  subroutine calculate_ILS_stretch_jacobian(my_instrument, coeff_idx, SV, &
       band, i_win, i_fp, N_spec, &
       ILS_delta_lambda, ILS_relative_response, &
       ILS_stretch, this_dispersion, &
       l1b_wl_idx_min, l1b_wl_idx_max, &
       hires_grid, center_pixel, radiance_calc_work, &
       radiance_calc_work_hi, ILS_stretch_jacobian)

    implicit none
    class(generic_instrument), intent(in) :: my_instrument
    integer, intent(in) :: coeff_idx
    type(statevector) :: SV
    integer, intent(in) :: band
    integer, intent(in) :: i_win
    integer, intent(in) :: i_fp
    integer, intent(in) :: N_spec
    double precision, intent(in) :: ILS_delta_lambda(:,:)
    double precision, intent(in) :: ILS_relative_response(:,:)
    double precision, intent(in) :: ILS_stretch(:)
    double precision, intent(in) :: this_dispersion(:)
    integer, intent(in) :: l1b_wl_idx_min
    integer, intent(in) :: l1b_wl_idx_max
    double precision, intent(in) :: hires_grid(:)
    integer, intent(in) :: center_pixel
    double precision, intent(in) :: radiance_calc_work(:)
    double precision, intent(in) :: radiance_calc_work_hi(:)
    double precision, intent(inout) :: ILS_stretch_jacobian(:)

    ! Function name for logging
    character(len=99), parameter :: fname = "calculate_ILS_stretch_jacobian"
    double precision, allocatable :: ILS_delta_lambda_pert(:,:)
    double precision, allocatable :: ILS_stretch_pert(:)
    double precision, allocatable :: radiance_tmp_work(:)
    logical :: ILS_success

    integer :: i, l

    allocate(ILS_delta_lambda_pert, mold=ILS_delta_lambda)
    allocate(ILS_stretch_pert, mold=ILS_stretch)
    allocate(radiance_tmp_work, mold=radiance_calc_work)

    ILS_stretch_pert(:) = 0.0d0
    ILS_delta_lambda_pert(:,:) = ILS_delta_lambda(:,:)

    do i=1, N_spec
       do l=1, SV%num_ils_stretch
          if (l == coeff_idx) then
             ILS_stretch_pert(i) = ILS_stretch_pert(i) &
                  + (dble(i - center_pixel) ** (l-1) * (MCS%window(i_win)%ils_stretch_pert(l) &
                  + SV%svsv(SV%idx_ils_stretch(l))))
          else
             ILS_stretch_pert(i) = ILS_stretch_pert(i) &
                  + (dble(i - center_pixel) ** (l-1) * SV%svsv(SV%idx_ils_stretch(l)))
          end if
       end do
    end do

    do i=1, N_spec
       ILS_delta_lambda_pert(:,i) = ILS_delta_lambda_pert(:,i) * ILS_stretch_pert(i)
    end do

    select type(my_instrument)
    type is (oco2_instrument)

       call oco_type_convolution(hires_grid, radiance_calc_work_hi, &
            ILS_delta_lambda_pert(:,:), &
            ILS_relative_response(:,:), &
            this_dispersion, radiance_tmp_work, &
            ILS_success)

       if (.not. ILS_success) then
          call logger%error(fname, "ILS convolution error.")
          return
       end if

    end select

    ILS_stretch_jacobian(:) = (radiance_tmp_work - radiance_calc_work) &
         / MCS%window(i_win)%ils_stretch_pert(coeff_idx)

  end subroutine calculate_ILS_stretch_jacobian



end module jacobians_mod
