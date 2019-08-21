#define xvec_zero		XCAT(TYPE_PREFIX, vec_zero)
#define xvec_copy		XCAT(TYPE_PREFIX, vec_copy)
#define xvec_print		XCAT(TYPE_PREFIX, vec_print)
#define xvec_scale		XCAT(TYPE_PREFIX, vec_scale)
#define xvec_add 		XCAT(TYPE_PREFIX, vec_add)
#define xvec_sub 		XCAT(TYPE_PREFIX, vec_sub)
#define xvec_lincmb		XCAT(TYPE_PREFIX, vec_lincmb)
#define xvec_dot		XCAT(TYPE_PREFIX, vec_dot)
#define xvec_mag 		XCAT(TYPE_PREFIX, vec_mag)
#define xvec_unit		XCAT(TYPE_PREFIX, vec_unit)
#define xvec_cross		XCAT(TYPE_PREFIX, vec_cross)
#define xvec_dyadic		XCAT(TYPE_PREFIX, vec_dyadic)
#define xvec_inv		XCAT(TYPE_PREFIX, vec_inv

#define xmat_zero		XCAT(TYPE_PREFIX, mat_zero)
#define xmat_init		XCAT(TYPE_PREFIX, mat_init)
#define xmat_copy		XCAT(TYPE_PREFIX, mat_copy)
#define xmat_print		XCAT(TYPE_PREFIX, mat_print)
#define xmat_trans		XCAT(TYPE_PREFIX, mat_trans)
#define xmat_p_one_norm		XCAT(TYPE_PREFIX, mat_p_one_norm)
#define xmat_p_inf_norm		XCAT(TYPE_PREFIX, mat_p_inf_norm)
#define xmat_frob_norm		XCAT(TYPE_PREFIX, mat_frob_norm)
#define xmat_add		XCAT(TYPE_PREFIX, mat_add)
#define xmat_add_diag		XCAT(TYPE_PREFIX, mat_add_diag)
#define xmat_add_trans		XCAT(TYPE_PREFIX, mat_add_trans)
#define xmat_sub		XCAT(TYPE_PREFIX, mat_sub)
#define xmat_sub_diag		XCAT(TYPE_PREFIX, mat_sub_diag)
#define xmat_diag_sub		XCAT(TYPE_PREFIX, mat_diag_sub)
#define xmat_sub_trans		XCAT(TYPE_PREFIX, mat_sub_trans)
#define xmat_i_sub		XCAT(TYPE_PREFIX, mat_i_sub)
#define xmat_scale		XCAT(TYPE_PREFIX, mat_scale)
#define xm_v_mul		XCAT(TYPE_PREFIX, m_v_mul)
#define xv_m_mul		XCAT(TYPE_PREFIX, v_m_mul)
#define xm_v_diag_mul		XCAT(TYPE_PREFIX, m_v_diag_mul)
#define xm_v_dinv_mul		XCAT(TYPE_PREFIX, m_v_dinv_mul)
#define xmat_mul		XCAT(TYPE_PREFIX, mat_mul)
#define xmat_gxgxmx		XCAT(TYPE_PREFIX, mat_gxgxmx)
#define xmat_mul_diag		XCAT(TYPE_PREFIX, mat_mul_diag)
#define xmat_gxdxmx		XCAT(TYPE_PREFIX, mat_gxdxmx)
#define xmat_mul_dinv		XCAT(TYPE_PREFIX, mat_mul_dinv)
#define xmat_diag_mul		XCAT(TYPE_PREFIX, mat_diag_mul)
#define xmat_dinv_mul		XCAT(TYPE_PREFIX, mat_dinv_mul)
#define xmat_pow_count		XCAT(TYPE_PREFIX, mat_pow_count)
#define xmat_pow		XCAT(TYPE_PREFIX, mat_pow)

#define xmat_pocon		XCAT(TYPE_PREFIX, mat_pocon)
#define xmat_gecon		XCAT(TYPE_PREFIX, mat_gecon)
#define xmat_potrf		XCAT(TYPE_PREFIX, mat_potrf)
#define xmat_potri		XCAT(TYPE_PREFIX, mat_potri)
#define xmat_potrs		XCAT(TYPE_PREFIX, mat_potrs)
#define xmat_getrf		XCAT(TYPE_PREFIX, mat_getrf)
#define xmat_getri		XCAT(TYPE_PREFIX, mat_getri)
#define xmat_getrs		XCAT(TYPE_PREFIX, mat_getrs)
#define xmat_getrs2		XCAT(TYPE_PREFIX, mat_getrs2)
#define xmat_gttrf		XCAT(TYPE_PREFIX, mat_gttrf)
#define xmat_gttrs		XCAT(TYPE_PREFIX, mat_gttrs)

#define xgbtrf_			XCAT(TYPE_PREFIX, gbtrf_)
#define xgbtrs_			XCAT(TYPE_PREFIX, gbtrs_)

#define get_work_x1		XCAT(XCAT(get_work_, TYPE_POSTFIX), 1)
#define get_work_x2		XCAT(XCAT(get_work_, TYPE_POSTFIX), 2)
#define get_work_x3		XCAT(XCAT(get_work_, TYPE_POSTFIX), 3)

#define copy_to_band_storage_x	XCAT(copy_to_band_storage_,  TYPE_POSTFIX)
#define copy_to_band_storage2_x	XCAT(copy_to_band_storage2_, TYPE_POSTFIX)
