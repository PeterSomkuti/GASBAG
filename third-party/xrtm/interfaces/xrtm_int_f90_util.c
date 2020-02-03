/*******************************************************************************
 * This file was generated by bindx version 0.01.  Edit at your own risk.
 ******************************************************************************/

#include <gutil.h>

#include <xrtm_interface.h>

#ifdef __cplusplus
extern "C" {
#endif


int xrtm_set_levels_b_l_n_bindx_f90(xrtm_data *d, double *levels_b_l)
{
     double **levels_b_l2;
     levels_b_l2 = array_from_mem2_d(levels_b_l, d->n_layers + 1, d->n_derivs);
     if (xrtm_set_levels_b_l_n(d, levels_b_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_levels_b_l_n()\n");
          return -1;
     }
     free_array2_d(levels_b_l2);
     return 0;
}


int xrtm_set_g_l_nn_bindx_f90(xrtm_data *d, double *g_l)
{
     double **g_l2;
     g_l2 = array_from_mem2_d(g_l, d->n_layers, d->n_derivs);
     if (xrtm_set_g_l_nn(d, g_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_g_l_nn()\n");
          return -1;
     }
     free_array2_d(g_l2);
     return 0;
}


int xrtm_set_coef_1_bindx_f90(xrtm_data *d, int i_layer, int n_coef_layer, double *coef)
{
     double **coef2;
     coef2 = array_from_mem2_d(coef, d->n_elem, n_coef_layer);
     if (xrtm_set_coef_1(d, i_layer, n_coef_layer, coef2)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_1()\n");
          return -1;
     }
     free_array2_d(coef2);
     return 0;
}


int xrtm_set_coef_n_bindx_f90(xrtm_data *d, int *n_coef_layer, double *coef)
{
     double ***coef2;
     coef2 = array_from_mem3_d(coef, d->n_layers, d->n_elem, d->n_coef);
     if (xrtm_set_coef_n(d, n_coef_layer, coef2)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_n()\n");
          return -1;
     }
     free_array3_d(coef2);
     return 0;
}


int xrtm_set_coef_l_11_bindx_f90(xrtm_data *d, int i_layer, int i_deriv, double *coef_l)
{
     double **coef_l2;
     coef_l2 = array_from_mem2_d(coef_l, d->n_elem, d->n_coef_layer[i_layer]);
     if (xrtm_set_coef_l_11(d, i_layer, i_deriv, coef_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_l_11()\n");
          return -1;
     }
     free_array2_d(coef_l2);
     return 0;
}


int xrtm_set_coef_l_n1_bindx_f90(xrtm_data *d, int i_deriv, double *coef_l)
{
     double ***coef_l2;
     coef_l2 = array_from_mem3_d(coef_l, d->n_layers, d->n_elem, d->n_coef);
     if (xrtm_set_coef_l_n1(d, i_deriv, coef_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_l_n1()\n");
          return -1;
     }
     free_array3_d(coef_l2);
     return 0;
}


int xrtm_set_coef_l_1n_bindx_f90(xrtm_data *d, int i_layer, double *coef_l)
{
     double ***coef_l2;
     coef_l2 = array_from_mem3_d(coef_l, d->n_derivs, d->n_elem, d->n_coef_layer[i_layer]);
     if (xrtm_set_coef_l_1n(d, i_layer, coef_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_l_1n()\n");
          return -1;
     }
     free_array3_d(coef_l2);
     return 0;
}


int xrtm_set_coef_l_nn_bindx_f90(xrtm_data *d, double *coef_l)
{
     double ****coef_l2;
     coef_l2 = array_from_mem4_d(coef_l, d->n_layers, d->n_derivs, d->n_elem, d->n_coef);
     if (xrtm_set_coef_l_nn(d, coef_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_coef_l_nn()\n");
          return -1;
     }
     free_array4_d(coef_l2);
     return 0;
}


int xrtm_set_omega_l_nn_bindx_f90(xrtm_data *d, double *omega_l)
{
     double **omega_l2;
     omega_l2 = array_from_mem2_d(omega_l, d->n_layers, d->n_derivs);
     if (xrtm_set_omega_l_nn(d, omega_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_omega_l_nn()\n");
          return -1;
     }
     free_array2_d(omega_l2);
     return 0;
}


int xrtm_set_ltau_l_nn_bindx_f90(xrtm_data *d, double *ltau_l)
{
     double **ltau_l2;
     ltau_l2 = array_from_mem2_d(ltau_l, d->n_layers, d->n_derivs);
     if (xrtm_set_ltau_l_nn(d, ltau_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_ltau_l_nn()\n");
          return -1;
     }
     free_array2_d(ltau_l2);
     return 0;
}


int xrtm_set_kernel_params_l_nn_bindx_f90(xrtm_data *d, int i_kernel, double *params_l)
{
     double **params_l2;
     params_l2 = array_from_mem2_d(params_l, d->n_derivs, kernel_n_params((enum xrtm_kernel_type) d->kernels[i_kernel]));
     if (xrtm_set_kernel_params_l_nn(d, i_kernel, params_l2)) {
          fprintf(stderr, "ERROR: xrtm_set_kernel_params_l_nn()\n");
          return -1;
     }
     free_array2_d(params_l2);
     return 0;
}


int xrtm_solution_bindx_f90(xrtm_data *d, enum xrtm_solver_mask solver, int solutions, int n_out_phis, double *out_phis, double *I_p, double *I_m, double *K_p, double *K_m, double *mean_p, double *mean_m, double *mean_p_l, double *mean_m_l, double *flux_p, double *flux_m, double *flux_p_l, double *flux_m_l, double *flux_div, double *flux_div_l)
{
     double **out_phis2;
     double ****I_p2;
     double ****I_m2;
     double *****K_p2;
     double *****K_m2;
     double **mean_p_l2;
     double **mean_m_l2;
     double **flux_p_l2;
     double **flux_m_l2;
     double **flux_div_l2;
     out_phis2 = array_from_mem2_d(out_phis, d->n_umus, n_out_phis);
     I_p2 = array_from_mem4_d(I_p, d->n_ulevels, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     I_m2 = array_from_mem4_d(I_m, d->n_ulevels, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     K_p2 = array_from_mem5_d(K_p, d->n_ulevels, d->n_derivs, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     K_m2 = array_from_mem5_d(K_m, d->n_ulevels, d->n_derivs, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     mean_p_l2 = array_from_mem2_d(mean_p_l, d->n_ulevels, d->n_derivs);
     mean_m_l2 = array_from_mem2_d(mean_m_l, d->n_ulevels, d->n_derivs);
     flux_p_l2 = array_from_mem2_d(flux_p_l, d->n_ulevels, d->n_derivs);
     flux_m_l2 = array_from_mem2_d(flux_m_l, d->n_ulevels, d->n_derivs);
     flux_div_l2 = array_from_mem2_d(flux_div_l, d->n_ulevels, d->n_derivs);
     if (xrtm_solution(d, solver, solutions, n_out_phis, out_phis2, I_p2, I_m2, K_p2, K_m2, mean_p, mean_m, mean_p_l2, mean_m_l2, flux_p, flux_m, flux_p_l2, flux_m_l2, flux_div, flux_div_l2)) {
          fprintf(stderr, "ERROR: xrtm_solution()\n");
          return -1;
     }
     free_array2_d(out_phis2);
     free_array4_d(I_p2);
     free_array4_d(I_m2);
     free_array5_d(K_p2);
     free_array5_d(K_m2);
     free_array2_d(mean_p_l2);
     free_array2_d(mean_m_l2);
     free_array2_d(flux_p_l2);
     free_array2_d(flux_m_l2);
     free_array2_d(flux_div_l2);
     return 0;
}


int xrtm_radiance_bindx_f90(xrtm_data *d, enum xrtm_solver_mask solver, int n_out_phis, double *out_phis, double *I_p, double *I_m, double *K_p, double *K_m)
{
     double **out_phis2;
     double ****I_p2;
     double ****I_m2;
     double *****K_p2;
     double *****K_m2;
     out_phis2 = array_from_mem2_d(out_phis, d->n_umus, n_out_phis);
     I_p2 = array_from_mem4_d(I_p, d->n_ulevels, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     I_m2 = array_from_mem4_d(I_m, d->n_ulevels, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     K_p2 = array_from_mem5_d(K_p, d->n_ulevels, d->n_derivs, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     K_m2 = array_from_mem5_d(K_m, d->n_ulevels, d->n_derivs, d->n_umus == 0 ? d->n_quad : d->n_umus, n_out_phis, d->n_stokes);
     if (xrtm_radiance(d, solver, n_out_phis, out_phis2, I_p2, I_m2, K_p2, K_m2)) {
          fprintf(stderr, "ERROR: xrtm_radiance()\n");
          return -1;
     }
     free_array2_d(out_phis2);
     free_array4_d(I_p2);
     free_array4_d(I_m2);
     free_array5_d(K_p2);
     free_array5_d(K_m2);
     return 0;
}


int xrtm_mean_radiance_bindx_f90(xrtm_data *d, enum xrtm_solver_mask solver, double *mean_p, double *mean_m, double *mean_p_l, double *mean_m_l)
{
     double **mean_p_l2;
     double **mean_m_l2;
     mean_p_l2 = array_from_mem2_d(mean_p_l, d->n_ulevels, d->n_derivs);
     mean_m_l2 = array_from_mem2_d(mean_m_l, d->n_ulevels, d->n_derivs);
     if (xrtm_mean_radiance(d, solver, mean_p, mean_m, mean_p_l2, mean_m_l2)) {
          fprintf(stderr, "ERROR: xrtm_mean_radiance()\n");
          return -1;
     }
     free_array2_d(mean_p_l2);
     free_array2_d(mean_m_l2);
     return 0;
}


int xrtm_flux_bindx_f90(xrtm_data *d, enum xrtm_solver_mask solver, double *flux_p, double *flux_m, double *flux_p_l, double *flux_m_l)
{
     double **flux_p_l2;
     double **flux_m_l2;
     flux_p_l2 = array_from_mem2_d(flux_p_l, d->n_ulevels, d->n_derivs);
     flux_m_l2 = array_from_mem2_d(flux_m_l, d->n_ulevels, d->n_derivs);
     if (xrtm_flux(d, solver, flux_p, flux_m, flux_p_l2, flux_m_l2)) {
          fprintf(stderr, "ERROR: xrtm_flux()\n");
          return -1;
     }
     free_array2_d(flux_p_l2);
     free_array2_d(flux_m_l2);
     return 0;
}


int xrtm_flux_divergence_bindx_f90(xrtm_data *d, enum xrtm_solver_mask solver, double *flux_div, double *flux_div_l)
{
     double **flux_div_l2;
     flux_div_l2 = array_from_mem2_d(flux_div_l, d->n_ulevels, d->n_derivs);
     if (xrtm_flux_divergence(d, solver, flux_div, flux_div_l2)) {
          fprintf(stderr, "ERROR: xrtm_flux_divergence()\n");
          return -1;
     }
     free_array2_d(flux_div_l2);
     return 0;
}


#ifdef __cplusplus
}
#endif
