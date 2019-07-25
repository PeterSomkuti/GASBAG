/* xrtm_matrix.c */
void dvec_scale_add(double *a, double *b, double *c, double alpha, long n);
void dmat_scale_trans(double **a, double **b, double alpha, long m, long n);
void dmat_scale_add(double **a, double **b, double **c, double alpha, long m, long n);
double dmat_p_one_norm_A(double **r, double **t, long n);
double dmat_p_inf_norm_A(double **r, double **t, long n);
double dmat_frob_norm_A(double **r, double **t, long n);
int positive_definite(double **a, int n);
void dmat_sym(double **a, long m, double alpha);
void dsym_add(double **a, double **b, double **c, long m, int flag);
void dsym_sub(double **a, double **b, double **c, long m, int flag);
void dsym_mul_kernel(double **a, double **b, double alpha, long m, long o, double **c, double beta, int flag, int thresh);
void dsym_mul(double **a, double **b, long m, long o, double **c, int flag);
void dsym_mma(double **a1, double **b1, double **a2, double **b2, long m, long o, double **c, int flag);
void dsym_mms(double **a1, double **b1, double **a2, double **b2, long m, long o, double **c, int flag);
void dsym_m_a(double **a1, double **b1, double **a2, long m, long o, double **c, int flag);
void dsym_m_s(double **a1, double **b1, double **a2, long m, long o, double **c, int flag);
void dsym_mta(double **a, double **b, long m, long o, double **c);
void dsym_mts(double **a, double **b, long m, long o, double **c);
void zmat_real(dcomplex **z, double **d, long m, long n);
void zmat_imag(dcomplex **z, double **d, long m, long n);
void dmat_complex(double **dr, double **di, dcomplex **z, long m, long n);
void dzmat_diag_mul(double *a, dcomplex **b, dcomplex **c, long m, long n);
void dzmat_mul(double **a, dcomplex **b, long m, long n, long o, dcomplex **c);
void dzmat_mul2(double **a, dcomplex **b, long m, long n, long o, dcomplex **c, double **don, double **dmn1, double **dmn2);
void copy_to_band_storage_d(double **ab, double **a, double alpha, int m, int n, int kl, int ku, int i1, int j1);
void copy_to_band_storage2_d(double **ab, double **a, double alpha, double *beta, int m, int n, int kl, int ku, int i1, int j1);
void copy_to_band_storage_dc(dcomplex **ab, dcomplex **a, dcomplex alpha, int m, int n, int kl, int ku, int i1, int j1);
void copy_to_band_storage2_dc(dcomplex **ab, dcomplex **a, dcomplex alpha, dcomplex *beta, int m, int n, int kl, int ku, int i1, int j1);