/* header file for clapack 3.2.1 */
#pragma once

#include "matcl-blas-lapack/blas_loader/clapack_types.h"

// based on clapack.h but only lapack declarations & different namespace

namespace raw_blas_lapack
{

#ifdef __cplusplus 	
extern "C" {	
#endif		

/* Subroutine */ int cbdsqr_(char *uplo, i_type_wr *n, i_type_wr *ncvt, i_type_wr *
	nru, i_type_wr *ncc, s_type_wr *d__, s_type_wr *e, c_type_wr *vt, i_type_wr *ldvt, 
	c_type_wr *u, i_type_wr *ldu, c_type_wr *c__, i_type_wr *ldc, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cgbbrd_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *ncc, 
	 i_type_wr *kl, i_type_wr *ku, c_type_wr *ab, i_type_wr *ldab, s_type_wr *d__, 
	s_type_wr *e, c_type_wr *q, i_type_wr *ldq, c_type_wr *pt, i_type_wr *ldpt, 
	c_type_wr *c__, i_type_wr *ldc, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgbcon_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 c_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, s_type_wr *anorm, s_type_wr *rcond, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgbequ_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 c_type_wr *ab, i_type_wr *ldab, s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr 
	*colcnd, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int cgbequb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, c_type_wr *ab, i_type_wr *ldab, s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, 
	s_type_wr *colcnd, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int cgbrfs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, i_type_wr *
	ldafb, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *
	ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int cgbrfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *
	afb, i_type_wr *ldafb, i_type_wr *ipiv, s_type_wr *r__, s_type_wr *c__, c_type_wr *b, 
	i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *berr, 
	i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *err_bnds_comp__, 
	i_type_wr *nparams, s_type_wr *params, c_type_wr *work, s_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int cgbsv_(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	nrhs, c_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, c_type_wr *b, i_type_wr *
	ldb, i_type_wr *info);

/* Subroutine */ int cgbsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *kl, 
	 i_type_wr *ku, i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, 
	 i_type_wr *ldafb, i_type_wr *ipiv, char *equed, s_type_wr *r__, s_type_wr *c__, 
	c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr 
	*ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgbsvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *
	afb, i_type_wr *ldafb, i_type_wr *ipiv, char *equed, s_type_wr *r__, s_type_wr *c__, 
	 c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, 
	s_type_wr *rpvgrw, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *
	err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *
	params, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgbtf2_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 c_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int cgbtrf_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 c_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int cgbtrs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, c_type_wr 
	*b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cgebak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, s_type_wr *scale, i_type_wr *m, c_type_wr *v, i_type_wr *ldv, 
	i_type_wr *info);

/* Subroutine */ int cgebal_(char *job, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	i_type_wr *ilo, i_type_wr *ihi, s_type_wr *scale, i_type_wr *info);

/* Subroutine */ int cgebd2_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *d__, s_type_wr *e, c_type_wr *tauq, c_type_wr *taup, c_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int cgebrd_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *d__, s_type_wr *e, c_type_wr *tauq, c_type_wr *taup, c_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgecon_(char *norm, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *anorm, s_type_wr *rcond, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgeequ_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *colcnd, s_type_wr *amax, 
	i_type_wr *info);

/* Subroutine */ int cgeequb_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *
	lda, s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *colcnd, s_type_wr *amax, 
	i_type_wr *info);

/* Subroutine */ int cgees_(char *jobvs, char *sort, sel_fun_wr select, i_type_wr *n, 
	c_type_wr *a, i_type_wr *lda, i_type_wr *sdim, c_type_wr *w, c_type_wr *vs, 
	i_type_wr *ldvs, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, l_type_wr *
	bwork, i_type_wr *info);

/* Subroutine */ int cgeesx_(char *jobvs, char *sort, sel_fun_wr select, char *
	sense, i_type_wr *n, c_type_wr *a, i_type_wr *lda, i_type_wr *sdim, c_type_wr *
	w, c_type_wr *vs, i_type_wr *ldvs, s_type_wr *rconde, s_type_wr *rcondv, c_type_wr *
	work, i_type_wr *lwork, s_type_wr *rwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int cgeev_(char *jobvl, char *jobvr, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *w, c_type_wr *vl, i_type_wr *ldvl, c_type_wr *vr, 
	i_type_wr *ldvr, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int cgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *w, c_type_wr *vl, 
	i_type_wr *ldvl, c_type_wr *vr, i_type_wr *ldvr, i_type_wr *ilo, i_type_wr *ihi, 
	 s_type_wr *scale, s_type_wr *abnrm, s_type_wr *rconde, s_type_wr *rcondv, c_type_wr *work, 
	i_type_wr *lwork, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgegs_(char *jobvsl, char *jobvsr, i_type_wr *n, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *alpha, c_type_wr *
	beta, c_type_wr *vsl, i_type_wr *ldvsl, c_type_wr *vsr, i_type_wr *ldvsr, 
	c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgegv_(char *jobvl, char *jobvr, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *alpha, c_type_wr *beta, 
	 c_type_wr *vl, i_type_wr *ldvl, c_type_wr *vr, i_type_wr *ldvr, c_type_wr *
	work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgehd2_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, c_type_wr *
	a, i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cgehrd_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, c_type_wr *
	a, i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr 
	*info);

/* Subroutine */ int cgelq2_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cgelqf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgels_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgelsd_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *s, s_type_wr *rcond, 
	i_type_wr *rank, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int cgelss_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *s, s_type_wr *rcond, 
	i_type_wr *rank, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int cgelsx_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, i_type_wr *jpvt, s_type_wr *rcond, 
	 i_type_wr *rank, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgelsy_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, i_type_wr *jpvt, s_type_wr *rcond, 
	 i_type_wr *rank, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int cgeql2_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cgeqlf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgeqp3_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *jpvt, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, s_type_wr *
	rwork, i_type_wr *info);

/* Subroutine */ int cgeqpf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *jpvt, c_type_wr *tau, c_type_wr *work, s_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int cgeqr2_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cgeqrf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgerfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, c_type_wr *
	b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgerfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, s_type_wr *r__, s_type_wr *c__, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, 
	i_type_wr *ldx, s_type_wr *rcond, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *
	err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *
	params, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgerq2_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cgerqf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgesc2_(i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *
	rhs, i_type_wr *ipiv, i_type_wr *jpiv, s_type_wr *scale);

/* Subroutine */ int cgesdd_(char *jobz, i_type_wr *m, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, s_type_wr *s, c_type_wr *u, i_type_wr *ldu, c_type_wr *vt, i_type_wr 
	*ldvt, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int cgesv_(i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, i_type_wr *
	lda, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cgesvd_(char *jobu, char *jobvt, i_type_wr *m, i_type_wr *n, 
	c_type_wr *a, i_type_wr *lda, s_type_wr *s, c_type_wr *u, i_type_wr *ldu, c_type_wr *
	vt, i_type_wr *ldvt, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cgesvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, char *equed, s_type_wr *r__, s_type_wr *c__, c_type_wr *b, i_type_wr *ldb, 
	c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgesvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, char *equed, s_type_wr *r__, s_type_wr *c__, c_type_wr *b, i_type_wr *ldb, 
	c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *rpvgrw, s_type_wr *berr, 
	i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *err_bnds_comp__, 
	i_type_wr *nparams, s_type_wr *params, c_type_wr *work, s_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int cgetc2_(i_type_wr *n, c_type_wr *a, i_type_wr *lda, i_type_wr *
	ipiv, i_type_wr *jpiv, i_type_wr *info);

/* Subroutine */ int cgetf2_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int cgetrf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int cgetri_(i_type_wr *n, c_type_wr *a, i_type_wr *lda, i_type_wr *
	ipiv, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgetrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int cggbak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, s_type_wr *lscale, s_type_wr *rscale, i_type_wr *m, c_type_wr *v, 
	i_type_wr *ldv, i_type_wr *info);

/* Subroutine */ int cggbal_(char *job, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *b, i_type_wr *ldb, i_type_wr *ilo, i_type_wr *ihi, s_type_wr *lscale, 
	s_type_wr *rscale, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int cgges_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *
	ldb, i_type_wr *sdim, c_type_wr *alpha, c_type_wr *beta, c_type_wr *vsl, 
	i_type_wr *ldvsl, c_type_wr *vsr, i_type_wr *ldvsr, c_type_wr *work, i_type_wr *
	lwork, s_type_wr *rwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int cggesx_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, char *sense, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, 
	 i_type_wr *ldb, i_type_wr *sdim, c_type_wr *alpha, c_type_wr *beta, c_type_wr *
	vsl, i_type_wr *ldvsl, c_type_wr *vsr, i_type_wr *ldvsr, s_type_wr *rconde, s_type_wr 
	*rcondv, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *iwork, 
	i_type_wr *liwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int cggev_(char *jobvl, char *jobvr, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *alpha, c_type_wr *beta, 
	 c_type_wr *vl, i_type_wr *ldvl, c_type_wr *vr, i_type_wr *ldvr, c_type_wr *
	work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
	 c_type_wr *alpha, c_type_wr *beta, c_type_wr *vl, i_type_wr *ldvl, c_type_wr *
	vr, i_type_wr *ldvr, i_type_wr *ilo, i_type_wr *ihi, s_type_wr *lscale, s_type_wr *
	rscale, s_type_wr *abnrm, s_type_wr *bbnrm, s_type_wr *rconde, s_type_wr *rcondv, c_type_wr 
	*work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *iwork, l_type_wr *bwork, 
	i_type_wr *info);

/* Subroutine */ int cggglm_(i_type_wr *n, i_type_wr *m, i_type_wr *p, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *d__, c_type_wr *x, 
	c_type_wr *y, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cgghrd_(char *compq, char *compz, i_type_wr *n, i_type_wr *
	ilo, i_type_wr *ihi, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
	 c_type_wr *q, i_type_wr *ldq, c_type_wr *z__, i_type_wr *ldz, i_type_wr *info);

/* Subroutine */ int cgglse_(i_type_wr *m, i_type_wr *n, i_type_wr *p, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *c__, c_type_wr *d__, 
	c_type_wr *x, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cggqrf_(i_type_wr *n, i_type_wr *m, i_type_wr *p, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *taua, c_type_wr *b, i_type_wr *ldb, c_type_wr *taub, 
	c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cggrqf_(i_type_wr *m, i_type_wr *p, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *taua, c_type_wr *b, i_type_wr *ldb, c_type_wr *taub, 
	c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cggsvd_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *n, i_type_wr *p, i_type_wr *k, i_type_wr *l, c_type_wr *a, i_type_wr *
	lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *alpha, s_type_wr *beta, c_type_wr *u, 
	i_type_wr *ldu, c_type_wr *v, i_type_wr *ldv, c_type_wr *q, i_type_wr *ldq, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int cggsvp_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr 
	*ldb, s_type_wr *tola, s_type_wr *tolb, i_type_wr *k, i_type_wr *l, c_type_wr *u, 
	i_type_wr *ldu, c_type_wr *v, i_type_wr *ldv, c_type_wr *q, i_type_wr *ldq, 
	i_type_wr *iwork, s_type_wr *rwork, c_type_wr *tau, c_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int cgtcon_(char *norm, i_type_wr *n, c_type_wr *dl, c_type_wr *
	d__, c_type_wr *du, c_type_wr *du2, i_type_wr *ipiv, s_type_wr *anorm, s_type_wr *
	rcond, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cgtrfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	dl, c_type_wr *d__, c_type_wr *du, c_type_wr *dlf, c_type_wr *df, c_type_wr *
	duf, c_type_wr *du2, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *
	x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cgtsv_(i_type_wr *n, i_type_wr *nrhs, c_type_wr *dl, c_type_wr *
	d__, c_type_wr *du, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cgtsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *dl, c_type_wr *d__, c_type_wr *du, c_type_wr *dlf, c_type_wr *
	df, c_type_wr *duf, c_type_wr *du2, i_type_wr *ipiv, c_type_wr *b, i_type_wr *
	ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cgttrf_(i_type_wr *n, c_type_wr *dl, c_type_wr *d__, c_type_wr *
	du, c_type_wr *du2, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int cgttrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	dl, c_type_wr *d__, c_type_wr *du, c_type_wr *du2, i_type_wr *ipiv, c_type_wr *
	b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cgtts2_(i_type_wr *itrans, i_type_wr *n, i_type_wr *nrhs, 
	c_type_wr *dl, c_type_wr *d__, c_type_wr *du, c_type_wr *du2, i_type_wr *ipiv, 
	c_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int chbev_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	c_type_wr *ab, i_type_wr *ldab, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int chbevd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	c_type_wr *ab, i_type_wr *ldab, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, 
	c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *lrwork, i_type_wr *
	iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int chbevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *kd, c_type_wr *ab, i_type_wr *ldab, c_type_wr *q, i_type_wr *ldq, 
	s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, i_type_wr *
	m, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int chbgst_(char *vect, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, c_type_wr *ab, i_type_wr *ldab, c_type_wr *bb, i_type_wr *ldbb, 
	c_type_wr *x, i_type_wr *ldx, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int chbgv_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, c_type_wr *ab, i_type_wr *ldab, c_type_wr *bb, i_type_wr *ldbb, 
	s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int chbgvd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, c_type_wr *ab, i_type_wr *ldab, c_type_wr *bb, i_type_wr *ldbb, 
	s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, i_type_wr *lwork, 
	s_type_wr *rwork, i_type_wr *lrwork, i_type_wr *iwork, i_type_wr *liwork, 
	i_type_wr *info);

/* Subroutine */ int chbgvx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *ka, i_type_wr *kb, c_type_wr *ab, i_type_wr *ldab, c_type_wr *bb, 
	i_type_wr *ldbb, c_type_wr *q, i_type_wr *ldq, s_type_wr *vl, s_type_wr *vu, i_type_wr *
	il, i_type_wr *iu, s_type_wr *abstol, i_type_wr *m, s_type_wr *w, c_type_wr *z__, 
	i_type_wr *ldz, c_type_wr *work, s_type_wr *rwork, i_type_wr *iwork, i_type_wr *
	ifail, i_type_wr *info);

/* Subroutine */ int chbtrd_(char *vect, char *uplo, i_type_wr *n, i_type_wr *kd, 
	c_type_wr *ab, i_type_wr *ldab, s_type_wr *d__, s_type_wr *e, c_type_wr *q, i_type_wr *
	ldq, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int checon_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, s_type_wr *anorm, s_type_wr *rcond, c_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int cheequb_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *
	lda, s_type_wr *s, s_type_wr *scond, s_type_wr *amax, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cheev_(char *jobz, char *uplo, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, s_type_wr *w, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cheevd_(char *jobz, char *uplo, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, s_type_wr *w, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, 
	i_type_wr *lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int cheevr_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	c_type_wr *a, i_type_wr *lda, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *
	iu, s_type_wr *abstol, i_type_wr *m, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, 
	i_type_wr *isuppz, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *
	lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int cheevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	c_type_wr *a, i_type_wr *lda, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *
	iu, s_type_wr *abstol, i_type_wr *m, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, 
	c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *iwork, i_type_wr *
	ifail, i_type_wr *info);

/* Subroutine */ int chegs2_(i_type_wr *itype, char *uplo, i_type_wr *n, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int chegst_(i_type_wr *itype, char *uplo, i_type_wr *n, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int chegv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *w, 
	c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int chegvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *w, 
	c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *lrwork, i_type_wr *
	iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int chegvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
	s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, i_type_wr *
	m, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, i_type_wr *lwork, 
	 s_type_wr *rwork, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int cherfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, c_type_wr *
	b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cherfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, s_type_wr *s, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, 
	s_type_wr *rcond, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, 
	 s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *params, c_type_wr *work, 
	 s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int chesv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, 
	 i_type_wr *lda, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *work, 
	 i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int chesvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, 
	 s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int chesvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, char *equed, s_type_wr *s, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, 
	i_type_wr *ldx, s_type_wr *rcond, s_type_wr *rpvgrw, s_type_wr *berr, i_type_wr *
	n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *
	nparams, s_type_wr *params, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int chetd2_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *d__, s_type_wr *e, c_type_wr *tau, i_type_wr *info);

/* Subroutine */ int chetf2_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int chetrd_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *d__, s_type_wr *e, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int chetrf_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int chetri_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int chetrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int chfrk_(char *transr, char *uplo, char *trans, i_type_wr *n, 
	 i_type_wr *k, s_type_wr *alpha, c_type_wr *a, i_type_wr *lda, s_type_wr *beta, 
	c_type_wr *c__);

/* Subroutine */ int chgeqz_(char *job, char *compq, char *compz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, c_type_wr *h__, i_type_wr *ldh, c_type_wr *t, 
	i_type_wr *ldt, c_type_wr *alpha, c_type_wr *beta, c_type_wr *q, i_type_wr *ldq, 
	 c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, i_type_wr *lwork, s_type_wr *
	rwork, i_type_wr *info);

/* Character */ void chla_transtype__(char *ret_val, ftn_len_wr ret_val_len, 
	i_type_wr *trans);

/* Subroutine */ int chpcon_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	ipiv, s_type_wr *anorm, s_type_wr *rcond, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int chpev_(char *jobz, char *uplo, i_type_wr *n, c_type_wr *ap, 
	s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int chpevd_(char *jobz, char *uplo, i_type_wr *n, c_type_wr *ap, 
	s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, i_type_wr *lwork, 
	s_type_wr *rwork, i_type_wr *lrwork, i_type_wr *iwork, i_type_wr *liwork, 
	i_type_wr *info);

/* Subroutine */ int chpevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	c_type_wr *ap, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *
	abstol, i_type_wr *m, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, c_type_wr *
	work, s_type_wr *rwork, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int chpgst_(i_type_wr *itype, char *uplo, i_type_wr *n, c_type_wr *
	ap, c_type_wr *bp, i_type_wr *info);

/* Subroutine */ int chpgv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, c_type_wr *ap, c_type_wr *bp, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int chpgvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, c_type_wr *ap, c_type_wr *bp, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, 
	c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, i_type_wr *lrwork, i_type_wr *
	iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int chpgvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, c_type_wr *ap, c_type_wr *bp, s_type_wr *vl, s_type_wr *vu, 
	i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, i_type_wr *m, s_type_wr *w, c_type_wr *
	z__, i_type_wr *ldz, c_type_wr *work, s_type_wr *rwork, i_type_wr *iwork, 
	i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int chprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, c_type_wr *afp, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, 
	 i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int chpsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int chpsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *ap, c_type_wr *afp, i_type_wr *ipiv, c_type_wr *b, i_type_wr *
	ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int chptrd_(char *uplo, i_type_wr *n, c_type_wr *ap, s_type_wr *d__, 
	s_type_wr *e, c_type_wr *tau, i_type_wr *info);

/* Subroutine */ int chptrf_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	ipiv, i_type_wr *info);

/* Subroutine */ int chptri_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	ipiv, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int chptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int chsein_(char *side, char *eigsrc, char *initv, l_type_wr *
	select, i_type_wr *n, c_type_wr *h__, i_type_wr *ldh, c_type_wr *w, c_type_wr *
	vl, i_type_wr *ldvl, c_type_wr *vr, i_type_wr *ldvr, i_type_wr *mm, i_type_wr *
	m, c_type_wr *work, s_type_wr *rwork, i_type_wr *ifaill, i_type_wr *ifailr, 
	i_type_wr *info);

/* Subroutine */ int chseqr_(char *job, char *compz, i_type_wr *n, i_type_wr *ilo, 
	 i_type_wr *ihi, c_type_wr *h__, i_type_wr *ldh, c_type_wr *w, c_type_wr *z__, 
	i_type_wr *ldz, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cla_gbamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *kl, i_type_wr *ku, s_type_wr *alpha, c_type_wr *ab, i_type_wr *ldab, 
	c_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy);

d_type_wr cla_gbrcond_c__(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, i_type_wr *ldafb, i_type_wr *
	ipiv, s_type_wr *c__, l_type_wr *capply, i_type_wr *info, c_type_wr *work, s_type_wr *
	rwork, ftn_len_wr trans_len);

d_type_wr cla_gbrcond_x__(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, i_type_wr *ldafb, i_type_wr *
	ipiv, c_type_wr *x, i_type_wr *info, c_type_wr *work, s_type_wr *rwork, ftn_len_wr 
	trans_len);

/* Subroutine */ int cla_gbrfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *nrhs, 
	c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, i_type_wr *ldafb, i_type_wr *
	ipiv, l_type_wr *colequ, s_type_wr *c__, c_type_wr *b, i_type_wr *ldb, c_type_wr *
	y, i_type_wr *ldy, s_type_wr *berr_out__, i_type_wr *n_norms__, s_type_wr *errs_n__,
	 s_type_wr *errs_c__, c_type_wr *res, s_type_wr *ayb, c_type_wr *dy, c_type_wr *
	y_tail__, s_type_wr *rcond, i_type_wr *ithresh, s_type_wr *rthresh, s_type_wr *dz_ub__,
	 l_type_wr *ignore_cwise__, i_type_wr *info);

d_type_wr cla_gbrpvgrw__(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	ncols, c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, i_type_wr *ldafb);

/* Subroutine */ int cla_geamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, s_type_wr 
	*alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, s_type_wr *
	beta, s_type_wr *y, i_type_wr *incy);

d_type_wr cla_gercond_c__(char *trans, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *c__, l_type_wr *capply,
	 i_type_wr *info, c_type_wr *work, s_type_wr *rwork, ftn_len_wr trans_len);

d_type_wr cla_gercond_x__(char *trans, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, c_type_wr *x, i_type_wr *info, 
	c_type_wr *work, s_type_wr *rwork, ftn_len_wr trans_len);

/* Subroutine */ int cla_gerfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, s_type_wr *c__,
	 c_type_wr *b, i_type_wr *ldb, c_type_wr *y, i_type_wr *ldy, s_type_wr *berr_out__,
	 i_type_wr *n_norms__, s_type_wr *errs_n__, s_type_wr *errs_c__, c_type_wr *res, 
	s_type_wr *ayb, c_type_wr *dy, c_type_wr *y_tail__, s_type_wr *rcond, i_type_wr *
	ithresh, s_type_wr *rthresh, s_type_wr *dz_ub__, l_type_wr *ignore_cwise__, 
	i_type_wr *info);

/* Subroutine */ int cla_heamv__(i_type_wr *uplo, i_type_wr *n, s_type_wr *alpha, 
	c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr 
	*y, i_type_wr *incy);

d_type_wr cla_hercond_c__(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *c__, l_type_wr *capply,
	 i_type_wr *info, c_type_wr *work, s_type_wr *rwork, ftn_len_wr uplo_len);

d_type_wr cla_hercond_x__(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, c_type_wr *x, i_type_wr *info, 
	c_type_wr *work, s_type_wr *rwork, ftn_len_wr uplo_len);

/* Subroutine */ int cla_herfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, 
	i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, s_type_wr *c__, c_type_wr *b, 
	i_type_wr *ldb, c_type_wr *y, i_type_wr *ldy, s_type_wr *berr_out__, i_type_wr *
	n_norms__, s_type_wr *errs_n__, s_type_wr *errs_c__, c_type_wr *res, s_type_wr *ayb, 
	c_type_wr *dy, c_type_wr *y_tail__, s_type_wr *rcond, i_type_wr *ithresh, s_type_wr *
	rthresh, s_type_wr *dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info, 
	ftn_len_wr uplo_len);

d_type_wr cla_herpvgrw__(char *uplo, i_type_wr *n, i_type_wr *info, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *work, 
	ftn_len_wr uplo_len);

/* Subroutine */ int cla_lin_berr__(i_type_wr *n, i_type_wr *nz, i_type_wr *nrhs, 
	c_type_wr *res, s_type_wr *ayb, s_type_wr *berr);

d_type_wr cla_porcond_c__(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, s_type_wr *c__, l_type_wr *capply, i_type_wr *info,
	 c_type_wr *work, s_type_wr *rwork, ftn_len_wr uplo_len);

d_type_wr cla_porcond_x__(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, c_type_wr *x, i_type_wr *info, c_type_wr *work, 
	s_type_wr *rwork, ftn_len_wr uplo_len);

/* Subroutine */ int cla_porfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, 
	i_type_wr *ldaf, l_type_wr *colequ, s_type_wr *c__, c_type_wr *b, i_type_wr *ldb, 
	c_type_wr *y, i_type_wr *ldy, s_type_wr *berr_out__, i_type_wr *n_norms__, s_type_wr *
	errs_n__, s_type_wr *errs_c__, c_type_wr *res, s_type_wr *ayb, c_type_wr *dy, 
	c_type_wr *y_tail__, s_type_wr *rcond, i_type_wr *ithresh, s_type_wr *rthresh, s_type_wr 
	*dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info, ftn_len_wr uplo_len);

d_type_wr cla_porpvgrw__(char *uplo, i_type_wr *ncols, c_type_wr *a, i_type_wr *
	lda, c_type_wr *af, i_type_wr *ldaf, s_type_wr *work, ftn_len_wr uplo_len);

d_type_wr cla_rpvgrw__(i_type_wr *n, i_type_wr *ncols, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf);

/* Subroutine */ int cla_syamv__(i_type_wr *uplo, i_type_wr *n, s_type_wr *alpha, 
	c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr 
	*y, i_type_wr *incy);

d_type_wr cla_syrcond_c__(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *c__, l_type_wr *capply,
	 i_type_wr *info, c_type_wr *work, s_type_wr *rwork, ftn_len_wr uplo_len);

d_type_wr cla_syrcond_x__(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, c_type_wr *x, i_type_wr *info, 
	c_type_wr *work, s_type_wr *rwork, ftn_len_wr uplo_len);

/* Subroutine */ int cla_syrfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, 
	i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, s_type_wr *c__, c_type_wr *b, 
	i_type_wr *ldb, c_type_wr *y, i_type_wr *ldy, s_type_wr *berr_out__, i_type_wr *
	n_norms__, s_type_wr *errs_n__, s_type_wr *errs_c__, c_type_wr *res, s_type_wr *ayb, 
	c_type_wr *dy, c_type_wr *y_tail__, s_type_wr *rcond, i_type_wr *ithresh, s_type_wr *
	rthresh, s_type_wr *dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info, 
	ftn_len_wr uplo_len);

d_type_wr cla_syrpvgrw__(char *uplo, i_type_wr *n, i_type_wr *info, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *work, 
	ftn_len_wr uplo_len);

/* Subroutine */ int cla_wwaddw__(i_type_wr *n, c_type_wr *x, c_type_wr *y, c_type_wr 
	*w);

/* Subroutine */ int clabrd_(i_type_wr *m, i_type_wr *n, i_type_wr *nb, c_type_wr *a, 
	i_type_wr *lda, s_type_wr *d__, s_type_wr *e, c_type_wr *tauq, c_type_wr *taup, 
	c_type_wr *x, i_type_wr *ldx, c_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int clacgv_(i_type_wr *n, c_type_wr *x, i_type_wr *incx);

/* Subroutine */ int clacn2_(i_type_wr *n, c_type_wr *v, c_type_wr *x, s_type_wr *est, 
	i_type_wr *kase, i_type_wr *isave);

/* Subroutine */ int clacon_(i_type_wr *n, c_type_wr *v, c_type_wr *x, s_type_wr *est, 
	i_type_wr *kase);

/* Subroutine */ int clacp2_(char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, c_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int clacpy_(char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int clacrm_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *b, i_type_wr *ldb, c_type_wr *c__, i_type_wr *ldc, s_type_wr *rwork);

/* Subroutine */ int clacrt_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *
	cy, i_type_wr *incy, c_type_wr *c__, c_type_wr *s);

/* c_type_wr */ void cladiv_(c_type_wr * ret_val, c_type_wr *x, c_type_wr *y);

/* Subroutine */ int claed0_(i_type_wr *qsiz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, 
	c_type_wr *q, i_type_wr *ldq, c_type_wr *qstore, i_type_wr *ldqs, s_type_wr *rwork, 
	 i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int claed7_(i_type_wr *n, i_type_wr *cutpnt, i_type_wr *qsiz, 
	i_type_wr *tlvls, i_type_wr *curlvl, i_type_wr *curpbm, s_type_wr *d__, c_type_wr *
	q, i_type_wr *ldq, s_type_wr *rho, i_type_wr *indxq, s_type_wr *qstore, i_type_wr *
	qptr, i_type_wr *prmptr, i_type_wr *perm, i_type_wr *givptr, i_type_wr *
	givcol, s_type_wr *givnum, c_type_wr *work, s_type_wr *rwork, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int claed8_(i_type_wr *k, i_type_wr *n, i_type_wr *qsiz, c_type_wr *
	q, i_type_wr *ldq, s_type_wr *d__, s_type_wr *rho, i_type_wr *cutpnt, s_type_wr *z__, 
	s_type_wr *dlamda, c_type_wr *q2, i_type_wr *ldq2, s_type_wr *w, i_type_wr *indxp, 
	i_type_wr *indx, i_type_wr *indxq, i_type_wr *perm, i_type_wr *givptr, 
	i_type_wr *givcol, s_type_wr *givnum, i_type_wr *info);

/* Subroutine */ int claein_(l_type_wr *rightv, l_type_wr *noinit, i_type_wr *n, 
	c_type_wr *h__, i_type_wr *ldh, c_type_wr *w, c_type_wr *v, c_type_wr *b, 
	i_type_wr *ldb, s_type_wr *rwork, s_type_wr *eps3, s_type_wr *smlnum, i_type_wr *info);

/* Subroutine */ int claesy_(c_type_wr *a, c_type_wr *b, c_type_wr *c__, c_type_wr *
	rt1, c_type_wr *rt2, c_type_wr *evscal, c_type_wr *cs1, c_type_wr *sn1);

/* Subroutine */ int claev2_(c_type_wr *a, c_type_wr *b, c_type_wr *c__, s_type_wr *rt1, 
	s_type_wr *rt2, s_type_wr *cs1, c_type_wr *sn1);

/* Subroutine */ int clag2z_(i_type_wr *m, i_type_wr *n, c_type_wr *sa, i_type_wr *
	ldsa, z_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int clags2_(l_type_wr *upper, s_type_wr *a1, c_type_wr *a2, s_type_wr *a3, 
	s_type_wr *b1, c_type_wr *b2, s_type_wr *b3, s_type_wr *csu, c_type_wr *snu, s_type_wr *csv, 
	c_type_wr *snv, s_type_wr *csq, c_type_wr *snq);

/* Subroutine */ int clagtm_(char *trans, i_type_wr *n, i_type_wr *nrhs, s_type_wr *
	alpha, c_type_wr *dl, c_type_wr *d__, c_type_wr *du, c_type_wr *x, i_type_wr *
	ldx, s_type_wr *beta, c_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int clahef_(char *uplo, i_type_wr *n, i_type_wr *nb, i_type_wr *kb, 
	 c_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, c_type_wr *w, i_type_wr *ldw, 
	i_type_wr *info);

/* Subroutine */ int clahqr_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, c_type_wr *h__, i_type_wr *ldh, c_type_wr *w, 
	i_type_wr *iloz, i_type_wr *ihiz, c_type_wr *z__, i_type_wr *ldz, i_type_wr *
	info);

/* Subroutine */ int clahr2_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *t, i_type_wr *ldt, c_type_wr *y, 
	i_type_wr *ldy);

/* Subroutine */ int clahrd_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *t, i_type_wr *ldt, c_type_wr *y, 
	i_type_wr *ldy);

/* Subroutine */ int claic1_(i_type_wr *job, i_type_wr *j, c_type_wr *x, s_type_wr *sest, 
	 c_type_wr *w, c_type_wr *gamma, s_type_wr *sestpr, c_type_wr *s, c_type_wr *c__);

/* Subroutine */ int clals0_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, i_type_wr *nrhs, c_type_wr *b, i_type_wr *ldb, c_type_wr *bx, 
	i_type_wr *ldbx, i_type_wr *perm, i_type_wr *givptr, i_type_wr *givcol, 
	i_type_wr *ldgcol, s_type_wr *givnum, i_type_wr *ldgnum, s_type_wr *poles, s_type_wr *
	difl, s_type_wr *difr, s_type_wr *z__, i_type_wr *k, s_type_wr *c__, s_type_wr *s, s_type_wr *
	rwork, i_type_wr *info);

/* Subroutine */ int clalsa_(i_type_wr *icompq, i_type_wr *smlsiz, i_type_wr *n, 
	i_type_wr *nrhs, c_type_wr *b, i_type_wr *ldb, c_type_wr *bx, i_type_wr *ldbx, 
	s_type_wr *u, i_type_wr *ldu, s_type_wr *vt, i_type_wr *k, s_type_wr *difl, s_type_wr *difr, 
	s_type_wr *z__, s_type_wr *poles, i_type_wr *givptr, i_type_wr *givcol, i_type_wr *
	ldgcol, i_type_wr *perm, s_type_wr *givnum, s_type_wr *c__, s_type_wr *s, s_type_wr *rwork, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int clalsd_(char *uplo, i_type_wr *smlsiz, i_type_wr *n, i_type_wr 
	*nrhs, s_type_wr *d__, s_type_wr *e, c_type_wr *b, i_type_wr *ldb, s_type_wr *rcond, 
	i_type_wr *rank, c_type_wr *work, s_type_wr *rwork, i_type_wr *iwork, i_type_wr *
	info);

s_type_wr clangb_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, c_type_wr *
	ab, i_type_wr *ldab, s_type_wr *work);

s_type_wr clange_(char *norm, i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *
	lda, s_type_wr *work);

s_type_wr clangt_(char *norm, i_type_wr *n, c_type_wr *dl, c_type_wr *d__, c_type_wr 
	*du);

s_type_wr clanhb_(char *norm, char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *
	ab, i_type_wr *ldab, s_type_wr *work);

s_type_wr clanhe_(char *norm, char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *
	lda, s_type_wr *work);

s_type_wr clanhf_(char *norm, char *transr, char *uplo, i_type_wr *n, c_type_wr *
	a, s_type_wr *work);

s_type_wr clanhp_(char *norm, char *uplo, i_type_wr *n, c_type_wr *ap, s_type_wr *
	work);

s_type_wr clanhs_(char *norm, i_type_wr *n, c_type_wr *a, i_type_wr *lda, s_type_wr *
	work);

s_type_wr clanht_(char *norm, i_type_wr *n, s_type_wr *d__, c_type_wr *e);

s_type_wr clansb_(char *norm, char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *
	ab, i_type_wr *ldab, s_type_wr *work);

s_type_wr clansp_(char *norm, char *uplo, i_type_wr *n, c_type_wr *ap, s_type_wr *
	work);

s_type_wr clansy_(char *norm, char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *
	lda, s_type_wr *work);

s_type_wr clantb_(char *norm, char *uplo, char *diag, i_type_wr *n, i_type_wr *k, 
	 c_type_wr *ab, i_type_wr *ldab, s_type_wr *work);

s_type_wr clantp_(char *norm, char *uplo, char *diag, i_type_wr *n, c_type_wr *
	ap, s_type_wr *work);

s_type_wr clantr_(char *norm, char *uplo, char *diag, i_type_wr *m, i_type_wr *n, 
	 c_type_wr *a, i_type_wr *lda, s_type_wr *work);

/* Subroutine */ int clapll_(i_type_wr *n, c_type_wr *x, i_type_wr *incx, c_type_wr *
	y, i_type_wr *incy, s_type_wr *ssmin);

/* Subroutine */ int clapmt_(l_type_wr *forwrd, i_type_wr *m, i_type_wr *n, c_type_wr 
	*x, i_type_wr *ldx, i_type_wr *k);

/* Subroutine */ int claqgb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 c_type_wr *ab, i_type_wr *ldab, s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr 
	*colcnd, s_type_wr *amax, char *equed);

/* Subroutine */ int claqge_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *colcnd, s_type_wr *amax, char *
	equed);

/* Subroutine */ int claqhb_(char *uplo, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, 
	 i_type_wr *ldab, s_type_wr *s, s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int claqhe_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *s, s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int claqhp_(char *uplo, i_type_wr *n, c_type_wr *ap, s_type_wr *s, 
	s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int claqp2_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, c_type_wr 
	*a, i_type_wr *lda, i_type_wr *jpvt, c_type_wr *tau, s_type_wr *vn1, s_type_wr *vn2, 
	c_type_wr *work);

/* Subroutine */ int claqps_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, i_type_wr 
	*nb, i_type_wr *kb, c_type_wr *a, i_type_wr *lda, i_type_wr *jpvt, c_type_wr *
	tau, s_type_wr *vn1, s_type_wr *vn2, c_type_wr *auxv, c_type_wr *f, i_type_wr *ldf);

/* Subroutine */ int claqr0_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, c_type_wr *h__, i_type_wr *ldh, c_type_wr *w, 
	i_type_wr *iloz, i_type_wr *ihiz, c_type_wr *z__, i_type_wr *ldz, c_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int claqr1_(i_type_wr *n, c_type_wr *h__, i_type_wr *ldh, c_type_wr *
	s1, c_type_wr *s2, c_type_wr *v);

/* Subroutine */ int claqr2_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, c_type_wr *h__, i_type_wr *ldh, 
	 i_type_wr *iloz, i_type_wr *ihiz, c_type_wr *z__, i_type_wr *ldz, i_type_wr *
	ns, i_type_wr *nd, c_type_wr *sh, c_type_wr *v, i_type_wr *ldv, i_type_wr *nh, 
	c_type_wr *t, i_type_wr *ldt, i_type_wr *nv, c_type_wr *wv, i_type_wr *ldwv, 
	c_type_wr *work, i_type_wr *lwork);

/* Subroutine */ int claqr3_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, c_type_wr *h__, i_type_wr *ldh, 
	 i_type_wr *iloz, i_type_wr *ihiz, c_type_wr *z__, i_type_wr *ldz, i_type_wr *
	ns, i_type_wr *nd, c_type_wr *sh, c_type_wr *v, i_type_wr *ldv, i_type_wr *nh, 
	c_type_wr *t, i_type_wr *ldt, i_type_wr *nv, c_type_wr *wv, i_type_wr *ldwv, 
	c_type_wr *work, i_type_wr *lwork);

/* Subroutine */ int claqr4_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, c_type_wr *h__, i_type_wr *ldh, c_type_wr *w, 
	i_type_wr *iloz, i_type_wr *ihiz, c_type_wr *z__, i_type_wr *ldz, c_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int claqr5_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *kacc22, 
	i_type_wr *n, i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nshfts, c_type_wr *s, 
	 c_type_wr *h__, i_type_wr *ldh, i_type_wr *iloz, i_type_wr *ihiz, c_type_wr *
	z__, i_type_wr *ldz, c_type_wr *v, i_type_wr *ldv, c_type_wr *u, i_type_wr *ldu, 
	 i_type_wr *nv, c_type_wr *wv, i_type_wr *ldwv, i_type_wr *nh, c_type_wr *wh, 
	i_type_wr *ldwh);

/* Subroutine */ int claqsb_(char *uplo, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, 
	 i_type_wr *ldab, s_type_wr *s, s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int claqsp_(char *uplo, i_type_wr *n, c_type_wr *ap, s_type_wr *s, 
	s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int claqsy_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *s, s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int clar1v_(i_type_wr *n, i_type_wr *b1, i_type_wr *bn, s_type_wr *
	lambda, s_type_wr *d__, s_type_wr *l, s_type_wr *ld, s_type_wr *lld, s_type_wr *pivmin, s_type_wr *
	gaptol, c_type_wr *z__, l_type_wr *wantnc, i_type_wr *negcnt, s_type_wr *ztz, 
	s_type_wr *mingma, i_type_wr *r__, i_type_wr *isuppz, s_type_wr *nrminv, s_type_wr *
	resid, s_type_wr *rqcorr, s_type_wr *work);

/* Subroutine */ int clar2v_(i_type_wr *n, c_type_wr *x, c_type_wr *y, c_type_wr *z__, 
	 i_type_wr *incx, s_type_wr *c__, c_type_wr *s, i_type_wr *incc);

/* Subroutine */ int clarcm_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	c_type_wr *b, i_type_wr *ldb, c_type_wr *c__, i_type_wr *ldc, s_type_wr *rwork);

/* Subroutine */ int clarf_(char *side, i_type_wr *m, i_type_wr *n, c_type_wr *v, 
	i_type_wr *incv, c_type_wr *tau, c_type_wr *c__, i_type_wr *ldc, c_type_wr *
	work);

/* Subroutine */ int clarfb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *v, i_type_wr *ldv, 
	c_type_wr *t, i_type_wr *ldt, c_type_wr *c__, i_type_wr *ldc, c_type_wr *work, 
	i_type_wr *ldwork);

/* Subroutine */ int clarfg_(i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *
	incx, c_type_wr *tau);

/* Subroutine */ int clarfp_(i_type_wr *n, c_type_wr *alpha, c_type_wr *x, i_type_wr *
	incx, c_type_wr *tau);

/* Subroutine */ int clarft_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, c_type_wr *v, i_type_wr *ldv, c_type_wr *tau, c_type_wr *t, i_type_wr *ldt);

/* Subroutine */ int clarfx_(char *side, i_type_wr *m, i_type_wr *n, c_type_wr *v, 
	c_type_wr *tau, c_type_wr *c__, i_type_wr *ldc, c_type_wr *work);

/* Subroutine */ int clargv_(i_type_wr *n, c_type_wr *x, i_type_wr *incx, c_type_wr *
	y, i_type_wr *incy, s_type_wr *c__, i_type_wr *incc);

/* Subroutine */ int clarnv_(i_type_wr *idist, i_type_wr *iseed, i_type_wr *n, 
	c_type_wr *x);

/* Subroutine */ int clarrv_(i_type_wr *n, s_type_wr *vl, s_type_wr *vu, s_type_wr *d__, s_type_wr *
	l, s_type_wr *pivmin, i_type_wr *isplit, i_type_wr *m, i_type_wr *dol, i_type_wr *
	dou, s_type_wr *minrgp, s_type_wr *rtol1, s_type_wr *rtol2, s_type_wr *w, s_type_wr *werr, 
	s_type_wr *wgap, i_type_wr *iblock, i_type_wr *indexw, s_type_wr *gers, c_type_wr *
	z__, i_type_wr *ldz, i_type_wr *isuppz, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int clarscl2_(i_type_wr *m, i_type_wr *n, s_type_wr *d__, c_type_wr *x, 
	i_type_wr *ldx);

/* Subroutine */ int clartg_(c_type_wr *f, c_type_wr *g, s_type_wr *cs, c_type_wr *sn, 
	c_type_wr *r__);

/* Subroutine */ int clartv_(i_type_wr *n, c_type_wr *x, i_type_wr *incx, c_type_wr *
	y, i_type_wr *incy, s_type_wr *c__, c_type_wr *s, i_type_wr *incc);

/* Subroutine */ int clarz_(char *side, i_type_wr *m, i_type_wr *n, i_type_wr *l, 
	c_type_wr *v, i_type_wr *incv, c_type_wr *tau, c_type_wr *c__, i_type_wr *ldc, 
	c_type_wr *work);

/* Subroutine */ int clarzb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, i_type_wr *l, c_type_wr *v, 
	i_type_wr *ldv, c_type_wr *t, i_type_wr *ldt, c_type_wr *c__, i_type_wr *ldc, 
	c_type_wr *work, i_type_wr *ldwork);

/* Subroutine */ int clarzt_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, c_type_wr *v, i_type_wr *ldv, c_type_wr *tau, c_type_wr *t, i_type_wr *ldt);

/* Subroutine */ int clascl_(char *type__, i_type_wr *kl, i_type_wr *ku, s_type_wr *
	cfrom, s_type_wr *cto, i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	i_type_wr *info);

/* Subroutine */ int clascl2_(i_type_wr *m, i_type_wr *n, s_type_wr *d__, c_type_wr *x, 
	i_type_wr *ldx);

/* Subroutine */ int claset_(char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *
	alpha, c_type_wr *beta, c_type_wr *a, i_type_wr *lda);

/* Subroutine */ int clasr_(char *side, char *pivot, char *direct, i_type_wr *m, 
	 i_type_wr *n, s_type_wr *c__, s_type_wr *s, c_type_wr *a, i_type_wr *lda);

/* Subroutine */ int classq_(i_type_wr *n, c_type_wr *x, i_type_wr *incx, s_type_wr *
	scale, s_type_wr *sumsq);

/* Subroutine */ int claswp_(i_type_wr *n, c_type_wr *a, i_type_wr *lda, i_type_wr *
	k1, i_type_wr *k2, i_type_wr *ipiv, i_type_wr *incx);

/* Subroutine */ int clasyf_(char *uplo, i_type_wr *n, i_type_wr *nb, i_type_wr *kb, 
	 c_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, c_type_wr *w, i_type_wr *ldw, 
	i_type_wr *info);

/* Subroutine */ int clatbs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, i_type_wr *ldab, c_type_wr *
	x, s_type_wr *scale, s_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int clatdf_(i_type_wr *ijob, i_type_wr *n, c_type_wr *z__, i_type_wr 
	*ldz, c_type_wr *rhs, s_type_wr *rdsum, s_type_wr *rdscal, i_type_wr *ipiv, i_type_wr 
	*jpiv);

/* Subroutine */ int clatps_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, c_type_wr *ap, c_type_wr *x, s_type_wr *scale, s_type_wr *cnorm, 
	 i_type_wr *info);

/* Subroutine */ int clatrd_(char *uplo, i_type_wr *n, i_type_wr *nb, c_type_wr *a, 
	i_type_wr *lda, s_type_wr *e, c_type_wr *tau, c_type_wr *w, i_type_wr *ldw);

/* Subroutine */ int clatrs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *x, s_type_wr *scale, 
	 s_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int clatrz_(i_type_wr *m, i_type_wr *n, i_type_wr *l, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work);

/* Subroutine */ int clatzm_(char *side, i_type_wr *m, i_type_wr *n, c_type_wr *v, 
	i_type_wr *incv, c_type_wr *tau, c_type_wr *c1, c_type_wr *c2, i_type_wr *ldc, 
	c_type_wr *work);

/* Subroutine */ int clauu2_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *info);

/* Subroutine */ int clauum_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *info);

/* Subroutine */ int cpbcon_(char *uplo, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, 
	 i_type_wr *ldab, s_type_wr *anorm, s_type_wr *rcond, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cpbequ_(char *uplo, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, 
	 i_type_wr *ldab, s_type_wr *s, s_type_wr *scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int cpbrfs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, i_type_wr *ldafb, 
	c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *
	berr, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cpbstf_(char *uplo, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, 
	 i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int cpbsv_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int cpbsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *kd, 
	i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *afb, i_type_wr *
	ldafb, char *equed, s_type_wr *s, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, 
	i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, 
	s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cpbtf2_(char *uplo, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, 
	 i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int cpbtrf_(char *uplo, i_type_wr *n, i_type_wr *kd, c_type_wr *ab, 
	 i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int cpbtrs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int cpftrf_(char *transr, char *uplo, i_type_wr *n, c_type_wr *a, 
	 i_type_wr *info);

/* Subroutine */ int cpftri_(char *transr, char *uplo, i_type_wr *n, c_type_wr *a, 
	 i_type_wr *info);

/* Subroutine */ int cpftrs_(char *transr, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cpocon_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 s_type_wr *anorm, s_type_wr *rcond, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cpoequ_(i_type_wr *n, c_type_wr *a, i_type_wr *lda, s_type_wr *s, 
	s_type_wr *scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int cpoequb_(i_type_wr *n, c_type_wr *a, i_type_wr *lda, s_type_wr *s, 
	s_type_wr *scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int cporfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, c_type_wr *b, i_type_wr *ldb, 
	 c_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, 
	s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cporfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, s_type_wr *s, 
	c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr 
	*berr, i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *
	err_bnds_comp__, i_type_wr *nparams, s_type_wr *params, c_type_wr *work, s_type_wr *
	rwork, i_type_wr *info);

/* Subroutine */ int cposv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, 
	 i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cposvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, char *
	equed, s_type_wr *s, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, 
	s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cposvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, char *
	equed, s_type_wr *s, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, 
	s_type_wr *rcond, s_type_wr *rpvgrw, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *
	err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *
	params, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cpotf2_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *info);

/* Subroutine */ int cpotrf_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *info);

/* Subroutine */ int cpotri_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *info);

/* Subroutine */ int cpotrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cppcon_(char *uplo, i_type_wr *n, c_type_wr *ap, s_type_wr *anorm, 
	 s_type_wr *rcond, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cppequ_(char *uplo, i_type_wr *n, c_type_wr *ap, s_type_wr *s, 
	s_type_wr *scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int cpprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, c_type_wr *afp, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, 
	s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cppsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cppsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *ap, c_type_wr *afp, char *equed, s_type_wr *s, c_type_wr *b, 
	i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr 
	*berr, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cpptrf_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	info);

/* Subroutine */ int cpptri_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	info);

/* Subroutine */ int cpptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cpstf2_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *piv, i_type_wr *rank, s_type_wr *tol, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int cpstrf_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *piv, i_type_wr *rank, s_type_wr *tol, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int cptcon_(i_type_wr *n, s_type_wr *d__, c_type_wr *e, s_type_wr *anorm, 
	s_type_wr *rcond, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cpteqr_(char *compz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, 
	c_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int cptrfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, 
	 c_type_wr *e, s_type_wr *df, c_type_wr *ef, c_type_wr *b, i_type_wr *ldb, c_type_wr 
	*x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cptsv_(i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, c_type_wr *e, 
	c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cptsvx_(char *fact, i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, 
	 c_type_wr *e, s_type_wr *df, c_type_wr *ef, c_type_wr *b, i_type_wr *ldb, c_type_wr 
	*x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, 
	s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int cpttrf_(i_type_wr *n, s_type_wr *d__, c_type_wr *e, i_type_wr *info);

/* Subroutine */ int cpttrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, 
	 c_type_wr *e, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cptts2_(i_type_wr *iuplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *
	d__, c_type_wr *e, c_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int crot_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *
	cy, i_type_wr *incy, s_type_wr *c__, c_type_wr *s);

/* Subroutine */ int cspcon_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	ipiv, s_type_wr *anorm, s_type_wr *rcond, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cspmv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *
	ap, c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *
	incy);

/* Subroutine */ int cspr_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, 
	 i_type_wr *incx, c_type_wr *ap);

/* Subroutine */ int csprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, c_type_wr *afp, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, 
	 i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int cspsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int cspsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *ap, c_type_wr *afp, i_type_wr *ipiv, c_type_wr *b, i_type_wr *
	ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int csptrf_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	ipiv, i_type_wr *info);

/* Subroutine */ int csptri_(char *uplo, i_type_wr *n, c_type_wr *ap, i_type_wr *
	ipiv, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int csptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	ap, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int csrscl_(i_type_wr *n, s_type_wr *sa, c_type_wr *sx, i_type_wr *incx);

/* Subroutine */ int cstedc_(char *compz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, 
	c_type_wr *z__, i_type_wr *ldz, c_type_wr *work, i_type_wr *lwork, s_type_wr *
	rwork, i_type_wr *lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *
	info);

/* Subroutine */ int cstegr_(char *jobz, char *range, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, 
	i_type_wr *m, s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, 
	s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *
	info);

/* Subroutine */ int cstein_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, i_type_wr *m, s_type_wr 
	*w, i_type_wr *iblock, i_type_wr *isplit, c_type_wr *z__, i_type_wr *ldz, 
	s_type_wr *work, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int cstemr_(char *jobz, char *range, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, i_type_wr *m, 
	s_type_wr *w, c_type_wr *z__, i_type_wr *ldz, i_type_wr *nzc, i_type_wr *isuppz, 
	l_type_wr *tryrac, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *
	liwork, i_type_wr *info);

/* Subroutine */ int csteqr_(char *compz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, 
	c_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int csycon_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, s_type_wr *anorm, s_type_wr *rcond, c_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int csyequb_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *
	lda, s_type_wr *s, s_type_wr *scond, s_type_wr *amax, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int csymv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *
	a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y, 
	 i_type_wr *incy);

/* Subroutine */ int csyr_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, 
	 i_type_wr *incx, c_type_wr *a, i_type_wr *lda);

/* Subroutine */ int csyrfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, c_type_wr *
	b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int csyrfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, s_type_wr *s, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, 
	s_type_wr *rcond, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, 
	 s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *params, c_type_wr *work, 
	 s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int csysv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *a, 
	 i_type_wr *lda, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *work, 
	 i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int csysvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, 
	 s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, i_type_wr *lwork, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int csysvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, char *equed, s_type_wr *s, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, 
	i_type_wr *ldx, s_type_wr *rcond, s_type_wr *rpvgrw, s_type_wr *berr, i_type_wr *
	n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *
	nparams, s_type_wr *params, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int csytf2_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int csytrf_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int csytri_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 i_type_wr *ipiv, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int csytrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, c_type_wr *
	a, i_type_wr *lda, i_type_wr *ipiv, c_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int ctbcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	i_type_wr *kd, c_type_wr *ab, i_type_wr *ldab, s_type_wr *rcond, c_type_wr *work, 
	s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ctbrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *b, 
	i_type_wr *ldb, c_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, 
	c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ctbtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, c_type_wr *ab, i_type_wr *ldab, c_type_wr *b, 
	i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int ctfsm_(char *transr, char *side, char *uplo, char *trans, 
	 char *diag, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, 
	c_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int ctftri_(char *transr, char *uplo, char *diag, i_type_wr *n, 
	 c_type_wr *a, i_type_wr *info);

/* Subroutine */ int ctfttp_(char *transr, char *uplo, i_type_wr *n, c_type_wr *
	arf, c_type_wr *ap, i_type_wr *info);

/* Subroutine */ int ctfttr_(char *transr, char *uplo, i_type_wr *n, c_type_wr *
	arf, c_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ctgevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, c_type_wr *s, i_type_wr *lds, c_type_wr *p, i_type_wr *ldp, 
	c_type_wr *vl, i_type_wr *ldvl, c_type_wr *vr, i_type_wr *ldvr, i_type_wr *mm, 
	i_type_wr *m, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ctgex2_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, 
	c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *q, 
	i_type_wr *ldq, c_type_wr *z__, i_type_wr *ldz, i_type_wr *j1, i_type_wr *info);

/* Subroutine */ int ctgexc_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, 
	c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *q, 
	i_type_wr *ldq, c_type_wr *z__, i_type_wr *ldz, i_type_wr *ifst, i_type_wr *
	ilst, i_type_wr *info);

/* Subroutine */ int ctgsen_(i_type_wr *ijob, l_type_wr *wantq, l_type_wr *wantz, 
	l_type_wr *select, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, 
	i_type_wr *ldb, c_type_wr *alpha, c_type_wr *beta, c_type_wr *q, i_type_wr *ldq, 
	 c_type_wr *z__, i_type_wr *ldz, i_type_wr *m, s_type_wr *pl, s_type_wr *pr, s_type_wr *
	dif, c_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, 
	i_type_wr *info);

/* Subroutine */ int ctgsja_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, i_type_wr *k, i_type_wr *l, c_type_wr *a, i_type_wr *
	lda, c_type_wr *b, i_type_wr *ldb, s_type_wr *tola, s_type_wr *tolb, s_type_wr *alpha, 
	s_type_wr *beta, c_type_wr *u, i_type_wr *ldu, c_type_wr *v, i_type_wr *ldv, 
	c_type_wr *q, i_type_wr *ldq, c_type_wr *work, i_type_wr *ncycle, i_type_wr *
	info);

/* Subroutine */ int ctgsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
	c_type_wr *vl, i_type_wr *ldvl, c_type_wr *vr, i_type_wr *ldvr, s_type_wr *s, s_type_wr 
	*dif, i_type_wr *mm, i_type_wr *m, c_type_wr *work, i_type_wr *lwork, i_type_wr 
	*iwork, i_type_wr *info);

/* Subroutine */ int ctgsy2_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *d__, i_type_wr *ldd, c_type_wr *e, i_type_wr *lde, 
	c_type_wr *f, i_type_wr *ldf, s_type_wr *scale, s_type_wr *rdsum, s_type_wr *rdscal, 
	i_type_wr *info);

/* Subroutine */ int ctgsyl_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *d__, i_type_wr *ldd, c_type_wr *e, i_type_wr *lde, 
	c_type_wr *f, i_type_wr *ldf, s_type_wr *scale, s_type_wr *dif, c_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int ctpcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	c_type_wr *ap, s_type_wr *rcond, c_type_wr *work, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ctprfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, c_type_wr *ap, c_type_wr *b, i_type_wr *ldb, c_type_wr *x, 
	i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int ctptri_(char *uplo, char *diag, i_type_wr *n, c_type_wr *ap, 
	i_type_wr *info);

/* Subroutine */ int ctptrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, c_type_wr *ap, c_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int ctpttf_(char *transr, char *uplo, i_type_wr *n, c_type_wr *
	ap, c_type_wr *arf, i_type_wr *info);

/* Subroutine */ int ctpttr_(char *uplo, i_type_wr *n, c_type_wr *ap, c_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ctrcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	c_type_wr *a, i_type_wr *lda, s_type_wr *rcond, c_type_wr *work, s_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int ctrevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, c_type_wr *t, i_type_wr *ldt, c_type_wr *vl, i_type_wr *ldvl, 
	c_type_wr *vr, i_type_wr *ldvr, i_type_wr *mm, i_type_wr *m, c_type_wr *work, 
	s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ctrexc_(char *compq, i_type_wr *n, c_type_wr *t, i_type_wr *
	ldt, c_type_wr *q, i_type_wr *ldq, i_type_wr *ifst, i_type_wr *ilst, i_type_wr *
	info);

/* Subroutine */ int ctrrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
	c_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, c_type_wr *work, s_type_wr 
	*rwork, i_type_wr *info);

/* Subroutine */ int ctrsen_(char *job, char *compq, l_type_wr *select, i_type_wr 
	*n, c_type_wr *t, i_type_wr *ldt, c_type_wr *q, i_type_wr *ldq, c_type_wr *w, 
	i_type_wr *m, s_type_wr *s, s_type_wr *sep, c_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int ctrsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, c_type_wr *t, i_type_wr *ldt, c_type_wr *vl, i_type_wr *ldvl, 
	c_type_wr *vr, i_type_wr *ldvr, s_type_wr *s, s_type_wr *sep, i_type_wr *mm, i_type_wr *
	m, c_type_wr *work, i_type_wr *ldwork, s_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ctrsyl_(char *trana, char *tranb, i_type_wr *isgn, i_type_wr 
	*m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
	c_type_wr *c__, i_type_wr *ldc, s_type_wr *scale, i_type_wr *info);

/* Subroutine */ int ctrti2_(char *uplo, char *diag, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ctrtri_(char *uplo, char *diag, i_type_wr *n, c_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ctrtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int ctrttf_(char *transr, char *uplo, i_type_wr *n, c_type_wr *a, 
	 i_type_wr *lda, c_type_wr *arf, i_type_wr *info);

/* Subroutine */ int ctrttp_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *ap, i_type_wr *info);

/* Subroutine */ int ctzrqf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, i_type_wr *info);

/* Subroutine */ int ctzrzf_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cung2l_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cung2r_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cungbr_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int cunghr_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, c_type_wr *
	a, i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr 
	*info);

/* Subroutine */ int cungl2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cunglq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int cungql_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int cungqr_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int cungr2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cungrq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, c_type_wr *a, 
	i_type_wr *lda, c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int cungtr_(char *uplo, i_type_wr *n, c_type_wr *a, i_type_wr *lda, 
	 c_type_wr *tau, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cunm2l_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cunm2r_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cunmbr_(char *vect, char *side, char *trans, i_type_wr *m, 
	i_type_wr *n, i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, 
	c_type_wr *c__, i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int cunmhr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, 
	c_type_wr *c__, i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int cunml2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cunmlq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cunmql_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cunmqr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cunmr2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cunmr3_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, 
	c_type_wr *c__, i_type_wr *ldc, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cunmrq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cunmrz_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, 
	c_type_wr *c__, i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int cunmtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, c_type_wr *a, i_type_wr *lda, c_type_wr *tau, c_type_wr *c__, 
	i_type_wr *ldc, c_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int cupgtr_(char *uplo, i_type_wr *n, c_type_wr *ap, c_type_wr *
	tau, c_type_wr *q, i_type_wr *ldq, c_type_wr *work, i_type_wr *info);

/* Subroutine */ int cupmtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, c_type_wr *ap, c_type_wr *tau, c_type_wr *c__, i_type_wr *ldc, 
	c_type_wr *work, i_type_wr *info);

/* Subroutine */ int dbdsdc_(char *uplo, char *compq, i_type_wr *n, d_type_wr *
	d__, d_type_wr *e, d_type_wr *u, i_type_wr *ldu, d_type_wr *vt, 
	i_type_wr *ldvt, d_type_wr *q, i_type_wr *iq, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dbdsqr_(char *uplo, i_type_wr *n, i_type_wr *ncvt, i_type_wr *
	nru, i_type_wr *ncc, d_type_wr *d__, d_type_wr *e, d_type_wr *vt, 
	i_type_wr *ldvt, d_type_wr *u, i_type_wr *ldu, d_type_wr *c__, i_type_wr *
	ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int ddisna_(char *job, i_type_wr *m, i_type_wr *n, d_type_wr *
	d__, d_type_wr *sep, i_type_wr *info);

/* Subroutine */ int dgbbrd_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *ncc, 
	 i_type_wr *kl, i_type_wr *ku, d_type_wr *ab, i_type_wr *ldab, d_type_wr *
	d__, d_type_wr *e, d_type_wr *q, i_type_wr *ldq, d_type_wr *pt, 
	i_type_wr *ldpt, d_type_wr *c__, i_type_wr *ldc, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dgbcon_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 d_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, d_type_wr *anorm, 
	d_type_wr *rcond, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgbequ_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 d_type_wr *ab, i_type_wr *ldab, d_type_wr *r__, d_type_wr *c__, 
	d_type_wr *rowcnd, d_type_wr *colcnd, d_type_wr *amax, i_type_wr *
	info);

/* Subroutine */ int dgbequb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, d_type_wr *ab, i_type_wr *ldab, d_type_wr *r__, d_type_wr *c__, 
	d_type_wr *rowcnd, d_type_wr *colcnd, d_type_wr *amax, i_type_wr *
	info);

/* Subroutine */ int dgbrfs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, d_type_wr *afb, 
	i_type_wr *ldafb, i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgbrfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, 
	d_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, d_type_wr *r__, 
	d_type_wr *c__, d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *
	ldx, d_type_wr *rcond, d_type_wr *berr, i_type_wr *n_err_bnds__, 
	d_type_wr *err_bnds_norm__, d_type_wr *err_bnds_comp__, i_type_wr *
	nparams, d_type_wr *params, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dgbsv_(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	nrhs, d_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, d_type_wr *b, 
	i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dgbsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *kl, 
	 i_type_wr *ku, i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, 
	d_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, char *equed, 
	d_type_wr *r__, d_type_wr *c__, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, 
	d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgbsvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, 
	d_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, char *equed, 
	d_type_wr *r__, d_type_wr *c__, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *rpvgrw, 
	d_type_wr *berr, i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, 
	d_type_wr *err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgbtf2_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 d_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int dgbtrf_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 d_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int dgbtrs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, 
	d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dgebak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, d_type_wr *scale, i_type_wr *m, d_type_wr *v, i_type_wr *
	ldv, i_type_wr *info);

/* Subroutine */ int dgebal_(char *job, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *ilo, i_type_wr *ihi, d_type_wr *scale, i_type_wr *info);

/* Subroutine */ int dgebd2_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *d__, d_type_wr *e, d_type_wr *tauq, d_type_wr *
	taup, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dgebrd_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *d__, d_type_wr *e, d_type_wr *tauq, d_type_wr *
	taup, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgecon_(char *norm, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *anorm, d_type_wr *rcond, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dgeequ_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *r__, d_type_wr *c__, d_type_wr *rowcnd, d_type_wr 
	*colcnd, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int dgeequb_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *r__, d_type_wr *c__, d_type_wr *rowcnd, d_type_wr 
	*colcnd, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int dgees_(char *jobvs, char *sort, sel_fun_wr select, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, i_type_wr *sdim, d_type_wr *wr, 
	d_type_wr *wi, d_type_wr *vs, i_type_wr *ldvs, d_type_wr *work, 
	i_type_wr *lwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int dgeesx_(char *jobvs, char *sort, sel_fun_wr select, char *
	sense, i_type_wr *n, d_type_wr *a, i_type_wr *lda, i_type_wr *sdim, 
	d_type_wr *wr, d_type_wr *wi, d_type_wr *vs, i_type_wr *ldvs, 
	d_type_wr *rconde, d_type_wr *rcondv, d_type_wr *work, i_type_wr *
	lwork, i_type_wr *iwork, i_type_wr *liwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int dgeev_(char *jobvl, char *jobvr, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *wr, d_type_wr *wi, d_type_wr *vl, 
	i_type_wr *ldvl, d_type_wr *vr, i_type_wr *ldvr, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *wr, 
	d_type_wr *wi, d_type_wr *vl, i_type_wr *ldvl, d_type_wr *vr, 
	i_type_wr *ldvr, i_type_wr *ilo, i_type_wr *ihi, d_type_wr *scale, 
	d_type_wr *abnrm, d_type_wr *rconde, d_type_wr *rcondv, d_type_wr 
	*work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgegs_(char *jobvsl, char *jobvsr, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	alphar, d_type_wr *alphai, d_type_wr *beta, d_type_wr *vsl, 
	i_type_wr *ldvsl, d_type_wr *vsr, i_type_wr *ldvsr, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgegv_(char *jobvl, char *jobvr, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *alphar, 
	d_type_wr *alphai, d_type_wr *beta, d_type_wr *vl, i_type_wr *ldvl, 
	d_type_wr *vr, i_type_wr *ldvr, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dgehd2_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dgehrd_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgejsv_(char *joba, char *jobu, char *jobv, char *jobr, 
	char *jobt, char *jobp, i_type_wr *m, i_type_wr *n, d_type_wr *a, 
	i_type_wr *lda, d_type_wr *sva, d_type_wr *u, i_type_wr *ldu, 
	d_type_wr *v, i_type_wr *ldv, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgelq2_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dgelqf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgels_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgelsd_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	s, d_type_wr *rcond, i_type_wr *rank, d_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgelss_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	s, d_type_wr *rcond, i_type_wr *rank, d_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int dgelsx_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, i_type_wr *
	jpvt, d_type_wr *rcond, i_type_wr *rank, d_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int dgelsy_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, i_type_wr *
	jpvt, d_type_wr *rcond, i_type_wr *rank, d_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int dgeql2_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dgeqlf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgeqp3_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *jpvt, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int dgeqpf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *jpvt, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dgeqr2_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dgeqrf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgerfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, 
	d_type_wr *ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dgerfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, d_type_wr *r__, d_type_wr *c__, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, 
	d_type_wr *berr, i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, 
	d_type_wr *err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgerq2_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dgerqf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgesc2_(i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *rhs, i_type_wr *ipiv, i_type_wr *jpiv, d_type_wr *scale);

/* Subroutine */ int dgesdd_(char *jobz, i_type_wr *m, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *s, d_type_wr *u, i_type_wr *ldu, 
	d_type_wr *vt, i_type_wr *ldvt, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dgesv_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *a, i_type_wr 
	*lda, i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dgesvd_(char *jobu, char *jobvt, i_type_wr *m, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *s, d_type_wr *u, i_type_wr *
	ldu, d_type_wr *vt, i_type_wr *ldvt, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dgesvj_(char *joba, char *jobu, char *jobv, i_type_wr *m, 
	i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *sva, i_type_wr *mv, 
	 d_type_wr *v, i_type_wr *ldv, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dgesvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, char *equed, d_type_wr *r__, d_type_wr *c__, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	rcond, d_type_wr *ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dgesvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, char *equed, d_type_wr *r__, d_type_wr *c__, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	rcond, d_type_wr *rpvgrw, d_type_wr *berr, i_type_wr *n_err_bnds__, 
	d_type_wr *err_bnds_norm__, d_type_wr *err_bnds_comp__, i_type_wr *
	nparams, d_type_wr *params, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dgetc2_(i_type_wr *n, d_type_wr *a, i_type_wr *lda, i_type_wr 
	*ipiv, i_type_wr *jpiv, i_type_wr *info);

/* Subroutine */ int dgetf2_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int dgetrf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int dgetri_(i_type_wr *n, d_type_wr *a, i_type_wr *lda, i_type_wr 
	*ipiv, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dgetrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, d_type_wr *b, i_type_wr *
	ldb, i_type_wr *info);

/* Subroutine */ int dggbak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, d_type_wr *lscale, d_type_wr *rscale, i_type_wr *m, 
	d_type_wr *v, i_type_wr *ldv, i_type_wr *info);

/* Subroutine */ int dggbal_(char *job, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *b, i_type_wr *ldb, i_type_wr *ilo, i_type_wr *ihi, 
	d_type_wr *lscale, d_type_wr *rscale, d_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int dgges_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
	i_type_wr *ldb, i_type_wr *sdim, d_type_wr *alphar, d_type_wr *alphai, 
	d_type_wr *beta, d_type_wr *vsl, i_type_wr *ldvsl, d_type_wr *vsr, 
	i_type_wr *ldvsr, d_type_wr *work, i_type_wr *lwork, l_type_wr *bwork, 
	i_type_wr *info);

/* Subroutine */ int dggesx_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, char *sense, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *b, i_type_wr *ldb, i_type_wr *sdim, d_type_wr *alphar, 
	d_type_wr *alphai, d_type_wr *beta, d_type_wr *vsl, i_type_wr *ldvsl, 
	 d_type_wr *vsr, i_type_wr *ldvsr, d_type_wr *rconde, d_type_wr *
	rcondv, d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *
	liwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int dggev_(char *jobvl, char *jobvr, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *alphar, 
	d_type_wr *alphai, d_type_wr *beta, d_type_wr *vl, i_type_wr *ldvl, 
	d_type_wr *vr, i_type_wr *ldvr, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *alphar, d_type_wr *alphai, d_type_wr *
	beta, d_type_wr *vl, i_type_wr *ldvl, d_type_wr *vr, i_type_wr *ldvr, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *lscale, d_type_wr *rscale, 
	d_type_wr *abnrm, d_type_wr *bbnrm, d_type_wr *rconde, d_type_wr *
	rcondv, d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, l_type_wr *
	bwork, i_type_wr *info);

/* Subroutine */ int dggglm_(i_type_wr *n, i_type_wr *m, i_type_wr *p, d_type_wr *
	a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *d__, 
	d_type_wr *x, d_type_wr *y, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dgghrd_(char *compq, char *compz, i_type_wr *n, i_type_wr *
	ilo, i_type_wr *ihi, d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *q, i_type_wr *ldq, d_type_wr *z__, i_type_wr *
	ldz, i_type_wr *info);

/* Subroutine */ int dgglse_(i_type_wr *m, i_type_wr *n, i_type_wr *p, d_type_wr *
	a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *c__, 
	d_type_wr *d__, d_type_wr *x, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dggqrf_(i_type_wr *n, i_type_wr *m, i_type_wr *p, d_type_wr *
	a, i_type_wr *lda, d_type_wr *taua, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *taub, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dggrqf_(i_type_wr *m, i_type_wr *p, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *taua, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *taub, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dggsvd_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *n, i_type_wr *p, i_type_wr *k, i_type_wr *l, d_type_wr *a, 
	i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *alpha, 
	d_type_wr *beta, d_type_wr *u, i_type_wr *ldu, d_type_wr *v, i_type_wr 
	*ldv, d_type_wr *q, i_type_wr *ldq, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dggsvp_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *tola, d_type_wr *tolb, i_type_wr *k, i_type_wr 
	*l, d_type_wr *u, i_type_wr *ldu, d_type_wr *v, i_type_wr *ldv, 
	d_type_wr *q, i_type_wr *ldq, i_type_wr *iwork, d_type_wr *tau, 
	d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dgsvj0_(char *jobv, i_type_wr *m, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *d__, d_type_wr *sva, i_type_wr *mv, 
	d_type_wr *v, i_type_wr *ldv, d_type_wr *eps, d_type_wr *sfmin, 
	d_type_wr *tol, i_type_wr *nsweep, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dgsvj1_(char *jobv, i_type_wr *m, i_type_wr *n, i_type_wr *n1, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *d__, d_type_wr *sva, 
	i_type_wr *mv, d_type_wr *v, i_type_wr *ldv, d_type_wr *eps, d_type_wr 
	*sfmin, d_type_wr *tol, i_type_wr *nsweep, d_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int dgtcon_(char *norm, i_type_wr *n, d_type_wr *dl, 
	d_type_wr *d__, d_type_wr *du, d_type_wr *du2, i_type_wr *ipiv, 
	d_type_wr *anorm, d_type_wr *rcond, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dgtrfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *dl, d_type_wr *d__, d_type_wr *du, d_type_wr *dlf, 
	d_type_wr *df, d_type_wr *duf, d_type_wr *du2, i_type_wr *ipiv, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int dgtsv_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *dl, 
	d_type_wr *d__, d_type_wr *du, d_type_wr *b, i_type_wr *ldb, i_type_wr 
	*info);

/* Subroutine */ int dgtsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *dl, d_type_wr *d__, d_type_wr *du, d_type_wr *
	dlf, d_type_wr *df, d_type_wr *duf, d_type_wr *du2, i_type_wr *ipiv, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	rcond, d_type_wr *ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dgttrf_(i_type_wr *n, d_type_wr *dl, d_type_wr *d__, 
	d_type_wr *du, d_type_wr *du2, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int dgttrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *dl, d_type_wr *d__, d_type_wr *du, d_type_wr *du2, 
	i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dgtts2_(i_type_wr *itrans, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *dl, d_type_wr *d__, d_type_wr *du, d_type_wr *du2, 
	i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int dhgeqz_(char *job, char *compq, char *compz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *h__, i_type_wr *ldh, d_type_wr 
	*t, i_type_wr *ldt, d_type_wr *alphar, d_type_wr *alphai, d_type_wr *
	beta, d_type_wr *q, i_type_wr *ldq, d_type_wr *z__, i_type_wr *ldz, 
	d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dhsein_(char *side, char *eigsrc, char *initv, l_type_wr *
	select, i_type_wr *n, d_type_wr *h__, i_type_wr *ldh, d_type_wr *wr, 
	d_type_wr *wi, d_type_wr *vl, i_type_wr *ldvl, d_type_wr *vr, 
	i_type_wr *ldvr, i_type_wr *mm, i_type_wr *m, d_type_wr *work, i_type_wr *
	ifaill, i_type_wr *ifailr, i_type_wr *info);

/* Subroutine */ int dhseqr_(char *job, char *compz, i_type_wr *n, i_type_wr *ilo, 
	 i_type_wr *ihi, d_type_wr *h__, i_type_wr *ldh, d_type_wr *wr, 
	d_type_wr *wi, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

l_type_wr disnan_(d_type_wr *din);

/* Subroutine */ int dla_gbamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *kl, i_type_wr *ku, d_type_wr *alpha, d_type_wr *ab, i_type_wr *
	ldab, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, 
	i_type_wr *incy);

d_type_wr dla_gbrcond__(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	d_type_wr *ab, i_type_wr *ldab, d_type_wr *afb, i_type_wr *ldafb, 
	i_type_wr *ipiv, i_type_wr *cmode, d_type_wr *c__, i_type_wr *info, 
	d_type_wr *work, i_type_wr *iwork, ftn_len_wr trans_len);

/* Subroutine */ int dla_gbrfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *nrhs, 
	d_type_wr *ab, i_type_wr *ldab, d_type_wr *afb, i_type_wr *ldafb, 
	i_type_wr *ipiv, l_type_wr *colequ, d_type_wr *c__, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *y, i_type_wr *ldy, d_type_wr *berr_out__, 
	i_type_wr *n_norms__, d_type_wr *errs_n__, d_type_wr *errs_c__, 
	d_type_wr *res, d_type_wr *ayb, d_type_wr *dy, d_type_wr *
	y_tail__, d_type_wr *rcond, i_type_wr *ithresh, d_type_wr *rthresh, 
	d_type_wr *dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info);

d_type_wr dla_gbrpvgrw__(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	ncols, d_type_wr *ab, i_type_wr *ldab, d_type_wr *afb, i_type_wr *ldafb);

/* Subroutine */ int dla_geamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, 
	d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *x, 
	i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy);

d_type_wr dla_gercond__(char *trans, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
	 d_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, i_type_wr *cmode, 
	d_type_wr *c__, i_type_wr *info, d_type_wr *work, i_type_wr *iwork, 
	ftn_len_wr trans_len);

/* Subroutine */ int dla_gerfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *nrhs, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, 
	d_type_wr *c__, d_type_wr *b, i_type_wr *ldb, d_type_wr *y, i_type_wr *
	ldy, d_type_wr *berr_out__, i_type_wr *n_norms__, d_type_wr *errs_n__,
	 d_type_wr *errs_c__, d_type_wr *res, d_type_wr *ayb, d_type_wr *
	dy, d_type_wr *y_tail__, d_type_wr *rcond, i_type_wr *ithresh, 
	d_type_wr *rthresh, d_type_wr *dz_ub__, l_type_wr *ignore_cwise__, 
	i_type_wr *info);

/* Subroutine */ int dla_lin_berr__(i_type_wr *n, i_type_wr *nz, i_type_wr *nrhs, 
	d_type_wr *res, d_type_wr *ayb, d_type_wr *berr);

d_type_wr dla_porcond__(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *af, i_type_wr *ldaf, i_type_wr *cmode, d_type_wr *c__, 
	i_type_wr *info, d_type_wr *work, i_type_wr *iwork, ftn_len_wr uplo_len);

/* Subroutine */ int dla_porfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *
	af, i_type_wr *ldaf, l_type_wr *colequ, d_type_wr *c__, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *y, i_type_wr *ldy, d_type_wr *berr_out__, 
	i_type_wr *n_norms__, d_type_wr *errs_n__, d_type_wr *errs_c__, 
	d_type_wr *res, d_type_wr *ayb, d_type_wr *dy, d_type_wr *
	y_tail__, d_type_wr *rcond, i_type_wr *ithresh, d_type_wr *rthresh, 
	d_type_wr *dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info, ftn_len_wr 
	uplo_len);

d_type_wr dla_porpvgrw__(char *uplo, i_type_wr *ncols, d_type_wr *a, i_type_wr *
	lda, d_type_wr *af, i_type_wr *ldaf, d_type_wr *work, ftn_len_wr uplo_len);

d_type_wr dla_rpvgrw__(i_type_wr *n, i_type_wr *ncols, d_type_wr *a, i_type_wr *
	lda, d_type_wr *af, i_type_wr *ldaf);

/* Subroutine */ int dla_syamv__(i_type_wr *uplo, i_type_wr *n, d_type_wr *alpha,
	 d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx, 
	d_type_wr *beta, d_type_wr *y, i_type_wr *incy);

d_type_wr dla_syrcond__(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, i_type_wr *cmode, 
	d_type_wr *c__, i_type_wr *info, d_type_wr *work, i_type_wr *iwork, 
	ftn_len_wr uplo_len);

/* Subroutine */ int dla_syrfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *
	af, i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, d_type_wr *c__, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *y, i_type_wr *ldy, d_type_wr *
	berr_out__, i_type_wr *n_norms__, d_type_wr *errs_n__, d_type_wr *
	errs_c__, d_type_wr *res, d_type_wr *ayb, d_type_wr *dy, 
	d_type_wr *y_tail__, d_type_wr *rcond, i_type_wr *ithresh, d_type_wr 
	*rthresh, d_type_wr *dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info,
	 ftn_len_wr uplo_len);

d_type_wr dla_syrpvgrw__(char *uplo, i_type_wr *n, i_type_wr *info, d_type_wr *
	a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, 
	d_type_wr *work, ftn_len_wr uplo_len);

/* Subroutine */ int dla_wwaddw__(i_type_wr *n, d_type_wr *x, d_type_wr *y, 
	d_type_wr *w);

/* Subroutine */ int dlabad_(d_type_wr *small, d_type_wr *large);

/* Subroutine */ int dlabrd_(i_type_wr *m, i_type_wr *n, i_type_wr *nb, d_type_wr *
	a, i_type_wr *lda, d_type_wr *d__, d_type_wr *e, d_type_wr *tauq, 
	d_type_wr *taup, d_type_wr *x, i_type_wr *ldx, d_type_wr *y, i_type_wr 
	*ldy);

/* Subroutine */ int dlacn2_(i_type_wr *n, d_type_wr *v, d_type_wr *x, 
	i_type_wr *isgn, d_type_wr *est, i_type_wr *kase, i_type_wr *isave);

/* Subroutine */ int dlacon_(i_type_wr *n, d_type_wr *v, d_type_wr *x, 
	i_type_wr *isgn, d_type_wr *est, i_type_wr *kase);

/* Subroutine */ int dlacpy_(char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int dladiv_(d_type_wr *a, d_type_wr *b, d_type_wr *c__, 
	d_type_wr *d__, d_type_wr *p, d_type_wr *q);

/* Subroutine */ int dlae2_(d_type_wr *a, d_type_wr *b, d_type_wr *c__, 
	d_type_wr *rt1, d_type_wr *rt2);

/* Subroutine */ int dlaebz_(i_type_wr *ijob, i_type_wr *nitmax, i_type_wr *n, 
	i_type_wr *mmax, i_type_wr *minp, i_type_wr *nbmin, d_type_wr *abstol, 
	d_type_wr *reltol, d_type_wr *pivmin, d_type_wr *d__, d_type_wr *
	e, d_type_wr *e2, i_type_wr *nval, d_type_wr *ab, d_type_wr *c__, 
	i_type_wr *mout, i_type_wr *nab, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dlaed0_(i_type_wr *icompq, i_type_wr *qsiz, i_type_wr *n, 
	d_type_wr *d__, d_type_wr *e, d_type_wr *q, i_type_wr *ldq, 
	d_type_wr *qstore, i_type_wr *ldqs, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dlaed1_(i_type_wr *n, d_type_wr *d__, d_type_wr *q, 
	i_type_wr *ldq, i_type_wr *indxq, d_type_wr *rho, i_type_wr *cutpnt, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dlaed2_(i_type_wr *k, i_type_wr *n, i_type_wr *n1, d_type_wr *
	d__, d_type_wr *q, i_type_wr *ldq, i_type_wr *indxq, d_type_wr *rho, 
	d_type_wr *z__, d_type_wr *dlamda, d_type_wr *w, d_type_wr *q2, 
	i_type_wr *indx, i_type_wr *indxc, i_type_wr *indxp, i_type_wr *coltyp, 
	i_type_wr *info);

/* Subroutine */ int dlaed3_(i_type_wr *k, i_type_wr *n, i_type_wr *n1, d_type_wr *
	d__, d_type_wr *q, i_type_wr *ldq, d_type_wr *rho, d_type_wr *dlamda, 
	 d_type_wr *q2, i_type_wr *indx, i_type_wr *ctot, d_type_wr *w, 
	d_type_wr *s, i_type_wr *info);

/* Subroutine */ int dlaed4_(i_type_wr *n, i_type_wr *i__, d_type_wr *d__, 
	d_type_wr *z__, d_type_wr *delta, d_type_wr *rho, d_type_wr *dlam, 
	 i_type_wr *info);

/* Subroutine */ int dlaed5_(i_type_wr *i__, d_type_wr *d__, d_type_wr *z__, 
	d_type_wr *delta, d_type_wr *rho, d_type_wr *dlam);

/* Subroutine */ int dlaed6_(i_type_wr *kniter, l_type_wr *orgati, d_type_wr *
	rho, d_type_wr *d__, d_type_wr *z__, d_type_wr *finit, d_type_wr *
	tau, i_type_wr *info);

/* Subroutine */ int dlaed7_(i_type_wr *icompq, i_type_wr *n, i_type_wr *qsiz, 
	i_type_wr *tlvls, i_type_wr *curlvl, i_type_wr *curpbm, d_type_wr *d__, 
	d_type_wr *q, i_type_wr *ldq, i_type_wr *indxq, d_type_wr *rho, i_type_wr 
	*cutpnt, d_type_wr *qstore, i_type_wr *qptr, i_type_wr *prmptr, i_type_wr *
	perm, i_type_wr *givptr, i_type_wr *givcol, d_type_wr *givnum, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dlaed8_(i_type_wr *icompq, i_type_wr *k, i_type_wr *n, i_type_wr 
	*qsiz, d_type_wr *d__, d_type_wr *q, i_type_wr *ldq, i_type_wr *indxq, 
	d_type_wr *rho, i_type_wr *cutpnt, d_type_wr *z__, d_type_wr *dlamda, 
	 d_type_wr *q2, i_type_wr *ldq2, d_type_wr *w, i_type_wr *perm, i_type_wr 
	*givptr, i_type_wr *givcol, d_type_wr *givnum, i_type_wr *indxp, i_type_wr 
	*indx, i_type_wr *info);

/* Subroutine */ int dlaed9_(i_type_wr *k, i_type_wr *kstart, i_type_wr *kstop, 
	i_type_wr *n, d_type_wr *d__, d_type_wr *q, i_type_wr *ldq, d_type_wr *
	rho, d_type_wr *dlamda, d_type_wr *w, d_type_wr *s, i_type_wr *lds, 
	i_type_wr *info);

/* Subroutine */ int dlaeda_(i_type_wr *n, i_type_wr *tlvls, i_type_wr *curlvl, 
	i_type_wr *curpbm, i_type_wr *prmptr, i_type_wr *perm, i_type_wr *givptr, 
	i_type_wr *givcol, d_type_wr *givnum, d_type_wr *q, i_type_wr *qptr, 
	d_type_wr *z__, d_type_wr *ztemp, i_type_wr *info);

/* Subroutine */ int dlaein_(l_type_wr *rightv, l_type_wr *noinit, i_type_wr *n, 
	d_type_wr *h__, i_type_wr *ldh, d_type_wr *wr, d_type_wr *wi, 
	d_type_wr *vr, d_type_wr *vi, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *work, d_type_wr *eps3, d_type_wr *smlnum, d_type_wr *
	bignum, i_type_wr *info);

/* Subroutine */ int dlaev2_(d_type_wr *a, d_type_wr *b, d_type_wr *c__, 
	d_type_wr *rt1, d_type_wr *rt2, d_type_wr *cs1, d_type_wr *sn1);

/* Subroutine */ int dlaexc_(l_type_wr *wantq, i_type_wr *n, d_type_wr *t, 
	i_type_wr *ldt, d_type_wr *q, i_type_wr *ldq, i_type_wr *j1, i_type_wr *n1, 
	i_type_wr *n2, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlag2_(d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *safmin, d_type_wr *scale1, d_type_wr *
	scale2, d_type_wr *wr1, d_type_wr *wr2, d_type_wr *wi);

/* Subroutine */ int dlag2s_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, s_type_wr *sa, i_type_wr *ldsa, i_type_wr *info);

/* Subroutine */ int dlags2_(l_type_wr *upper, d_type_wr *a1, d_type_wr *a2, 
	d_type_wr *a3, d_type_wr *b1, d_type_wr *b2, d_type_wr *b3, 
	d_type_wr *csu, d_type_wr *snu, d_type_wr *csv, d_type_wr *snv, 
	d_type_wr *csq, d_type_wr *snq);

/* Subroutine */ int dlagtf_(i_type_wr *n, d_type_wr *a, d_type_wr *lambda, 
	d_type_wr *b, d_type_wr *c__, d_type_wr *tol, d_type_wr *d__, 
	i_type_wr *in, i_type_wr *info);

/* Subroutine */ int dlagtm_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *alpha, d_type_wr *dl, d_type_wr *d__, d_type_wr *du, 
	d_type_wr *x, i_type_wr *ldx, d_type_wr *beta, d_type_wr *b, i_type_wr 
	*ldb);

/* Subroutine */ int dlagts_(i_type_wr *job, i_type_wr *n, d_type_wr *a, 
	d_type_wr *b, d_type_wr *c__, d_type_wr *d__, i_type_wr *in, 
	d_type_wr *y, d_type_wr *tol, i_type_wr *info);

/* Subroutine */ int dlagv2_(d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *alphar, d_type_wr *alphai, d_type_wr *
	beta, d_type_wr *csl, d_type_wr *snl, d_type_wr *csr, d_type_wr *
	snr);

/* Subroutine */ int dlahqr_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *h__, i_type_wr *ldh, d_type_wr 
	*wr, d_type_wr *wi, i_type_wr *iloz, i_type_wr *ihiz, d_type_wr *z__, 
	i_type_wr *ldz, i_type_wr *info);

/* Subroutine */ int dlahr2_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *t, i_type_wr *ldt, 
	d_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int dlahrd_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *t, i_type_wr *ldt, 
	d_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int dlaic1_(i_type_wr *job, i_type_wr *j, d_type_wr *x, 
	d_type_wr *sest, d_type_wr *w, d_type_wr *gamma, d_type_wr *
	sestpr, d_type_wr *s, d_type_wr *c__);

l_type_wr dlaisnan_(d_type_wr *din1, d_type_wr *din2);

/* Subroutine */ int dlaln2_(l_type_wr *ltrans, i_type_wr *na, i_type_wr *nw, 
	d_type_wr *smin, d_type_wr *ca, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *d1, d_type_wr *d2, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *wr, d_type_wr *wi, d_type_wr *x, i_type_wr *ldx, 
	d_type_wr *scale, d_type_wr *xnorm, i_type_wr *info);

/* Subroutine */ int dlals0_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, i_type_wr *nrhs, d_type_wr *b, i_type_wr *ldb, d_type_wr 
	*bx, i_type_wr *ldbx, i_type_wr *perm, i_type_wr *givptr, i_type_wr *givcol, 
	i_type_wr *ldgcol, d_type_wr *givnum, i_type_wr *ldgnum, d_type_wr *
	poles, d_type_wr *difl, d_type_wr *difr, d_type_wr *z__, i_type_wr *
	k, d_type_wr *c__, d_type_wr *s, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlalsa_(i_type_wr *icompq, i_type_wr *smlsiz, i_type_wr *n, 
	i_type_wr *nrhs, d_type_wr *b, i_type_wr *ldb, d_type_wr *bx, i_type_wr *
	ldbx, d_type_wr *u, i_type_wr *ldu, d_type_wr *vt, i_type_wr *k, 
	d_type_wr *difl, d_type_wr *difr, d_type_wr *z__, d_type_wr *
	poles, i_type_wr *givptr, i_type_wr *givcol, i_type_wr *ldgcol, i_type_wr *
	perm, d_type_wr *givnum, d_type_wr *c__, d_type_wr *s, d_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dlalsd_(char *uplo, i_type_wr *smlsiz, i_type_wr *n, i_type_wr 
	*nrhs, d_type_wr *d__, d_type_wr *e, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *rcond, i_type_wr *rank, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dlamrg_(i_type_wr *n1, i_type_wr *n2, d_type_wr *a, i_type_wr 
	*dtrd1, i_type_wr *dtrd2, i_type_wr *index);

i_type_wr dlaneg_(i_type_wr *n, d_type_wr *d__, d_type_wr *lld, d_type_wr *
	sigma, d_type_wr *pivmin, i_type_wr *r__);

d_type_wr dlangb_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	d_type_wr *ab, i_type_wr *ldab, d_type_wr *work);

d_type_wr dlange_(char *norm, i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr 
	*lda, d_type_wr *work);

d_type_wr dlangt_(char *norm, i_type_wr *n, d_type_wr *dl, d_type_wr *d__, 
	d_type_wr *du);

d_type_wr dlanhs_(char *norm, i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *work);

d_type_wr dlansb_(char *norm, char *uplo, i_type_wr *n, i_type_wr *k, d_type_wr 
	*ab, i_type_wr *ldab, d_type_wr *work);

d_type_wr dlansf_(char *norm, char *transr, char *uplo, i_type_wr *n, 
	d_type_wr *a, d_type_wr *work);

d_type_wr dlansp_(char *norm, char *uplo, i_type_wr *n, d_type_wr *ap, 
	d_type_wr *work);

d_type_wr dlanst_(char *norm, i_type_wr *n, d_type_wr *d__, d_type_wr *e);

d_type_wr dlansy_(char *norm, char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr 
	*lda, d_type_wr *work);

d_type_wr dlantb_(char *norm, char *uplo, char *diag, i_type_wr *n, i_type_wr *k, 
	 d_type_wr *ab, i_type_wr *ldab, d_type_wr *work);

d_type_wr dlantp_(char *norm, char *uplo, char *diag, i_type_wr *n, d_type_wr 
	*ap, d_type_wr *work);

d_type_wr dlantr_(char *norm, char *uplo, char *diag, i_type_wr *m, i_type_wr *n, 
	 d_type_wr *a, i_type_wr *lda, d_type_wr *work);

/* Subroutine */ int dlanv2_(d_type_wr *a, d_type_wr *b, d_type_wr *c__, 
	d_type_wr *d__, d_type_wr *rt1r, d_type_wr *rt1i, d_type_wr *rt2r, 
	 d_type_wr *rt2i, d_type_wr *cs, d_type_wr *sn);

/* Subroutine */ int dlapll_(i_type_wr *n, d_type_wr *x, i_type_wr *incx, 
	d_type_wr *y, i_type_wr *incy, d_type_wr *ssmin);

/* Subroutine */ int dlapmt_(l_type_wr *forwrd, i_type_wr *m, i_type_wr *n, 
	d_type_wr *x, i_type_wr *ldx, i_type_wr *k);

d_type_wr dlapy2_(d_type_wr *x, d_type_wr *y);

d_type_wr dlapy3_(d_type_wr *x, d_type_wr *y, d_type_wr *z__);

/* Subroutine */ int dlaqgb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 d_type_wr *ab, i_type_wr *ldab, d_type_wr *r__, d_type_wr *c__, 
	d_type_wr *rowcnd, d_type_wr *colcnd, d_type_wr *amax, char *equed);

/* Subroutine */ int dlaqge_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *r__, d_type_wr *c__, d_type_wr *rowcnd, d_type_wr 
	*colcnd, d_type_wr *amax, char *equed);

/* Subroutine */ int dlaqp2_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, 
	d_type_wr *a, i_type_wr *lda, i_type_wr *jpvt, d_type_wr *tau, 
	d_type_wr *vn1, d_type_wr *vn2, d_type_wr *work);

/* Subroutine */ int dlaqps_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, i_type_wr 
	*nb, i_type_wr *kb, d_type_wr *a, i_type_wr *lda, i_type_wr *jpvt, 
	d_type_wr *tau, d_type_wr *vn1, d_type_wr *vn2, d_type_wr *auxv, 
	d_type_wr *f, i_type_wr *ldf);

/* Subroutine */ int dlaqr0_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *h__, i_type_wr *ldh, d_type_wr 
	*wr, d_type_wr *wi, i_type_wr *iloz, i_type_wr *ihiz, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dlaqr1_(i_type_wr *n, d_type_wr *h__, i_type_wr *ldh, 
	d_type_wr *sr1, d_type_wr *si1, d_type_wr *sr2, d_type_wr *si2, 
	d_type_wr *v);

/* Subroutine */ int dlaqr2_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, d_type_wr *h__, i_type_wr *
	ldh, i_type_wr *iloz, i_type_wr *ihiz, d_type_wr *z__, i_type_wr *ldz, 
	i_type_wr *ns, i_type_wr *nd, d_type_wr *sr, d_type_wr *si, d_type_wr *
	v, i_type_wr *ldv, i_type_wr *nh, d_type_wr *t, i_type_wr *ldt, i_type_wr *
	nv, d_type_wr *wv, i_type_wr *ldwv, d_type_wr *work, i_type_wr *lwork);

/* Subroutine */ int dlaqr3_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, d_type_wr *h__, i_type_wr *
	ldh, i_type_wr *iloz, i_type_wr *ihiz, d_type_wr *z__, i_type_wr *ldz, 
	i_type_wr *ns, i_type_wr *nd, d_type_wr *sr, d_type_wr *si, d_type_wr *
	v, i_type_wr *ldv, i_type_wr *nh, d_type_wr *t, i_type_wr *ldt, i_type_wr *
	nv, d_type_wr *wv, i_type_wr *ldwv, d_type_wr *work, i_type_wr *lwork);

/* Subroutine */ int dlaqr4_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *h__, i_type_wr *ldh, d_type_wr 
	*wr, d_type_wr *wi, i_type_wr *iloz, i_type_wr *ihiz, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dlaqr5_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *kacc22, 
	i_type_wr *n, i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nshfts, d_type_wr 
	*sr, d_type_wr *si, d_type_wr *h__, i_type_wr *ldh, i_type_wr *iloz, 
	i_type_wr *ihiz, d_type_wr *z__, i_type_wr *ldz, d_type_wr *v, i_type_wr *
	ldv, d_type_wr *u, i_type_wr *ldu, i_type_wr *nv, d_type_wr *wv, 
	i_type_wr *ldwv, i_type_wr *nh, d_type_wr *wh, i_type_wr *ldwh);

/* Subroutine */ int dlaqsb_(char *uplo, i_type_wr *n, i_type_wr *kd, d_type_wr *
	ab, i_type_wr *ldab, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, 
	 char *equed);

/* Subroutine */ int dlaqsp_(char *uplo, i_type_wr *n, d_type_wr *ap, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, char *equed);

/* Subroutine */ int dlaqsy_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, char *equed);

/* Subroutine */ int dlaqtr_(l_type_wr *ltran, l_type_wr *ls_type_wr, i_type_wr *n, 
	d_type_wr *t, i_type_wr *ldt, d_type_wr *b, d_type_wr *w, d_type_wr 
	*scale, d_type_wr *x, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlar1v_(i_type_wr *n, i_type_wr *b1, i_type_wr *bn, d_type_wr 
	*lambda, d_type_wr *d__, d_type_wr *l, d_type_wr *ld, d_type_wr *
	lld, d_type_wr *pivmin, d_type_wr *gaptol, d_type_wr *z__, l_type_wr 
	*wantnc, i_type_wr *negcnt, d_type_wr *ztz, d_type_wr *mingma, 
	i_type_wr *r__, i_type_wr *isuppz, d_type_wr *nrminv, d_type_wr *resid, 
	d_type_wr *rqcorr, d_type_wr *work);

/* Subroutine */ int dlar2v_(i_type_wr *n, d_type_wr *x, d_type_wr *y, 
	d_type_wr *z__, i_type_wr *incx, d_type_wr *c__, d_type_wr *s, 
	i_type_wr *incc);

/* Subroutine */ int dlarf_(char *side, i_type_wr *m, i_type_wr *n, d_type_wr *v, 
	 i_type_wr *incv, d_type_wr *tau, d_type_wr *c__, i_type_wr *ldc, 
	d_type_wr *work);

/* Subroutine */ int dlarfb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *v, i_type_wr *
	ldv, d_type_wr *t, i_type_wr *ldt, d_type_wr *c__, i_type_wr *ldc, 
	d_type_wr *work, i_type_wr *ldwork);

/* Subroutine */ int dlarfg_(i_type_wr *n, d_type_wr *alpha, d_type_wr *x, 
	i_type_wr *incx, d_type_wr *tau);

/* Subroutine */ int dlarfp_(i_type_wr *n, d_type_wr *alpha, d_type_wr *x, 
	i_type_wr *incx, d_type_wr *tau);

/* Subroutine */ int dlarft_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, d_type_wr *v, i_type_wr *ldv, d_type_wr *tau, d_type_wr *t, 
	i_type_wr *ldt);

/* Subroutine */ int dlarfx_(char *side, i_type_wr *m, i_type_wr *n, d_type_wr *
	v, d_type_wr *tau, d_type_wr *c__, i_type_wr *ldc, d_type_wr *work);

/* Subroutine */ int dlargv_(i_type_wr *n, d_type_wr *x, i_type_wr *incx, 
	d_type_wr *y, i_type_wr *incy, d_type_wr *c__, i_type_wr *incc);

/* Subroutine */ int dlarnv_(i_type_wr *idist, i_type_wr *iseed, i_type_wr *n, 
	d_type_wr *x);

/* Subroutine */ int dlarra_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	d_type_wr *e2, d_type_wr *spltol, d_type_wr *tnrm, i_type_wr *nsplit, 
	 i_type_wr *isplit, i_type_wr *info);

/* Subroutine */ int dlarrb_(i_type_wr *n, d_type_wr *d__, d_type_wr *lld, 
	i_type_wr *ifirst, i_type_wr *ilast, d_type_wr *rtol1, d_type_wr *rtol2, 
	 i_type_wr *offset, d_type_wr *w, d_type_wr *wgap, d_type_wr *werr, 
	d_type_wr *work, i_type_wr *iwork, d_type_wr *pivmin, d_type_wr *
	spdiam, i_type_wr *twist, i_type_wr *info);

/* Subroutine */ int dlarrc_(char *jobt, i_type_wr *n, d_type_wr *vl, 
	d_type_wr *vu, d_type_wr *d__, d_type_wr *e, d_type_wr *pivmin, 
	i_type_wr *eigcnt, i_type_wr *lcnt, i_type_wr *rcnt, i_type_wr *info);

/* Subroutine */ int dlarrd_(char *range, char *order, i_type_wr *n, d_type_wr 
	*vl, d_type_wr *vu, i_type_wr *il, i_type_wr *iu, d_type_wr *gers, 
	d_type_wr *reltol, d_type_wr *d__, d_type_wr *e, d_type_wr *e2, 
	d_type_wr *pivmin, i_type_wr *nsplit, i_type_wr *isplit, i_type_wr *m, 
	d_type_wr *w, d_type_wr *werr, d_type_wr *wl, d_type_wr *wu, 
	i_type_wr *iblock, i_type_wr *indexw, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dlarre_(char *range, i_type_wr *n, d_type_wr *vl, 
	d_type_wr *vu, i_type_wr *il, i_type_wr *iu, d_type_wr *d__, d_type_wr 
	*e, d_type_wr *e2, d_type_wr *rtol1, d_type_wr *rtol2, d_type_wr *
	spltol, i_type_wr *nsplit, i_type_wr *isplit, i_type_wr *m, d_type_wr *w, 
	d_type_wr *werr, d_type_wr *wgap, i_type_wr *iblock, i_type_wr *indexw, 
	d_type_wr *gers, d_type_wr *pivmin, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dlarrf_(i_type_wr *n, d_type_wr *d__, d_type_wr *l, 
	d_type_wr *ld, i_type_wr *clstrt, i_type_wr *clend, d_type_wr *w, 
	d_type_wr *wgap, d_type_wr *werr, d_type_wr *spdiam, d_type_wr *
	clgapl, d_type_wr *clgapr, d_type_wr *pivmin, d_type_wr *sigma, 
	d_type_wr *dplus, d_type_wr *lplus, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlarrj_(i_type_wr *n, d_type_wr *d__, d_type_wr *e2, 
	i_type_wr *ifirst, i_type_wr *ilast, d_type_wr *rtol, i_type_wr *offset, 
	d_type_wr *w, d_type_wr *werr, d_type_wr *work, i_type_wr *iwork, 
	d_type_wr *pivmin, d_type_wr *spdiam, i_type_wr *info);

/* Subroutine */ int dlarrk_(i_type_wr *n, i_type_wr *iw, d_type_wr *gl, 
	d_type_wr *gu, d_type_wr *d__, d_type_wr *e2, d_type_wr *pivmin, 
	d_type_wr *reltol, d_type_wr *w, d_type_wr *werr, i_type_wr *info);

/* Subroutine */ int dlarrr_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	i_type_wr *info);

/* Subroutine */ int dlarrv_(i_type_wr *n, d_type_wr *vl, d_type_wr *vu, 
	d_type_wr *d__, d_type_wr *l, d_type_wr *pivmin, i_type_wr *isplit, 
	i_type_wr *m, i_type_wr *dol, i_type_wr *dou, d_type_wr *minrgp, 
	d_type_wr *rtol1, d_type_wr *rtol2, d_type_wr *w, d_type_wr *werr, 
	 d_type_wr *wgap, i_type_wr *iblock, i_type_wr *indexw, d_type_wr *gers, 
	 d_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, d_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dlarscl2_(i_type_wr *m, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *x, i_type_wr *ldx);

/* Subroutine */ int dlartg_(d_type_wr *f, d_type_wr *g, d_type_wr *cs, 
	d_type_wr *sn, d_type_wr *r__);

/* Subroutine */ int dlartv_(i_type_wr *n, d_type_wr *x, i_type_wr *incx, 
	d_type_wr *y, i_type_wr *incy, d_type_wr *c__, d_type_wr *s, i_type_wr 
	*incc);

/* Subroutine */ int dlaruv_(i_type_wr *iseed, i_type_wr *n, d_type_wr *x);

/* Subroutine */ int dlarz_(char *side, i_type_wr *m, i_type_wr *n, i_type_wr *l, 
	d_type_wr *v, i_type_wr *incv, d_type_wr *tau, d_type_wr *c__, 
	i_type_wr *ldc, d_type_wr *work);

/* Subroutine */ int dlarzb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, i_type_wr *l, d_type_wr *v, 
	 i_type_wr *ldv, d_type_wr *t, i_type_wr *ldt, d_type_wr *c__, i_type_wr *
	ldc, d_type_wr *work, i_type_wr *ldwork);

/* Subroutine */ int dlarzt_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, d_type_wr *v, i_type_wr *ldv, d_type_wr *tau, d_type_wr *t, 
	i_type_wr *ldt);

/* Subroutine */ int dlas2_(d_type_wr *f, d_type_wr *g, d_type_wr *h__, 
	d_type_wr *ssmin, d_type_wr *ssmax);

/* Subroutine */ int dlascl_(char *type__, i_type_wr *kl, i_type_wr *ku, 
	d_type_wr *cfrom, d_type_wr *cto, i_type_wr *m, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int dlascl2_(i_type_wr *m, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *x, i_type_wr *ldx);

/* Subroutine */ int dlasd0_(i_type_wr *n, i_type_wr *sqre, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *u, i_type_wr *ldu, d_type_wr *vt, i_type_wr *
	ldvt, i_type_wr *smlsiz, i_type_wr *iwork, d_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int dlasd1_(i_type_wr *nl, i_type_wr *nr, i_type_wr *sqre, 
	d_type_wr *d__, d_type_wr *alpha, d_type_wr *beta, d_type_wr *u, 
	i_type_wr *ldu, d_type_wr *vt, i_type_wr *ldvt, i_type_wr *idxq, i_type_wr *
	iwork, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlasd2_(i_type_wr *nl, i_type_wr *nr, i_type_wr *sqre, i_type_wr 
	*k, d_type_wr *d__, d_type_wr *z__, d_type_wr *alpha, d_type_wr *
	beta, d_type_wr *u, i_type_wr *ldu, d_type_wr *vt, i_type_wr *ldvt, 
	d_type_wr *dsigma, d_type_wr *u2, i_type_wr *ldu2, d_type_wr *vt2, 
	i_type_wr *ldvt2, i_type_wr *idxp, i_type_wr *idx, i_type_wr *idxc, i_type_wr *
	idxq, i_type_wr *coltyp, i_type_wr *info);

/* Subroutine */ int dlasd3_(i_type_wr *nl, i_type_wr *nr, i_type_wr *sqre, i_type_wr 
	*k, d_type_wr *d__, d_type_wr *q, i_type_wr *ldq, d_type_wr *dsigma, 
	d_type_wr *u, i_type_wr *ldu, d_type_wr *u2, i_type_wr *ldu2, 
	d_type_wr *vt, i_type_wr *ldvt, d_type_wr *vt2, i_type_wr *ldvt2, 
	i_type_wr *idxc, i_type_wr *ctot, d_type_wr *z__, i_type_wr *info);

/* Subroutine */ int dlasd4_(i_type_wr *n, i_type_wr *i__, d_type_wr *d__, 
	d_type_wr *z__, d_type_wr *delta, d_type_wr *rho, d_type_wr *
	sigma, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlasd5_(i_type_wr *i__, d_type_wr *d__, d_type_wr *z__, 
	d_type_wr *delta, d_type_wr *rho, d_type_wr *dsigma, d_type_wr *
	work);

/* Subroutine */ int dlasd6_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, d_type_wr *d__, d_type_wr *vf, d_type_wr *vl, 
	d_type_wr *alpha, d_type_wr *beta, i_type_wr *idxq, i_type_wr *perm, 
	i_type_wr *givptr, i_type_wr *givcol, i_type_wr *ldgcol, d_type_wr *givnum, 
	 i_type_wr *ldgnum, d_type_wr *poles, d_type_wr *difl, d_type_wr *
	difr, d_type_wr *z__, i_type_wr *k, d_type_wr *c__, d_type_wr *s, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dlasd7_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, i_type_wr *k, d_type_wr *d__, d_type_wr *z__, 
	d_type_wr *zw, d_type_wr *vf, d_type_wr *vfw, d_type_wr *vl, 
	d_type_wr *vlw, d_type_wr *alpha, d_type_wr *beta, d_type_wr *
	dsigma, i_type_wr *idx, i_type_wr *idxp, i_type_wr *idxq, i_type_wr *perm, 
	i_type_wr *givptr, i_type_wr *givcol, i_type_wr *ldgcol, d_type_wr *givnum, 
	 i_type_wr *ldgnum, d_type_wr *c__, d_type_wr *s, i_type_wr *info);

/* Subroutine */ int dlasd8_(i_type_wr *icompq, i_type_wr *k, d_type_wr *d__, 
	d_type_wr *z__, d_type_wr *vf, d_type_wr *vl, d_type_wr *difl, 
	d_type_wr *difr, i_type_wr *lddifr, d_type_wr *dsigma, d_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int dlasda_(i_type_wr *icompq, i_type_wr *smlsiz, i_type_wr *n, 
	i_type_wr *sqre, d_type_wr *d__, d_type_wr *e, d_type_wr *u, i_type_wr 
	*ldu, d_type_wr *vt, i_type_wr *k, d_type_wr *difl, d_type_wr *difr, 
	d_type_wr *z__, d_type_wr *poles, i_type_wr *givptr, i_type_wr *givcol, 
	i_type_wr *ldgcol, i_type_wr *perm, d_type_wr *givnum, d_type_wr *c__, 
	d_type_wr *s, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dlasdq_(char *uplo, i_type_wr *sqre, i_type_wr *n, i_type_wr *
	ncvt, i_type_wr *nru, i_type_wr *ncc, d_type_wr *d__, d_type_wr *e, 
	d_type_wr *vt, i_type_wr *ldvt, d_type_wr *u, i_type_wr *ldu, 
	d_type_wr *c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlasdt_(i_type_wr *n, i_type_wr *lvl, i_type_wr *nd, i_type_wr *
	inode, i_type_wr *ndiml, i_type_wr *ndimr, i_type_wr *msub);

/* Subroutine */ int dlaset_(char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *
	alpha, d_type_wr *beta, d_type_wr *a, i_type_wr *lda);

/* Subroutine */ int dlasq1_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlasq2_(i_type_wr *n, d_type_wr *z__, i_type_wr *info);

/* Subroutine */ int dlasq3_(i_type_wr *i0, i_type_wr *n0, d_type_wr *z__, 
	i_type_wr *pp, d_type_wr *dmin__, d_type_wr *sigma, d_type_wr *desig, 
	 d_type_wr *qmax, i_type_wr *nfail, i_type_wr *iter, i_type_wr *ndiv, 
	l_type_wr *ieee, i_type_wr *ttype, d_type_wr *dmin1, d_type_wr *dmin2, 
	d_type_wr *dn, d_type_wr *dn1, d_type_wr *dn2, d_type_wr *g, 
	d_type_wr *tau);

/* Subroutine */ int dlasq4_(i_type_wr *i0, i_type_wr *n0, d_type_wr *z__, 
	i_type_wr *pp, i_type_wr *n0in, d_type_wr *dmin__, d_type_wr *dmin1, 
	d_type_wr *dmin2, d_type_wr *dn, d_type_wr *dn1, d_type_wr *dn2, 
	d_type_wr *tau, i_type_wr *ttype, d_type_wr *g);

/* Subroutine */ int dlasq5_(i_type_wr *i0, i_type_wr *n0, d_type_wr *z__, 
	i_type_wr *pp, d_type_wr *tau, d_type_wr *dmin__, d_type_wr *dmin1, 
	d_type_wr *dmin2, d_type_wr *dn, d_type_wr *dnm1, d_type_wr *dnm2, 
	 l_type_wr *ieee);

/* Subroutine */ int dlasq6_(i_type_wr *i0, i_type_wr *n0, d_type_wr *z__, 
	i_type_wr *pp, d_type_wr *dmin__, d_type_wr *dmin1, d_type_wr *dmin2, 
	 d_type_wr *dn, d_type_wr *dnm1, d_type_wr *dnm2);

/* Subroutine */ int dlasr_(char *side, char *pivot, char *direct, i_type_wr *m, 
	 i_type_wr *n, d_type_wr *c__, d_type_wr *s, d_type_wr *a, i_type_wr *
	lda);

/* Subroutine */ int dlasrt_(char *id, i_type_wr *n, d_type_wr *d__, i_type_wr *
	info);

/* Subroutine */ int dlassq_(i_type_wr *n, d_type_wr *x, i_type_wr *incx, 
	d_type_wr *scale, d_type_wr *sumsq);

/* Subroutine */ int dlasv2_(d_type_wr *f, d_type_wr *g, d_type_wr *h__, 
	d_type_wr *ssmin, d_type_wr *ssmax, d_type_wr *snr, d_type_wr *
	csr, d_type_wr *snl, d_type_wr *csl);

/* Subroutine */ int dlaswp_(i_type_wr *n, d_type_wr *a, i_type_wr *lda, i_type_wr 
	*k1, i_type_wr *k2, i_type_wr *ipiv, i_type_wr *incx);

/* Subroutine */ int dlasy2_(l_type_wr *ltranl, l_type_wr *ltranr, i_type_wr *isgn, 
	i_type_wr *n1, i_type_wr *n2, d_type_wr *tl, i_type_wr *ldtl, d_type_wr *
	tr, i_type_wr *ldtr, d_type_wr *b, i_type_wr *ldb, d_type_wr *scale, 
	d_type_wr *x, i_type_wr *ldx, d_type_wr *xnorm, i_type_wr *info);

/* Subroutine */ int dlasyf_(char *uplo, i_type_wr *n, i_type_wr *nb, i_type_wr *kb, 
	 d_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, d_type_wr *w, i_type_wr *
	ldw, i_type_wr *info);

/* Subroutine */ int dlat2s_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, s_type_wr *sa, i_type_wr *ldsa, i_type_wr *info);

/* Subroutine */ int dlatbs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, i_type_wr *kd, d_type_wr *ab, i_type_wr *ldab, 
	d_type_wr *x, d_type_wr *scale, d_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int dlatdf_(i_type_wr *ijob, i_type_wr *n, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *rhs, d_type_wr *rdsum, d_type_wr *rdscal, 
	i_type_wr *ipiv, i_type_wr *jpiv);

/* Subroutine */ int dlatps_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, d_type_wr *ap, d_type_wr *x, d_type_wr *scale, 
	d_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int dlatrd_(char *uplo, i_type_wr *n, i_type_wr *nb, d_type_wr *
	a, i_type_wr *lda, d_type_wr *e, d_type_wr *tau, d_type_wr *w, 
	i_type_wr *ldw);

/* Subroutine */ int dlatrs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *x, 
	d_type_wr *scale, d_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int dlatrz_(i_type_wr *m, i_type_wr *n, i_type_wr *l, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work);

/* Subroutine */ int dlatzm_(char *side, i_type_wr *m, i_type_wr *n, d_type_wr *
	v, i_type_wr *incv, d_type_wr *tau, d_type_wr *c1, d_type_wr *c2, 
	i_type_wr *ldc, d_type_wr *work);

/* Subroutine */ int dlauu2_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *info);

/* Subroutine */ int dlauum_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *info);

/* Subroutine */ int dopgtr_(char *uplo, i_type_wr *n, d_type_wr *ap, 
	d_type_wr *tau, d_type_wr *q, i_type_wr *ldq, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dopmtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, d_type_wr *ap, d_type_wr *tau, d_type_wr *c__, i_type_wr 
	*ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dorg2l_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dorg2r_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dorgbr_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dorghr_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dorgl2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dorglq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dorgql_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dorgqr_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dorgr2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dorgrq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, d_type_wr *
	a, i_type_wr *lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dorgtr_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dorm2l_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dorm2r_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dormbr_(char *vect, char *side, char *trans, i_type_wr *m, 
	i_type_wr *n, i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, 
	d_type_wr *c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dormhr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *a, i_type_wr *lda, d_type_wr *
	tau, d_type_wr *c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dorml2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dormlq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dormql_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dormqr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dormr2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dormr3_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, 
	d_type_wr *c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dormrq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dormrz_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, 
	d_type_wr *c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dormtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *tau, d_type_wr *
	c__, i_type_wr *ldc, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dpbcon_(char *uplo, i_type_wr *n, i_type_wr *kd, d_type_wr *
	ab, i_type_wr *ldab, d_type_wr *anorm, d_type_wr *rcond, d_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dpbequ_(char *uplo, i_type_wr *n, i_type_wr *kd, d_type_wr *
	ab, i_type_wr *ldab, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, 
	 i_type_wr *info);

/* Subroutine */ int dpbrfs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, d_type_wr *ab, i_type_wr *ldab, d_type_wr *afb, i_type_wr *ldafb, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int dpbstf_(char *uplo, i_type_wr *n, i_type_wr *kd, d_type_wr *
	ab, i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int dpbsv_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, d_type_wr *ab, i_type_wr *ldab, d_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int dpbsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *kd, 
	i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, d_type_wr *afb, 
	i_type_wr *ldafb, char *equed, d_type_wr *s, d_type_wr *b, i_type_wr *
	ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, 
	 d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dpbtf2_(char *uplo, i_type_wr *n, i_type_wr *kd, d_type_wr *
	ab, i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int dpbtrf_(char *uplo, i_type_wr *n, i_type_wr *kd, d_type_wr *
	ab, i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int dpbtrs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, d_type_wr *ab, i_type_wr *ldab, d_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int dpftrf_(char *transr, char *uplo, i_type_wr *n, d_type_wr 
	*a, i_type_wr *info);

/* Subroutine */ int dpftri_(char *transr, char *uplo, i_type_wr *n, d_type_wr 
	*a, i_type_wr *info);

/* Subroutine */ int dpftrs_(char *transr, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dpocon_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *anorm, d_type_wr *rcond, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dpoequ_(i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int dpoequb_(i_type_wr *n, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int dporfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int dporfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	d_type_wr *s, d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *
	ldx, d_type_wr *rcond, d_type_wr *berr, i_type_wr *n_err_bnds__, 
	d_type_wr *err_bnds_norm__, d_type_wr *err_bnds_comp__, i_type_wr *
	nparams, d_type_wr *params, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dposv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, d_type_wr 
	*a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dposvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	char *equed, d_type_wr *s, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, d_type_wr *
	berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dposvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	char *equed, d_type_wr *s, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *rpvgrw, d_type_wr *
	berr, i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, d_type_wr *
	err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, d_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dpotf2_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *info);

/* Subroutine */ int dpotrf_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *info);

/* Subroutine */ int dpotri_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *info);

/* Subroutine */ int dpotrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int dppcon_(char *uplo, i_type_wr *n, d_type_wr *ap, 
	d_type_wr *anorm, d_type_wr *rcond, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dppequ_(char *uplo, i_type_wr *n, d_type_wr *ap, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int dpprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *ap, d_type_wr *afp, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dppsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, d_type_wr 
	*ap, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dppsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *ap, d_type_wr *afp, char *equed, d_type_wr *s, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	rcond, d_type_wr *ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dpptrf_(char *uplo, i_type_wr *n, d_type_wr *ap, i_type_wr *
	info);

/* Subroutine */ int dpptri_(char *uplo, i_type_wr *n, d_type_wr *ap, i_type_wr *
	info);

/* Subroutine */ int dpptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *ap, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dpstf2_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *piv, i_type_wr *rank, d_type_wr *tol, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dpstrf_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *piv, i_type_wr *rank, d_type_wr *tol, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dptcon_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	d_type_wr *anorm, d_type_wr *rcond, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dpteqr_(char *compz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dptrfs_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *df, d_type_wr *ef, d_type_wr *b, i_type_wr 
	*ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	 d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dptsv_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dptsvx_(char *fact, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *d__, d_type_wr *e, d_type_wr *df, d_type_wr *ef, 
	d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *
	rcond, d_type_wr *ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int dpttrf_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	i_type_wr *info);

/* Subroutine */ int dpttrs_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dptts2_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int drscl_(i_type_wr *n, d_type_wr *sa, d_type_wr *sx, 
	i_type_wr *incx);

/* Subroutine */ int dsbev_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	d_type_wr *ab, i_type_wr *ldab, d_type_wr *w, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dsbevd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	d_type_wr *ab, i_type_wr *ldab, d_type_wr *w, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dsbevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *kd, d_type_wr *ab, i_type_wr *ldab, d_type_wr *q, i_type_wr *
	ldq, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, i_type_wr *iu, 
	d_type_wr *abstol, i_type_wr *m, d_type_wr *w, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *iwork, i_type_wr *ifail, 
	i_type_wr *info);

/* Subroutine */ int dsbgst_(char *vect, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, d_type_wr *ab, i_type_wr *ldab, d_type_wr *bb, i_type_wr *
	ldbb, d_type_wr *x, i_type_wr *ldx, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dsbgv_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, d_type_wr *ab, i_type_wr *ldab, d_type_wr *bb, i_type_wr *
	ldbb, d_type_wr *w, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dsbgvd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, d_type_wr *ab, i_type_wr *ldab, d_type_wr *bb, i_type_wr *
	ldbb, d_type_wr *w, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dsbgvx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *ka, i_type_wr *kb, d_type_wr *ab, i_type_wr *ldab, d_type_wr *
	bb, i_type_wr *ldbb, d_type_wr *q, i_type_wr *ldq, d_type_wr *vl, 
	d_type_wr *vu, i_type_wr *il, i_type_wr *iu, d_type_wr *abstol, i_type_wr 
	*m, d_type_wr *w, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int dsbtrd_(char *vect, char *uplo, i_type_wr *n, i_type_wr *kd, 
	d_type_wr *ab, i_type_wr *ldab, d_type_wr *d__, d_type_wr *e, 
	d_type_wr *q, i_type_wr *ldq, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dsfrk_(char *transr, char *uplo, char *trans, i_type_wr *n, 
	 i_type_wr *k, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, 
	d_type_wr *beta, d_type_wr *c__);

/* Subroutine */ int dsgesv_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	x, i_type_wr *ldx, d_type_wr *work, s_type_wr *swork, i_type_wr *iter, 
	i_type_wr *info);

/* Subroutine */ int dspcon_(char *uplo, i_type_wr *n, d_type_wr *ap, i_type_wr *
	ipiv, d_type_wr *anorm, d_type_wr *rcond, d_type_wr *work, i_type_wr 
	*iwork, i_type_wr *info);

/* Subroutine */ int dspev_(char *jobz, char *uplo, i_type_wr *n, d_type_wr *
	ap, d_type_wr *w, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dspevd_(char *jobz, char *uplo, i_type_wr *n, d_type_wr *
	ap, d_type_wr *w, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dspevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	d_type_wr *ap, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, i_type_wr *
	iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *iwork, i_type_wr *ifail, 
	i_type_wr *info);

/* Subroutine */ int dspgst_(i_type_wr *itype, char *uplo, i_type_wr *n, 
	d_type_wr *ap, d_type_wr *bp, i_type_wr *info);

/* Subroutine */ int dspgv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, d_type_wr *ap, d_type_wr *bp, d_type_wr *w, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dspgvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, d_type_wr *ap, d_type_wr *bp, d_type_wr *w, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dspgvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, d_type_wr *ap, d_type_wr *bp, d_type_wr *vl, 
	d_type_wr *vu, i_type_wr *il, i_type_wr *iu, d_type_wr *abstol, i_type_wr 
	*m, d_type_wr *w, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int dsposv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	x, i_type_wr *ldx, d_type_wr *work, s_type_wr *swork, i_type_wr *iter, 
	i_type_wr *info);

/* Subroutine */ int dsprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *ap, d_type_wr *afp, i_type_wr *ipiv, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, 
	d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dspsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, d_type_wr 
	*ap, i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dspsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *ap, d_type_wr *afp, i_type_wr *ipiv, d_type_wr *b, 
	i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, 
	d_type_wr *ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dsptrd_(char *uplo, i_type_wr *n, d_type_wr *ap, 
	d_type_wr *d__, d_type_wr *e, d_type_wr *tau, i_type_wr *info);

/* Subroutine */ int dsptrf_(char *uplo, i_type_wr *n, d_type_wr *ap, i_type_wr *
	ipiv, i_type_wr *info);

/* Subroutine */ int dsptri_(char *uplo, i_type_wr *n, d_type_wr *ap, i_type_wr *
	ipiv, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dsptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *ap, i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int dstebz_(char *range, char *order, i_type_wr *n, d_type_wr 
	*vl, d_type_wr *vu, i_type_wr *il, i_type_wr *iu, d_type_wr *abstol, 
	d_type_wr *d__, d_type_wr *e, i_type_wr *m, i_type_wr *nsplit, 
	d_type_wr *w, i_type_wr *iblock, i_type_wr *isplit, d_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dstedc_(char *compz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dstegr_(char *jobz, char *range, i_type_wr *n, d_type_wr *
	d__, d_type_wr *e, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, 
	i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, 
	d_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dstein_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	i_type_wr *m, d_type_wr *w, i_type_wr *iblock, i_type_wr *isplit, 
	d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int dstemr_(char *jobz, char *range, i_type_wr *n, d_type_wr *
	d__, d_type_wr *e, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, 
	i_type_wr *iu, i_type_wr *m, d_type_wr *w, d_type_wr *z__, i_type_wr *ldz, 
	 i_type_wr *nzc, i_type_wr *isuppz, l_type_wr *tryrac, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dsteqr_(char *compz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dsterf_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	i_type_wr *info);

/* Subroutine */ int dstev_(char *jobz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int dstevd_(char *jobz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dstevr_(char *jobz, char *range, i_type_wr *n, d_type_wr *
	d__, d_type_wr *e, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, 
	i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, 
	d_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dstevx_(char *jobz, char *range, i_type_wr *n, d_type_wr *
	d__, d_type_wr *e, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, 
	i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, 
	d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int dsycon_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *ipiv, d_type_wr *anorm, d_type_wr *rcond, d_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dsyequb_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, d_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int dsyev_(char *jobz, char *uplo, i_type_wr *n, d_type_wr *a, 
	 i_type_wr *lda, d_type_wr *w, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int dsyevd_(char *jobz, char *uplo, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, d_type_wr *w, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dsyevr_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *vl, d_type_wr *vu, i_type_wr *
	il, i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, 
	d_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dsyevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *vl, d_type_wr *vu, i_type_wr *
	il, i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, 
	d_type_wr *z__, i_type_wr *ldz, d_type_wr *work, i_type_wr *lwork, 
	i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int dsygs2_(i_type_wr *itype, char *uplo, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int dsygst_(i_type_wr *itype, char *uplo, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int dsygv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *w, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dsygvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *w, d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int dsygvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr 
	*ldb, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, i_type_wr *iu, 
	d_type_wr *abstol, i_type_wr *m, d_type_wr *w, d_type_wr *z__, 
	i_type_wr *ldz, d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, 
	i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int dsyrfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, i_type_wr *
	ipiv, d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, 
	d_type_wr *ferr, d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dsyrfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, d_type_wr *s, d_type_wr *b, i_type_wr *ldb, d_type_wr 
	*x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *berr, i_type_wr *
	n_err_bnds__, d_type_wr *err_bnds_norm__, d_type_wr *
	err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, d_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dsysv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, d_type_wr 
	*a, i_type_wr *lda, i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dsysvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, d_type_wr *b, i_type_wr *ldb, d_type_wr *x, i_type_wr *
	ldx, d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, 
	d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dsysvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, char *equed, d_type_wr *s, d_type_wr *b, i_type_wr *
	ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *
	rpvgrw, d_type_wr *berr, i_type_wr *n_err_bnds__, d_type_wr *
	err_bnds_norm__, d_type_wr *err_bnds_comp__, i_type_wr *nparams, 
	d_type_wr *params, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dsytd2_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *d__, d_type_wr *e, d_type_wr *tau, i_type_wr *info);

/* Subroutine */ int dsytf2_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int dsytrd_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *d__, d_type_wr *e, d_type_wr *tau, d_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dsytrf_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *ipiv, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dsytri_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, i_type_wr *ipiv, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dsytrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, d_type_wr *b, i_type_wr *
	ldb, i_type_wr *info);

/* Subroutine */ int dtbcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	i_type_wr *kd, d_type_wr *ab, i_type_wr *ldab, d_type_wr *rcond, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dtbrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, d_type_wr 
	*b, i_type_wr *ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, 
	d_type_wr *berr, d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dtbtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, d_type_wr *ab, i_type_wr *ldab, d_type_wr 
	*b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int dtfsm_(char *transr, char *side, char *uplo, char *trans, 
	 char *diag, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, 
	 d_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int dtftri_(char *transr, char *uplo, char *diag, i_type_wr *n, 
	 d_type_wr *a, i_type_wr *info);

/* Subroutine */ int dtfttp_(char *transr, char *uplo, i_type_wr *n, d_type_wr 
	*arf, d_type_wr *ap, i_type_wr *info);

/* Subroutine */ int dtfttr_(char *transr, char *uplo, i_type_wr *n, d_type_wr 
	*arf, d_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int dtgevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, d_type_wr *s, i_type_wr *lds, d_type_wr *p, i_type_wr *ldp, 
	d_type_wr *vl, i_type_wr *ldvl, d_type_wr *vr, i_type_wr *ldvr, i_type_wr 
	*mm, i_type_wr *m, d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dtgex2_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	q, i_type_wr *ldq, d_type_wr *z__, i_type_wr *ldz, i_type_wr *j1, i_type_wr *
	n1, i_type_wr *n2, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dtgexc_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *
	q, i_type_wr *ldq, d_type_wr *z__, i_type_wr *ldz, i_type_wr *ifst, 
	i_type_wr *ilst, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int dtgsen_(i_type_wr *ijob, l_type_wr *wantq, l_type_wr *wantz, 
	l_type_wr *select, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *
	b, i_type_wr *ldb, d_type_wr *alphar, d_type_wr *alphai, d_type_wr *
	beta, d_type_wr *q, i_type_wr *ldq, d_type_wr *z__, i_type_wr *ldz, 
	i_type_wr *m, d_type_wr *pl, d_type_wr *pr, d_type_wr *dif, 
	d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, 
	i_type_wr *info);

/* Subroutine */ int dtgsja_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, i_type_wr *k, i_type_wr *l, d_type_wr *a, 
	i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *tola, 
	d_type_wr *tolb, d_type_wr *alpha, d_type_wr *beta, d_type_wr *u, 
	i_type_wr *ldu, d_type_wr *v, i_type_wr *ldv, d_type_wr *q, i_type_wr *
	ldq, d_type_wr *work, i_type_wr *ncycle, i_type_wr *info);

/* Subroutine */ int dtgsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *vl, i_type_wr *ldvl, d_type_wr *vr, i_type_wr *ldvr, 
	d_type_wr *s, d_type_wr *dif, i_type_wr *mm, i_type_wr *m, d_type_wr *
	work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dtgsy2_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *c__, i_type_wr *ldc, d_type_wr *d__, i_type_wr *ldd, 
	d_type_wr *e, i_type_wr *lde, d_type_wr *f, i_type_wr *ldf, d_type_wr *
	scale, d_type_wr *rdsum, d_type_wr *rdscal, i_type_wr *iwork, i_type_wr 
	*pq, i_type_wr *info);

/* Subroutine */ int dtgsyl_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *c__, i_type_wr *ldc, d_type_wr *d__, i_type_wr *ldd, 
	d_type_wr *e, i_type_wr *lde, d_type_wr *f, i_type_wr *ldf, d_type_wr *
	scale, d_type_wr *dif, d_type_wr *work, i_type_wr *lwork, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dtpcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	d_type_wr *ap, d_type_wr *rcond, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int dtprfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, d_type_wr *ap, d_type_wr *b, i_type_wr *ldb, 
	d_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dtptri_(char *uplo, char *diag, i_type_wr *n, d_type_wr *
	ap, i_type_wr *info);

/* Subroutine */ int dtptrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, d_type_wr *ap, d_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int dtpttf_(char *transr, char *uplo, i_type_wr *n, d_type_wr 
	*ap, d_type_wr *arf, i_type_wr *info);

/* Subroutine */ int dtpttr_(char *uplo, i_type_wr *n, d_type_wr *ap, 
	d_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int dtrcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	d_type_wr *a, i_type_wr *lda, d_type_wr *rcond, d_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dtrevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, d_type_wr *t, i_type_wr *ldt, d_type_wr *vl, i_type_wr *
	ldvl, d_type_wr *vr, i_type_wr *ldvr, i_type_wr *mm, i_type_wr *m, 
	d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dtrexc_(char *compq, i_type_wr *n, d_type_wr *t, i_type_wr *
	ldt, d_type_wr *q, i_type_wr *ldq, i_type_wr *ifst, i_type_wr *ilst, 
	d_type_wr *work, i_type_wr *info);

/* Subroutine */ int dtrrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *
	ldb, d_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	d_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int dtrsen_(char *job, char *compq, l_type_wr *select, i_type_wr 
	*n, d_type_wr *t, i_type_wr *ldt, d_type_wr *q, i_type_wr *ldq, 
	d_type_wr *wr, d_type_wr *wi, i_type_wr *m, d_type_wr *s, d_type_wr 
	*sep, d_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *
	liwork, i_type_wr *info);

/* Subroutine */ int dtrsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, d_type_wr *t, i_type_wr *ldt, d_type_wr *vl, i_type_wr *
	ldvl, d_type_wr *vr, i_type_wr *ldvr, d_type_wr *s, d_type_wr *sep, 
	i_type_wr *mm, i_type_wr *m, d_type_wr *work, i_type_wr *ldwork, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int dtrsyl_(char *trana, char *tranb, i_type_wr *isgn, i_type_wr 
	*m, i_type_wr *n, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *
	ldb, d_type_wr *c__, i_type_wr *ldc, d_type_wr *scale, i_type_wr *info);

/* Subroutine */ int dtrti2_(char *uplo, char *diag, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int dtrtri_(char *uplo, char *diag, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int dtrtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *
	ldb, i_type_wr *info);

/* Subroutine */ int dtrttf_(char *transr, char *uplo, i_type_wr *n, d_type_wr 
	*a, i_type_wr *lda, d_type_wr *arf, i_type_wr *info);

/* Subroutine */ int dtrttp_(char *uplo, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *ap, i_type_wr *info);

/* Subroutine */ int dtzrqf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, i_type_wr *info);

/* Subroutine */ int dtzrzf_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, d_type_wr *tau, d_type_wr *work, i_type_wr *lwork, i_type_wr *info);

d_type_wr dzsum1_(i_type_wr *n, z_type_wr *cx, i_type_wr *incx);

i_type_wr icmax1_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx);

i_type_wr ieeeck_(i_type_wr *ispec, s_type_wr *zero, s_type_wr *one);

i_type_wr ilaclc_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda);

i_type_wr ilaclr_(i_type_wr *m, i_type_wr *n, c_type_wr *a, i_type_wr *lda);

i_type_wr iladiag_(char *diag);

i_type_wr iladlc_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *lda);

i_type_wr iladlr_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *lda);

i_type_wr ilaenv_(i_type_wr *ispec, char *name__, char *opts, i_type_wr *n1, 
	i_type_wr *n2, i_type_wr *n3, i_type_wr *n4);

i_type_wr ilaprec_(char *prec);

i_type_wr ilaslc_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda);

i_type_wr ilaslr_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda);

i_type_wr ilatrans_(char *trans);

i_type_wr ilauplo_(char *uplo);

/* Subroutine */ int ilaver_(i_type_wr *vers_major__, i_type_wr *vers_minor__, 
	i_type_wr *vers_patch__);

i_type_wr ilazlc_(i_type_wr *m, i_type_wr *n, z_type_wr *a, i_type_wr *lda);

i_type_wr ilazlr_(i_type_wr *m, i_type_wr *n, z_type_wr *a, i_type_wr *lda);

i_type_wr iparmq_(i_type_wr *ispec, char *name__, char *opts, i_type_wr *n, i_type_wr 
	*ilo, i_type_wr *ihi, i_type_wr *lwork);

i_type_wr izmax1_(i_type_wr *n, z_type_wr *cx, i_type_wr *incx);

l_type_wr lsamen_(i_type_wr *n, char *ca, char *cb);

i_type_wr smaxloc_(s_type_wr *a, i_type_wr *dimm);

/* Subroutine */ int sbdsdc_(char *uplo, char *compq, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *u, i_type_wr *ldu, s_type_wr *vt, i_type_wr *ldvt, s_type_wr *q, 
	i_type_wr *iq, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sbdsqr_(char *uplo, i_type_wr *n, i_type_wr *ncvt, i_type_wr *
	nru, i_type_wr *ncc, s_type_wr *d__, s_type_wr *e, s_type_wr *vt, i_type_wr *ldvt, s_type_wr *
	u, i_type_wr *ldu, s_type_wr *c__, i_type_wr *ldc, s_type_wr *work, i_type_wr *info);

d_type_wr scsum1_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx);

/* Subroutine */ int sdisna_(char *job, i_type_wr *m, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *sep, i_type_wr *info);

/* Subroutine */ int sgbbrd_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *ncc, 
	 i_type_wr *kl, i_type_wr *ku, s_type_wr *ab, i_type_wr *ldab, s_type_wr *d__, s_type_wr *
	e, s_type_wr *q, i_type_wr *ldq, s_type_wr *pt, i_type_wr *ldpt, s_type_wr *c__, i_type_wr 
	*ldc, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgbcon_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 s_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, s_type_wr *anorm, s_type_wr *rcond, 
	s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgbequ_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 s_type_wr *ab, i_type_wr *ldab, s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *
	colcnd, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int sgbequb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, s_type_wr *ab, i_type_wr *ldab, s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr 
	*colcnd, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int sgbrfs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, i_type_wr *ldafb, 
	 i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *
	ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgbrfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, 
	i_type_wr *ldafb, i_type_wr *ipiv, s_type_wr *r__, s_type_wr *c__, s_type_wr *b, i_type_wr 
	*ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *berr, i_type_wr *
	n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *
	nparams, s_type_wr *params, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgbsv_(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	nrhs, s_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int sgbsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *kl, 
	 i_type_wr *ku, i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, 
	i_type_wr *ldafb, i_type_wr *ipiv, char *equed, s_type_wr *r__, s_type_wr *c__, 
	s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, 
	 s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgbsvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, 
	i_type_wr *ldafb, i_type_wr *ipiv, char *equed, s_type_wr *r__, s_type_wr *c__, 
	s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *
	rpvgrw, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, 
	s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *params, s_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgbtf2_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 s_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int sgbtrf_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 s_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int sgbtrs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, s_type_wr *b, 
	i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sgebak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, s_type_wr *scale, i_type_wr *m, s_type_wr *v, i_type_wr *ldv, i_type_wr 
	*info);

/* Subroutine */ int sgebal_(char *job, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ilo, i_type_wr *ihi, s_type_wr *scale, i_type_wr *info);

/* Subroutine */ int sgebd2_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *d__, s_type_wr *e, s_type_wr *tauq, s_type_wr *taup, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgebrd_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *d__, s_type_wr *e, s_type_wr *tauq, s_type_wr *taup, s_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int sgecon_(char *norm, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *anorm, s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgeequ_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *colcnd, s_type_wr *amax, i_type_wr 
	*info);

/* Subroutine */ int sgeequb_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *colcnd, s_type_wr *amax, i_type_wr 
	*info);

/* Subroutine */ int sgees_(char *jobvs, char *sort, sel_fun_wr select, i_type_wr *n, 
	s_type_wr *a, i_type_wr *lda, i_type_wr *sdim, s_type_wr *wr, s_type_wr *wi, s_type_wr *vs, 
	i_type_wr *ldvs, s_type_wr *work, i_type_wr *lwork, l_type_wr *bwork, i_type_wr *
	info);

/* Subroutine */ int sgeesx_(char *jobvs, char *sort, sel_fun_wr select, char *
	sense, i_type_wr *n, s_type_wr *a, i_type_wr *lda, i_type_wr *sdim, s_type_wr *wr, 
	s_type_wr *wi, s_type_wr *vs, i_type_wr *ldvs, s_type_wr *rconde, s_type_wr *rcondv, s_type_wr *
	work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, l_type_wr *bwork, 
	 i_type_wr *info);

/* Subroutine */ int sgeev_(char *jobvl, char *jobvr, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *wr, s_type_wr *wi, s_type_wr *vl, i_type_wr *ldvl, s_type_wr *vr, 
	i_type_wr *ldvr, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *wr, s_type_wr *wi, s_type_wr *
	vl, i_type_wr *ldvl, s_type_wr *vr, i_type_wr *ldvr, i_type_wr *ilo, i_type_wr *
	ihi, s_type_wr *scale, s_type_wr *abnrm, s_type_wr *rconde, s_type_wr *rcondv, s_type_wr *work, 
	 i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgegs_(char *jobvsl, char *jobvsr, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *alphar, s_type_wr *alphai, s_type_wr 
	*beta, s_type_wr *vsl, i_type_wr *ldvsl, s_type_wr *vsr, i_type_wr *ldvsr, s_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgegv_(char *jobvl, char *jobvr, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *alphar, s_type_wr *alphai, s_type_wr 
	*beta, s_type_wr *vl, i_type_wr *ldvl, s_type_wr *vr, i_type_wr *ldvr, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgehd2_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgehrd_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgejsv_(char *joba, char *jobu, char *jobv, char *jobr, 
	char *jobt, char *jobp, i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	 s_type_wr *sva, s_type_wr *u, i_type_wr *ldu, s_type_wr *v, i_type_wr *ldv, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgelq2_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgelqf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgels_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgelsd_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *s, s_type_wr *rcond, i_type_wr *
	rank, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgelss_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *s, s_type_wr *rcond, i_type_wr *
	rank, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgelsx_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, i_type_wr *jpvt, s_type_wr *rcond, 
	i_type_wr *rank, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgelsy_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, i_type_wr *jpvt, s_type_wr *rcond, 
	i_type_wr *rank, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgeql2_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgeqlf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgeqp3_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *jpvt, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgeqpf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *jpvt, s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgeqr2_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgeqrf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgerfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *b, 
	i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, s_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgerfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, 
	s_type_wr *r__, s_type_wr *c__, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, 
	s_type_wr *rcond, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, 
	 s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *params, s_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgerq2_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgerqf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgesc2_(i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *rhs, 
	i_type_wr *ipiv, i_type_wr *jpiv, s_type_wr *scale);

/* Subroutine */ int sgesdd_(char *jobz, i_type_wr *m, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *s, s_type_wr *u, i_type_wr *ldu, s_type_wr *vt, i_type_wr *ldvt, 
	 s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgesv_(i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sgesvd_(char *jobu, char *jobvt, i_type_wr *m, i_type_wr *n, 
	s_type_wr *a, i_type_wr *lda, s_type_wr *s, s_type_wr *u, i_type_wr *ldu, s_type_wr *vt, 
	i_type_wr *ldvt, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgesvj_(char *joba, char *jobu, char *jobv, i_type_wr *m, 
	i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *sva, i_type_wr *mv, s_type_wr *v, 
	i_type_wr *ldv, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgesvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, 
	char *equed, s_type_wr *r__, s_type_wr *c__, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, 
	i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgesvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, 
	char *equed, s_type_wr *r__, s_type_wr *c__, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, 
	i_type_wr *ldx, s_type_wr *rcond, s_type_wr *rpvgrw, s_type_wr *berr, i_type_wr *
	n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *
	nparams, s_type_wr *params, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgetc2_(i_type_wr *n, s_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, 
	 i_type_wr *jpiv, i_type_wr *info);

/* Subroutine */ int sgetf2_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int sgetrf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int sgetri_(i_type_wr *n, s_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, 
	 s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgetrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sggbak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, s_type_wr *lscale, s_type_wr *rscale, i_type_wr *m, s_type_wr *v, 
	i_type_wr *ldv, i_type_wr *info);

/* Subroutine */ int sggbal_(char *job, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *b, i_type_wr *ldb, i_type_wr *ilo, i_type_wr *ihi, s_type_wr *lscale, s_type_wr 
	*rscale, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgges_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, 
	i_type_wr *sdim, s_type_wr *alphar, s_type_wr *alphai, s_type_wr *beta, s_type_wr *vsl, 
	i_type_wr *ldvsl, s_type_wr *vsr, i_type_wr *ldvsr, s_type_wr *work, i_type_wr *lwork, 
	 l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int sggesx_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, char *sense, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, 
	i_type_wr *ldb, i_type_wr *sdim, s_type_wr *alphar, s_type_wr *alphai, s_type_wr *beta, 
	s_type_wr *vsl, i_type_wr *ldvsl, s_type_wr *vsr, i_type_wr *ldvsr, s_type_wr *rconde, 
	s_type_wr *rcondv, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *
	liwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int sggev_(char *jobvl, char *jobvr, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *alphar, s_type_wr *alphai, s_type_wr 
	*beta, s_type_wr *vl, i_type_wr *ldvl, s_type_wr *vr, i_type_wr *ldvr, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr 
	*alphar, s_type_wr *alphai, s_type_wr *beta, s_type_wr *vl, i_type_wr *ldvl, s_type_wr *vr, 
	i_type_wr *ldvr, i_type_wr *ilo, i_type_wr *ihi, s_type_wr *lscale, s_type_wr *rscale, 
	 s_type_wr *abnrm, s_type_wr *bbnrm, s_type_wr *rconde, s_type_wr *rcondv, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int sggglm_(i_type_wr *n, i_type_wr *m, i_type_wr *p, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *d__, s_type_wr *x, s_type_wr *y, 
	s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgghrd_(char *compq, char *compz, i_type_wr *n, i_type_wr *
	ilo, i_type_wr *ihi, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr 
	*q, i_type_wr *ldq, s_type_wr *z__, i_type_wr *ldz, i_type_wr *info);

/* Subroutine */ int sgglse_(i_type_wr *m, i_type_wr *n, i_type_wr *p, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *c__, s_type_wr *d__, s_type_wr *x, 
	s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sggqrf_(i_type_wr *n, i_type_wr *m, i_type_wr *p, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *taua, s_type_wr *b, i_type_wr *ldb, s_type_wr *taub, s_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sggrqf_(i_type_wr *m, i_type_wr *p, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *taua, s_type_wr *b, i_type_wr *ldb, s_type_wr *taub, s_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sggsvd_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *n, i_type_wr *p, i_type_wr *k, i_type_wr *l, s_type_wr *a, i_type_wr *lda, 
	 s_type_wr *b, i_type_wr *ldb, s_type_wr *alpha, s_type_wr *beta, s_type_wr *u, i_type_wr *
	ldu, s_type_wr *v, i_type_wr *ldv, s_type_wr *q, i_type_wr *ldq, s_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sggsvp_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, 
	s_type_wr *tola, s_type_wr *tolb, i_type_wr *k, i_type_wr *l, s_type_wr *u, i_type_wr *ldu, 
	 s_type_wr *v, i_type_wr *ldv, s_type_wr *q, i_type_wr *ldq, i_type_wr *iwork, s_type_wr *
	tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sgsvj0_(char *jobv, i_type_wr *m, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *d__, s_type_wr *sva, i_type_wr *mv, s_type_wr *v, i_type_wr *
	ldv, s_type_wr *eps, s_type_wr *sfmin, s_type_wr *tol, i_type_wr *nsweep, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgsvj1_(char *jobv, i_type_wr *m, i_type_wr *n, i_type_wr *n1, 
	s_type_wr *a, i_type_wr *lda, s_type_wr *d__, s_type_wr *sva, i_type_wr *mv, s_type_wr *v, 
	i_type_wr *ldv, s_type_wr *eps, s_type_wr *sfmin, s_type_wr *tol, i_type_wr *nsweep, 
	s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sgtcon_(char *norm, i_type_wr *n, s_type_wr *dl, s_type_wr *d__, 
	s_type_wr *du, s_type_wr *du2, i_type_wr *ipiv, s_type_wr *anorm, s_type_wr *rcond, s_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgtrfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, s_type_wr *dl, 
	 s_type_wr *d__, s_type_wr *du, s_type_wr *dlf, s_type_wr *df, s_type_wr *duf, s_type_wr *du2, 
	i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *
	ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sgtsv_(i_type_wr *n, i_type_wr *nrhs, s_type_wr *dl, s_type_wr *d__, 
	s_type_wr *du, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sgtsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *dl, s_type_wr *d__, s_type_wr *du, s_type_wr *dlf, s_type_wr *df, s_type_wr *duf, 
	s_type_wr *du2, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *
	ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int sgttrf_(i_type_wr *n, s_type_wr *dl, s_type_wr *d__, s_type_wr *du, s_type_wr *
	du2, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int sgttrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, s_type_wr *dl, 
	 s_type_wr *d__, s_type_wr *du, s_type_wr *du2, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, 
	 i_type_wr *info);

/* Subroutine */ int sgtts2_(i_type_wr *itrans, i_type_wr *n, i_type_wr *nrhs, s_type_wr 
	*dl, s_type_wr *d__, s_type_wr *du, s_type_wr *du2, i_type_wr *ipiv, s_type_wr *b, i_type_wr *
	ldb);

/* Subroutine */ int shgeqz_(char *job, char *compq, char *compz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, s_type_wr *h__, i_type_wr *ldh, s_type_wr *t, i_type_wr 
	*ldt, s_type_wr *alphar, s_type_wr *alphai, s_type_wr *beta, s_type_wr *q, i_type_wr *ldq, 
	s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int shsein_(char *side, char *eigsrc, char *initv, l_type_wr *
	select, i_type_wr *n, s_type_wr *h__, i_type_wr *ldh, s_type_wr *wr, s_type_wr *wi, s_type_wr 
	*vl, i_type_wr *ldvl, s_type_wr *vr, i_type_wr *ldvr, i_type_wr *mm, i_type_wr *m, 
	s_type_wr *work, i_type_wr *ifaill, i_type_wr *ifailr, i_type_wr *info);

/* Subroutine */ int shseqr_(char *job, char *compz, i_type_wr *n, i_type_wr *ilo, 
	 i_type_wr *ihi, s_type_wr *h__, i_type_wr *ldh, s_type_wr *wr, s_type_wr *wi, s_type_wr *z__, 
	 i_type_wr *ldz, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

l_type_wr sisnan_(s_type_wr *sin__);

/* Subroutine */ int sla_gbamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *kl, i_type_wr *ku, s_type_wr *alpha, s_type_wr *ab, i_type_wr *ldab, s_type_wr *
	x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy);

d_type_wr sla_gbrcond__(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, 
	i_type_wr *cmode, s_type_wr *c__, i_type_wr *info, s_type_wr *work, i_type_wr *iwork, 
	ftn_len_wr trans_len);

/* Subroutine */ int sla_gbrfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *nrhs, 
	s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, 
	l_type_wr *colequ, s_type_wr *c__, s_type_wr *b, i_type_wr *ldb, s_type_wr *y, i_type_wr *
	ldy, s_type_wr *berr_out__, i_type_wr *n_norms__, s_type_wr *errs_n__, s_type_wr *
	errs_c__, s_type_wr *res, s_type_wr *ayb, s_type_wr *dy, s_type_wr *y_tail__, s_type_wr *rcond,
	 i_type_wr *ithresh, s_type_wr *rthresh, s_type_wr *dz_ub__, l_type_wr *
	ignore_cwise__, i_type_wr *info);

d_type_wr sla_gbrpvgrw__(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	ncols, s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, i_type_wr *ldafb);

/* Subroutine */ int sla_geamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, s_type_wr 
	*alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, 
	s_type_wr *y, i_type_wr *incy);

d_type_wr sla_gercond__(char *trans, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr 
	*af, i_type_wr *ldaf, i_type_wr *ipiv, i_type_wr *cmode, s_type_wr *c__, i_type_wr 
	*info, s_type_wr *work, i_type_wr *iwork, ftn_len_wr trans_len);

/* Subroutine */ int sla_gerfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *
	af, i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, s_type_wr *c__, s_type_wr *b,
	 i_type_wr *ldb, s_type_wr *y, i_type_wr *ldy, s_type_wr *berr_out__, i_type_wr *
	n_norms__, s_type_wr *errs_n__, s_type_wr *errs_c__, s_type_wr *res, s_type_wr *ayb, s_type_wr 
	*dy, s_type_wr *y_tail__, s_type_wr *rcond, i_type_wr *ithresh, s_type_wr *rthresh, 
	s_type_wr *dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info);

/* Subroutine */ int sla_lin_berr__(i_type_wr *n, i_type_wr *nz, i_type_wr *nrhs, 
	s_type_wr *res, s_type_wr *ayb, s_type_wr *berr);

d_type_wr sla_porcond__(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *
	af, i_type_wr *ldaf, i_type_wr *cmode, s_type_wr *c__, i_type_wr *info, s_type_wr *
	work, i_type_wr *iwork, ftn_len_wr uplo_len);

/* Subroutine */ int sla_porfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *
	ldaf, l_type_wr *colequ, s_type_wr *c__, s_type_wr *b, i_type_wr *ldb, s_type_wr *y, 
	i_type_wr *ldy, s_type_wr *berr_out__, i_type_wr *n_norms__, s_type_wr *errs_n__, 
	s_type_wr *errs_c__, s_type_wr *res, s_type_wr *ayb, s_type_wr *dy, s_type_wr *y_tail__, s_type_wr *
	rcond, i_type_wr *ithresh, s_type_wr *rthresh, s_type_wr *dz_ub__, l_type_wr *
	ignore_cwise__, i_type_wr *info, ftn_len_wr uplo_len);

d_type_wr sla_porpvgrw__(char *uplo, i_type_wr *ncols, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *af, i_type_wr *ldaf, s_type_wr *work, ftn_len_wr uplo_len);

d_type_wr sla_rpvgrw__(i_type_wr *n, i_type_wr *ncols, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *af, i_type_wr *ldaf);

/* Subroutine */ int sla_syamv__(i_type_wr *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr 
	*a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, 
	i_type_wr *incy);

d_type_wr sla_syrcond__(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *
	af, i_type_wr *ldaf, i_type_wr *ipiv, i_type_wr *cmode, s_type_wr *c__, i_type_wr *
	info, s_type_wr *work, i_type_wr *iwork, ftn_len_wr uplo_len);

/* Subroutine */ int sla_syrfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, l_type_wr *colequ, s_type_wr *c__, s_type_wr *b, i_type_wr *
	ldb, s_type_wr *y, i_type_wr *ldy, s_type_wr *berr_out__, i_type_wr *n_norms__, 
	s_type_wr *errs_n__, s_type_wr *errs_c__, s_type_wr *res, s_type_wr *ayb, s_type_wr *dy, s_type_wr *
	y_tail__, s_type_wr *rcond, i_type_wr *ithresh, s_type_wr *rthresh, s_type_wr *dz_ub__,
	 l_type_wr *ignore_cwise__, i_type_wr *info, ftn_len_wr uplo_len);

d_type_wr sla_syrpvgrw__(char *uplo, i_type_wr *n, i_type_wr *info, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *work, 
	ftn_len_wr uplo_len);

/* Subroutine */ int sla_wwaddw__(i_type_wr *n, s_type_wr *x, s_type_wr *y, s_type_wr *w);

/* Subroutine */ int slabad_(s_type_wr *small, s_type_wr *large);

/* Subroutine */ int slabrd_(i_type_wr *m, i_type_wr *n, i_type_wr *nb, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *d__, s_type_wr *e, s_type_wr *tauq, s_type_wr *taup, s_type_wr *x, 
	i_type_wr *ldx, s_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int slacn2_(i_type_wr *n, s_type_wr *v, s_type_wr *x, i_type_wr *isgn, 
	s_type_wr *est, i_type_wr *kase, i_type_wr *isave);

/* Subroutine */ int slacon_(i_type_wr *n, s_type_wr *v, s_type_wr *x, i_type_wr *isgn, 
	s_type_wr *est, i_type_wr *kase);

/* Subroutine */ int slacpy_(char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int sladiv_(s_type_wr *a, s_type_wr *b, s_type_wr *c__, s_type_wr *d__, s_type_wr *p, 
	s_type_wr *q);

/* Subroutine */ int slae2_(s_type_wr *a, s_type_wr *b, s_type_wr *c__, s_type_wr *rt1, s_type_wr *rt2);

/* Subroutine */ int slaebz_(i_type_wr *ijob, i_type_wr *nitmax, i_type_wr *n, 
	i_type_wr *mmax, i_type_wr *minp, i_type_wr *nbmin, s_type_wr *abstol, s_type_wr *
	reltol, s_type_wr *pivmin, s_type_wr *d__, s_type_wr *e, s_type_wr *e2, i_type_wr *nval, 
	s_type_wr *ab, s_type_wr *c__, i_type_wr *mout, i_type_wr *nab, s_type_wr *work, i_type_wr 
	*iwork, i_type_wr *info);

/* Subroutine */ int slaed0_(i_type_wr *icompq, i_type_wr *qsiz, i_type_wr *n, s_type_wr 
	*d__, s_type_wr *e, s_type_wr *q, i_type_wr *ldq, s_type_wr *qstore, i_type_wr *ldqs, 
	s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int slaed1_(i_type_wr *n, s_type_wr *d__, s_type_wr *q, i_type_wr *ldq, 
	i_type_wr *indxq, s_type_wr *rho, i_type_wr *cutpnt, s_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int slaed2_(i_type_wr *k, i_type_wr *n, i_type_wr *n1, s_type_wr *d__, 
	s_type_wr *q, i_type_wr *ldq, i_type_wr *indxq, s_type_wr *rho, s_type_wr *z__, s_type_wr *
	dlamda, s_type_wr *w, s_type_wr *q2, i_type_wr *indx, i_type_wr *indxc, i_type_wr *
	indxp, i_type_wr *coltyp, i_type_wr *info);

/* Subroutine */ int slaed3_(i_type_wr *k, i_type_wr *n, i_type_wr *n1, s_type_wr *d__, 
	s_type_wr *q, i_type_wr *ldq, s_type_wr *rho, s_type_wr *dlamda, s_type_wr *q2, i_type_wr *
	indx, i_type_wr *ctot, s_type_wr *w, s_type_wr *s, i_type_wr *info);

/* Subroutine */ int slaed4_(i_type_wr *n, i_type_wr *i__, s_type_wr *d__, s_type_wr *z__, 
	s_type_wr *delta, s_type_wr *rho, s_type_wr *dlam, i_type_wr *info);

/* Subroutine */ int slaed5_(i_type_wr *i__, s_type_wr *d__, s_type_wr *z__, s_type_wr *delta, 
	s_type_wr *rho, s_type_wr *dlam);

/* Subroutine */ int slaed6_(i_type_wr *kniter, l_type_wr *orgati, s_type_wr *rho, 
	s_type_wr *d__, s_type_wr *z__, s_type_wr *finit, s_type_wr *tau, i_type_wr *info);

/* Subroutine */ int slaed7_(i_type_wr *icompq, i_type_wr *n, i_type_wr *qsiz, 
	i_type_wr *tlvls, i_type_wr *curlvl, i_type_wr *curpbm, s_type_wr *d__, s_type_wr *q, 
	i_type_wr *ldq, i_type_wr *indxq, s_type_wr *rho, i_type_wr *cutpnt, s_type_wr *
	qstore, i_type_wr *qptr, i_type_wr *prmptr, i_type_wr *perm, i_type_wr *
	givptr, i_type_wr *givcol, s_type_wr *givnum, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int slaed8_(i_type_wr *icompq, i_type_wr *k, i_type_wr *n, i_type_wr 
	*qsiz, s_type_wr *d__, s_type_wr *q, i_type_wr *ldq, i_type_wr *indxq, s_type_wr *rho, 
	i_type_wr *cutpnt, s_type_wr *z__, s_type_wr *dlamda, s_type_wr *q2, i_type_wr *ldq2, 
	s_type_wr *w, i_type_wr *perm, i_type_wr *givptr, i_type_wr *givcol, s_type_wr *
	givnum, i_type_wr *indxp, i_type_wr *indx, i_type_wr *info);

/* Subroutine */ int slaed9_(i_type_wr *k, i_type_wr *kstart, i_type_wr *kstop, 
	i_type_wr *n, s_type_wr *d__, s_type_wr *q, i_type_wr *ldq, s_type_wr *rho, s_type_wr *dlamda, 
	 s_type_wr *w, s_type_wr *s, i_type_wr *lds, i_type_wr *info);

/* Subroutine */ int slaeda_(i_type_wr *n, i_type_wr *tlvls, i_type_wr *curlvl, 
	i_type_wr *curpbm, i_type_wr *prmptr, i_type_wr *perm, i_type_wr *givptr, 
	i_type_wr *givcol, s_type_wr *givnum, s_type_wr *q, i_type_wr *qptr, s_type_wr *z__, 
	s_type_wr *ztemp, i_type_wr *info);

/* Subroutine */ int slaein_(l_type_wr *rightv, l_type_wr *noinit, i_type_wr *n, 
	s_type_wr *h__, i_type_wr *ldh, s_type_wr *wr, s_type_wr *wi, s_type_wr *vr, s_type_wr *vi, s_type_wr 
	*b, i_type_wr *ldb, s_type_wr *work, s_type_wr *eps3, s_type_wr *smlnum, s_type_wr *bignum, 
	i_type_wr *info);

/* Subroutine */ int slaev2_(s_type_wr *a, s_type_wr *b, s_type_wr *c__, s_type_wr *rt1, s_type_wr *
	rt2, s_type_wr *cs1, s_type_wr *sn1);

/* Subroutine */ int slaexc_(l_type_wr *wantq, i_type_wr *n, s_type_wr *t, i_type_wr *
	ldt, s_type_wr *q, i_type_wr *ldq, i_type_wr *j1, i_type_wr *n1, i_type_wr *n2, 
	s_type_wr *work, i_type_wr *info);

/* Subroutine */ int slag2_(s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, 
	s_type_wr *safmin, s_type_wr *scale1, s_type_wr *scale2, s_type_wr *wr1, s_type_wr *wr2, s_type_wr *
	wi);

/* Subroutine */ int slag2d_(i_type_wr *m, i_type_wr *n, s_type_wr *sa, i_type_wr *ldsa, 
	d_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int slags2_(l_type_wr *upper, s_type_wr *a1, s_type_wr *a2, s_type_wr *a3, 
	s_type_wr *b1, s_type_wr *b2, s_type_wr *b3, s_type_wr *csu, s_type_wr *snu, s_type_wr *csv, s_type_wr *
	snv, s_type_wr *csq, s_type_wr *snq);

/* Subroutine */ int slagtf_(i_type_wr *n, s_type_wr *a, s_type_wr *lambda, s_type_wr *b, s_type_wr 
	*c__, s_type_wr *tol, s_type_wr *d__, i_type_wr *in, i_type_wr *info);

/* Subroutine */ int slagtm_(char *trans, i_type_wr *n, i_type_wr *nrhs, s_type_wr *
	alpha, s_type_wr *dl, s_type_wr *d__, s_type_wr *du, s_type_wr *x, i_type_wr *ldx, s_type_wr *
	beta, s_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int slagts_(i_type_wr *job, i_type_wr *n, s_type_wr *a, s_type_wr *b, s_type_wr 
	*c__, s_type_wr *d__, i_type_wr *in, s_type_wr *y, s_type_wr *tol, i_type_wr *info);

/* Subroutine */ int slagv2_(s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, 
	s_type_wr *alphar, s_type_wr *alphai, s_type_wr *beta, s_type_wr *csl, s_type_wr *snl, s_type_wr *
	csr, s_type_wr *snr);

/* Subroutine */ int slahqr_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, s_type_wr *h__, i_type_wr *ldh, s_type_wr *wr, s_type_wr *
	wi, i_type_wr *iloz, i_type_wr *ihiz, s_type_wr *z__, i_type_wr *ldz, i_type_wr *
	info);

/* Subroutine */ int slahr2_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *t, i_type_wr *ldt, s_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int slahrd_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *t, i_type_wr *ldt, s_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int slaic1_(i_type_wr *job, i_type_wr *j, s_type_wr *x, s_type_wr *sest, 
	s_type_wr *w, s_type_wr *gamma, s_type_wr *sestpr, s_type_wr *s, s_type_wr *c__);

l_type_wr slaisnan_(s_type_wr *sin1, s_type_wr *sin2);

/* Subroutine */ int slaln2_(l_type_wr *ltrans, i_type_wr *na, i_type_wr *nw, s_type_wr *
	smin, s_type_wr *ca, s_type_wr *a, i_type_wr *lda, s_type_wr *d1, s_type_wr *d2, s_type_wr *b, 
	i_type_wr *ldb, s_type_wr *wr, s_type_wr *wi, s_type_wr *x, i_type_wr *ldx, s_type_wr *scale, 
	s_type_wr *xnorm, i_type_wr *info);

/* Subroutine */ int slals0_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, i_type_wr *nrhs, s_type_wr *b, i_type_wr *ldb, s_type_wr *bx, 
	i_type_wr *ldbx, i_type_wr *perm, i_type_wr *givptr, i_type_wr *givcol, 
	i_type_wr *ldgcol, s_type_wr *givnum, i_type_wr *ldgnum, s_type_wr *poles, s_type_wr *
	difl, s_type_wr *difr, s_type_wr *z__, i_type_wr *k, s_type_wr *c__, s_type_wr *s, s_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int slalsa_(i_type_wr *icompq, i_type_wr *smlsiz, i_type_wr *n, 
	i_type_wr *nrhs, s_type_wr *b, i_type_wr *ldb, s_type_wr *bx, i_type_wr *ldbx, s_type_wr *
	u, i_type_wr *ldu, s_type_wr *vt, i_type_wr *k, s_type_wr *difl, s_type_wr *difr, s_type_wr *
	z__, s_type_wr *poles, i_type_wr *givptr, i_type_wr *givcol, i_type_wr *ldgcol, 
	i_type_wr *perm, s_type_wr *givnum, s_type_wr *c__, s_type_wr *s, s_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int slalsd_(char *uplo, i_type_wr *smlsiz, i_type_wr *n, i_type_wr 
	*nrhs, s_type_wr *d__, s_type_wr *e, s_type_wr *b, i_type_wr *ldb, s_type_wr *rcond, 
	i_type_wr *rank, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int slamrg_(i_type_wr *n1, i_type_wr *n2, s_type_wr *a, i_type_wr *
	strd1, i_type_wr *strd2, i_type_wr *index);

i_type_wr slaneg_(i_type_wr *n, s_type_wr *d__, s_type_wr *lld, s_type_wr *sigma, s_type_wr *pivmin, 
	i_type_wr *r__);

s_type_wr slangb_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, s_type_wr *ab, 
	 i_type_wr *ldab, s_type_wr *work);

s_type_wr slange_(char *norm, i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *work);

s_type_wr slangt_(char *norm, i_type_wr *n, s_type_wr *dl, s_type_wr *d__, s_type_wr *du);

s_type_wr slanhs_(char *norm, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *work);

s_type_wr slansb_(char *norm, char *uplo, i_type_wr *n, i_type_wr *k, s_type_wr *ab, 
	i_type_wr *ldab, s_type_wr *work);

s_type_wr slansf_(char *norm, char *transr, char *uplo, i_type_wr *n, s_type_wr *a, 
	s_type_wr *work);

s_type_wr slansp_(char *norm, char *uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *work);

s_type_wr slanst_(char *norm, i_type_wr *n, s_type_wr *d__, s_type_wr *e);

s_type_wr slansy_(char *norm, char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *work);

s_type_wr slantb_(char *norm, char *uplo, char *diag, i_type_wr *n, i_type_wr *k, 
	 s_type_wr *ab, i_type_wr *ldab, s_type_wr *work);

s_type_wr slantp_(char *norm, char *uplo, char *diag, i_type_wr *n, s_type_wr *ap, 
	s_type_wr *work);

s_type_wr slantr_(char *norm, char *uplo, char *diag, i_type_wr *m, i_type_wr *n, 
	 s_type_wr *a, i_type_wr *lda, s_type_wr *work);

/* Subroutine */ int slanv2_(s_type_wr *a, s_type_wr *b, s_type_wr *c__, s_type_wr *d__, s_type_wr *
	rt1r, s_type_wr *rt1i, s_type_wr *rt2r, s_type_wr *rt2i, s_type_wr *cs, s_type_wr *sn);

/* Subroutine */ int slapll_(i_type_wr *n, s_type_wr *x, i_type_wr *incx, s_type_wr *y, 
	i_type_wr *incy, s_type_wr *ssmin);

/* Subroutine */ int slapmt_(l_type_wr *forwrd, i_type_wr *m, i_type_wr *n, s_type_wr *x, 
	 i_type_wr *ldx, i_type_wr *k);

s_type_wr slapy2_(s_type_wr *x, s_type_wr *y);

s_type_wr slapy3_(s_type_wr *x, s_type_wr *y, s_type_wr *z__);

/* Subroutine */ int slaqgb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 s_type_wr *ab, i_type_wr *ldab, s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *
	colcnd, s_type_wr *amax, char *equed);

/* Subroutine */ int slaqge_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *r__, s_type_wr *c__, s_type_wr *rowcnd, s_type_wr *colcnd, s_type_wr *amax, char *
	equed);

/* Subroutine */ int slaqp2_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, s_type_wr *a, 
	 i_type_wr *lda, i_type_wr *jpvt, s_type_wr *tau, s_type_wr *vn1, s_type_wr *vn2, s_type_wr *
	work);

/* Subroutine */ int slaqps_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, i_type_wr 
	*nb, i_type_wr *kb, s_type_wr *a, i_type_wr *lda, i_type_wr *jpvt, s_type_wr *tau, 
	s_type_wr *vn1, s_type_wr *vn2, s_type_wr *auxv, s_type_wr *f, i_type_wr *ldf);

/* Subroutine */ int slaqr0_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, s_type_wr *h__, i_type_wr *ldh, s_type_wr *wr, s_type_wr *
	wi, i_type_wr *iloz, i_type_wr *ihiz, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, 
	 i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int slaqr1_(i_type_wr *n, s_type_wr *h__, i_type_wr *ldh, s_type_wr *sr1, 
	s_type_wr *si1, s_type_wr *sr2, s_type_wr *si2, s_type_wr *v);

/* Subroutine */ int slaqr2_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, s_type_wr *h__, i_type_wr *ldh, 
	i_type_wr *iloz, i_type_wr *ihiz, s_type_wr *z__, i_type_wr *ldz, i_type_wr *ns, 
	i_type_wr *nd, s_type_wr *sr, s_type_wr *si, s_type_wr *v, i_type_wr *ldv, i_type_wr *nh, 
	s_type_wr *t, i_type_wr *ldt, i_type_wr *nv, s_type_wr *wv, i_type_wr *ldwv, s_type_wr *
	work, i_type_wr *lwork);

/* Subroutine */ int slaqr3_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, s_type_wr *h__, i_type_wr *ldh, 
	i_type_wr *iloz, i_type_wr *ihiz, s_type_wr *z__, i_type_wr *ldz, i_type_wr *ns, 
	i_type_wr *nd, s_type_wr *sr, s_type_wr *si, s_type_wr *v, i_type_wr *ldv, i_type_wr *nh, 
	s_type_wr *t, i_type_wr *ldt, i_type_wr *nv, s_type_wr *wv, i_type_wr *ldwv, s_type_wr *
	work, i_type_wr *lwork);

/* Subroutine */ int slaqr4_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, s_type_wr *h__, i_type_wr *ldh, s_type_wr *wr, s_type_wr *
	wi, i_type_wr *iloz, i_type_wr *ihiz, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, 
	 i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int slaqr5_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *kacc22, 
	i_type_wr *n, i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nshfts, s_type_wr *sr, 
	s_type_wr *si, s_type_wr *h__, i_type_wr *ldh, i_type_wr *iloz, i_type_wr *ihiz, s_type_wr 
	*z__, i_type_wr *ldz, s_type_wr *v, i_type_wr *ldv, s_type_wr *u, i_type_wr *ldu, 
	i_type_wr *nv, s_type_wr *wv, i_type_wr *ldwv, i_type_wr *nh, s_type_wr *wh, i_type_wr *
	ldwh);

/* Subroutine */ int slaqsb_(char *uplo, i_type_wr *n, i_type_wr *kd, s_type_wr *ab, 
	i_type_wr *ldab, s_type_wr *s, s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int slaqsp_(char *uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *s, s_type_wr *
	scond, s_type_wr *amax, char *equed);

/* Subroutine */ int slaqsy_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *s, s_type_wr *scond, s_type_wr *amax, char *equed);

/* Subroutine */ int slaqtr_(l_type_wr *ltran, l_type_wr *ls_type_wr, i_type_wr *n, s_type_wr 
	*t, i_type_wr *ldt, s_type_wr *b, s_type_wr *w, s_type_wr *scale, s_type_wr *x, s_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int slar1v_(i_type_wr *n, i_type_wr *b1, i_type_wr *bn, s_type_wr *
	lambda, s_type_wr *d__, s_type_wr *l, s_type_wr *ld, s_type_wr *lld, s_type_wr *pivmin, s_type_wr *
	gaptol, s_type_wr *z__, l_type_wr *wantnc, i_type_wr *negcnt, s_type_wr *ztz, s_type_wr *
	mingma, i_type_wr *r__, i_type_wr *isuppz, s_type_wr *nrminv, s_type_wr *resid, 
	s_type_wr *rqcorr, s_type_wr *work);

/* Subroutine */ int slar2v_(i_type_wr *n, s_type_wr *x, s_type_wr *y, s_type_wr *z__, i_type_wr 
	*incx, s_type_wr *c__, s_type_wr *s, i_type_wr *incc);

/* Subroutine */ int slarf_(char *side, i_type_wr *m, i_type_wr *n, s_type_wr *v, 
	i_type_wr *incv, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, s_type_wr *work);

/* Subroutine */ int slarfb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *v, i_type_wr *ldv, 
	s_type_wr *t, i_type_wr *ldt, s_type_wr *c__, i_type_wr *ldc, s_type_wr *work, i_type_wr *
	ldwork);

/* Subroutine */ int slarfg_(i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, 
	s_type_wr *tau);

/* Subroutine */ int slarfp_(i_type_wr *n, s_type_wr *alpha, s_type_wr *x, i_type_wr *incx, 
	s_type_wr *tau);

/* Subroutine */ int slarft_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, s_type_wr *v, i_type_wr *ldv, s_type_wr *tau, s_type_wr *t, i_type_wr *ldt);

/* Subroutine */ int slarfx_(char *side, i_type_wr *m, i_type_wr *n, s_type_wr *v, 
	s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, s_type_wr *work);

/* Subroutine */ int slargv_(i_type_wr *n, s_type_wr *x, i_type_wr *incx, s_type_wr *y, 
	i_type_wr *incy, s_type_wr *c__, i_type_wr *incc);

/* Subroutine */ int slarnv_(i_type_wr *idist, i_type_wr *iseed, i_type_wr *n, s_type_wr 
	*x);

/* Subroutine */ int slarra_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, s_type_wr *e2, s_type_wr *
	spltol, s_type_wr *tnrm, i_type_wr *nsplit, i_type_wr *isplit, i_type_wr *info);

/* Subroutine */ int slarrb_(i_type_wr *n, s_type_wr *d__, s_type_wr *lld, i_type_wr *
	ifirst, i_type_wr *ilast, s_type_wr *rtol1, s_type_wr *rtol2, i_type_wr *offset, 
	s_type_wr *w, s_type_wr *wgap, s_type_wr *werr, s_type_wr *work, i_type_wr *iwork, s_type_wr *
	pivmin, s_type_wr *spdiam, i_type_wr *twist, i_type_wr *info);

/* Subroutine */ int slarrc_(char *jobt, i_type_wr *n, s_type_wr *vl, s_type_wr *vu, s_type_wr 
	*d__, s_type_wr *e, s_type_wr *pivmin, i_type_wr *eigcnt, i_type_wr *lcnt, i_type_wr *
	rcnt, i_type_wr *info);

/* Subroutine */ int slarrd_(char *range, char *order, i_type_wr *n, s_type_wr *vl, 
	s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *gers, s_type_wr *reltol, s_type_wr *
	d__, s_type_wr *e, s_type_wr *e2, s_type_wr *pivmin, i_type_wr *nsplit, i_type_wr *
	isplit, i_type_wr *m, s_type_wr *w, s_type_wr *werr, s_type_wr *wl, s_type_wr *wu, i_type_wr *
	iblock, i_type_wr *indexw, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int slarre_(char *range, i_type_wr *n, s_type_wr *vl, s_type_wr *vu, 
	i_type_wr *il, i_type_wr *iu, s_type_wr *d__, s_type_wr *e, s_type_wr *e2, s_type_wr *rtol1, 
	s_type_wr *rtol2, s_type_wr *spltol, i_type_wr *nsplit, i_type_wr *isplit, i_type_wr *
	m, s_type_wr *w, s_type_wr *werr, s_type_wr *wgap, i_type_wr *iblock, i_type_wr *indexw, 
	s_type_wr *gers, s_type_wr *pivmin, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int slarrf_(i_type_wr *n, s_type_wr *d__, s_type_wr *l, s_type_wr *ld, 
	i_type_wr *clstrt, i_type_wr *clend, s_type_wr *w, s_type_wr *wgap, s_type_wr *werr, 
	s_type_wr *spdiam, s_type_wr *clgapl, s_type_wr *clgapr, s_type_wr *pivmin, s_type_wr *sigma, 
	s_type_wr *dplus, s_type_wr *lplus, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int slarrj_(i_type_wr *n, s_type_wr *d__, s_type_wr *e2, i_type_wr *ifirst, 
	 i_type_wr *ilast, s_type_wr *rtol, i_type_wr *offset, s_type_wr *w, s_type_wr *werr, 
	s_type_wr *work, i_type_wr *iwork, s_type_wr *pivmin, s_type_wr *spdiam, i_type_wr *info);

/* Subroutine */ int slarrk_(i_type_wr *n, i_type_wr *iw, s_type_wr *gl, s_type_wr *gu, 
	s_type_wr *d__, s_type_wr *e2, s_type_wr *pivmin, s_type_wr *reltol, s_type_wr *w, s_type_wr *werr, 
	i_type_wr *info);

/* Subroutine */ int slarrr_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, i_type_wr *info);

/* Subroutine */ int slarrv_(i_type_wr *n, s_type_wr *vl, s_type_wr *vu, s_type_wr *d__, s_type_wr *
	l, s_type_wr *pivmin, i_type_wr *isplit, i_type_wr *m, i_type_wr *dol, i_type_wr *
	dou, s_type_wr *minrgp, s_type_wr *rtol1, s_type_wr *rtol2, s_type_wr *w, s_type_wr *werr, 
	s_type_wr *wgap, i_type_wr *iblock, i_type_wr *indexw, s_type_wr *gers, s_type_wr *z__, 
	i_type_wr *ldz, i_type_wr *isuppz, s_type_wr *work, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int slarscl2_(i_type_wr *m, i_type_wr *n, s_type_wr *d__, s_type_wr *x, 
	i_type_wr *ldx);

/* Subroutine */ int slartg_(s_type_wr *f, s_type_wr *g, s_type_wr *cs, s_type_wr *sn, s_type_wr *r__);

/* Subroutine */ int slartv_(i_type_wr *n, s_type_wr *x, i_type_wr *incx, s_type_wr *y, 
	i_type_wr *incy, s_type_wr *c__, s_type_wr *s, i_type_wr *incc);

/* Subroutine */ int slaruv_(i_type_wr *iseed, i_type_wr *n, s_type_wr *x);

/* Subroutine */ int slarz_(char *side, i_type_wr *m, i_type_wr *n, i_type_wr *l, 
	s_type_wr *v, i_type_wr *incv, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, s_type_wr *
	work);

/* Subroutine */ int slarzb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, i_type_wr *l, s_type_wr *v, 
	i_type_wr *ldv, s_type_wr *t, i_type_wr *ldt, s_type_wr *c__, i_type_wr *ldc, s_type_wr *
	work, i_type_wr *ldwork);

/* Subroutine */ int slarzt_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, s_type_wr *v, i_type_wr *ldv, s_type_wr *tau, s_type_wr *t, i_type_wr *ldt);

/* Subroutine */ int slas2_(s_type_wr *f, s_type_wr *g, s_type_wr *h__, s_type_wr *ssmin, s_type_wr *
	ssmax);

/* Subroutine */ int slascl_(char *type__, i_type_wr *kl, i_type_wr *ku, s_type_wr *
	cfrom, s_type_wr *cto, i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *info);

/* Subroutine */ int slascl2_(i_type_wr *m, i_type_wr *n, s_type_wr *d__, s_type_wr *x, 
	i_type_wr *ldx);

/* Subroutine */ int slasd0_(i_type_wr *n, i_type_wr *sqre, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *u, i_type_wr *ldu, s_type_wr *vt, i_type_wr *ldvt, i_type_wr *smlsiz, 
	i_type_wr *iwork, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int slasd1_(i_type_wr *nl, i_type_wr *nr, i_type_wr *sqre, s_type_wr *
	d__, s_type_wr *alpha, s_type_wr *beta, s_type_wr *u, i_type_wr *ldu, s_type_wr *vt, 
	i_type_wr *ldvt, i_type_wr *idxq, i_type_wr *iwork, s_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int slasd2_(i_type_wr *nl, i_type_wr *nr, i_type_wr *sqre, i_type_wr 
	*k, s_type_wr *d__, s_type_wr *z__, s_type_wr *alpha, s_type_wr *beta, s_type_wr *u, i_type_wr *
	ldu, s_type_wr *vt, i_type_wr *ldvt, s_type_wr *dsigma, s_type_wr *u2, i_type_wr *ldu2, 
	s_type_wr *vt2, i_type_wr *ldvt2, i_type_wr *idxp, i_type_wr *idx, i_type_wr *idxc, 
	 i_type_wr *idxq, i_type_wr *coltyp, i_type_wr *info);

/* Subroutine */ int slasd3_(i_type_wr *nl, i_type_wr *nr, i_type_wr *sqre, i_type_wr 
	*k, s_type_wr *d__, s_type_wr *q, i_type_wr *ldq, s_type_wr *dsigma, s_type_wr *u, i_type_wr *
	ldu, s_type_wr *u2, i_type_wr *ldu2, s_type_wr *vt, i_type_wr *ldvt, s_type_wr *vt2, 
	i_type_wr *ldvt2, i_type_wr *idxc, i_type_wr *ctot, s_type_wr *z__, i_type_wr *
	info);

/* Subroutine */ int slasd4_(i_type_wr *n, i_type_wr *i__, s_type_wr *d__, s_type_wr *z__, 
	s_type_wr *delta, s_type_wr *rho, s_type_wr *sigma, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int slasd5_(i_type_wr *i__, s_type_wr *d__, s_type_wr *z__, s_type_wr *delta, 
	s_type_wr *rho, s_type_wr *dsigma, s_type_wr *work);

/* Subroutine */ int slasd6_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, s_type_wr *d__, s_type_wr *vf, s_type_wr *vl, s_type_wr *alpha, s_type_wr *beta, 
	 i_type_wr *idxq, i_type_wr *perm, i_type_wr *givptr, i_type_wr *givcol, 
	i_type_wr *ldgcol, s_type_wr *givnum, i_type_wr *ldgnum, s_type_wr *poles, s_type_wr *
	difl, s_type_wr *difr, s_type_wr *z__, i_type_wr *k, s_type_wr *c__, s_type_wr *s, s_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int slasd7_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, i_type_wr *k, s_type_wr *d__, s_type_wr *z__, s_type_wr *zw, s_type_wr *vf, 
	s_type_wr *vfw, s_type_wr *vl, s_type_wr *vlw, s_type_wr *alpha, s_type_wr *beta, s_type_wr *dsigma, 
	 i_type_wr *idx, i_type_wr *idxp, i_type_wr *idxq, i_type_wr *perm, i_type_wr *
	givptr, i_type_wr *givcol, i_type_wr *ldgcol, s_type_wr *givnum, i_type_wr *
	ldgnum, s_type_wr *c__, s_type_wr *s, i_type_wr *info);

/* Subroutine */ int slasd8_(i_type_wr *icompq, i_type_wr *k, s_type_wr *d__, s_type_wr *
	z__, s_type_wr *vf, s_type_wr *vl, s_type_wr *difl, s_type_wr *difr, i_type_wr *lddifr, 
	s_type_wr *dsigma, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int slasda_(i_type_wr *icompq, i_type_wr *smlsiz, i_type_wr *n, 
	i_type_wr *sqre, s_type_wr *d__, s_type_wr *e, s_type_wr *u, i_type_wr *ldu, s_type_wr *vt, 
	i_type_wr *k, s_type_wr *difl, s_type_wr *difr, s_type_wr *z__, s_type_wr *poles, i_type_wr *
	givptr, i_type_wr *givcol, i_type_wr *ldgcol, i_type_wr *perm, s_type_wr *givnum, 
	 s_type_wr *c__, s_type_wr *s, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int slasdq_(char *uplo, i_type_wr *sqre, i_type_wr *n, i_type_wr *
	ncvt, i_type_wr *nru, i_type_wr *ncc, s_type_wr *d__, s_type_wr *e, s_type_wr *vt, 
	i_type_wr *ldvt, s_type_wr *u, i_type_wr *ldu, s_type_wr *c__, i_type_wr *ldc, s_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int slasdt_(i_type_wr *n, i_type_wr *lvl, i_type_wr *nd, i_type_wr *
	inode, i_type_wr *ndiml, i_type_wr *ndimr, i_type_wr *msub);

/* Subroutine */ int slaset_(char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, 
	s_type_wr *beta, s_type_wr *a, i_type_wr *lda);

/* Subroutine */ int slasq1_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, s_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int slasq2_(i_type_wr *n, s_type_wr *z__, i_type_wr *info);

/* Subroutine */ int slasq3_(i_type_wr *i0, i_type_wr *n0, s_type_wr *z__, i_type_wr *pp, 
	 s_type_wr *dmin__, s_type_wr *sigma, s_type_wr *desig, s_type_wr *qmax, i_type_wr *nfail, 
	i_type_wr *iter, i_type_wr *ndiv, l_type_wr *ieee, i_type_wr *ttype, s_type_wr *
	dmin1, s_type_wr *dmin2, s_type_wr *dn, s_type_wr *dn1, s_type_wr *dn2, s_type_wr *g, s_type_wr *
	tau);

/* Subroutine */ int slasq4_(i_type_wr *i0, i_type_wr *n0, s_type_wr *z__, i_type_wr *pp, 
	 i_type_wr *n0in, s_type_wr *dmin__, s_type_wr *dmin1, s_type_wr *dmin2, s_type_wr *dn, 
	s_type_wr *dn1, s_type_wr *dn2, s_type_wr *tau, i_type_wr *ttype, s_type_wr *g);

/* Subroutine */ int slasq5_(i_type_wr *i0, i_type_wr *n0, s_type_wr *z__, i_type_wr *pp, 
	 s_type_wr *tau, s_type_wr *dmin__, s_type_wr *dmin1, s_type_wr *dmin2, s_type_wr *dn, s_type_wr *
	dnm1, s_type_wr *dnm2, l_type_wr *ieee);

/* Subroutine */ int slasq6_(i_type_wr *i0, i_type_wr *n0, s_type_wr *z__, i_type_wr *pp, 
	 s_type_wr *dmin__, s_type_wr *dmin1, s_type_wr *dmin2, s_type_wr *dn, s_type_wr *dnm1, s_type_wr *
	dnm2);

/* Subroutine */ int slasr_(char *side, char *pivot, char *direct, i_type_wr *m, 
	 i_type_wr *n, s_type_wr *c__, s_type_wr *s, s_type_wr *a, i_type_wr *lda);

/* Subroutine */ int slasrt_(char *id, i_type_wr *n, s_type_wr *d__, i_type_wr *info);

/* Subroutine */ int slassq_(i_type_wr *n, s_type_wr *x, i_type_wr *incx, s_type_wr *scale, 
	s_type_wr *sumsq);

/* Subroutine */ int slasv2_(s_type_wr *f, s_type_wr *g, s_type_wr *h__, s_type_wr *ssmin, s_type_wr *
	ssmax, s_type_wr *snr, s_type_wr *csr, s_type_wr *snl, s_type_wr *csl);

/* Subroutine */ int slaswp_(i_type_wr *n, s_type_wr *a, i_type_wr *lda, i_type_wr *k1, 
	i_type_wr *k2, i_type_wr *ipiv, i_type_wr *incx);

/* Subroutine */ int slasy2_(l_type_wr *ltranl, l_type_wr *ltranr, i_type_wr *isgn, 
	i_type_wr *n1, i_type_wr *n2, s_type_wr *tl, i_type_wr *ldtl, s_type_wr *tr, i_type_wr *
	ldtr, s_type_wr *b, i_type_wr *ldb, s_type_wr *scale, s_type_wr *x, i_type_wr *ldx, s_type_wr 
	*xnorm, i_type_wr *info);

/* Subroutine */ int slasyf_(char *uplo, i_type_wr *n, i_type_wr *nb, i_type_wr *kb, 
	 s_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, s_type_wr *w, i_type_wr *ldw, i_type_wr 
	*info);

/* Subroutine */ int slatbs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, i_type_wr *kd, s_type_wr *ab, i_type_wr *ldab, s_type_wr *x, 
	s_type_wr *scale, s_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int slatdf_(i_type_wr *ijob, i_type_wr *n, s_type_wr *z__, i_type_wr *
	ldz, s_type_wr *rhs, s_type_wr *rdsum, s_type_wr *rdscal, i_type_wr *ipiv, i_type_wr *
	jpiv);

/* Subroutine */ int slatps_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, s_type_wr *ap, s_type_wr *x, s_type_wr *scale, s_type_wr *cnorm, 
	i_type_wr *info);

/* Subroutine */ int slatrd_(char *uplo, i_type_wr *n, i_type_wr *nb, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *e, s_type_wr *tau, s_type_wr *w, i_type_wr *ldw);

/* Subroutine */ int slatrs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *x, s_type_wr *scale, s_type_wr 
	*cnorm, i_type_wr *info);

/* Subroutine */ int slatrz_(i_type_wr *m, i_type_wr *n, i_type_wr *l, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work);

/* Subroutine */ int slatzm_(char *side, i_type_wr *m, i_type_wr *n, s_type_wr *v, 
	i_type_wr *incv, s_type_wr *tau, s_type_wr *c1, s_type_wr *c2, i_type_wr *ldc, s_type_wr *
	work);

/* Subroutine */ int slauu2_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *info);

/* Subroutine */ int slauum_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *info);

/* Subroutine */ int sopgtr_(char *uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *tau, 
	s_type_wr *q, i_type_wr *ldq, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sopmtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, s_type_wr *ap, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, s_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int sorg2l_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sorg2r_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sorgbr_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr 
	*info);

/* Subroutine */ int sorghr_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sorgl2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sorglq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sorgql_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sorgqr_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sorgr2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sorgrq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sorgtr_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sorm2l_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sorm2r_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sormbr_(char *vect, char *side, char *trans, i_type_wr *m, 
	i_type_wr *n, i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, 
	i_type_wr *ldc, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sormhr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *
	c__, i_type_wr *ldc, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sorml2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sormlq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sormql_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sormqr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sormr2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sormr3_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, 
	i_type_wr *ldc, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sormrq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sormrz_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, 
	i_type_wr *ldc, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int sormtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *tau, s_type_wr *c__, i_type_wr *ldc, 
	 s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int spbcon_(char *uplo, i_type_wr *n, i_type_wr *kd, s_type_wr *ab, 
	i_type_wr *ldab, s_type_wr *anorm, s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int spbequ_(char *uplo, i_type_wr *n, i_type_wr *kd, s_type_wr *ab, 
	i_type_wr *ldab, s_type_wr *s, s_type_wr *scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int spbrfs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, i_type_wr *ldafb, s_type_wr *b, 
	i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, s_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int spbstf_(char *uplo, i_type_wr *n, i_type_wr *kd, s_type_wr *ab, 
	i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int spbsv_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int spbsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *kd, 
	i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *afb, i_type_wr *ldafb, 
	char *equed, s_type_wr *s, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, 
	s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int spbtf2_(char *uplo, i_type_wr *n, i_type_wr *kd, s_type_wr *ab, 
	i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int spbtrf_(char *uplo, i_type_wr *n, i_type_wr *kd, s_type_wr *ab, 
	i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int spbtrs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int spftrf_(char *transr, char *uplo, i_type_wr *n, s_type_wr *a, 
	i_type_wr *info);

/* Subroutine */ int spftri_(char *transr, char *uplo, i_type_wr *n, s_type_wr *a, 
	i_type_wr *info);

/* Subroutine */ int spftrs_(char *transr, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int spocon_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *anorm, s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int spoequ_(i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *s, s_type_wr 
	*scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int spoequb_(i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *s, 
	s_type_wr *scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int sporfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, 
	 i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int sporfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, s_type_wr *s, s_type_wr *
	b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *berr, 
	i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *err_bnds_comp__, 
	i_type_wr *nparams, s_type_wr *params, s_type_wr *work, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int sposv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sposvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, char *equed, 
	s_type_wr *s, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, 
	s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sposvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, char *equed, 
	s_type_wr *s, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, 
	s_type_wr *rpvgrw, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *
	err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *
	params, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int spotf2_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *info);

/* Subroutine */ int spotrf_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *info);

/* Subroutine */ int spotri_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *info);

/* Subroutine */ int spotrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sppcon_(char *uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *anorm, 
	s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sppequ_(char *uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *s, s_type_wr *
	scond, s_type_wr *amax, i_type_wr *info);

/* Subroutine */ int spprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *ap, 
	s_type_wr *afp, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, 
	s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sppsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *ap, 
	s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sppsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *ap, s_type_wr *afp, char *equed, s_type_wr *s, s_type_wr *b, i_type_wr *
	ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, s_type_wr 
	*work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int spptrf_(char *uplo, i_type_wr *n, s_type_wr *ap, i_type_wr *info);

/* Subroutine */ int spptri_(char *uplo, i_type_wr *n, s_type_wr *ap, i_type_wr *info);

/* Subroutine */ int spptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *ap, 
	s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int spstf2_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *piv, i_type_wr *rank, s_type_wr *tol, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int spstrf_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *piv, i_type_wr *rank, s_type_wr *tol, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sptcon_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, s_type_wr *anorm, 
	s_type_wr *rcond, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int spteqr_(char *compz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sptrfs_(i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *df, s_type_wr *ef, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, 
	s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sptsv_(i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sptsvx_(char *fact, i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, 
	 s_type_wr *e, s_type_wr *df, s_type_wr *ef, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr 
	*ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int spttrf_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, i_type_wr *info);

/* Subroutine */ int spttrs_(i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sptts2_(i_type_wr *n, i_type_wr *nrhs, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int srscl_(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx);

/* Subroutine */ int ssbev_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	s_type_wr *ab, i_type_wr *ldab, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, 
	 i_type_wr *info);

/* Subroutine */ int ssbevd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	s_type_wr *ab, i_type_wr *ldab, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, 
	 i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int ssbevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *kd, s_type_wr *ab, i_type_wr *ldab, s_type_wr *q, i_type_wr *ldq, s_type_wr *vl, 
	 s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, i_type_wr *m, s_type_wr *
	w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *iwork, i_type_wr *
	ifail, i_type_wr *info);

/* Subroutine */ int ssbgst_(char *vect, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, s_type_wr *ab, i_type_wr *ldab, s_type_wr *bb, i_type_wr *ldbb, s_type_wr *
	x, i_type_wr *ldx, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int ssbgv_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, s_type_wr *ab, i_type_wr *ldab, s_type_wr *bb, i_type_wr *ldbb, s_type_wr *
	w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int ssbgvd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, s_type_wr *ab, i_type_wr *ldab, s_type_wr *bb, i_type_wr *ldbb, s_type_wr *
	w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *lwork, i_type_wr *
	iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int ssbgvx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *ka, i_type_wr *kb, s_type_wr *ab, i_type_wr *ldab, s_type_wr *bb, i_type_wr *
	ldbb, s_type_wr *q, i_type_wr *ldq, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr 
	*iu, s_type_wr *abstol, i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr 
	*work, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int ssbtrd_(char *vect, char *uplo, i_type_wr *n, i_type_wr *kd, 
	s_type_wr *ab, i_type_wr *ldab, s_type_wr *d__, s_type_wr *e, s_type_wr *q, i_type_wr *ldq, 
	s_type_wr *work, i_type_wr *info);

/* Subroutine */ int ssfrk_(char *transr, char *uplo, char *trans, i_type_wr *n, 
	 i_type_wr *k, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *beta, s_type_wr *
	c__);

/* Subroutine */ int sspcon_(char *uplo, i_type_wr *n, s_type_wr *ap, i_type_wr *ipiv, 
	s_type_wr *anorm, s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sspev_(char *jobz, char *uplo, i_type_wr *n, s_type_wr *ap, 
	s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sspevd_(char *jobz, char *uplo, i_type_wr *n, s_type_wr *ap, 
	s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *lwork, i_type_wr 
	*iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int sspevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	s_type_wr *ap, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, 
	i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *
	iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int sspgst_(i_type_wr *itype, char *uplo, i_type_wr *n, s_type_wr *ap, 
	 s_type_wr *bp, i_type_wr *info);

/* Subroutine */ int sspgv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, s_type_wr *ap, s_type_wr *bp, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int sspgvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, s_type_wr *ap, s_type_wr *bp, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int sspgvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *bp, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, 
	 i_type_wr *iu, s_type_wr *abstol, i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *
	ldz, s_type_wr *work, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int ssprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *ap, 
	s_type_wr *afp, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *
	ldx, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int sspsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *ap, 
	i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sspsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *ap, s_type_wr *afp, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, s_type_wr 
	*x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int ssptrd_(char *uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *tau, i_type_wr *info);

/* Subroutine */ int ssptrf_(char *uplo, i_type_wr *n, s_type_wr *ap, i_type_wr *ipiv, 
	i_type_wr *info);

/* Subroutine */ int ssptri_(char *uplo, i_type_wr *n, s_type_wr *ap, i_type_wr *ipiv, 
	s_type_wr *work, i_type_wr *info);

/* Subroutine */ int ssptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *ap, 
	i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int sstebz_(char *range, char *order, i_type_wr *n, s_type_wr *vl, 
	s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, s_type_wr *d__, s_type_wr *e, 
	i_type_wr *m, i_type_wr *nsplit, s_type_wr *w, i_type_wr *iblock, i_type_wr *
	isplit, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int sstedc_(char *compz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int sstegr_(char *jobz, char *range, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, 
	i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, s_type_wr *
	work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int sstein_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, i_type_wr *m, s_type_wr 
	*w, i_type_wr *iblock, i_type_wr *isplit, s_type_wr *z__, i_type_wr *ldz, s_type_wr *
	work, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int sstemr_(char *jobz, char *range, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, i_type_wr *m, 
	s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, i_type_wr *nzc, i_type_wr *isuppz, 
	l_type_wr *tryrac, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *
	liwork, i_type_wr *info);

/* Subroutine */ int ssteqr_(char *compz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, 
	s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int ssterf_(i_type_wr *n, s_type_wr *d__, s_type_wr *e, i_type_wr *info);

/* Subroutine */ int sstev_(char *jobz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, s_type_wr *
	z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int sstevd_(char *jobz, i_type_wr *n, s_type_wr *d__, s_type_wr *e, s_type_wr 
	*z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int sstevr_(char *jobz, char *range, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, 
	i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, s_type_wr *
	work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int sstevx_(char *jobz, char *range, i_type_wr *n, s_type_wr *d__, 
	s_type_wr *e, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, 
	i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *
	iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int ssycon_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, s_type_wr *anorm, s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int ssyequb_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *s, s_type_wr *scond, s_type_wr *amax, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int ssyev_(char *jobz, char *uplo, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *w, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int ssyevd_(char *jobz, char *uplo, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *w, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int ssyevr_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	s_type_wr *a, i_type_wr *lda, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, 
	s_type_wr *abstol, i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, i_type_wr *
	isuppz, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, 
	i_type_wr *info);

/* Subroutine */ int ssyevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	s_type_wr *a, i_type_wr *lda, s_type_wr *vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, 
	s_type_wr *abstol, i_type_wr *m, s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *
	work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int ssygs2_(i_type_wr *itype, char *uplo, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int ssygst_(i_type_wr *itype, char *uplo, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int ssygv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *w, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int ssygvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *w, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int ssygvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *
	vl, s_type_wr *vu, i_type_wr *il, i_type_wr *iu, s_type_wr *abstol, i_type_wr *m, 
	s_type_wr *w, s_type_wr *z__, i_type_wr *ldz, s_type_wr *work, i_type_wr *lwork, i_type_wr 
	*iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int ssyrfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, s_type_wr *b, 
	i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, s_type_wr *
	work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int ssyrfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, 
	s_type_wr *s, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, 
	s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *err_bnds_norm__, s_type_wr *
	err_bnds_comp__, i_type_wr *nparams, s_type_wr *params, s_type_wr *work, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int ssysv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int ssysvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, 
	s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *rcond, s_type_wr *ferr, 
	 s_type_wr *berr, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int ssysvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, 
	char *equed, s_type_wr *s, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, 
	s_type_wr *rcond, s_type_wr *rpvgrw, s_type_wr *berr, i_type_wr *n_err_bnds__, s_type_wr *
	err_bnds_norm__, s_type_wr *err_bnds_comp__, i_type_wr *nparams, s_type_wr *
	params, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int ssytd2_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *d__, s_type_wr *e, s_type_wr *tau, i_type_wr *info);

/* Subroutine */ int ssytf2_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int ssytrd_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *d__, s_type_wr *e, s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int ssytrf_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int ssytri_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int ssytrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, s_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int stbcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	i_type_wr *kd, s_type_wr *ab, i_type_wr *ldab, s_type_wr *rcond, s_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int stbrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *b, i_type_wr 
	*ldb, s_type_wr *x, i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int stbtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, s_type_wr *ab, i_type_wr *ldab, s_type_wr *b, i_type_wr 
	*ldb, i_type_wr *info);

/* Subroutine */ int stfsm_(char *transr, char *side, char *uplo, char *trans, 
	 char *diag, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, s_type_wr *b, 
	i_type_wr *ldb);

/* Subroutine */ int stftri_(char *transr, char *uplo, char *diag, i_type_wr *n, 
	 s_type_wr *a, i_type_wr *info);

/* Subroutine */ int stfttp_(char *transr, char *uplo, i_type_wr *n, s_type_wr *arf, 
	s_type_wr *ap, i_type_wr *info);

/* Subroutine */ int stfttr_(char *transr, char *uplo, i_type_wr *n, s_type_wr *arf, 
	s_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int stgevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, s_type_wr *s, i_type_wr *lds, s_type_wr *p, i_type_wr *ldp, s_type_wr *vl, 
	i_type_wr *ldvl, s_type_wr *vr, i_type_wr *ldvr, i_type_wr *mm, i_type_wr *m, s_type_wr 
	*work, i_type_wr *info);

/* Subroutine */ int stgex2_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, s_type_wr 
	*a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *q, i_type_wr *ldq, s_type_wr *
	z__, i_type_wr *ldz, i_type_wr *j1, i_type_wr *n1, i_type_wr *n2, s_type_wr *work, 
	i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int stgexc_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, s_type_wr 
	*a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *q, i_type_wr *ldq, s_type_wr *
	z__, i_type_wr *ldz, i_type_wr *ifst, i_type_wr *ilst, s_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int stgsen_(i_type_wr *ijob, l_type_wr *wantq, l_type_wr *wantz, 
	l_type_wr *select, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *
	ldb, s_type_wr *alphar, s_type_wr *alphai, s_type_wr *beta, s_type_wr *q, i_type_wr *ldq, 
	s_type_wr *z__, i_type_wr *ldz, i_type_wr *m, s_type_wr *pl, s_type_wr *pr, s_type_wr *dif, 
	s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *
	info);

/* Subroutine */ int stgsja_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, i_type_wr *k, i_type_wr *l, s_type_wr *a, i_type_wr *lda, 
	 s_type_wr *b, i_type_wr *ldb, s_type_wr *tola, s_type_wr *tolb, s_type_wr *alpha, s_type_wr *
	beta, s_type_wr *u, i_type_wr *ldu, s_type_wr *v, i_type_wr *ldv, s_type_wr *q, i_type_wr *
	ldq, s_type_wr *work, i_type_wr *ncycle, i_type_wr *info);

/* Subroutine */ int stgsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *vl, 
	i_type_wr *ldvl, s_type_wr *vr, i_type_wr *ldvr, s_type_wr *s, s_type_wr *dif, i_type_wr *
	mm, i_type_wr *m, s_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *
	info);

/* Subroutine */ int stgsy2_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *c__, i_type_wr *
	ldc, s_type_wr *d__, i_type_wr *ldd, s_type_wr *e, i_type_wr *lde, s_type_wr *f, i_type_wr 
	*ldf, s_type_wr *scale, s_type_wr *rdsum, s_type_wr *rdscal, i_type_wr *iwork, i_type_wr 
	*pq, i_type_wr *info);

/* Subroutine */ int stgsyl_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *c__, i_type_wr *
	ldc, s_type_wr *d__, i_type_wr *ldd, s_type_wr *e, i_type_wr *lde, s_type_wr *f, i_type_wr 
	*ldf, s_type_wr *scale, s_type_wr *dif, s_type_wr *work, i_type_wr *lwork, i_type_wr *
	iwork, i_type_wr *info);

/* Subroutine */ int stpcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	s_type_wr *ap, s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int stprfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, s_type_wr *ap, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, i_type_wr *ldx, 
	 s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int stptri_(char *uplo, char *diag, i_type_wr *n, s_type_wr *ap, 
	i_type_wr *info);

/* Subroutine */ int stptrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, s_type_wr *ap, s_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int stpttf_(char *transr, char *uplo, i_type_wr *n, s_type_wr *ap, 
	s_type_wr *arf, i_type_wr *info);

/* Subroutine */ int stpttr_(char *uplo, i_type_wr *n, s_type_wr *ap, s_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int strcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	s_type_wr *a, i_type_wr *lda, s_type_wr *rcond, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int strevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, s_type_wr *t, i_type_wr *ldt, s_type_wr *vl, i_type_wr *ldvl, s_type_wr *vr, 
	i_type_wr *ldvr, i_type_wr *mm, i_type_wr *m, s_type_wr *work, i_type_wr *info);

/* Subroutine */ int strexc_(char *compq, i_type_wr *n, s_type_wr *t, i_type_wr *ldt, 
	s_type_wr *q, i_type_wr *ldq, i_type_wr *ifst, i_type_wr *ilst, s_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int strrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *x, 
	i_type_wr *ldx, s_type_wr *ferr, s_type_wr *berr, s_type_wr *work, i_type_wr *iwork, 
	i_type_wr *info);

/* Subroutine */ int strsen_(char *job, char *compq, l_type_wr *select, i_type_wr 
	*n, s_type_wr *t, i_type_wr *ldt, s_type_wr *q, i_type_wr *ldq, s_type_wr *wr, s_type_wr *wi, 
	i_type_wr *m, s_type_wr *s, s_type_wr *sep, s_type_wr *work, i_type_wr *lwork, i_type_wr *
	iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int strsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, s_type_wr *t, i_type_wr *ldt, s_type_wr *vl, i_type_wr *ldvl, s_type_wr *vr, 
	i_type_wr *ldvr, s_type_wr *s, s_type_wr *sep, i_type_wr *mm, i_type_wr *m, s_type_wr *
	work, i_type_wr *ldwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int strsyl_(char *trana, char *tranb, i_type_wr *isgn, i_type_wr 
	*m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, s_type_wr *
	c__, i_type_wr *ldc, s_type_wr *scale, i_type_wr *info);

/* Subroutine */ int strti2_(char *uplo, char *diag, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int strtri_(char *uplo, char *diag, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int strtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int strttf_(char *transr, char *uplo, i_type_wr *n, s_type_wr *a, 
	i_type_wr *lda, s_type_wr *arf, i_type_wr *info);

/* Subroutine */ int strttp_(char *uplo, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *ap, i_type_wr *info);

/* Subroutine */ int stzrqf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, i_type_wr *info);

/* Subroutine */ int stzrzf_(i_type_wr *m, i_type_wr *n, s_type_wr *a, i_type_wr *lda, 
	s_type_wr *tau, s_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int xerbla_array__(char *srname_array__, i_type_wr *
	srname_len__, i_type_wr *info, ftn_len_wr srname_array_len);

/* Subroutine */ int zbdsqr_(char *uplo, i_type_wr *n, i_type_wr *ncvt, i_type_wr *
	nru, i_type_wr *ncc, d_type_wr *d__, d_type_wr *e, z_type_wr *vt, 
	i_type_wr *ldvt, z_type_wr *u, i_type_wr *ldu, z_type_wr *c__, 
	i_type_wr *ldc, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zcgesv_(i_type_wr *n, i_type_wr *nrhs, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, z_type_wr *work, c_type_wr *swork, 
	d_type_wr *rwork, i_type_wr *iter, i_type_wr *info);

/* Subroutine */ int zcposv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, z_type_wr *work, c_type_wr *swork, 
	d_type_wr *rwork, i_type_wr *iter, i_type_wr *info);

/* Subroutine */ int zdrscl_(i_type_wr *n, d_type_wr *sa, z_type_wr *sx, 
	i_type_wr *incx);

/* Subroutine */ int zgbbrd_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *ncc, 
	 i_type_wr *kl, i_type_wr *ku, z_type_wr *ab, i_type_wr *ldab, 
	d_type_wr *d__, d_type_wr *e, z_type_wr *q, i_type_wr *ldq, 
	z_type_wr *pt, i_type_wr *ldpt, z_type_wr *c__, i_type_wr *ldc, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgbcon_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 z_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, d_type_wr *anorm, 
	d_type_wr *rcond, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zgbequ_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 z_type_wr *ab, i_type_wr *ldab, d_type_wr *r__, d_type_wr *c__, 
	d_type_wr *rowcnd, d_type_wr *colcnd, d_type_wr *amax, i_type_wr *
	info);

/* Subroutine */ int zgbequb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, z_type_wr *ab, i_type_wr *ldab, d_type_wr *r__, d_type_wr *
	c__, d_type_wr *rowcnd, d_type_wr *colcnd, d_type_wr *amax, 
	i_type_wr *info);

/* Subroutine */ int zgbrfs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, z_type_wr *
	afb, i_type_wr *ldafb, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgbrfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, 
	z_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, d_type_wr *r__, 
	d_type_wr *c__, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, 
	i_type_wr *ldx, d_type_wr *rcond, d_type_wr *berr, i_type_wr *
	n_err_bnds__, d_type_wr *err_bnds_norm__, d_type_wr *
	err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgbsv_(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	nrhs, z_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, z_type_wr *
	b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zgbsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *kl, 
	 i_type_wr *ku, i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, 
	z_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, char *equed, 
	d_type_wr *r__, d_type_wr *c__, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, 
	d_type_wr *berr, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zgbsvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	kl, i_type_wr *ku, i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, 
	z_type_wr *afb, i_type_wr *ldafb, i_type_wr *ipiv, char *equed, 
	d_type_wr *r__, d_type_wr *c__, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *rpvgrw, 
	 d_type_wr *berr, i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, 
	 d_type_wr *err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgbtf2_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 z_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zgbtrf_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 z_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zgbtrs_(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *
	ku, i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, i_type_wr *ipiv, 
	z_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zgebak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, d_type_wr *scale, i_type_wr *m, z_type_wr *v, 
	i_type_wr *ldv, i_type_wr *info);

/* Subroutine */ int zgebal_(char *job, i_type_wr *n, z_type_wr *a, i_type_wr 
	*lda, i_type_wr *ilo, i_type_wr *ihi, d_type_wr *scale, i_type_wr *info);

/* Subroutine */ int zgebd2_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *d__, d_type_wr *e, z_type_wr *tauq, 
	z_type_wr *taup, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zgebrd_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *d__, d_type_wr *e, z_type_wr *tauq, 
	z_type_wr *taup, z_type_wr *work, i_type_wr *lwork, i_type_wr *
	info);

/* Subroutine */ int zgecon_(char *norm, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *anorm, d_type_wr *rcond, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgeequ_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *r__, d_type_wr *c__, d_type_wr *rowcnd, 
	d_type_wr *colcnd, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int zgeequb_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *r__, d_type_wr *c__, d_type_wr *rowcnd, 
	d_type_wr *colcnd, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int zgees_(char *jobvs, char *sort, sel_fun_wr select, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *sdim, z_type_wr *w, 
	z_type_wr *vs, i_type_wr *ldvs, z_type_wr *work, i_type_wr *lwork, 
	 d_type_wr *rwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int zgeesx_(char *jobvs, char *sort, sel_fun_wr select, char *
	sense, i_type_wr *n, z_type_wr *a, i_type_wr *lda, i_type_wr *sdim, 
	z_type_wr *w, z_type_wr *vs, i_type_wr *ldvs, d_type_wr *
	rconde, d_type_wr *rcondv, z_type_wr *work, i_type_wr *lwork, 
	d_type_wr *rwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int zgeev_(char *jobvl, char *jobvr, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *w, z_type_wr *vl, 
	i_type_wr *ldvl, z_type_wr *vr, i_type_wr *ldvr, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgeevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *w, 
	z_type_wr *vl, i_type_wr *ldvl, z_type_wr *vr, i_type_wr *ldvr, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *scale, d_type_wr *abnrm, 
	d_type_wr *rconde, d_type_wr *rcondv, z_type_wr *work, i_type_wr *
	lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgegs_(char *jobvsl, char *jobvsr, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *alpha, z_type_wr *beta, z_type_wr *vsl, 
	i_type_wr *ldvsl, z_type_wr *vsr, i_type_wr *ldvsr, z_type_wr *
	work, i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgegv_(char *jobvl, char *jobvr, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *alpha, z_type_wr *beta, z_type_wr *vl, i_type_wr 
	*ldvl, z_type_wr *vr, i_type_wr *ldvr, z_type_wr *work, i_type_wr 
	*lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgehd2_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zgehrd_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zgelq2_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zgelqf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zgels_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zgelsd_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	d_type_wr *s, d_type_wr *rcond, i_type_wr *rank, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int zgelss_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	d_type_wr *s, d_type_wr *rcond, i_type_wr *rank, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgelsx_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *jpvt, d_type_wr *rcond, i_type_wr *rank, z_type_wr *work, 
	d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgelsy_(i_type_wr *m, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *jpvt, d_type_wr *rcond, i_type_wr *rank, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgeql2_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zgeqlf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zgeqp3_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *jpvt, z_type_wr *tau, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgeqpf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *jpvt, z_type_wr *tau, z_type_wr *work, 
	d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgeqr2_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zgeqrf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zgerfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, 
	i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, z_type_wr *work, 
	 d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgerfsx_(char *trans, char *equed, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, d_type_wr *r__, d_type_wr *c__, z_type_wr *
	b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, 
	d_type_wr *berr, i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, 
	d_type_wr *err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgerq2_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zgerqf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zgesc2_(i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	z_type_wr *rhs, i_type_wr *ipiv, i_type_wr *jpiv, d_type_wr *scale);

/* Subroutine */ int zgesdd_(char *jobz, i_type_wr *m, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *s, z_type_wr *u, 
	i_type_wr *ldu, z_type_wr *vt, i_type_wr *ldvt, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int zgesv_(i_type_wr *n, i_type_wr *nrhs, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, i_type_wr *
	info);

/* Subroutine */ int zgesvd_(char *jobu, char *jobvt, i_type_wr *m, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *s, z_type_wr *u, 
	i_type_wr *ldu, z_type_wr *vt, i_type_wr *ldvt, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgesvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, char *equed, d_type_wr *r__, d_type_wr *c__, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgesvxx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, char *equed, d_type_wr *r__, d_type_wr *c__, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *rcond, d_type_wr *rpvgrw, d_type_wr *berr, i_type_wr *
	n_err_bnds__, d_type_wr *err_bnds_norm__, d_type_wr *
	err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgetc2_(i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, i_type_wr *jpiv, i_type_wr *info);

/* Subroutine */ int zgetf2_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zgetrf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zgetri_(i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	i_type_wr *ipiv, z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zgetrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, z_type_wr *b, 
	i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zggbak_(char *job, char *side, i_type_wr *n, i_type_wr *ilo, 
	i_type_wr *ihi, d_type_wr *lscale, d_type_wr *rscale, i_type_wr *m, 
	z_type_wr *v, i_type_wr *ldv, i_type_wr *info);

/* Subroutine */ int zggbal_(char *job, i_type_wr *n, z_type_wr *a, i_type_wr 
	*lda, z_type_wr *b, i_type_wr *ldb, i_type_wr *ilo, i_type_wr *ihi, 
	d_type_wr *lscale, d_type_wr *rscale, d_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int zgges_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
	i_type_wr *ldb, i_type_wr *sdim, z_type_wr *alpha, z_type_wr *
	beta, z_type_wr *vsl, i_type_wr *ldvsl, z_type_wr *vsr, i_type_wr 
	*ldvsr, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, 
	l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int zggesx_(char *jobvsl, char *jobvsr, char *sort, sel_fun_wr 
	selctg, char *sense, i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	z_type_wr *b, i_type_wr *ldb, i_type_wr *sdim, z_type_wr *alpha, 
	z_type_wr *beta, z_type_wr *vsl, i_type_wr *ldvsl, 
	z_type_wr *vsr, i_type_wr *ldvsr, d_type_wr *rconde, d_type_wr *
	rcondv, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, 
	i_type_wr *iwork, i_type_wr *liwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int zggev_(char *jobvl, char *jobvr, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *alpha, z_type_wr *beta, z_type_wr *vl, i_type_wr 
	*ldvl, z_type_wr *vr, i_type_wr *ldvr, z_type_wr *work, i_type_wr 
	*lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zggevx_(char *balanc, char *jobvl, char *jobvr, char *
	sense, i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *alpha, z_type_wr *beta, 
	z_type_wr *vl, i_type_wr *ldvl, z_type_wr *vr, i_type_wr *ldvr, 
	i_type_wr *ilo, i_type_wr *ihi, d_type_wr *lscale, d_type_wr *rscale, 
	d_type_wr *abnrm, d_type_wr *bbnrm, d_type_wr *rconde, d_type_wr *
	rcondv, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, 
	i_type_wr *iwork, l_type_wr *bwork, i_type_wr *info);

/* Subroutine */ int zggglm_(i_type_wr *n, i_type_wr *m, i_type_wr *p, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *d__, z_type_wr *x, z_type_wr *y, z_type_wr 
	*work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zgghrd_(char *compq, char *compz, i_type_wr *n, i_type_wr *
	ilo, i_type_wr *ihi, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *q, i_type_wr *ldq, z_type_wr *z__, 
	i_type_wr *ldz, i_type_wr *info);

/* Subroutine */ int zgglse_(i_type_wr *m, i_type_wr *n, i_type_wr *p, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *c__, z_type_wr *d__, z_type_wr *x, 
	z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zggqrf_(i_type_wr *n, i_type_wr *m, i_type_wr *p, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *taua, z_type_wr *b, 
	 i_type_wr *ldb, z_type_wr *taub, z_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int zggrqf_(i_type_wr *m, i_type_wr *p, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *taua, z_type_wr *b, 
	 i_type_wr *ldb, z_type_wr *taub, z_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int zggsvd_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *n, i_type_wr *p, i_type_wr *k, i_type_wr *l, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, d_type_wr *alpha, 
	d_type_wr *beta, z_type_wr *u, i_type_wr *ldu, z_type_wr *v, 
	i_type_wr *ldv, z_type_wr *q, i_type_wr *ldq, z_type_wr *work, 
	d_type_wr *rwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int zggsvp_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr 
	*b, i_type_wr *ldb, d_type_wr *tola, d_type_wr *tolb, i_type_wr *k, 
	i_type_wr *l, z_type_wr *u, i_type_wr *ldu, z_type_wr *v, i_type_wr 
	*ldv, z_type_wr *q, i_type_wr *ldq, i_type_wr *iwork, d_type_wr *
	rwork, z_type_wr *tau, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zgtcon_(char *norm, i_type_wr *n, z_type_wr *dl, 
	z_type_wr *d__, z_type_wr *du, z_type_wr *du2, i_type_wr *
	ipiv, d_type_wr *anorm, d_type_wr *rcond, z_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int zgtrfs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *dl, z_type_wr *d__, z_type_wr *du, 
	z_type_wr *dlf, z_type_wr *df, z_type_wr *duf, 
	z_type_wr *du2, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zgtsv_(i_type_wr *n, i_type_wr *nrhs, z_type_wr *dl, 
	z_type_wr *d__, z_type_wr *du, z_type_wr *b, i_type_wr *ldb, 
	 i_type_wr *info);

/* Subroutine */ int zgtsvx_(char *fact, char *trans, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *dl, z_type_wr *d__, z_type_wr *du, 
	z_type_wr *dlf, z_type_wr *df, z_type_wr *duf, 
	z_type_wr *du2, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, 
	d_type_wr *berr, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zgttrf_(i_type_wr *n, z_type_wr *dl, z_type_wr *
	d__, z_type_wr *du, z_type_wr *du2, i_type_wr *ipiv, i_type_wr *
	info);

/* Subroutine */ int zgttrs_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *dl, z_type_wr *d__, z_type_wr *du, 
	z_type_wr *du2, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zgtts2_(i_type_wr *itrans, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *dl, z_type_wr *d__, z_type_wr *du, 
	z_type_wr *du2, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int zhbev_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *w, z_type_wr *z__, 
	i_type_wr *ldz, z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zhbevd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *w, z_type_wr *z__, 
	i_type_wr *ldz, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, 
	i_type_wr *lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zhbevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *kd, z_type_wr *ab, i_type_wr *ldab, z_type_wr *q, 
	i_type_wr *ldq, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, i_type_wr *
	iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, z_type_wr *z__, 
	 i_type_wr *ldz, z_type_wr *work, d_type_wr *rwork, i_type_wr *iwork, 
	 i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int zhbgst_(char *vect, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, z_type_wr *ab, i_type_wr *ldab, z_type_wr *bb, 
	i_type_wr *ldbb, z_type_wr *x, i_type_wr *ldx, z_type_wr *work, 
	d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zhbgv_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, z_type_wr *ab, i_type_wr *ldab, z_type_wr *bb, 
	i_type_wr *ldbb, d_type_wr *w, z_type_wr *z__, i_type_wr *ldz, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zhbgvd_(char *jobz, char *uplo, i_type_wr *n, i_type_wr *ka, 
	i_type_wr *kb, z_type_wr *ab, i_type_wr *ldab, z_type_wr *bb, 
	i_type_wr *ldbb, d_type_wr *w, z_type_wr *z__, i_type_wr *ldz, 
	z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, i_type_wr *
	lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zhbgvx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	i_type_wr *ka, i_type_wr *kb, z_type_wr *ab, i_type_wr *ldab, 
	z_type_wr *bb, i_type_wr *ldbb, z_type_wr *q, i_type_wr *ldq, 
	d_type_wr *vl, d_type_wr *vu, i_type_wr *il, i_type_wr *iu, d_type_wr *
	abstol, i_type_wr *m, d_type_wr *w, z_type_wr *z__, i_type_wr *ldz, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *iwork, i_type_wr *
	ifail, i_type_wr *info);

/* Subroutine */ int zhbtrd_(char *vect, char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *d__, d_type_wr *e, 
	z_type_wr *q, i_type_wr *ldq, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zhecon_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, d_type_wr *anorm, d_type_wr *rcond, 
	z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zheequb_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, 
	z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zheev_(char *jobz, char *uplo, i_type_wr *n, z_type_wr 
	*a, i_type_wr *lda, d_type_wr *w, z_type_wr *work, i_type_wr *lwork, 
	d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zheevd_(char *jobz, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *w, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *lrwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zheevr_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *vl, d_type_wr *vu, 
	i_type_wr *il, i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *
	w, z_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, z_type_wr *
	work, i_type_wr *lwork, d_type_wr *rwork, i_type_wr *lrwork, i_type_wr *
	iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zheevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *vl, d_type_wr *vu, 
	i_type_wr *il, i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *
	w, z_type_wr *z__, i_type_wr *ldz, z_type_wr *work, i_type_wr *
	lwork, d_type_wr *rwork, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *
	info);

/* Subroutine */ int zhegs2_(i_type_wr *itype, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zhegst_(i_type_wr *itype, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zhegv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	d_type_wr *w, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, 
	 i_type_wr *info);

/* Subroutine */ int zhegvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	d_type_wr *w, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, 
	 i_type_wr *lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zhegvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
	i_type_wr *ldb, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, i_type_wr *
	iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, z_type_wr *z__, 
	 i_type_wr *ldz, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, 
	 i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int zherfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, 
	i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, z_type_wr *work, 
	 d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zherfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, d_type_wr *s, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *berr, 
	i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, d_type_wr *
	err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zhesv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zhesvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, 
	 i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, 
	z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zhesvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, char *equed, d_type_wr *s, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, 
	d_type_wr *rpvgrw, d_type_wr *berr, i_type_wr *n_err_bnds__, 
	d_type_wr *err_bnds_norm__, d_type_wr *err_bnds_comp__, i_type_wr *
	nparams, d_type_wr *params, z_type_wr *work, d_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int zhetd2_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *d__, d_type_wr *e, z_type_wr *tau, 
	i_type_wr *info);

/* Subroutine */ int zhetf2_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zhetrd_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *d__, d_type_wr *e, z_type_wr *tau, 
	z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zhetrf_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, z_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int zhetri_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zhetrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, z_type_wr *b, 
	i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zhfrk_(char *transr, char *uplo, char *trans, i_type_wr *n, 
	 i_type_wr *k, d_type_wr *alpha, z_type_wr *a, i_type_wr *lda, 
	d_type_wr *beta, z_type_wr *c__);

/* Subroutine */ int zhgeqz_(char *job, char *compq, char *compz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, z_type_wr *h__, i_type_wr *ldh, 
	z_type_wr *t, i_type_wr *ldt, z_type_wr *alpha, z_type_wr *
	beta, z_type_wr *q, i_type_wr *ldq, z_type_wr *z__, i_type_wr *
	ldz, z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zhpcon_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *ipiv, d_type_wr *anorm, d_type_wr *rcond, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zhpev_(char *jobz, char *uplo, i_type_wr *n, z_type_wr 
	*ap, d_type_wr *w, z_type_wr *z__, i_type_wr *ldz, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zhpevd_(char *jobz, char *uplo, i_type_wr *n, 
	z_type_wr *ap, d_type_wr *w, z_type_wr *z__, i_type_wr *ldz, 
	z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, i_type_wr *
	lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zhpevx_(char *jobz, char *range, char *uplo, i_type_wr *n, 
	z_type_wr *ap, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, 
	i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, 
	z_type_wr *z__, i_type_wr *ldz, z_type_wr *work, d_type_wr *
	rwork, i_type_wr *iwork, i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int zhpgst_(i_type_wr *itype, char *uplo, i_type_wr *n, 
	z_type_wr *ap, z_type_wr *bp, i_type_wr *info);

/* Subroutine */ int zhpgv_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, z_type_wr *ap, z_type_wr *bp, d_type_wr *w, z_type_wr 
	*z__, i_type_wr *ldz, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zhpgvd_(i_type_wr *itype, char *jobz, char *uplo, i_type_wr *
	n, z_type_wr *ap, z_type_wr *bp, d_type_wr *w, z_type_wr 
	*z__, i_type_wr *ldz, z_type_wr *work, i_type_wr *lwork, d_type_wr *
	rwork, i_type_wr *lrwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *
	info);

/* Subroutine */ int zhpgvx_(i_type_wr *itype, char *jobz, char *range, char *
	uplo, i_type_wr *n, z_type_wr *ap, z_type_wr *bp, d_type_wr *
	vl, d_type_wr *vu, i_type_wr *il, i_type_wr *iu, d_type_wr *abstol, 
	i_type_wr *m, d_type_wr *w, z_type_wr *z__, i_type_wr *ldz, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *iwork, i_type_wr *
	ifail, i_type_wr *info);

/* Subroutine */ int zhprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, z_type_wr *afp, i_type_wr *ipiv, z_type_wr *
	b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, 
	d_type_wr *berr, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zhpsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zhpsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *ap, z_type_wr *afp, i_type_wr *ipiv, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zhptrd_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	d_type_wr *d__, d_type_wr *e, z_type_wr *tau, i_type_wr *info);

/* Subroutine */ int zhptrf_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zhptri_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *ipiv, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zhptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zhsein_(char *side, char *eigsrc, char *initv, l_type_wr *
	select, i_type_wr *n, z_type_wr *h__, i_type_wr *ldh, z_type_wr *
	w, z_type_wr *vl, i_type_wr *ldvl, z_type_wr *vr, i_type_wr *ldvr, 
	 i_type_wr *mm, i_type_wr *m, z_type_wr *work, d_type_wr *rwork, 
	i_type_wr *ifaill, i_type_wr *ifailr, i_type_wr *info);

/* Subroutine */ int zhseqr_(char *job, char *compz, i_type_wr *n, i_type_wr *ilo, 
	 i_type_wr *ihi, z_type_wr *h__, i_type_wr *ldh, z_type_wr *w, 
	z_type_wr *z__, i_type_wr *ldz, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zla_gbamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *kl, i_type_wr *ku, d_type_wr *alpha, z_type_wr *ab, 
	i_type_wr *ldab, z_type_wr *x, i_type_wr *incx, d_type_wr *beta, 
	d_type_wr *y, i_type_wr *incy);

d_type_wr zla_gbrcond_c__(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	z_type_wr *ab, i_type_wr *ldab, z_type_wr *afb, i_type_wr *ldafb, 
	i_type_wr *ipiv, d_type_wr *c__, l_type_wr *capply, i_type_wr *info, 
	z_type_wr *work, d_type_wr *rwork, ftn_len_wr trans_len);

d_type_wr zla_gbrcond_x__(char *trans, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	z_type_wr *ab, i_type_wr *ldab, z_type_wr *afb, i_type_wr *ldafb, 
	i_type_wr *ipiv, z_type_wr *x, i_type_wr *info, z_type_wr *work, 
	d_type_wr *rwork, ftn_len_wr trans_len);

/* Subroutine */ int zla_gbrfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *nrhs, 
	z_type_wr *ab, i_type_wr *ldab, z_type_wr *afb, i_type_wr *ldafb, 
	i_type_wr *ipiv, l_type_wr *colequ, d_type_wr *c__, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *y, i_type_wr *ldy, d_type_wr *berr_out__, 
	i_type_wr *n_norms__, d_type_wr *errs_n__, d_type_wr *errs_c__, 
	z_type_wr *res, d_type_wr *ayb, z_type_wr *dy, z_type_wr 
	*y_tail__, d_type_wr *rcond, i_type_wr *ithresh, d_type_wr *rthresh, 
	d_type_wr *dz_ub__, l_type_wr *ignore_cwise__, i_type_wr *info);

d_type_wr zla_gbrpvgrw__(i_type_wr *n, i_type_wr *kl, i_type_wr *ku, i_type_wr *
	ncols, z_type_wr *ab, i_type_wr *ldab, z_type_wr *afb, i_type_wr *
	ldafb);

/* Subroutine */ int zla_geamv__(i_type_wr *trans, i_type_wr *m, i_type_wr *n, 
	d_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *x, 
	i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy);

d_type_wr zla_gercond_c__(char *trans, i_type_wr *n, z_type_wr *a, i_type_wr 
	*lda, z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, d_type_wr *
	c__, l_type_wr *capply, i_type_wr *info, z_type_wr *work, d_type_wr *
	rwork, ftn_len_wr trans_len);

d_type_wr zla_gercond_x__(char *trans, i_type_wr *n, z_type_wr *a, i_type_wr 
	*lda, z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, z_type_wr *
	x, i_type_wr *info, z_type_wr *work, d_type_wr *rwork, ftn_len_wr 
	trans_len);

/* Subroutine */ int zla_gerfsx_extended__(i_type_wr *prec_type__, i_type_wr *
	trans_type__, i_type_wr *n, i_type_wr *nrhs, z_type_wr *a, i_type_wr *
	lda, z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ,
	 d_type_wr *c__, z_type_wr *b, i_type_wr *ldb, z_type_wr *y, 
	i_type_wr *ldy, d_type_wr *berr_out__, i_type_wr *n_norms__, d_type_wr *
	errs_n__, d_type_wr *errs_c__, z_type_wr *res, d_type_wr *ayb, 
	z_type_wr *dy, z_type_wr *y_tail__, d_type_wr *rcond, 
	i_type_wr *ithresh, d_type_wr *rthresh, d_type_wr *dz_ub__, l_type_wr *
	ignore_cwise__, i_type_wr *info);

/* Subroutine */ int zla_heamv__(i_type_wr *uplo, i_type_wr *n, d_type_wr *alpha,
	 z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx, 
	d_type_wr *beta, d_type_wr *y, i_type_wr *incy);

d_type_wr zla_hercond_c__(char *uplo, i_type_wr *n, z_type_wr *a, i_type_wr *
	lda, z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, d_type_wr *c__,
	 l_type_wr *capply, i_type_wr *info, z_type_wr *work, d_type_wr *
	rwork, ftn_len_wr uplo_len);

d_type_wr zla_hercond_x__(char *uplo, i_type_wr *n, z_type_wr *a, i_type_wr *
	lda, z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, z_type_wr *
	x, i_type_wr *info, z_type_wr *work, d_type_wr *rwork, ftn_len_wr 
	uplo_len);

/* Subroutine */ int zla_herfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, z_type_wr *a, i_type_wr *lda, 
	z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, 
	d_type_wr *c__, z_type_wr *b, i_type_wr *ldb, z_type_wr *y, 
	i_type_wr *ldy, d_type_wr *berr_out__, i_type_wr *n_norms__, d_type_wr *
	errs_n__, d_type_wr *errs_c__, z_type_wr *res, d_type_wr *ayb, 
	z_type_wr *dy, z_type_wr *y_tail__, d_type_wr *rcond, 
	i_type_wr *ithresh, d_type_wr *rthresh, d_type_wr *dz_ub__, l_type_wr *
	ignore_cwise__, i_type_wr *info, ftn_len_wr uplo_len);

d_type_wr zla_herpvgrw__(char *uplo, i_type_wr *n, i_type_wr *info, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, d_type_wr *work, ftn_len_wr uplo_len);

/* Subroutine */ int zla_lin_berr__(i_type_wr *n, i_type_wr *nz, i_type_wr *nrhs, 
	z_type_wr *res, d_type_wr *ayb, d_type_wr *berr);

d_type_wr zla_porcond_c__(char *uplo, i_type_wr *n, z_type_wr *a, i_type_wr *
	lda, z_type_wr *af, i_type_wr *ldaf, d_type_wr *c__, l_type_wr *
	capply, i_type_wr *info, z_type_wr *work, d_type_wr *rwork, ftn_len_wr 
	uplo_len);

d_type_wr zla_porcond_x__(char *uplo, i_type_wr *n, z_type_wr *a, i_type_wr *
	lda, z_type_wr *af, i_type_wr *ldaf, z_type_wr *x, i_type_wr *
	info, z_type_wr *work, d_type_wr *rwork, ftn_len_wr uplo_len);

/* Subroutine */ int zla_porfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, z_type_wr *a, i_type_wr *lda, 
	z_type_wr *af, i_type_wr *ldaf, l_type_wr *colequ, d_type_wr *c__, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *y, i_type_wr *ldy, 
	d_type_wr *berr_out__, i_type_wr *n_norms__, d_type_wr *errs_n__, 
	d_type_wr *errs_c__, z_type_wr *res, d_type_wr *ayb, 
	z_type_wr *dy, z_type_wr *y_tail__, d_type_wr *rcond, 
	i_type_wr *ithresh, d_type_wr *rthresh, d_type_wr *dz_ub__, l_type_wr *
	ignore_cwise__, i_type_wr *info, ftn_len_wr uplo_len);

d_type_wr zla_porpvgrw__(char *uplo, i_type_wr *ncols, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *af, i_type_wr *ldaf, d_type_wr *work, 
	ftn_len_wr uplo_len);

d_type_wr zla_rpvgrw__(i_type_wr *n, i_type_wr *ncols, z_type_wr *a, i_type_wr 
	*lda, z_type_wr *af, i_type_wr *ldaf);

/* Subroutine */ int zla_syamv__(i_type_wr *uplo, i_type_wr *n, d_type_wr *alpha,
	 z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx, 
	d_type_wr *beta, d_type_wr *y, i_type_wr *incy);

d_type_wr zla_syrcond_c__(char *uplo, i_type_wr *n, z_type_wr *a, i_type_wr *
	lda, z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, d_type_wr *c__,
	 l_type_wr *capply, i_type_wr *info, z_type_wr *work, d_type_wr *
	rwork, ftn_len_wr uplo_len);

d_type_wr zla_syrcond_x__(char *uplo, i_type_wr *n, z_type_wr *a, i_type_wr *
	lda, z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, z_type_wr *
	x, i_type_wr *info, z_type_wr *work, d_type_wr *rwork, ftn_len_wr 
	uplo_len);

/* Subroutine */ int zla_syrfsx_extended__(i_type_wr *prec_type__, char *uplo, 
	i_type_wr *n, i_type_wr *nrhs, z_type_wr *a, i_type_wr *lda, 
	z_type_wr *af, i_type_wr *ldaf, i_type_wr *ipiv, l_type_wr *colequ, 
	d_type_wr *c__, z_type_wr *b, i_type_wr *ldb, z_type_wr *y, 
	i_type_wr *ldy, d_type_wr *berr_out__, i_type_wr *n_norms__, d_type_wr *
	errs_n__, d_type_wr *errs_c__, z_type_wr *res, d_type_wr *ayb, 
	z_type_wr *dy, z_type_wr *y_tail__, d_type_wr *rcond, 
	i_type_wr *ithresh, d_type_wr *rthresh, d_type_wr *dz_ub__, l_type_wr *
	ignore_cwise__, i_type_wr *info, ftn_len_wr uplo_len);

d_type_wr zla_syrpvgrw__(char *uplo, i_type_wr *n, i_type_wr *info, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, d_type_wr *work, ftn_len_wr uplo_len);

/* Subroutine */ int zla_wwaddw__(i_type_wr *n, z_type_wr *x, z_type_wr 
	*y, z_type_wr *w);

/* Subroutine */ int zlabrd_(i_type_wr *m, i_type_wr *n, i_type_wr *nb, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *d__, d_type_wr *e, 
	z_type_wr *tauq, z_type_wr *taup, z_type_wr *x, i_type_wr *
	ldx, z_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int zlacgv_(i_type_wr *n, z_type_wr *x, i_type_wr *incx);

/* Subroutine */ int zlacn2_(i_type_wr *n, z_type_wr *v, z_type_wr *x, 
	d_type_wr *est, i_type_wr *kase, i_type_wr *isave);

/* Subroutine */ int zlacon_(i_type_wr *n, z_type_wr *v, z_type_wr *x, 
	d_type_wr *est, i_type_wr *kase);

/* Subroutine */ int zlacp2_(char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *
	a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int zlacpy_(char *uplo, i_type_wr *m, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int zlacrm_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, z_type_wr *c__, 
	i_type_wr *ldc, d_type_wr *rwork);

/* Subroutine */ int zlacrt_(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, 
	z_type_wr *cy, i_type_wr *incy, z_type_wr *c__, z_type_wr *
	s);

/* Double c_type_wr */ void zladiv_(z_type_wr * ret_val, z_type_wr *x, 
	z_type_wr *y);

/* Subroutine */ int zlaed0_(i_type_wr *qsiz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, z_type_wr *q, i_type_wr *ldq, z_type_wr *qstore, 
	i_type_wr *ldqs, d_type_wr *rwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int zlaed7_(i_type_wr *n, i_type_wr *cutpnt, i_type_wr *qsiz, 
	i_type_wr *tlvls, i_type_wr *curlvl, i_type_wr *curpbm, d_type_wr *d__, 
	z_type_wr *q, i_type_wr *ldq, d_type_wr *rho, i_type_wr *indxq, 
	d_type_wr *qstore, i_type_wr *qptr, i_type_wr *prmptr, i_type_wr *perm, 
	i_type_wr *givptr, i_type_wr *givcol, d_type_wr *givnum, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int zlaed8_(i_type_wr *k, i_type_wr *n, i_type_wr *qsiz, 
	z_type_wr *q, i_type_wr *ldq, d_type_wr *d__, d_type_wr *rho, 
	i_type_wr *cutpnt, d_type_wr *z__, d_type_wr *dlamda, z_type_wr *
	q2, i_type_wr *ldq2, d_type_wr *w, i_type_wr *indxp, i_type_wr *indx, 
	i_type_wr *indxq, i_type_wr *perm, i_type_wr *givptr, i_type_wr *givcol, 
	d_type_wr *givnum, i_type_wr *info);

/* Subroutine */ int zlaein_(l_type_wr *rightv, l_type_wr *noinit, i_type_wr *n, 
	z_type_wr *h__, i_type_wr *ldh, z_type_wr *w, z_type_wr *v, 
	z_type_wr *b, i_type_wr *ldb, d_type_wr *rwork, d_type_wr *eps3, 
	d_type_wr *smlnum, i_type_wr *info);

/* Subroutine */ int zlaesy_(z_type_wr *a, z_type_wr *b, 
	z_type_wr *c__, z_type_wr *rt1, z_type_wr *rt2, 
	z_type_wr *evscal, z_type_wr *cs1, z_type_wr *sn1);

/* Subroutine */ int zlaev2_(z_type_wr *a, z_type_wr *b, 
	z_type_wr *c__, d_type_wr *rt1, d_type_wr *rt2, d_type_wr *cs1, 
	 z_type_wr *sn1);

/* Subroutine */ int zlag2c_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, c_type_wr *sa, i_type_wr *ldsa, i_type_wr *info);

/* Subroutine */ int zlags2_(l_type_wr *upper, d_type_wr *a1, z_type_wr *
	a2, d_type_wr *a3, d_type_wr *b1, z_type_wr *b2, d_type_wr *b3, 
	 d_type_wr *csu, z_type_wr *snu, d_type_wr *csv, z_type_wr *
	snv, d_type_wr *csq, z_type_wr *snq);

/* Subroutine */ int zlagtm_(char *trans, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *alpha, z_type_wr *dl, z_type_wr *d__, 
	z_type_wr *du, z_type_wr *x, i_type_wr *ldx, d_type_wr *beta, 
	z_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int zlahef_(char *uplo, i_type_wr *n, i_type_wr *nb, i_type_wr *kb, 
	 z_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, z_type_wr *w, 
	i_type_wr *ldw, i_type_wr *info);

/* Subroutine */ int zlahqr_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, z_type_wr *h__, i_type_wr *ldh, 
	z_type_wr *w, i_type_wr *iloz, i_type_wr *ihiz, z_type_wr *z__, 
	i_type_wr *ldz, i_type_wr *info);

/* Subroutine */ int zlahr2_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *t, 
	i_type_wr *ldt, z_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int zlahrd_(i_type_wr *n, i_type_wr *k, i_type_wr *nb, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *t, 
	i_type_wr *ldt, z_type_wr *y, i_type_wr *ldy);

/* Subroutine */ int zlaic1_(i_type_wr *job, i_type_wr *j, z_type_wr *x, 
	d_type_wr *sest, z_type_wr *w, z_type_wr *gamma, d_type_wr *
	sestpr, z_type_wr *s, z_type_wr *c__);

/* Subroutine */ int zlals0_(i_type_wr *icompq, i_type_wr *nl, i_type_wr *nr, 
	i_type_wr *sqre, i_type_wr *nrhs, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *bx, i_type_wr *ldbx, i_type_wr *perm, i_type_wr *givptr, 
	i_type_wr *givcol, i_type_wr *ldgcol, d_type_wr *givnum, i_type_wr *ldgnum, 
	 d_type_wr *poles, d_type_wr *difl, d_type_wr *difr, d_type_wr *
	z__, i_type_wr *k, d_type_wr *c__, d_type_wr *s, d_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int zlalsa_(i_type_wr *icompq, i_type_wr *smlsiz, i_type_wr *n, 
	i_type_wr *nrhs, z_type_wr *b, i_type_wr *ldb, z_type_wr *bx, 
	i_type_wr *ldbx, d_type_wr *u, i_type_wr *ldu, d_type_wr *vt, i_type_wr *
	k, d_type_wr *difl, d_type_wr *difr, d_type_wr *z__, d_type_wr *
	poles, i_type_wr *givptr, i_type_wr *givcol, i_type_wr *ldgcol, i_type_wr *
	perm, d_type_wr *givnum, d_type_wr *c__, d_type_wr *s, d_type_wr *
	rwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int zlalsd_(char *uplo, i_type_wr *smlsiz, i_type_wr *n, i_type_wr 
	*nrhs, d_type_wr *d__, d_type_wr *e, z_type_wr *b, i_type_wr *ldb, 
	 d_type_wr *rcond, i_type_wr *rank, z_type_wr *work, d_type_wr *
	rwork, i_type_wr *iwork, i_type_wr *info);

d_type_wr zlangb_(char *norm, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *work);

d_type_wr zlange_(char *norm, i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *work);

d_type_wr zlangt_(char *norm, i_type_wr *n, z_type_wr *dl, z_type_wr *
	d__, z_type_wr *du);

d_type_wr zlanhb_(char *norm, char *uplo, i_type_wr *n, i_type_wr *k, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *work);

d_type_wr zlanhe_(char *norm, char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *work);

d_type_wr zlanhf_(char *norm, char *transr, char *uplo, i_type_wr *n, 
	z_type_wr *a, d_type_wr *work);

d_type_wr zlanhp_(char *norm, char *uplo, i_type_wr *n, z_type_wr *ap, 
	d_type_wr *work);

d_type_wr zlanhs_(char *norm, i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	d_type_wr *work);

d_type_wr zlanht_(char *norm, i_type_wr *n, d_type_wr *d__, z_type_wr *e);

d_type_wr zlansb_(char *norm, char *uplo, i_type_wr *n, i_type_wr *k, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *work);

d_type_wr zlansp_(char *norm, char *uplo, i_type_wr *n, z_type_wr *ap, 
	d_type_wr *work);

d_type_wr zlansy_(char *norm, char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *work);

d_type_wr zlantb_(char *norm, char *uplo, char *diag, i_type_wr *n, i_type_wr *k, 
	 z_type_wr *ab, i_type_wr *ldab, d_type_wr *work);

d_type_wr zlantp_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	z_type_wr *ap, d_type_wr *work);

d_type_wr zlantr_(char *norm, char *uplo, char *diag, i_type_wr *m, i_type_wr *n, 
	 z_type_wr *a, i_type_wr *lda, d_type_wr *work);

/* Subroutine */ int zlapll_(i_type_wr *n, z_type_wr *x, i_type_wr *incx, 
	z_type_wr *y, i_type_wr *incy, d_type_wr *ssmin);

/* Subroutine */ int zlapmt_(l_type_wr *forwrd, i_type_wr *m, i_type_wr *n, 
	z_type_wr *x, i_type_wr *ldx, i_type_wr *k);

/* Subroutine */ int zlaqgb_(i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku, 
	 z_type_wr *ab, i_type_wr *ldab, d_type_wr *r__, d_type_wr *c__, 
	d_type_wr *rowcnd, d_type_wr *colcnd, d_type_wr *amax, char *equed);

/* Subroutine */ int zlaqge_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *r__, d_type_wr *c__, d_type_wr *rowcnd, 
	d_type_wr *colcnd, d_type_wr *amax, char *equed);

/* Subroutine */ int zlaqhb_(char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *s, d_type_wr *scond, 
	d_type_wr *amax, char *equed);

/* Subroutine */ int zlaqhe_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, 
	char *equed);

/* Subroutine */ int zlaqhp_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, char *equed);

/* Subroutine */ int zlaqp2_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *jpvt, z_type_wr *tau, 
	d_type_wr *vn1, d_type_wr *vn2, z_type_wr *work);

/* Subroutine */ int zlaqps_(i_type_wr *m, i_type_wr *n, i_type_wr *offset, i_type_wr 
	*nb, i_type_wr *kb, z_type_wr *a, i_type_wr *lda, i_type_wr *jpvt, 
	z_type_wr *tau, d_type_wr *vn1, d_type_wr *vn2, z_type_wr *
	auxv, z_type_wr *f, i_type_wr *ldf);

/* Subroutine */ int zlaqr0_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, z_type_wr *h__, i_type_wr *ldh, 
	z_type_wr *w, i_type_wr *iloz, i_type_wr *ihiz, z_type_wr *z__, 
	i_type_wr *ldz, z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zlaqr1_(i_type_wr *n, z_type_wr *h__, i_type_wr *ldh, 
	z_type_wr *s1, z_type_wr *s2, z_type_wr *v);

/* Subroutine */ int zlaqr2_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, z_type_wr *h__, 
	i_type_wr *ldh, i_type_wr *iloz, i_type_wr *ihiz, z_type_wr *z__, 
	i_type_wr *ldz, i_type_wr *ns, i_type_wr *nd, z_type_wr *sh, 
	z_type_wr *v, i_type_wr *ldv, i_type_wr *nh, z_type_wr *t, 
	i_type_wr *ldt, i_type_wr *nv, z_type_wr *wv, i_type_wr *ldwv, 
	z_type_wr *work, i_type_wr *lwork);

/* Subroutine */ int zlaqr3_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nw, z_type_wr *h__, 
	i_type_wr *ldh, i_type_wr *iloz, i_type_wr *ihiz, z_type_wr *z__, 
	i_type_wr *ldz, i_type_wr *ns, i_type_wr *nd, z_type_wr *sh, 
	z_type_wr *v, i_type_wr *ldv, i_type_wr *nh, z_type_wr *t, 
	i_type_wr *ldt, i_type_wr *nv, z_type_wr *wv, i_type_wr *ldwv, 
	z_type_wr *work, i_type_wr *lwork);

/* Subroutine */ int zlaqr4_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, z_type_wr *h__, i_type_wr *ldh, 
	z_type_wr *w, i_type_wr *iloz, i_type_wr *ihiz, z_type_wr *z__, 
	i_type_wr *ldz, z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zlaqr5_(l_type_wr *wantt, l_type_wr *wantz, i_type_wr *kacc22, 
	i_type_wr *n, i_type_wr *ktop, i_type_wr *kbot, i_type_wr *nshfts, 
	z_type_wr *s, z_type_wr *h__, i_type_wr *ldh, i_type_wr *iloz, 
	i_type_wr *ihiz, z_type_wr *z__, i_type_wr *ldz, z_type_wr *v, 
	i_type_wr *ldv, z_type_wr *u, i_type_wr *ldu, i_type_wr *nv, 
	z_type_wr *wv, i_type_wr *ldwv, i_type_wr *nh, z_type_wr *wh, 
	i_type_wr *ldwh);

/* Subroutine */ int zlaqsb_(char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *s, d_type_wr *scond, 
	d_type_wr *amax, char *equed);

/* Subroutine */ int zlaqsp_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, char *equed);

/* Subroutine */ int zlaqsy_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, 
	char *equed);

/* Subroutine */ int zlar1v_(i_type_wr *n, i_type_wr *b1, i_type_wr *bn, d_type_wr 
	*lambda, d_type_wr *d__, d_type_wr *l, d_type_wr *ld, d_type_wr *
	lld, d_type_wr *pivmin, d_type_wr *gaptol, z_type_wr *z__, 
	l_type_wr *wantnc, i_type_wr *negcnt, d_type_wr *ztz, d_type_wr *mingma, 
	 i_type_wr *r__, i_type_wr *isuppz, d_type_wr *nrminv, d_type_wr *resid, 
	 d_type_wr *rqcorr, d_type_wr *work);

/* Subroutine */ int zlar2v_(i_type_wr *n, z_type_wr *x, z_type_wr *y, 
	z_type_wr *z__, i_type_wr *incx, d_type_wr *c__, z_type_wr *s, 
	i_type_wr *incc);

/* Subroutine */ int zlarcm_(i_type_wr *m, i_type_wr *n, d_type_wr *a, i_type_wr *
	lda, z_type_wr *b, i_type_wr *ldb, z_type_wr *c__, i_type_wr *ldc, 
	 d_type_wr *rwork);

/* Subroutine */ int zlarf_(char *side, i_type_wr *m, i_type_wr *n, z_type_wr 
	*v, i_type_wr *incv, z_type_wr *tau, z_type_wr *c__, i_type_wr *
	ldc, z_type_wr *work);

/* Subroutine */ int zlarfb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, z_type_wr *v, i_type_wr 
	*ldv, z_type_wr *t, i_type_wr *ldt, z_type_wr *c__, i_type_wr *
	ldc, z_type_wr *work, i_type_wr *ldwork);

/* Subroutine */ int zlarfg_(i_type_wr *n, z_type_wr *alpha, z_type_wr *
	x, i_type_wr *incx, z_type_wr *tau);

/* Subroutine */ int zlarfp_(i_type_wr *n, z_type_wr *alpha, z_type_wr *
	x, i_type_wr *incx, z_type_wr *tau);

/* Subroutine */ int zlarft_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, z_type_wr *v, i_type_wr *ldv, z_type_wr *tau, z_type_wr *
	t, i_type_wr *ldt);

/* Subroutine */ int zlarfx_(char *side, i_type_wr *m, i_type_wr *n, 
	z_type_wr *v, z_type_wr *tau, z_type_wr *c__, i_type_wr *
	ldc, z_type_wr *work);

/* Subroutine */ int zlargv_(i_type_wr *n, z_type_wr *x, i_type_wr *incx, 
	z_type_wr *y, i_type_wr *incy, d_type_wr *c__, i_type_wr *incc);

/* Subroutine */ int zlarnv_(i_type_wr *idist, i_type_wr *iseed, i_type_wr *n, 
	z_type_wr *x);

/* Subroutine */ int zlarrv_(i_type_wr *n, d_type_wr *vl, d_type_wr *vu, 
	d_type_wr *d__, d_type_wr *l, d_type_wr *pivmin, i_type_wr *isplit, 
	i_type_wr *m, i_type_wr *dol, i_type_wr *dou, d_type_wr *minrgp, 
	d_type_wr *rtol1, d_type_wr *rtol2, d_type_wr *w, d_type_wr *werr, 
	 d_type_wr *wgap, i_type_wr *iblock, i_type_wr *indexw, d_type_wr *gers, 
	 z_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, d_type_wr *work, 
	i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int zlarscl2_(i_type_wr *m, i_type_wr *n, d_type_wr *d__, 
	z_type_wr *x, i_type_wr *ldx);

/* Subroutine */ int zlartg_(z_type_wr *f, z_type_wr *g, d_type_wr *
	cs, z_type_wr *sn, z_type_wr *r__);

/* Subroutine */ int zlartv_(i_type_wr *n, z_type_wr *x, i_type_wr *incx, 
	z_type_wr *y, i_type_wr *incy, d_type_wr *c__, z_type_wr *s, 
	i_type_wr *incc);

/* Subroutine */ int zlarz_(char *side, i_type_wr *m, i_type_wr *n, i_type_wr *l, 
	z_type_wr *v, i_type_wr *incv, z_type_wr *tau, z_type_wr *
	c__, i_type_wr *ldc, z_type_wr *work);

/* Subroutine */ int zlarzb_(char *side, char *trans, char *direct, char *
	storev, i_type_wr *m, i_type_wr *n, i_type_wr *k, i_type_wr *l, z_type_wr 
	*v, i_type_wr *ldv, z_type_wr *t, i_type_wr *ldt, z_type_wr *c__, 
	i_type_wr *ldc, z_type_wr *work, i_type_wr *ldwork);

/* Subroutine */ int zlarzt_(char *direct, char *storev, i_type_wr *n, i_type_wr *
	k, z_type_wr *v, i_type_wr *ldv, z_type_wr *tau, z_type_wr *
	t, i_type_wr *ldt);

/* Subroutine */ int zlascl_(char *type__, i_type_wr *kl, i_type_wr *ku, 
	d_type_wr *cfrom, d_type_wr *cto, i_type_wr *m, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int zlascl2_(i_type_wr *m, i_type_wr *n, d_type_wr *d__, 
	z_type_wr *x, i_type_wr *ldx);

/* Subroutine */ int zlaset_(char *uplo, i_type_wr *m, i_type_wr *n, 
	z_type_wr *alpha, z_type_wr *beta, z_type_wr *a, i_type_wr *
	lda);

/* Subroutine */ int zlasr_(char *side, char *pivot, char *direct, i_type_wr *m, 
	 i_type_wr *n, d_type_wr *c__, d_type_wr *s, z_type_wr *a, 
	i_type_wr *lda);

/* Subroutine */ int zlassq_(i_type_wr *n, z_type_wr *x, i_type_wr *incx, 
	d_type_wr *scale, d_type_wr *sumsq);

/* Subroutine */ int zlaswp_(i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	i_type_wr *k1, i_type_wr *k2, i_type_wr *ipiv, i_type_wr *incx);

/* Subroutine */ int zlasyf_(char *uplo, i_type_wr *n, i_type_wr *nb, i_type_wr *kb, 
	 z_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, z_type_wr *w, 
	i_type_wr *ldw, i_type_wr *info);

/* Subroutine */ int zlat2c_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, c_type_wr *sa, i_type_wr *ldsa, i_type_wr *info);

/* Subroutine */ int zlatbs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, i_type_wr *kd, z_type_wr *ab, i_type_wr *ldab, 
	z_type_wr *x, d_type_wr *scale, d_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int zlatdf_(i_type_wr *ijob, i_type_wr *n, z_type_wr *z__, 
	i_type_wr *ldz, z_type_wr *rhs, d_type_wr *rdsum, d_type_wr *
	rdscal, i_type_wr *ipiv, i_type_wr *jpiv);

/* Subroutine */ int zlatps_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, z_type_wr *ap, z_type_wr *x, d_type_wr *
	scale, d_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int zlatrd_(char *uplo, i_type_wr *n, i_type_wr *nb, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *e, z_type_wr *tau, 
	z_type_wr *w, i_type_wr *ldw);

/* Subroutine */ int zlatrs_(char *uplo, char *trans, char *diag, char *
	normin, i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *x, 
	d_type_wr *scale, d_type_wr *cnorm, i_type_wr *info);

/* Subroutine */ int zlatrz_(i_type_wr *m, i_type_wr *n, i_type_wr *l, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work);

/* Subroutine */ int zlatzm_(char *side, i_type_wr *m, i_type_wr *n, 
	z_type_wr *v, i_type_wr *incv, z_type_wr *tau, z_type_wr *
	c1, z_type_wr *c2, i_type_wr *ldc, z_type_wr *work);

/* Subroutine */ int zlauu2_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int zlauum_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int zpbcon_(char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *anorm, d_type_wr *
	rcond, z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zpbequ_(char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, d_type_wr *s, d_type_wr *scond, 
	d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int zpbrfs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, z_type_wr *ab, i_type_wr *ldab, z_type_wr *afb, i_type_wr *
	ldafb, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	 d_type_wr *ferr, d_type_wr *berr, z_type_wr *work, d_type_wr *
	rwork, i_type_wr *info);

/* Subroutine */ int zpbstf_(char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int zpbsv_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, z_type_wr *ab, i_type_wr *ldab, z_type_wr *b, i_type_wr *
	ldb, i_type_wr *info);

/* Subroutine */ int zpbsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *kd, 
	i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, z_type_wr *afb, 
	i_type_wr *ldafb, char *equed, d_type_wr *s, z_type_wr *b, i_type_wr 
	*ldb, z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *
	ferr, d_type_wr *berr, z_type_wr *work, d_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int zpbtf2_(char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int zpbtrf_(char *uplo, i_type_wr *n, i_type_wr *kd, 
	z_type_wr *ab, i_type_wr *ldab, i_type_wr *info);

/* Subroutine */ int zpbtrs_(char *uplo, i_type_wr *n, i_type_wr *kd, i_type_wr *
	nrhs, z_type_wr *ab, i_type_wr *ldab, z_type_wr *b, i_type_wr *
	ldb, i_type_wr *info);

/* Subroutine */ int zpftrf_(char *transr, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *info);

/* Subroutine */ int zpftri_(char *transr, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *info);

/* Subroutine */ int zpftrs_(char *transr, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, z_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zpocon_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *anorm, d_type_wr *rcond, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zpoequ_(i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int zpoequb_(i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int zporfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *ldaf, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *ferr, d_type_wr *berr, z_type_wr *work, d_type_wr *
	rwork, i_type_wr *info);

/* Subroutine */ int zporfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, d_type_wr *s, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, 
	 i_type_wr *ldx, d_type_wr *rcond, d_type_wr *berr, i_type_wr *
	n_err_bnds__, d_type_wr *err_bnds_norm__, d_type_wr *
	err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zposv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zposvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, char *equed, d_type_wr *s, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, 
	d_type_wr *berr, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zposvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, char *equed, d_type_wr *s, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *rpvgrw, 
	 d_type_wr *berr, i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, 
	 d_type_wr *err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zpotf2_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int zpotrf_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int zpotri_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int zpotrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zppcon_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	d_type_wr *anorm, d_type_wr *rcond, z_type_wr *work, d_type_wr 
	*rwork, i_type_wr *info);

/* Subroutine */ int zppequ_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	d_type_wr *s, d_type_wr *scond, d_type_wr *amax, i_type_wr *info);

/* Subroutine */ int zpprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, z_type_wr *afp, z_type_wr *b, i_type_wr *ldb, 
	 z_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zppsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, z_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zppsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *ap, z_type_wr *afp, char *equed, d_type_wr *
	s, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zpptrf_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *info);

/* Subroutine */ int zpptri_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *info);

/* Subroutine */ int zpptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, z_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zpstf2_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *piv, i_type_wr *rank, d_type_wr *tol, 
	d_type_wr *work, i_type_wr *info);

/* Subroutine */ int zpstrf_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *piv, i_type_wr *rank, d_type_wr *tol, 
	d_type_wr *work, i_type_wr *info);

/* Subroutine */ int zptcon_(i_type_wr *n, d_type_wr *d__, z_type_wr *e, 
	d_type_wr *anorm, d_type_wr *rcond, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zpteqr_(char *compz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, z_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int zptrfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *d__, z_type_wr *e, d_type_wr *df, z_type_wr *ef, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *ferr, d_type_wr *berr, z_type_wr *work, d_type_wr *
	rwork, i_type_wr *info);

/* Subroutine */ int zptsv_(i_type_wr *n, i_type_wr *nrhs, d_type_wr *d__, 
	z_type_wr *e, z_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int zptsvx_(char *fact, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *d__, z_type_wr *e, d_type_wr *df, z_type_wr *ef, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zpttrf_(i_type_wr *n, d_type_wr *d__, z_type_wr *e, 
	i_type_wr *info);

/* Subroutine */ int zpttrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *d__, z_type_wr *e, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zptts2_(i_type_wr *iuplo, i_type_wr *n, i_type_wr *nrhs, 
	d_type_wr *d__, z_type_wr *e, z_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int zrot_(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, 
	z_type_wr *cy, i_type_wr *incy, d_type_wr *c__, z_type_wr *s);

/* Subroutine */ int zspcon_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *ipiv, d_type_wr *anorm, d_type_wr *rcond, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zspmv_(char *uplo, i_type_wr *n, z_type_wr *alpha, 
	z_type_wr *ap, z_type_wr *x, i_type_wr *incx, z_type_wr *
	beta, z_type_wr *y, i_type_wr *incy);

/* Subroutine */ int zspr_(char *uplo, i_type_wr *n, z_type_wr *alpha, 
	z_type_wr *x, i_type_wr *incx, z_type_wr *ap);

/* Subroutine */ int zsprfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, z_type_wr *afp, i_type_wr *ipiv, z_type_wr *
	b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, 
	d_type_wr *berr, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int zspsv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zspsvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *ap, z_type_wr *afp, i_type_wr *ipiv, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zsptrf_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zsptri_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	i_type_wr *ipiv, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zsptrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *ap, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int zstedc_(char *compz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, z_type_wr *z__, i_type_wr *ldz, z_type_wr *work, 
	i_type_wr *lwork, d_type_wr *rwork, i_type_wr *lrwork, i_type_wr *iwork, 
	i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zstegr_(char *jobz, char *range, i_type_wr *n, d_type_wr *
	d__, d_type_wr *e, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, 
	i_type_wr *iu, d_type_wr *abstol, i_type_wr *m, d_type_wr *w, 
	z_type_wr *z__, i_type_wr *ldz, i_type_wr *isuppz, d_type_wr *work, 
	i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zstein_(i_type_wr *n, d_type_wr *d__, d_type_wr *e, 
	i_type_wr *m, d_type_wr *w, i_type_wr *iblock, i_type_wr *isplit, 
	z_type_wr *z__, i_type_wr *ldz, d_type_wr *work, i_type_wr *iwork, 
	i_type_wr *ifail, i_type_wr *info);

/* Subroutine */ int zstemr_(char *jobz, char *range, i_type_wr *n, d_type_wr *
	d__, d_type_wr *e, d_type_wr *vl, d_type_wr *vu, i_type_wr *il, 
	i_type_wr *iu, i_type_wr *m, d_type_wr *w, z_type_wr *z__, i_type_wr *
	ldz, i_type_wr *nzc, i_type_wr *isuppz, l_type_wr *tryrac, d_type_wr *work, 
	 i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, i_type_wr *info);

/* Subroutine */ int zsteqr_(char *compz, i_type_wr *n, d_type_wr *d__, 
	d_type_wr *e, z_type_wr *z__, i_type_wr *ldz, d_type_wr *work, 
	i_type_wr *info);

/* Subroutine */ int zsycon_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, d_type_wr *anorm, d_type_wr *rcond, 
	z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zsyequb_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, d_type_wr *s, d_type_wr *scond, d_type_wr *amax, 
	z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zsymv_(char *uplo, i_type_wr *n, z_type_wr *alpha, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx, 
	z_type_wr *beta, z_type_wr *y, i_type_wr *incy);

/* Subroutine */ int zsyr_(char *uplo, i_type_wr *n, z_type_wr *alpha, 
	z_type_wr *x, i_type_wr *incx, z_type_wr *a, i_type_wr *lda);

/* Subroutine */ int zsyrfs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *ldaf, 
	i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, 
	i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, z_type_wr *work, 
	 d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zsyrfsx_(char *uplo, char *equed, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, d_type_wr *s, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, d_type_wr *berr, 
	i_type_wr *n_err_bnds__, d_type_wr *err_bnds_norm__, d_type_wr *
	err_bnds_comp__, i_type_wr *nparams, d_type_wr *params, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zsysv_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zsysvx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, z_type_wr *b, i_type_wr *ldb, z_type_wr *x, 
	 i_type_wr *ldx, d_type_wr *rcond, d_type_wr *ferr, d_type_wr *berr, 
	z_type_wr *work, i_type_wr *lwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int zsysvxx_(char *fact, char *uplo, i_type_wr *n, i_type_wr *
	nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *af, i_type_wr *
	ldaf, i_type_wr *ipiv, char *equed, d_type_wr *s, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, d_type_wr *rcond, 
	d_type_wr *rpvgrw, d_type_wr *berr, i_type_wr *n_err_bnds__, 
	d_type_wr *err_bnds_norm__, d_type_wr *err_bnds_comp__, i_type_wr *
	nparams, d_type_wr *params, z_type_wr *work, d_type_wr *rwork, 
	i_type_wr *info);

/* Subroutine */ int zsytf2_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, i_type_wr *info);

/* Subroutine */ int zsytrf_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, z_type_wr *work, i_type_wr *lwork, 
	i_type_wr *info);

/* Subroutine */ int zsytri_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, i_type_wr *ipiv, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zsytrs_(char *uplo, i_type_wr *n, i_type_wr *nrhs, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *ipiv, z_type_wr *b, 
	i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int ztbcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	i_type_wr *kd, z_type_wr *ab, i_type_wr *ldab, d_type_wr *rcond, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ztbrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, 
	d_type_wr *ferr, d_type_wr *berr, z_type_wr *work, d_type_wr *
	rwork, i_type_wr *info);

/* Subroutine */ int ztbtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *kd, i_type_wr *nrhs, z_type_wr *ab, i_type_wr *ldab, 
	z_type_wr *b, i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int ztfsm_(char *transr, char *side, char *uplo, char *trans, 
	 char *diag, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, 
	z_type_wr *a, z_type_wr *b, i_type_wr *ldb);

/* Subroutine */ int ztftri_(char *transr, char *uplo, char *diag, i_type_wr *n, 
	 z_type_wr *a, i_type_wr *info);

/* Subroutine */ int ztfttp_(char *transr, char *uplo, i_type_wr *n, 
	z_type_wr *arf, z_type_wr *ap, i_type_wr *info);

/* Subroutine */ int ztfttr_(char *transr, char *uplo, i_type_wr *n, 
	z_type_wr *arf, z_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ztgevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, z_type_wr *s, i_type_wr *lds, z_type_wr *p, i_type_wr 
	*ldp, z_type_wr *vl, i_type_wr *ldvl, z_type_wr *vr, i_type_wr *
	ldvr, i_type_wr *mm, i_type_wr *m, z_type_wr *work, d_type_wr *rwork, 
	 i_type_wr *info);

/* Subroutine */ int ztgex2_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *q, i_type_wr *ldq, z_type_wr *z__, i_type_wr *ldz, 
	i_type_wr *j1, i_type_wr *info);

/* Subroutine */ int ztgexc_(l_type_wr *wantq, l_type_wr *wantz, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *q, i_type_wr *ldq, z_type_wr *z__, i_type_wr *ldz, 
	i_type_wr *ifst, i_type_wr *ilst, i_type_wr *info);

/* Subroutine */ int ztgsen_(i_type_wr *ijob, l_type_wr *wantq, l_type_wr *wantz, 
	l_type_wr *select, i_type_wr *n, z_type_wr *a, i_type_wr *lda, 
	z_type_wr *b, i_type_wr *ldb, z_type_wr *alpha, z_type_wr *
	beta, z_type_wr *q, i_type_wr *ldq, z_type_wr *z__, i_type_wr *
	ldz, i_type_wr *m, d_type_wr *pl, d_type_wr *pr, d_type_wr *dif, 
	z_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *liwork, 
	i_type_wr *info);

/* Subroutine */ int ztgsja_(char *jobu, char *jobv, char *jobq, i_type_wr *m, 
	i_type_wr *p, i_type_wr *n, i_type_wr *k, i_type_wr *l, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, d_type_wr *tola, 
	d_type_wr *tolb, d_type_wr *alpha, d_type_wr *beta, z_type_wr *
	u, i_type_wr *ldu, z_type_wr *v, i_type_wr *ldv, z_type_wr *q, 
	i_type_wr *ldq, z_type_wr *work, i_type_wr *ncycle, i_type_wr *info);

/* Subroutine */ int ztgsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr 
	*ldb, z_type_wr *vl, i_type_wr *ldvl, z_type_wr *vr, i_type_wr *
	ldvr, d_type_wr *s, d_type_wr *dif, i_type_wr *mm, i_type_wr *m, 
	z_type_wr *work, i_type_wr *lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int ztgsy2_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *d__, i_type_wr *ldd, 
	z_type_wr *e, i_type_wr *lde, z_type_wr *f, i_type_wr *ldf, 
	d_type_wr *scale, d_type_wr *rdsum, d_type_wr *rdscal, i_type_wr *
	info);

/* Subroutine */ int ztgsyl_(char *trans, i_type_wr *ijob, i_type_wr *m, i_type_wr *
	n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *d__, i_type_wr *ldd, 
	z_type_wr *e, i_type_wr *lde, z_type_wr *f, i_type_wr *ldf, 
	d_type_wr *scale, d_type_wr *dif, z_type_wr *work, i_type_wr *
	lwork, i_type_wr *iwork, i_type_wr *info);

/* Subroutine */ int ztpcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	z_type_wr *ap, d_type_wr *rcond, z_type_wr *work, d_type_wr 
	*rwork, i_type_wr *info);

/* Subroutine */ int ztprfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, z_type_wr *ap, z_type_wr *b, i_type_wr *ldb, 
	z_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, d_type_wr *berr, 
	z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ztptri_(char *uplo, char *diag, i_type_wr *n, 
	z_type_wr *ap, i_type_wr *info);

/* Subroutine */ int ztptrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, z_type_wr *ap, z_type_wr *b, i_type_wr *ldb, 
	i_type_wr *info);

/* Subroutine */ int ztpttf_(char *transr, char *uplo, i_type_wr *n, 
	z_type_wr *ap, z_type_wr *arf, i_type_wr *info);

/* Subroutine */ int ztpttr_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ztrcon_(char *norm, char *uplo, char *diag, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, d_type_wr *rcond, z_type_wr *
	work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ztrevc_(char *side, char *howmny, l_type_wr *select, 
	i_type_wr *n, z_type_wr *t, i_type_wr *ldt, z_type_wr *vl, 
	i_type_wr *ldvl, z_type_wr *vr, i_type_wr *ldvr, i_type_wr *mm, i_type_wr 
	*m, z_type_wr *work, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ztrexc_(char *compq, i_type_wr *n, z_type_wr *t, 
	i_type_wr *ldt, z_type_wr *q, i_type_wr *ldq, i_type_wr *ifst, i_type_wr *
	ilst, i_type_wr *info);

/* Subroutine */ int ztrrfs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *x, i_type_wr *ldx, d_type_wr *ferr, 
	d_type_wr *berr, z_type_wr *work, d_type_wr *rwork, i_type_wr *
	info);

/* Subroutine */ int ztrsen_(char *job, char *compq, l_type_wr *select, i_type_wr 
	*n, z_type_wr *t, i_type_wr *ldt, z_type_wr *q, i_type_wr *ldq, 
	z_type_wr *w, i_type_wr *m, d_type_wr *s, d_type_wr *sep, 
	z_type_wr *work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int ztrsna_(char *job, char *howmny, l_type_wr *select, 
	i_type_wr *n, z_type_wr *t, i_type_wr *ldt, z_type_wr *vl, 
	i_type_wr *ldvl, z_type_wr *vr, i_type_wr *ldvr, d_type_wr *s, 
	d_type_wr *sep, i_type_wr *mm, i_type_wr *m, z_type_wr *work, 
	i_type_wr *ldwork, d_type_wr *rwork, i_type_wr *info);

/* Subroutine */ int ztrsyl_(char *trana, char *tranb, i_type_wr *isgn, i_type_wr 
	*m, i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
	i_type_wr *ldb, z_type_wr *c__, i_type_wr *ldc, d_type_wr *scale, 
	i_type_wr *info);

/* Subroutine */ int ztrti2_(char *uplo, char *diag, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ztrtri_(char *uplo, char *diag, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, i_type_wr *info);

/* Subroutine */ int ztrtrs_(char *uplo, char *trans, char *diag, i_type_wr *n, 
	i_type_wr *nrhs, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
	i_type_wr *ldb, i_type_wr *info);

/* Subroutine */ int ztrttf_(char *transr, char *uplo, i_type_wr *n, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *arf, i_type_wr *info);

/* Subroutine */ int ztrttp_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *ap, i_type_wr *info);

/* Subroutine */ int ztzrqf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, i_type_wr *info);

/* Subroutine */ int ztzrzf_(i_type_wr *m, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zung2l_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zung2r_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zungbr_(char *vect, i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zunghr_(i_type_wr *n, i_type_wr *ilo, i_type_wr *ihi, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zungl2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zunglq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zungql_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zungqr_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zungr2_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zungrq_(i_type_wr *m, i_type_wr *n, i_type_wr *k, 
	z_type_wr *a, i_type_wr *lda, z_type_wr *tau, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zungtr_(char *uplo, i_type_wr *n, z_type_wr *a, 
	i_type_wr *lda, z_type_wr *tau, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zunm2l_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zunm2r_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zunmbr_(char *vect, char *side, char *trans, i_type_wr *m, 
	i_type_wr *n, i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr 
	*tau, z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int zunmhr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *ilo, i_type_wr *ihi, z_type_wr *a, i_type_wr *lda, 
	z_type_wr *tau, z_type_wr *c__, i_type_wr *ldc, z_type_wr *
	work, i_type_wr *lwork, i_type_wr *info);

/* Subroutine */ int zunml2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zunmlq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zunmql_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zunmqr_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zunmr2_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int zunmr3_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, z_type_wr *a, i_type_wr *lda, z_type_wr 
	*tau, z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *
	info);

/* Subroutine */ int zunmrq_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zunmrz_(char *side, char *trans, i_type_wr *m, i_type_wr *n, 
	i_type_wr *k, i_type_wr *l, z_type_wr *a, i_type_wr *lda, z_type_wr 
	*tau, z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *
	lwork, i_type_wr *info);

/* Subroutine */ int zunmtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, z_type_wr *a, i_type_wr *lda, z_type_wr *tau, 
	z_type_wr *c__, i_type_wr *ldc, z_type_wr *work, i_type_wr *lwork, 
	 i_type_wr *info);

/* Subroutine */ int zupgtr_(char *uplo, i_type_wr *n, z_type_wr *ap, 
	z_type_wr *tau, z_type_wr *q, i_type_wr *ldq, z_type_wr *
	work, i_type_wr *info);

/* Subroutine */ int zupmtr_(char *side, char *uplo, char *trans, i_type_wr *m, 
	i_type_wr *n, z_type_wr *ap, z_type_wr *tau, z_type_wr *c__, 
	 i_type_wr *ldc, z_type_wr *work, i_type_wr *info);

/* Subroutine */ int dlamc1_(i_type_wr *beta, i_type_wr *t, l_type_wr *rnd, l_type_wr 
	*ieee1);

d_type_wr dsecnd_();

/* Subroutine */ int ilaver_(i_type_wr *vers_major__, i_type_wr *vers_minor__, 
	i_type_wr *vers_patch__);

d_type_wr second_();

s_type_wr slamch_(char *cmach);

/* Subroutine */ int slamc1_(i_type_wr *beta, i_type_wr *t, l_type_wr *rnd, l_type_wr 
	*ieee1);

/* Subroutine */ int slamc2_(i_type_wr *beta, i_type_wr *t, l_type_wr *rnd, s_type_wr *
		    eps, i_type_wr *emin, s_type_wr *rmin, i_type_wr *emax, s_type_wr *rmax);

s_type_wr slamc3_(s_type_wr *a, s_type_wr *b);

/* Subroutine */ int slamc4_(i_type_wr *emin, s_type_wr *start, i_type_wr *base);

/* Subroutine */ int slamc5_(i_type_wr *beta, i_type_wr *p, i_type_wr *emin,
		    l_type_wr *ieee, i_type_wr *emax, s_type_wr *rmax);


d_type_wr dlamch_(char *cmach);

/* Subroutine */ int dlamc1_(i_type_wr *beta, i_type_wr *t, l_type_wr *rnd, l_type_wr
		    *ieee1);

/* Subroutine */ int dlamc2_(i_type_wr *beta, i_type_wr *t, l_type_wr *rnd,
		    d_type_wr *eps, i_type_wr *emin, d_type_wr *rmin, i_type_wr *emax,
			    d_type_wr *rmax);

d_type_wr dlamc3_(d_type_wr *a, d_type_wr *b);

/* Subroutine */ int dlamc4_(i_type_wr *emin, d_type_wr *start, i_type_wr *base);

/* Subroutine */ int dlamc5_(i_type_wr *beta, i_type_wr *p, i_type_wr *emin,
		    l_type_wr *ieee, i_type_wr *emax, d_type_wr *rmax);

i_type_wr ilaenv_(i_type_wr *ispec, char *name__, char *opts, i_type_wr *n1, 
	i_type_wr *n2, i_type_wr *n3, i_type_wr *n4);

/* Subroutine */ int zsbmv_(char *uplo, i_type_wr *n, i_type_wr *k, z_type_wr 
	*alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *
	incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy);

/* Subroutine */ int csbmv_(char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *
	alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *
	beta, c_type_wr *y, i_type_wr *incy);

int dggsvd3_(const char* JOBU, const char* JOBV, const char* JOBQ, i_type_wr* M, i_type_wr* N, 
    i_type_wr* P, i_type_wr* K, i_type_wr* L, d_type_wr* ptr_A, i_type_wr* LDA, d_type_wr* ptr_B, i_type_wr* LDB, 
    d_type_wr* ptr_alpha, d_type_wr* ptr_beta, d_type_wr* ptr_U, i_type_wr* LDU, d_type_wr* ptr_V, i_type_wr* LDV, 
    d_type_wr* ptr_Q, i_type_wr* LDQ, d_type_wr* work, i_type_wr* lwork, i_type_wr* iwork, i_type_wr* info);

int sggsvd3_(const char* JOBU, const char* JOBV, const char* JOBQ, i_type_wr* M, i_type_wr* N, 
    i_type_wr* P, i_type_wr* K, i_type_wr* L, s_type_wr* ptr_A, i_type_wr* LDA, s_type_wr* ptr_B, i_type_wr* LDB, 
    s_type_wr* ptr_alpha, s_type_wr* ptr_beta, s_type_wr* ptr_U, i_type_wr* LDU, s_type_wr* ptr_V, i_type_wr* LDV, 
    s_type_wr* ptr_Q, i_type_wr* LDQ, s_type_wr* work, i_type_wr* lwork, i_type_wr* iwork, i_type_wr* info);

int cggsvd3_(const char* JOBU, const char* JOBV, const char* JOBQ, i_type_wr* M, i_type_wr* N, 
    i_type_wr* P, i_type_wr* K, i_type_wr* L, c_type_wr* ptr_A, i_type_wr* LDA, c_type_wr* ptr_B, i_type_wr* LDB, 
    s_type_wr* ptr_alpha, s_type_wr* ptr_beta, c_type_wr* ptr_U, i_type_wr* LDU, c_type_wr* ptr_V, i_type_wr* LDV, 
    c_type_wr* ptr_Q, i_type_wr* LDQ, c_type_wr* work, i_type_wr* lwork, s_type_wr* rwork, i_type_wr* iwork, i_type_wr* info);

int zggsvd3_(const char* JOBU, const char* JOBV, const char* JOBQ, i_type_wr* M, i_type_wr* N, 
    i_type_wr* P, i_type_wr* K, i_type_wr* L, z_type_wr* ptr_A, i_type_wr* LDA, z_type_wr* ptr_B, i_type_wr* LDB, 
    d_type_wr* ptr_alpha, d_type_wr* ptr_beta, z_type_wr* ptr_U, i_type_wr* LDU, z_type_wr* ptr_V, i_type_wr* LDV, 
    z_type_wr* ptr_Q, i_type_wr* LDQ, z_type_wr* work, i_type_wr* lwork, d_type_wr* rwork, i_type_wr* iwork, i_type_wr* info);

int dgesvdx_(const char* JOBU, const char* JOBVT, const char* RANGE, i_type_wr* M, i_type_wr* N, d_type_wr* A, 
    i_type_wr* LDA, d_type_wr* VL, d_type_wr* VU, i_type_wr* IL, i_type_wr* IU, i_type_wr* NS, d_type_wr* S, 
    d_type_wr* U, i_type_wr* LDU, d_type_wr* VT, i_type_wr* LDVT, d_type_wr* WORK, i_type_wr* LWORK, i_type_wr* IWORK, 
    i_type_wr* INFO);

int sgesvdx_(const char* JOBU, const char* JOBVT, const char* RANGE, i_type_wr* M, i_type_wr* N, s_type_wr* A, 
    i_type_wr* LDA, s_type_wr* VL, s_type_wr* VU, i_type_wr* IL, i_type_wr* IU, i_type_wr* NS, s_type_wr* S, 
    s_type_wr* U, i_type_wr* LDU, s_type_wr* VT, i_type_wr* LDVT, s_type_wr* WORK, i_type_wr* LWORK, i_type_wr* IWORK, 
    i_type_wr* INFO);

int cgesvdx_(const char* JOBU, const char* JOBVT, const char* RANGE, i_type_wr* M, i_type_wr* N, c_type_wr* A, 
    i_type_wr* LDA, s_type_wr* VL, s_type_wr* VU, i_type_wr* IL, i_type_wr* IU, i_type_wr* NS, s_type_wr* S, 
    c_type_wr* U, i_type_wr* LDU, c_type_wr* VT, i_type_wr* LDVT, c_type_wr* WORK, i_type_wr* LWORK, s_type_wr* RWORK,
    i_type_wr* IWORK, i_type_wr* INFO);

int zgesvdx_(const char* JOBU, const char* JOBVT, const char* RANGE, i_type_wr* M, i_type_wr* N, z_type_wr* A, 
    i_type_wr* LDA, d_type_wr* VL, d_type_wr* VU, i_type_wr* IL, i_type_wr* IU, i_type_wr* NS, d_type_wr* S, 
    z_type_wr* U, i_type_wr* LDU, z_type_wr* VT, i_type_wr* LDVT, z_type_wr* WORK, i_type_wr* LWORK, d_type_wr* RWORK,
    i_type_wr* IWORK, i_type_wr* INFO);

int dbdsvdx_(const char* UPLO, const char* JOBZ, const char* RANGE, i_type_wr* N, const d_type_wr* D, const d_type_wr* E, 
       d_type_wr* VL, d_type_wr* VU, i_type_wr* IL, i_type_wr* IU, i_type_wr* NS, d_type_wr* S, d_type_wr* Z, 
       i_type_wr* LDZ, d_type_wr* WORK, i_type_wr* IWORK, i_type_wr* INFO );

int sbdsvdx_(const char* UPLO, const char* JOBZ, const char* RANGE, i_type_wr* N, const s_type_wr* D, const s_type_wr* E, 
       s_type_wr* VL, s_type_wr* VU, i_type_wr* IL, i_type_wr* IU, i_type_wr* NS, s_type_wr* S, s_type_wr* Z, 
       i_type_wr* LDZ, s_type_wr* WORK, i_type_wr* IWORK, i_type_wr* INFO );

#ifdef __cplusplus
}
#endif

};
