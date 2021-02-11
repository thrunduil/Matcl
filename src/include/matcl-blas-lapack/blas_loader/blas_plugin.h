/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 *
 *  This program is distributed in the hope that it will be useful,
 *  but WITHOUT ANY WARRANTY; without even the implied warranty of
 *  MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 *  GNU General Public License for more details.
 *
 *  You should have received a copy of the GNU General Public License
 *  along with this program; if not, write to the Free Software
 *  Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#pragma once

#ifndef __unix__
	#define BLAS_PLUGIN_EXPORT __declspec(dllexport)
#else
	#define BLAS_PLUGIN_EXPORT
#endif

#include "matcl-blas-lapack/blas_loader/blas_types.h"

// Structure holding pointers to specific implementations of all
// blas functions
struct blas_plugin
{
    // Constructor which initializes pointers to blas functions implementations

    blas_plugin();    

    using get_blas_plugin_type          = blas_plugin* ();
    using get_blas_plugin_gpu_cpu_type  = blas_plugin* (const ::blas_plugin* plugin_cpu, 
                                            const ::blas_plugin* plugin_gpu);
    
    using get_num_threads_type      = i_type_wr (*)();
    using get_default_threads_type  = i_type_wr (*)();
    using set_num_threads_type      = void      (*)(i_type_wr*);
    using user_threads_allowed_type = bool      (*)();
    using get_name_type             = const char* (*)();
    using initialize_type           = void      (*)();

    using caxpy_func_type = i_type_wr (*)(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, 
                                          i_type_wr *incx, c_type_wr *cy, i_type_wr *incy);
    using ccopy_func_type = i_type_wr (*)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, 
                                          c_type_wr *cy, i_type_wr *incy);
    using cdotc_func_type = c_type_wr (*)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                                          i_type_wr *incy);
    using cdotu_func_type = c_type_wr (*)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                                          i_type_wr *incy);
    using cgbmv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl,
                                          i_type_wr *ku, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, 
                                          c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y,
                                          i_type_wr *incy);
    using cgemm_func_type = i_type_wr (*)(char *transa, char *transb, i_type_wr *m, i_type_wr *n, 
                                          i_type_wr *k, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, 
                                          c_type_wr *b, i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, 
                                          i_type_wr *ldc);
    using cgemv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                                          c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, 
                                          c_type_wr *beta, c_type_wr *y, i_type_wr *incy);
    using cgerc_func_type = i_type_wr (*)(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, 
                                          i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *a, 
                                          i_type_wr *lda);
    using cgeru_func_type = i_type_wr (*)(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x,
                                          i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *a,
                                          i_type_wr *lda);
    using chbmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, 
                                          c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx,
                                          c_type_wr *beta, c_type_wr *y, i_type_wr *incy);
    using chemm_func_type = i_type_wr (*)(char *side, char *uplo, i_type_wr *m, i_type_wr *n, 
                                          c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b,
                                          i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);
    using chemv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, 
                                          i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta,
                                          c_type_wr *y, i_type_wr *incy);
    using cher_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x,
                                          i_type_wr *incx, c_type_wr *a, i_type_wr *lda);
    using cher2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x,
                                          i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *a,
                                          i_type_wr *lda);
    using cher2k_func_type= i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b,
                                          i_type_wr *ldb, s_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);
    using cherk_func_type = i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k,
                                          s_type_wr *alpha, c_type_wr *a, i_type_wr *lda, s_type_wr *beta, 
                                          c_type_wr *c__, i_type_wr *ldc);
    using chpmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *ap,
                                          c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y,
                                          i_type_wr *incy);
    using chpr_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x,
                                          i_type_wr *incx, c_type_wr *ap);
    using chpr2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x,
                                          i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *ap);
    using crotg_func_type = i_type_wr (*)(c_type_wr *ca, c_type_wr *cb, s_type_wr *c__, c_type_wr *s);
    using cscal_func_type = i_type_wr (*)(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx);
    using csrot_func_type = i_type_wr (*)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                                          i_type_wr *incy, s_type_wr *c__, s_type_wr *s);
    using csscal_func_type= i_type_wr (*)(i_type_wr *n, s_type_wr *sa, c_type_wr *cx, i_type_wr *incx);
    using cswap_func_type = i_type_wr (*)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                                          i_type_wr *incy);
    using csymm_func_type = i_type_wr (*)(char *side, char *uplo, i_type_wr *m, i_type_wr *n, 
                                          c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b,
                                          i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);
    using csyr2k_func_type= i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k,
                                          c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, 
                                          i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);
    using csyrk_func_type = i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *beta, 
                                          c_type_wr *c__, i_type_wr *ldc);
    using ctbmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, 
                                          c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx);
    using ctbsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                                          c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx);
    using ctpmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap,
                                          c_type_wr *x, i_type_wr *incx);
    using ctpsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap,
                                          c_type_wr *x, i_type_wr *incx);
    using ctrmm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, 
                                          i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda,
                                          c_type_wr *b, i_type_wr *ldb);
    using ctrmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a,
                                          i_type_wr *lda, c_type_wr *x, i_type_wr *incx);
    using ctrsm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                                          i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda,
                                          c_type_wr *b, i_type_wr *ldb);
    using ctrsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a,
                                          i_type_wr *lda, c_type_wr *x, i_type_wr *incx);
    using dasum_func_type = d_type_wr (*)(i_type_wr *n, d_type_wr *dx, i_type_wr *incx);
    using daxpy_func_type = i_type_wr (*)(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx, 
                                          d_type_wr *dy, i_type_wr *incy);
    using dcabs1_func_type= d_type_wr (*)(z_type_wr *z__);
    using dcopy_func_type = i_type_wr (*)(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                                          i_type_wr *incy);
    using ddot_func_type  = d_type_wr (*)(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                                          i_type_wr *incy);
    using dgbmv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, 
                                          i_type_wr *ku, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda,
                                          d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, 
                                          i_type_wr *incy);
    using dgemm_func_type = i_type_wr (*)(char *transa, char *transb, i_type_wr *m, i_type_wr *n, 
                                          i_type_wr *k, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, 
                                          d_type_wr *b, i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, 
                                          i_type_wr *ldc);
    using dgemv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, 
                                          d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx,
                                          d_type_wr *beta, d_type_wr *y, i_type_wr *incy);
    using dger_func_type  = i_type_wr (*)(i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *x,
                                          i_type_wr *incx, d_type_wr *y, i_type_wr *incy, d_type_wr *a,
                                          i_type_wr *lda);
    using dnrm2_func_type = d_type_wr (*)(i_type_wr *n, d_type_wr *x, i_type_wr *incx);
    using drot_func_type  = i_type_wr (*)(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                                          i_type_wr *incy, d_type_wr *c__, d_type_wr *s);
    using drotg_func_type = i_type_wr (*)(d_type_wr *da, d_type_wr *db, d_type_wr *c__, d_type_wr *s);
    using drotm_func_type = i_type_wr (*)(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                                          i_type_wr *incy, d_type_wr *dparam);
    using drotmg_func_type= i_type_wr (*)(d_type_wr *dd1, d_type_wr *dd2, d_type_wr *dx1, d_type_wr *dy1,
                                          d_type_wr *dparam);
    using dsbmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, 
                                          d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx,
                                          d_type_wr *beta, d_type_wr *y, i_type_wr *incy);
    using dscal_func_type = i_type_wr (*)(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx);
    using dsdot_func_type = d_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                                          i_type_wr *incy);
    using dspmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *ap, 
                                          d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y,
                                          i_type_wr *incy);
    using dspr_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x,
                                          i_type_wr *incx, d_type_wr *ap);
    using dspr2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, 
                                          i_type_wr *incx, d_type_wr *y, i_type_wr *incy, d_type_wr *ap);
    using dswap_func_type = i_type_wr (*)(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                                          i_type_wr *incy);
    using dsymm_func_type = i_type_wr (*)(char *side, char *uplo, i_type_wr *m, i_type_wr *n, 
                                          d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b,
                                          i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc);
    using dsymv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *a,
                                          i_type_wr *lda, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, 
                                          d_type_wr *y, i_type_wr *incy);
    using dsyr_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x,
                                          i_type_wr *incx, d_type_wr *a, i_type_wr *lda);
    using dsyr2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, 
                                          i_type_wr *incx, d_type_wr *y, i_type_wr *incy, d_type_wr *a,
                                          i_type_wr *lda);
    using dsyr2k_func_type= i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
                                          i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc);
    using dsyrk_func_type = i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *beta,
                                          d_type_wr *c__, i_type_wr *ldc);
    using dtbmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                                          d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx);
    using dtbsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                                          d_type_wr *a, i_type_wr *lda, d_type_wr *x, i_type_wr *incx);
    using dtpmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap,
                                          d_type_wr *x, i_type_wr *incx);
    using dtpsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap, 
                                          d_type_wr *x, i_type_wr *incx);
    using dtrmm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                                          i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda,
                                          d_type_wr *b, i_type_wr *ldb);
    using dtrmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, 
                                          i_type_wr *lda, d_type_wr *x, i_type_wr *incx);
    using dtrsm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                                          i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, 
                                          d_type_wr *b, i_type_wr *ldb);
    using dtrsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a,
                                          i_type_wr *lda, d_type_wr *x, i_type_wr *incx);
    using dzasum_func_type= d_type_wr (*)(i_type_wr *n, z_type_wr *zx, i_type_wr *incx);
    using dznrm2_func_type= d_type_wr (*)(i_type_wr *n, z_type_wr *x, i_type_wr *incx);
    using icamax_func_type= i_type_wr (*)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx);
    using idamax_func_type= i_type_wr (*)(i_type_wr *n, d_type_wr *dx, i_type_wr *incx);
    using isamax_func_type= i_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx);
    using izamax_func_type= i_type_wr (*)(i_type_wr *n, z_type_wr *zx, i_type_wr *incx);
    using lsame_func_type = i_type_wr (*)(char *ca, char *cb);
    using sasum_func_type = s_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx);
    using saxpy_func_type = i_type_wr (*)(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx, 
                                          s_type_wr *sy, i_type_wr *incy);
    using scabs1_func_type= s_type_wr (*)(c_type_wr *z__);
    using scasum_func_type= s_type_wr (*)(i_type_wr *n, c_type_wr *cx, i_type_wr *incx);
    using scnrm2_func_type= s_type_wr (*)(i_type_wr *n, c_type_wr *x, i_type_wr *incx);
    using scopy_func_type = i_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                                          i_type_wr *incy);
    using sdot_func_type  = s_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                                          i_type_wr *incy);
    using sdsdot_func_type= s_type_wr (*)(i_type_wr *n, s_type_wr *sb, s_type_wr *sx, i_type_wr *incx,
                                          s_type_wr *sy, i_type_wr *incy);
    using sgbmv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl,
                                          i_type_wr *ku, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda,
                                          s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y, 
                                          i_type_wr *incy);
    using sgemm_func_type = i_type_wr (*)(char *transa, char *transb, i_type_wr *m, i_type_wr *n, 
                                          i_type_wr *k, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda,
                                          s_type_wr *b, i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, 
                                          i_type_wr *ldc);
    using sgemv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, 
                                          s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx, 
                                          s_type_wr *beta, s_type_wr *y, i_type_wr *incy);
    using sger_func_type  = i_type_wr (*)(i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, 
                                          i_type_wr *incx, s_type_wr *y, i_type_wr *incy, s_type_wr *a,
                                          i_type_wr *lda);
    using snrm2_func_type = s_type_wr (*)(i_type_wr *n, s_type_wr *x, i_type_wr *incx);
    using srot_func_type  = i_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy,
                                          i_type_wr *incy, s_type_wr *c__, s_type_wr *s);
    using srotg_func_type = i_type_wr (*)(s_type_wr *sa, s_type_wr *sb, s_type_wr *c__, s_type_wr *s);
    using srotm_func_type = i_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy,
                                          i_type_wr *incy, s_type_wr *sparam);
    using srotmg_func_type= i_type_wr (*)(s_type_wr *sd1, s_type_wr *sd2, s_type_wr *sx1, s_type_wr *sy1,
                                          s_type_wr *sparam);
    using ssbmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, 
                                          s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx,
                                          s_type_wr *beta, s_type_wr *y, i_type_wr *incy);
    using sscal_func_type = i_type_wr (*)(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx);
    using sspmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *ap, 
                                          s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                                          i_type_wr *incy);
    using sspr_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x,
                                          i_type_wr *incx, s_type_wr *ap);
    using sspr2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x,
                                          i_type_wr *incx, s_type_wr *y, i_type_wr *incy, s_type_wr *ap);
    using sswap_func_type = i_type_wr (*)(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy,
                                          i_type_wr *incy);
    using ssymm_func_type = i_type_wr (*)(char *side, char *uplo, i_type_wr *m, i_type_wr *n, 
                                          s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b,
                                          i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc);
    using ssymv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *a,
                                          i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta,
                                          s_type_wr *y, i_type_wr *incy);
    using ssyr_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, 
                                          i_type_wr *incx, s_type_wr *a, i_type_wr *lda);
    using ssyr2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, 
                                          i_type_wr *incx, s_type_wr *y, i_type_wr *incy, s_type_wr *a,
                                          i_type_wr *lda);
    using ssyr2k_func_type= i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, 
                                          i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc);
    using ssyrk_func_type = i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *beta,
                                          s_type_wr *c__, i_type_wr *ldc);
    using stbmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                                          s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx);
    using stbsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                                          s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx);
    using stpmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap,
                                          s_type_wr *x, i_type_wr *incx);
    using stpsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap, 
                                          s_type_wr *x, i_type_wr *incx);
    using strmm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                                          i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, 
                                          s_type_wr *b, i_type_wr *ldb);
    using strmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a,
                                          i_type_wr *lda, s_type_wr *x, i_type_wr *incx);
    using strsm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                                          i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda,
                                          s_type_wr *b, i_type_wr *ldb);
    using strsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a,
                                          i_type_wr *lda, s_type_wr *x, i_type_wr *incx);
    using xerbla_func_type= i_type_wr (*)(char *srname, i_type_wr *info);
    using zaxpy_func_type = i_type_wr (*)(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx,
                                          z_type_wr *zy, i_type_wr *incy);
    using zcopy_func_type = i_type_wr (*)(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                                          i_type_wr *incy);
    using zdotc_func_type = z_type_wr (*)(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                                          i_type_wr *incy);
    using zdotu_func_type = z_type_wr (*)(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                                          i_type_wr *incy);
    using zdrot_func_type = i_type_wr (*)(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy,
                                          i_type_wr *incy, d_type_wr *c__, d_type_wr *s);
    using zdscal_func_type= i_type_wr (*)(i_type_wr *n, d_type_wr *da, z_type_wr *zx, i_type_wr *incx);
    using zgbmv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, 
                                          i_type_wr *ku, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, 
                                          z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, 
                                          i_type_wr *incy);
    using zgemm_func_type = i_type_wr (*)(char *transa, char *transb, i_type_wr *m, i_type_wr *n,
                                          i_type_wr *k, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda,
                                          z_type_wr *b, i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, 
                                          i_type_wr *ldc);
    using zgemv_func_type = i_type_wr (*)(char *trans, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, 
                                          z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx, 
                                          z_type_wr *beta, z_type_wr *y, i_type_wr *incy);
    using zgerc_func_type = i_type_wr (*)(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, 
                                          i_type_wr *incx, z_type_wr *y, i_type_wr *incy, z_type_wr *a,
                                          i_type_wr *lda);
    using zgeru_func_type = i_type_wr (*)(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x,
                                          i_type_wr *incx, z_type_wr *y, i_type_wr *incy, z_type_wr *a, 
                                          i_type_wr *lda);
    using zhbmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, 
                                          z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx,
                                          z_type_wr *beta, z_type_wr *y, i_type_wr *incy);
    using zhemm_func_type = i_type_wr (*)(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                                          z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
                                          z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);
    using zhemv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, 
                                          i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta,
                                          z_type_wr *y, i_type_wr *incy);
    using zher_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x,
                                          i_type_wr *incx, z_type_wr *a, i_type_wr *lda);
    using zher2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, 
                                          i_type_wr *incx, z_type_wr *y, i_type_wr *incy, z_type_wr *a, 
                                          i_type_wr *lda);
    using zher2k_func_type= i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b,
                                          i_type_wr *ldb, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);
    using zherk_func_type = i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k,
                                          d_type_wr *alpha, z_type_wr *a, i_type_wr *lda, d_type_wr *beta,
                                          z_type_wr *c__, i_type_wr *ldc);
    using zhpmv_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *ap, 
                                          z_type_wr *x, i_type_wr *incx, z_type_wr *beta, z_type_wr *y, 
                                          i_type_wr *incy);
    using zhpr_func_type  = i_type_wr (*)(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, 
                                          i_type_wr *incx, z_type_wr *ap);
    using zhpr2_func_type = i_type_wr (*)(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, 
                                          i_type_wr *incx, z_type_wr *y, i_type_wr *incy, z_type_wr *ap);
    using zrotg_func_type = i_type_wr (*)(z_type_wr *ca, z_type_wr *cb, d_type_wr *c__, z_type_wr *s);
    using zscal_func_type = i_type_wr (*)(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx);
    using zswap_func_type = i_type_wr (*)(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                                          i_type_wr *incy);
    using zsymm_func_type = i_type_wr (*)(char *side, char *uplo, i_type_wr *m, i_type_wr *n,
                                          z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
                                          i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);
    using zsyr2k_func_type= i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b,
                                          i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);
    using zsyrk_func_type = i_type_wr (*)(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, 
                                          z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *beta,
                                          z_type_wr *c__, i_type_wr *ldc);
    using ztbmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                                          z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx);
    using ztbsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                                          z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx);
    using ztpmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap,
                                          z_type_wr *x, i_type_wr *incx);
    using ztpsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap,
                                          z_type_wr *x, i_type_wr *incx);
    using ztrmm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                                          i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda,
                                          z_type_wr *b, i_type_wr *ldb);
    using ztrmv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a,
                                          i_type_wr *lda, z_type_wr *x, i_type_wr *incx);
    using ztrsm_func_type = i_type_wr (*)(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, 
                                          i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, 
                                          z_type_wr *b, i_type_wr *ldb);
    using ztrsv_func_type = i_type_wr (*)(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a,
                                          i_type_wr *lda, z_type_wr *x, i_type_wr *incx);

    // return name of this plugin
    const get_name_type             get_name_fptr;

    // perform initialization; called after loading of this plugin
    const initialize_type           initialize_fptr;

    // perform initialization; called when blas::initialize_plugin() function
    // is called
    const initialize_type           force_initialization_fptr;

    const get_num_threads_type      get_num_threads_fptr;
    const get_default_threads_type  get_default_threads_fptr;
    const set_num_threads_type      set_num_threads_fptr;
    const user_threads_allowed_type are_user_threads_allowed_fptr;

    const caxpy_func_type   caxpy_fptr;
    const ccopy_func_type   ccopy_fptr;
    const cdotc_func_type   cdotc_fptr;
    const cdotu_func_type   cdotu_fptr;
    const cgbmv_func_type   cgbmv_fptr;
    const cgemm_func_type   cgemm_fptr;
    const cgemv_func_type   cgemv_fptr;
    const cgerc_func_type   cgerc_fptr;
    const cgeru_func_type   cgeru_fptr;
    const chbmv_func_type   chbmv_fptr;
    const chemm_func_type   chemm_fptr;
    const chemv_func_type   chemv_fptr;
    const cher_func_type    cher_fptr;
    const cher2_func_type   cher2_fptr;
    const cher2k_func_type  cher2k_fptr;
    const cherk_func_type   cherk_fptr;
    const chpmv_func_type   chpmv_fptr;
    const chpr_func_type    chpr_fptr;
    const chpr2_func_type   chpr2_fptr;
    const crotg_func_type   crotg_fptr;
    const cscal_func_type   cscal_fptr;
    const csrot_func_type   csrot_fptr;
    const csscal_func_type  csscal_fptr;
    const cswap_func_type   cswap_fptr;
    const csymm_func_type   csymm_fptr;
    const csyr2k_func_type  csyr2k_fptr;
    const csyrk_func_type   csyrk_fptr;
    const ctbmv_func_type   ctbmv_fptr;
    const ctbsv_func_type   ctbsv_fptr;
    const ctpmv_func_type   ctpmv_fptr;
    const ctpsv_func_type   ctpsv_fptr;
    const ctrmm_func_type   ctrmm_fptr;
    const ctrmv_func_type   ctrmv_fptr;
    const ctrsm_func_type   ctrsm_fptr;
    const ctrsv_func_type   ctrsv_fptr;
    const dasum_func_type   dasum_fptr;
    const daxpy_func_type   daxpy_fptr;
    const dcabs1_func_type  dcabs1_fptr;
    const dcopy_func_type   dcopy_fptr;
    const ddot_func_type    ddot_fptr;
    const dgbmv_func_type   dgbmv_fptr;
    const dgemm_func_type   dgemm_fptr;
    const dgemv_func_type   dgemv_fptr;
    const dger_func_type    dger_fptr;
    const dnrm2_func_type   dnrm2_fptr;
    const drot_func_type    drot_fptr;
    const drotg_func_type   drotg_fptr;
    const drotm_func_type   drotm_fptr;
    const drotmg_func_type  drotmg_fptr;
    const dsbmv_func_type   dsbmv_fptr;
    const dscal_func_type   dscal_fptr;
    const dsdot_func_type   dsdot_fptr;
    const dspmv_func_type   dspmv_fptr;
    const dspr_func_type    dspr_fptr;
    const dspr2_func_type   dspr2_fptr;
    const dswap_func_type   dswap_fptr;
    const dsymm_func_type   dsymm_fptr;
    const dsymv_func_type   dsymv_fptr;
    const dsyr_func_type    dsyr_fptr;
    const dsyr2_func_type   dsyr2_fptr;
    const dsyr2k_func_type  dsyr2k_fptr;
    const dsyrk_func_type   dsyrk_fptr;
    const dtbmv_func_type   dtbmv_fptr;
    const dtbsv_func_type   dtbsv_fptr;
    const dtpmv_func_type   dtpmv_fptr;
    const dtpsv_func_type   dtpsv_fptr;
    const dtrmm_func_type   dtrmm_fptr;
    const dtrmv_func_type   dtrmv_fptr;
    const dtrsm_func_type   dtrsm_fptr;
    const dtrsv_func_type   dtrsv_fptr;
    const dzasum_func_type  dzasum_fptr;
    const dznrm2_func_type  dznrm2_fptr;
    const icamax_func_type  icamax_fptr;
    const idamax_func_type  idamax_fptr;
    const isamax_func_type  isamax_fptr;
    const izamax_func_type  izamax_fptr;
    const sasum_func_type   sasum_fptr;
    const saxpy_func_type   saxpy_fptr;
    const scabs1_func_type  scabs1_fptr;
    const scasum_func_type  scasum_fptr;
    const scnrm2_func_type  scnrm2_fptr;
    const scopy_func_type   scopy_fptr;
    const sdot_func_type    sdot_fptr;
    const sdsdot_func_type  sdsdot_fptr;
    const sgbmv_func_type   sgbmv_fptr;
    const sgemm_func_type   sgemm_fptr;
    const sgemv_func_type   sgemv_fptr;
    const sger_func_type    sger_fptr;
    const snrm2_func_type   snrm2_fptr;
    const srot_func_type    srot_fptr;
    const srotg_func_type   srotg_fptr;
    const srotm_func_type   srotm_fptr;
    const srotmg_func_type  srotmg_fptr;
    const ssbmv_func_type   ssbmv_fptr;
    const sscal_func_type   sscal_fptr;
    const sspmv_func_type   sspmv_fptr;
    const sspr_func_type    sspr_fptr;
    const sspr2_func_type   sspr2_fptr;
    const sswap_func_type   sswap_fptr;
    const ssymm_func_type   ssymm_fptr;
    const ssymv_func_type   ssymv_fptr;
    const ssyr_func_type    ssyr_fptr;
    const ssyr2_func_type   ssyr2_fptr;
    const ssyr2k_func_type  ssyr2k_fptr;
    const ssyrk_func_type   ssyrk_fptr;
    const stbmv_func_type   stbmv_fptr;
    const stbsv_func_type   stbsv_fptr;
    const stpmv_func_type   stpmv_fptr;
    const stpsv_func_type   stpsv_fptr;
    const strmm_func_type   strmm_fptr;
    const strmv_func_type   strmv_fptr;
    const strsm_func_type   strsm_fptr;
    const strsv_func_type   strsv_fptr;
    const zaxpy_func_type   zaxpy_fptr;
    const zcopy_func_type   zcopy_fptr;
    const zdotc_func_type   zdotc_fptr;
    const zdotu_func_type   zdotu_fptr;
    const zdrot_func_type   zdrot_fptr;
    const zdscal_func_type  zdscal_fptr;
    const zgbmv_func_type   zgbmv_fptr;
    const zgemm_func_type   zgemm_fptr;
    const zgemv_func_type   zgemv_fptr;
    const zgerc_func_type   zgerc_fptr;
    const zgeru_func_type   zgeru_fptr;
    const zhbmv_func_type   zhbmv_fptr;
    const zhemm_func_type   zhemm_fptr;
    const zhemv_func_type   zhemv_fptr;
    const zher_func_type    zher_fptr;
    const zher2_func_type   zher2_fptr;
    const zher2k_func_type  zher2k_fptr;
    const zherk_func_type   zherk_fptr;
    const zhpmv_func_type   zhpmv_fptr;
    const zhpr_func_type    zhpr_fptr;
    const zhpr2_func_type   zhpr2_fptr;
    const zrotg_func_type   zrotg_fptr;
    const zscal_func_type   zscal_fptr;
    const zswap_func_type   zswap_fptr;
    const zsymm_func_type   zsymm_fptr;
    const zsyr2k_func_type  zsyr2k_fptr;
    const zsyrk_func_type   zsyrk_fptr;
    const ztbmv_func_type   ztbmv_fptr;
    const ztbsv_func_type   ztbsv_fptr;
    const ztpmv_func_type   ztpmv_fptr;
    const ztpsv_func_type   ztpsv_fptr;
    const ztrmm_func_type   ztrmm_fptr;
    const ztrmv_func_type   ztrmv_fptr;
    const ztrsm_func_type   ztrsm_fptr;
    const ztrsv_func_type   ztrsv_fptr;
};

extern "C"
{
    // Function returning a pointer to an implementation-specific
    // instance from the plugin
    //
    // Needs to be implemented in the corresponding plugin file

    BLAS_PLUGIN_EXPORT const blas_plugin* get_blas_plugin();
}
