/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2018
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

#include "matcl-blas-lapack/blas_loader/blas_plugin.h"

#ifndef __unix__
    #ifdef BLAS_LOADER_EXPORTS
	    #define BLAS_LOADER_EXPORT __declspec(dllexport)
    #else
        #define BLAS_LOADER_EXPORT __declspec(dllimport)
    #endif
#else
	#define BLAS_LOADER_EXPORT
#endif

// declaration of blas functions which are dynamically loaded on first blas
// use from an available blas implementation
//
// 
// In file blas_loader.cpp a global pointer to a blas_plugin instance
// is defined, which contains pointers to all blas functions from the specific
// blas implementation. Pointer to this specific blas_plugin is retrieved
// from the corresponding plugin dll by its blas_plugin::get_blas_plugin

namespace raw_blas_lapack
{

#ifdef __cplusplus 	
extern "C" {	
#endif

BLAS_LOADER_EXPORT 
i_type_wr get_num_threads_blas_kernel();

BLAS_LOADER_EXPORT 
i_type_wr get_default_threads_blas_kernel();

BLAS_LOADER_EXPORT 
void      set_num_threads_blas_kernel(i_type_wr n);

BLAS_LOADER_EXPORT 
bool      are_user_threads_allowed();

BLAS_LOADER_EXPORT 
i_type_wr caxpy_(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx,
                c_type_wr *cy, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr ccopy_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                 i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr cdotc_(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, 
                 c_type_wr *cy, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr cdotu_(c_type_wr *ret_val, i_type_wr *n, c_type_wr *cx, i_type_wr *incx, 
                 c_type_wr *cy, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr cgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                 c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *x, 
                 i_type_wr *incx, c_type_wr *beta, c_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr cgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                 c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, 
                 i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr cgemv_(char *trans, i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *a,
                 i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta, 
                 c_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr cgerc_(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, 
                 i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *a, 
                 i_type_wr *lda);

BLAS_LOADER_EXPORT 
i_type_wr cgeru_(i_type_wr *m, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, 
                 i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *a, 
                 i_type_wr *lda);

BLAS_LOADER_EXPORT 
i_type_wr chbmv_(char *uplo, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, 
                 c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx, 
                 c_type_wr *beta, c_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr chemm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, c_type_wr *alpha,
                 c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb,
                 c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr chemv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *a, 
                 i_type_wr *lda, c_type_wr *x, i_type_wr *incx, c_type_wr *beta,
                 c_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr cher_(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, 
                i_type_wr *incx, c_type_wr *a, i_type_wr *lda);

BLAS_LOADER_EXPORT 
i_type_wr cher2_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, 
                 i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *a,
                 i_type_wr *lda);

BLAS_LOADER_EXPORT 
i_type_wr cher2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha, 
                  c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
                  s_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr cherk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                 c_type_wr *a, i_type_wr *lda, s_type_wr *beta, c_type_wr *c__, 
                 i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr chpmv_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *ap, 
                 c_type_wr *x, i_type_wr *incx, c_type_wr *beta, c_type_wr *y,
                 i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr chpr_(char *uplo, i_type_wr *n, s_type_wr *alpha, c_type_wr *x, 
                i_type_wr *incx, c_type_wr *ap);

BLAS_LOADER_EXPORT 
i_type_wr chpr2_(char *uplo, i_type_wr *n, c_type_wr *alpha, c_type_wr *x, 
                 i_type_wr *incx, c_type_wr *y, i_type_wr *incy, c_type_wr *ap);

BLAS_LOADER_EXPORT 
i_type_wr crotg_(c_type_wr *ca, c_type_wr *cb, s_type_wr *c__, c_type_wr *s);

BLAS_LOADER_EXPORT 
i_type_wr cscal_(i_type_wr *n, c_type_wr *ca, c_type_wr *cx, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr csrot_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy, 
                 i_type_wr *incy, s_type_wr *c__, s_type_wr *s);

BLAS_LOADER_EXPORT 
i_type_wr csscal_(i_type_wr *n, s_type_wr *sa, c_type_wr *cx, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr cswap_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx, c_type_wr *cy,
                 i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr csymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, 
                 c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, c_type_wr *b, 
                 i_type_wr *ldb, c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr csyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                  c_type_wr *a, i_type_wr *lda, c_type_wr *b, i_type_wr *ldb, 
                  c_type_wr *beta, c_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr csyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, c_type_wr *alpha,
                 c_type_wr *a, i_type_wr *lda, c_type_wr *beta, c_type_wr *c__, 
                 i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr ctbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, 
                 c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr ctbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                 c_type_wr *a, i_type_wr *lda, c_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr ctpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap,
                 c_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr ctpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *ap,
                 c_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr ctrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                 i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, 
                 c_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT 
i_type_wr ctrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a,
                 i_type_wr *lda, c_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr ctrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, 
                 i_type_wr *n, c_type_wr *alpha, c_type_wr *a, i_type_wr *lda, 
                 c_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT 
i_type_wr ctrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, c_type_wr *a,
                 i_type_wr *lda, c_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
d_type_wr dasum_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr daxpy_(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx, 
                 d_type_wr *dy, i_type_wr *incy);

BLAS_LOADER_EXPORT 
d_type_wr dcabs1_(z_type_wr *z__);

BLAS_LOADER_EXPORT 
i_type_wr dcopy_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, 
                 i_type_wr *incy);

BLAS_LOADER_EXPORT 
d_type_wr ddot_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, 
                i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr dgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, 
                 i_type_wr *ku, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, 
                 d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, 
                 i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr dgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                 d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, 
                 i_type_wr *ldb, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr dgemv_(char *trans, i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *a,
                 i_type_wr *lda, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, 
                 d_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr dger_(i_type_wr *m, i_type_wr *n, d_type_wr *alpha, d_type_wr *x,
                i_type_wr *incx, d_type_wr *y, i_type_wr *incy, d_type_wr *a, 
                i_type_wr *lda);

BLAS_LOADER_EXPORT 
d_type_wr dnrm2_(i_type_wr *n, d_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr drot_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                i_type_wr *incy, d_type_wr *c__, d_type_wr *s);

BLAS_LOADER_EXPORT 
i_type_wr drotg_(d_type_wr *da, d_type_wr *db, d_type_wr *c__, d_type_wr *s);

BLAS_LOADER_EXPORT 
i_type_wr drotm_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy, 
                 i_type_wr *incy, d_type_wr *dparam);

BLAS_LOADER_EXPORT 
i_type_wr drotmg_(d_type_wr *dd1, d_type_wr *dd2, d_type_wr *dx1, d_type_wr *dy1,
                  d_type_wr *dparam);

BLAS_LOADER_EXPORT 
i_type_wr dsbmv_(char *uplo, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, d_type_wr *a, 
                 i_type_wr *lda, d_type_wr *x, i_type_wr *incx, d_type_wr *beta, 
                 d_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr dscal_(i_type_wr *n, d_type_wr *da, d_type_wr *dx, i_type_wr *incx);

BLAS_LOADER_EXPORT
d_type_wr dsdot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                 i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr dspmv_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *ap, 
                 d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y,
                 i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr dspr_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, 
                i_type_wr *incx, d_type_wr *ap);

BLAS_LOADER_EXPORT 
i_type_wr dspr2_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, 
                 i_type_wr *incx, d_type_wr *y, i_type_wr *incy, d_type_wr *ap);

BLAS_LOADER_EXPORT 
i_type_wr dswap_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx, d_type_wr *dy,
                 i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr dsymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, d_type_wr *alpha,
                 d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                 d_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr dsymv_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, 
                 d_type_wr *x, i_type_wr *incx, d_type_wr *beta, d_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT 
i_type_wr dsyr_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                d_type_wr *a, i_type_wr *lda);

BLAS_LOADER_EXPORT 
i_type_wr dsyr2_(char *uplo, i_type_wr *n, d_type_wr *alpha, d_type_wr *x, i_type_wr *incx,
                 d_type_wr *y, i_type_wr *incy, d_type_wr *a, i_type_wr *lda);

BLAS_LOADER_EXPORT 
i_type_wr dsyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha, 
                  d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb, d_type_wr *beta,
                  d_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr dsyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                 d_type_wr *a, i_type_wr *lda, d_type_wr *beta, d_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr dtbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                 i_type_wr *lda, d_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr dtbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k, d_type_wr *a,
                 i_type_wr *lda, d_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr dtpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap, d_type_wr *x,
                 i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr dtpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *ap, d_type_wr *x,
                 i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr dtrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                 d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT 
i_type_wr dtrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
                 d_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr dtrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, i_type_wr *n,
                 d_type_wr *alpha, d_type_wr *a, i_type_wr *lda, d_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT 
i_type_wr dtrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, d_type_wr *a, i_type_wr *lda,
                 d_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
d_type_wr dzasum_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx);

BLAS_LOADER_EXPORT
d_type_wr dznrm2_(i_type_wr *n, z_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr icamax_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr idamax_(i_type_wr *n, d_type_wr *dx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr isamax_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr izamax_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr lsame_(char *ca, char *cb);

BLAS_LOADER_EXPORT
d_type_wr sasum_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr saxpy_(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx, 
                 s_type_wr *sy, i_type_wr *incy);

BLAS_LOADER_EXPORT
d_type_wr scabs1_(c_type_wr *z__);

BLAS_LOADER_EXPORT
d_type_wr scasum_(i_type_wr *n, c_type_wr *cx, i_type_wr *incx);

BLAS_LOADER_EXPORT
d_type_wr scnrm2_(i_type_wr *n, c_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr scopy_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                 i_type_wr *incy);

BLAS_LOADER_EXPORT
s_type_wr sdot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                i_type_wr *incy);

BLAS_LOADER_EXPORT
d_type_wr sdsdot_(i_type_wr *n, s_type_wr *sb, s_type_wr *sx, i_type_wr *incx, 
                  s_type_wr *sy, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr sgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                 s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *x, 
                 i_type_wr *incx, s_type_wr *beta, s_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr sgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                 s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, s_type_wr *b, 
                 i_type_wr *ldb, s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr sgemv_(char *trans, i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *a,
                 i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta, 
                 s_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr sger_(i_type_wr *m, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, 
                i_type_wr *incx, s_type_wr *y, i_type_wr *incy, s_type_wr *a,
                i_type_wr *lda);

BLAS_LOADER_EXPORT
d_type_wr snrm2_(i_type_wr *n, s_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr srot_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                i_type_wr *incy, s_type_wr *c__, s_type_wr *s);

BLAS_LOADER_EXPORT
i_type_wr srotg_(s_type_wr *sa, s_type_wr *sb, s_type_wr *c__, s_type_wr *s);

BLAS_LOADER_EXPORT
i_type_wr srotm_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                 i_type_wr *incy, s_type_wr *sparam);

BLAS_LOADER_EXPORT 
i_type_wr srotmg_(s_type_wr *sd1, s_type_wr *sd2, s_type_wr *sx1, s_type_wr *sy1,
                  s_type_wr *sparam);

BLAS_LOADER_EXPORT
i_type_wr ssbmv_(char *uplo, i_type_wr *n, i_type_wr *k, s_type_wr *alpha, 
                 s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx, 
                 s_type_wr *beta, s_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr sscal_(i_type_wr *n, s_type_wr *sa, s_type_wr *sx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr sspmv_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *ap,
                 s_type_wr *x, i_type_wr *incx, s_type_wr *beta, s_type_wr *y,
                 i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr sspr_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x,
                i_type_wr *incx, s_type_wr *ap);

BLAS_LOADER_EXPORT 
i_type_wr sspr2_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, 
                 i_type_wr *incx, s_type_wr *y, i_type_wr *incy, s_type_wr *ap);

BLAS_LOADER_EXPORT
i_type_wr sswap_(i_type_wr *n, s_type_wr *sx, i_type_wr *incx, s_type_wr *sy, 
                 i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr ssymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, s_type_wr *alpha,
                 s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, 
                 s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr ssymv_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *a, 
                 i_type_wr *lda, s_type_wr *x, i_type_wr *incx, s_type_wr *beta,
                 s_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr ssyr_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x, 
                i_type_wr *incx, s_type_wr *a, i_type_wr *lda);

BLAS_LOADER_EXPORT
i_type_wr ssyr2_(char *uplo, i_type_wr *n, s_type_wr *alpha, s_type_wr *x,
                 i_type_wr *incx, s_type_wr *y, i_type_wr *incy, s_type_wr *a,
                 i_type_wr *lda);

BLAS_LOADER_EXPORT
i_type_wr ssyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                  s_type_wr *a, i_type_wr *lda, s_type_wr *b, i_type_wr *ldb, 
                  s_type_wr *beta, s_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr ssyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, s_type_wr *alpha,
                 s_type_wr *a, i_type_wr *lda, s_type_wr *beta, s_type_wr *c__, 
                 i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr stbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                 s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr stbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                 s_type_wr *a, i_type_wr *lda, s_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr stpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap,
                 s_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr stpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *ap,
                 s_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT 
i_type_wr strmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, 
                 i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, 
                 s_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT 
i_type_wr strmv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a, 
                 i_type_wr *lda, s_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr strsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m,
                 i_type_wr *n, s_type_wr *alpha, s_type_wr *a, i_type_wr *lda, 
                 s_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT
i_type_wr strsv_(char *uplo, char *trans, char *diag, i_type_wr *n, s_type_wr *a,
                 i_type_wr *lda, s_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr xerbla_(char *srname, i_type_wr *info);

BLAS_LOADER_EXPORT
i_type_wr zaxpy_(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx,
                 z_type_wr *zy, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zcopy_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy,
                 i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zdotc_(z_type_wr *ret_val, i_type_wr *n, z_type_wr *zx, i_type_wr *incx, 
                 z_type_wr *zy, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zdotu_(z_type_wr *ret_val, i_type_wr *n, z_type_wr *zx, i_type_wr *incx,
                 z_type_wr *zy, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zdrot_(i_type_wr *n, z_type_wr *cx, i_type_wr *incx, z_type_wr *cy, 
                 i_type_wr *incy, d_type_wr *c__, d_type_wr *s);

BLAS_LOADER_EXPORT
i_type_wr zdscal_(i_type_wr *n, d_type_wr *da, z_type_wr *zx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr zgbmv_(char *trans, i_type_wr *m, i_type_wr *n, i_type_wr *kl, i_type_wr *ku,
                 z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *x, 
                 i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zgemm_(char *transa, char *transb, i_type_wr *m, i_type_wr *n, i_type_wr *k,
                 z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, z_type_wr *b, 
                 i_type_wr *ldb, z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr zgemv_(char *trans, i_type_wr *m, i_type_wr *n, z_type_wr *alpha, 
                 z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx,
                 z_type_wr *beta, z_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zgerc_(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, 
                 i_type_wr *incx, z_type_wr *y, i_type_wr *incy, z_type_wr *a, 
                 i_type_wr *lda);

BLAS_LOADER_EXPORT
i_type_wr zgeru_(i_type_wr *m, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, 
                 i_type_wr *incx, z_type_wr *y, i_type_wr *incy, z_type_wr *a, 
                 i_type_wr *lda);

BLAS_LOADER_EXPORT
i_type_wr zhbmv_(char *uplo, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, 
                 z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx,
                 z_type_wr *beta, z_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zhemm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                 z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
                 z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr zhemv_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *a, 
                 i_type_wr *lda, z_type_wr *x, i_type_wr *incx, z_type_wr *beta,
                 z_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zher_(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, 
                i_type_wr *incx, z_type_wr *a, i_type_wr *lda);

BLAS_LOADER_EXPORT
i_type_wr zher2_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, 
                 i_type_wr *incx, z_type_wr *y, i_type_wr *incy, z_type_wr *a, 
                 i_type_wr *lda);

BLAS_LOADER_EXPORT 
i_type_wr zher2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, 
                  z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
                  d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT 
i_type_wr zherk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, d_type_wr *alpha,
                 z_type_wr *a, i_type_wr *lda, d_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr zhpmv_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *ap, z_type_wr *x,
                 i_type_wr *incx, z_type_wr *beta, z_type_wr *y, i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zhpr_(char *uplo, i_type_wr *n, d_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                z_type_wr *ap);

BLAS_LOADER_EXPORT
i_type_wr zhpr2_(char *uplo, i_type_wr *n, z_type_wr *alpha, z_type_wr *x, i_type_wr *incx,
                 z_type_wr *y, i_type_wr *incy, z_type_wr *ap);

BLAS_LOADER_EXPORT
i_type_wr zrotg_(z_type_wr *ca, z_type_wr *cb, d_type_wr *c__, z_type_wr *s);

BLAS_LOADER_EXPORT
i_type_wr zscal_(i_type_wr *n, z_type_wr *za, z_type_wr *zx, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr zswap_(i_type_wr *n, z_type_wr *zx, i_type_wr *incx, z_type_wr *zy, 
                 i_type_wr *incy);

BLAS_LOADER_EXPORT
i_type_wr zsymm_(char *side, char *uplo, i_type_wr *m, i_type_wr *n, z_type_wr *alpha,
                 z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
                 z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr zsyr2k_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha, 
                  z_type_wr *a, i_type_wr *lda, z_type_wr *b, i_type_wr *ldb, 
                  z_type_wr *beta, z_type_wr *c__, i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr zsyrk_(char *uplo, char *trans, i_type_wr *n, i_type_wr *k, z_type_wr *alpha,
                 z_type_wr *a, i_type_wr *lda, z_type_wr *beta, z_type_wr *c__, 
                 i_type_wr *ldc);

BLAS_LOADER_EXPORT
i_type_wr ztbmv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                 z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr ztbsv_(char *uplo, char *trans, char *diag, i_type_wr *n, i_type_wr *k,
                 z_type_wr *a, i_type_wr *lda, z_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr ztpmv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap,
                 z_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr ztpsv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *ap,
                 z_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr ztrmm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, 
                 i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, 
                 z_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT
i_type_wr ztrmv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a,
                 i_type_wr *lda, z_type_wr *x, i_type_wr *incx);

BLAS_LOADER_EXPORT
i_type_wr ztrsm_(char *side, char *uplo, char *transa, char *diag, i_type_wr *m, 
                 i_type_wr *n, z_type_wr *alpha, z_type_wr *a, i_type_wr *lda, 
                 z_type_wr *b, i_type_wr *ldb);

BLAS_LOADER_EXPORT
i_type_wr ztrsv_(char *uplo, char *trans, char *diag, i_type_wr *n, z_type_wr *a, 
                 i_type_wr *lda, z_type_wr *x, i_type_wr *incx);

#ifdef __cplusplus
}
#endif

};
