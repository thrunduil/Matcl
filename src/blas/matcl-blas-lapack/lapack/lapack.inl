/*
 * This file is a part of Matrix Computation Library (MATCL)
 *
 * Copyright (c) Pawe³ Kowal 2017 - 2021
 *
 * This program is free software; you can redistribute it and/or modify
 * it under the terms of the GNU General Public License as published by
 * the Free Software Foundation; either version 2 of the License, or
 * (at your option) any later version.
 *
 * This program is distributed in the hope that it will be useful,
 * but WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
 * GNU General Public License for more details.
 *
 * You should have received a copy of the GNU General Public License
 * along with this program; if not, write to the Free Software
 * Foundation, Inc., 59 Temple Place, Suite 330, Boston, MA 02111-1307 USA
 */

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-lapack/blas/blas.h"
#include "lacpy2.inl"

#include <vector>

#include "matcl-blas-lapack/blas/details/config_blas_lib.h"
#include "matcl-core/utils/workspace.h"

#ifndef INLINE_TYPE
#define INLINE_TYPE
#endif

namespace matcl { namespace lapack
{

//-----------------------------------------------------------------------
//                          ILAENV
//-----------------------------------------------------------------------

BLAS_EXPORT INLINE_TYPE 
i_type ilaenv(i_type ispec, const char *name, const char *opts, 
                           i_type n1, i_type n2, i_type n3, i_type n4)
{
    return LAPACK_NAME(ilaenv)(_rc(&ispec), _rc(name), _rc(opts), _rc(&n1), 
                               _rc(&n2), _rc(&n3), _rc(&n4));
};

//-----------------------------------------------------------------------
//                          laswp
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void laswp<s_type>(i_type n, s_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,
                   i_type incx)
{
    LAPACK_NAME(slaswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void laswp<d_type>(i_type n, d_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,
                   i_type incx)
{
    LAPACK_NAME(dlaswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void laswp<c_type>(i_type n, c_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,
                   i_type incx)
{
    LAPACK_NAME(claswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

template<> BLAS_EXPORT INLINE_TYPE 
void laswp<z_type>(i_type n, z_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,
                   i_type incx)
{
    LAPACK_NAME(zlaswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
void claswp(i_type n, c_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    LAPACK_NAME(claswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
void slaswp(i_type n, s_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    LAPACK_NAME(slaswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
void dlaswp(i_type n, d_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    LAPACK_NAME(dlaswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

BLAS_EXPORT INLINE_TYPE
void zlaswp(i_type n, z_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx)
{
    LAPACK_NAME(zlaswp)(_rc(&n), _rc(a),_rc(&lda),_rc(&k1),_rc(&k2),_rc(ipiv),_rc(&incx));
};

//-----------------------------------------------------------------------
//                          laset
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void laset<s_type>(const char *uplo, i_type m, i_type n, s_type alpha, s_type beta, 
                   s_type* a, i_type lda)
{
    LAPACK_NAME(slaset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void laset<d_type>(const char *uplo, i_type m, i_type n, d_type alpha, d_type beta, d_type* a, 
                   i_type lda)
{
    LAPACK_NAME(dlaset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void laset<c_type>(const char *uplo, i_type m, i_type n, c_type alpha, c_type beta, c_type* a, 
                   i_type lda)
{
    LAPACK_NAME(claset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

template<> BLAS_EXPORT INLINE_TYPE 
void laset<z_type>(const char *uplo, i_type m, i_type n, z_type alpha, z_type beta, z_type* a,
                   i_type lda)
{
    LAPACK_NAME(zlaset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void claset(const char *uplo, i_type m, i_type n, c_type alpha, c_type beta, c_type* a, i_type lda)
{
    LAPACK_NAME(claset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void slaset(const char *uplo, i_type m, i_type n, s_type alpha, s_type beta, s_type* a, i_type lda)
{
    LAPACK_NAME(slaset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void dlaset(const char *uplo, i_type m, i_type n, d_type alpha, d_type beta, d_type* a, i_type lda)
{
    LAPACK_NAME(dlaset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

BLAS_EXPORT INLINE_TYPE
void zlaset(const char *uplo, i_type m, i_type n, z_type alpha, z_type beta, z_type* a, i_type lda)
{
    LAPACK_NAME(zlaset)(_rc(uplo), _rc(&m), _rc(&n), _rc(&alpha), _rc(&beta), _rc(a), _rc(&lda));
};

//-----------------------------------------------------------------------
//                          lartg
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void lartg<s_type>(s_type f, s_type g, s_type *cs, s_type *sn,  s_type *r)
{
    LAPACK_NAME(slartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

template<> BLAS_EXPORT INLINE_TYPE 
void lartg<d_type>(d_type f, d_type g, d_type *cs, d_type *sn,  d_type *r)
{
    LAPACK_NAME(dlartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

template<> BLAS_EXPORT INLINE_TYPE 
void lartg<c_type>(c_type f, c_type g, s_type *cs, c_type *sn,  c_type *r)
{
    LAPACK_NAME(clartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

template<> BLAS_EXPORT INLINE_TYPE 
void lartg<z_type>(z_type f, z_type g, d_type *cs, z_type *sn,  z_type *r)
{
    LAPACK_NAME(zlartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

BLAS_EXPORT INLINE_TYPE
void clartg(c_type f, c_type g, s_type *cs, c_type *sn,  c_type *r)
{
    LAPACK_NAME(clartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

BLAS_EXPORT INLINE_TYPE
void dlartg(d_type f, d_type g, d_type *cs, d_type *sn, d_type *r)
{
    LAPACK_NAME(dlartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

BLAS_EXPORT INLINE_TYPE 
void slartg(s_type f, s_type g, s_type *cs, s_type *sn, s_type *r)
{
    LAPACK_NAME(slartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

BLAS_EXPORT INLINE_TYPE 
void zlartg(z_type f, z_type g, d_type * cs, z_type *sn, z_type *r)
{
    LAPACK_NAME(zlartg)(_rc(&f), _rc(&g), _rc(cs), _rc(sn),  _rc(r));
};

//-----------------------------------------------------------------------
//                          lange
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
s_type lange<s_type>(const char *norm,i_type m,i_type n,const s_type *a,i_type lda,
                     s_type *work)
{
    auto res = LAPACK_NAME(slange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda),
                                       _rc(work));
    return (s_type)res;
};

template<> BLAS_EXPORT INLINE_TYPE 
s_type lange<c_type>(const char *norm,i_type m,i_type n,const c_type *a,i_type lda,
                     s_type *work)
{
    auto res = LAPACK_NAME(clange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda),
                                       _rc(work));
    return (s_type)res;
};

template<> BLAS_EXPORT INLINE_TYPE 
d_type lange<d_type>(const char *norm,i_type m,i_type n,const d_type *a,i_type lda,
                     d_type *work)
{
    auto res = LAPACK_NAME(dlange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda), _rc(work));
    return res;
};

template<> BLAS_EXPORT INLINE_TYPE 
d_type lange<z_type>(const char *norm,i_type m,i_type n,const z_type *a,i_type lda,
                     d_type *work)
{
    auto res = LAPACK_NAME(zlange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(work));
    return res;
};

BLAS_EXPORT INLINE_TYPE
s_type clange(const char *norm,i_type m,i_type n,const c_type *a,i_type lda,s_type *work)
{
    return (s_type)LAPACK_NAME(clange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(work));
};

BLAS_EXPORT INLINE_TYPE
s_type slange(const char *norm,i_type m,i_type n,const s_type *a,i_type lda,s_type *work)
{
    return (s_type)LAPACK_NAME(slange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(work));
};

BLAS_EXPORT INLINE_TYPE
d_type dlange(const char *norm,i_type m,i_type n,const d_type *a,i_type lda,d_type *work)
{
    return LAPACK_NAME(dlange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(work));
};

BLAS_EXPORT INLINE_TYPE
d_type zlange(const char *norm,i_type m,i_type n,const z_type *a,i_type lda,d_type *work)
{
    return LAPACK_NAME(zlange)(_rc(norm),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(work));
};

//-----------------------------------------------------------------------
//                          lamch
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
s_type lamch<s_type>(const char *cmach)
{
    return (s_type)LAPACK_NAME(slamch)(_rc(cmach));
};

template<> BLAS_EXPORT INLINE_TYPE 
d_type lamch<d_type>(const char *cmach)
{
    return LAPACK_NAME(dlamch)(_rc(cmach));
};

template<> BLAS_EXPORT INLINE_TYPE 
s_type lamch<c_type>(const char *cmach)
{
    return (s_type)LAPACK_NAME(slamch)(_rc(cmach));
};

template<> BLAS_EXPORT INLINE_TYPE 
d_type lamch<z_type>(const char *cmach)
{
    return LAPACK_NAME(dlamch)(_rc(cmach));
};

BLAS_EXPORT INLINE_TYPE
s_type slamch(const char *cmach)
{
    return (s_type)LAPACK_NAME(slamch)(_rc(cmach));
};

BLAS_EXPORT INLINE_TYPE
d_type dlamch(const char *cmach)
{
    return LAPACK_NAME(dlamch)(_rc(cmach));
};

//-----------------------------------------------------------------------
//                          LACGV
//-----------------------------------------------------------------------

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
lacgv(i_type N, V* X, i_type INCX)
{};

template<> BLAS_EXPORT INLINE_TYPE 
void lacgv<c_type>(i_type N, c_type* X, i_type INCX)
{
    LAPACK_NAME(clacgv)(_rc(&N),_rc(X),_rc(&INCX));
};

template<> BLAS_EXPORT INLINE_TYPE 
void lacgv<d_type>(i_type N, d_type* X, i_type INCX)
{
    (void)N;
    (void)X;
    (void)INCX;
    return;
};

template<> BLAS_EXPORT INLINE_TYPE 
void lacgv<s_type>(i_type N, s_type* X, i_type INCX)
{
    (void)N;
    (void)X;
    (void)INCX;
    return;
};

template<> BLAS_EXPORT INLINE_TYPE 
void lacgv<z_type>(i_type N, z_type* X, i_type INCX)
{
    LAPACK_NAME(zlacgv)(_rc(&N),_rc(X),_rc(&INCX));
};

BLAS_EXPORT void clacgv(i_type N, c_type* X, i_type INCX)
{
    return lapack::lacgv<c_type>(N,X,INCX);
};

BLAS_EXPORT void dlacgv(i_type N, d_type* X, i_type INCX)
{
    return lapack::lacgv<d_type>(N,X,INCX);
};

BLAS_EXPORT void slacgv(i_type N, s_type* X, i_type INCX)
{
    return lapack::lacgv<s_type>(N,X,INCX);
};

BLAS_EXPORT void zlacgv(i_type N, z_type* X, i_type INCX)
{
    return lapack::lacgv<z_type>(N,X,INCX);
};

//-----------------------------------------------------------------------
//                          GBTRF
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrf<s_type>(i_type m,i_type n,i_type kl,i_type ku, s_type *ab,i_type ldab,
                   i_type *ipiv,i_type *info)
{
    LAPACK_NAME(sgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrf<d_type>(i_type m,i_type n,i_type kl,i_type ku, d_type *ab,i_type ldab,
                   i_type *ipiv,i_type *info)
{
    LAPACK_NAME(dgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrf<c_type>(i_type m,i_type n,i_type kl,i_type ku, c_type *ab,i_type ldab,
                   i_type *ipiv,i_type *info)
{
    LAPACK_NAME(cgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrf<z_type>(i_type m,i_type n,i_type kl,i_type ku, z_type *ab,i_type ldab,
                   i_type *ipiv,i_type *info)
{
    LAPACK_NAME(zgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void cgbtrf(i_type m,i_type n,i_type kl,i_type ku, c_type *ab,i_type ldab,i_type *ipiv,
            i_type *info)
{
    LAPACK_NAME(cgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void dgbtrf(i_type m,i_type n,i_type kl,i_type ku, d_type *ab,i_type ldab,i_type *ipiv,
            i_type *info)
{
    LAPACK_NAME(dgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void sgbtrf(i_type m,i_type n,i_type kl,i_type ku, s_type *ab,i_type ldab,i_type *ipiv,
            i_type *info)
{
    LAPACK_NAME(sgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void zgbtrf(i_type m,i_type n,i_type kl,i_type ku, z_type *ab,i_type ldab,i_type *ipiv,
            i_type *info)
{
    LAPACK_NAME(zgbtrf)(_rc(&m),_rc(&n),_rc(&kl),_rc(&ku), _rc(ab),_rc(&ldab),
                        _rc(ipiv),_rc(info));
};

//-----------------------------------------------------------------------
//                          GESV
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gesv<d_type>(i_type n,i_type nrhs, d_type *a,i_type lda,i_type *ipiv,d_type* b,
                  i_type ldb,i_type *info)
{
    LAPACK_NAME(dgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesv<c_type>(i_type n,i_type nrhs, c_type *a,i_type lda,i_type *ipiv,c_type* b,
                  i_type ldb,i_type *info)
{
    LAPACK_NAME(cgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesv<s_type>(i_type n,i_type nrhs, s_type *a,i_type lda,i_type *ipiv,s_type* b,
                  i_type ldb,i_type *info)
{
    LAPACK_NAME(sgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesv<z_type>(i_type n,i_type nrhs, z_type *a,i_type lda,i_type *ipiv,z_type* b,
                  i_type ldb,i_type *info)
{
    LAPACK_NAME(zgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dgesv(i_type n,i_type nrhs, d_type *a,i_type lda,i_type *ipiv,d_type* b,
           i_type ldb,i_type *info)
{
    LAPACK_NAME(dgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void cgesv(i_type n,i_type nrhs, c_type *a,i_type lda,i_type *ipiv,c_type* b,
           i_type ldb,i_type *info)
{
    LAPACK_NAME(cgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sgesv(i_type n,i_type nrhs, s_type *a,i_type lda,i_type *ipiv,s_type* b,
           i_type ldb,i_type *info)
{
    LAPACK_NAME(sgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zgesv(i_type n,i_type nrhs, z_type *a,i_type lda,i_type *ipiv,z_type* b,
           i_type ldb,i_type *info)
{
    LAPACK_NAME(zgesv)(_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                       _rc(&ldb),_rc(info));
};

//-----------------------------------------------------------------------
//                          GETRS
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void getrs<d_type>(const char *trans,i_type n,i_type nrhs,const d_type *a,i_type lda,
                   const i_type *ipiv,d_type *b,i_type ldb, i_type *info)
{
    LAPACK_NAME(dgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),
                        _rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void getrs<s_type>(const char *trans,i_type n,i_type nrhs,const s_type *a,i_type lda,
                   const i_type *ipiv,s_type *b,i_type ldb, i_type *info)
{
    LAPACK_NAME(sgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),
                        _rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void getrs<c_type>(const char *trans,i_type n,i_type nrhs,const c_type *a,i_type lda,
                   const i_type *ipiv,c_type *b,i_type ldb, i_type *info)
{
    LAPACK_NAME(cgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),_rc(b),
                        _rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void getrs<z_type>(const char *trans,i_type n,i_type nrhs,const z_type *a,i_type lda,
                   const i_type *ipiv,z_type *b,i_type ldb, i_type *info)
{
    LAPACK_NAME(zgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),
                        _rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void cgetrs(const char *trans,i_type n,i_type nrhs,const c_type *a,i_type lda,
            const i_type *ipiv,c_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(cgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),
                        _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sgetrs(const char *trans,i_type n,i_type nrhs,const s_type *a,i_type lda,
            const i_type *ipiv,s_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(sgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),
                        _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dgetrs(const char *trans,i_type n,i_type nrhs,const d_type *a,i_type lda,
            const i_type *ipiv,d_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(dgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),
                        _rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zgetrs(const char *trans,i_type n,i_type nrhs,const z_type *a,i_type lda,
            const i_type *ipiv,z_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(zgetrs)(_rc(trans),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),_rc(ipiv),
                        _rc(b),_rc(&ldb),_rc(info));
};

//-----------------------------------------------------------------------
//                          TRTRS
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void trtrs<d_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                   const d_type *a,i_type lda,d_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(dtrtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),
                        _rc(&lda),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trtrs<c_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                   const c_type *a,i_type lda,c_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(ctrtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),
                        _rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trtrs<s_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                   const s_type *a,i_type lda,s_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(strtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),
                        _rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void trtrs<z_type>(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                   const z_type *a,i_type lda,z_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(ztrtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),_rc(&lda),
                        _rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dtrtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
            const d_type *a,i_type lda,d_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(dtrtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),
                        _rc(&lda),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void ctrtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
            const c_type *a,i_type lda,c_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(ctrtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),
                        _rc(&lda),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void strtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
            const s_type *a,i_type lda,s_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(strtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),
                        _rc(&lda),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void ztrtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
            const z_type *a,i_type lda,z_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(ztrtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&nrhs),_rc(a),
                        _rc(&lda),_rc(b),_rc(&ldb),_rc(info));
};

//-----------------------------------------------------------------------
//                          GBSV
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gbsv<d_type>(i_type n,i_type kl,i_type ku,i_type nrhs,d_type *ab,i_type ldab, 
                  i_type *ipiv,d_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(dgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbsv<c_type>(i_type n,i_type kl,i_type ku,i_type nrhs,c_type *ab,i_type ldab, 
                  i_type *ipiv,c_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(cgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbsv<z_type>(i_type n,i_type kl,i_type ku,i_type nrhs,z_type *ab,i_type ldab, 
                  i_type *ipiv,z_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(zgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbsv<s_type>(i_type n,i_type kl,i_type ku,i_type nrhs,s_type *ab,i_type ldab, 
                  i_type *ipiv,s_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(sgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE  
void cgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,c_type *ab,i_type ldab,
           i_type *ipiv,c_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(cgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE  
void sgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,s_type *ab,i_type ldab,
           i_type *ipiv,s_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(sgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE  
void dgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,d_type *ab,i_type ldab,
           i_type *ipiv,d_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(dgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE  
void zgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,z_type *ab,i_type ldab,
           i_type *ipiv,z_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(zgbsv)(_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
                       _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

//-----------------------------------------------------------------------
//                          GBTRS
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrs<d_type>(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,
                   const d_type *ab, i_type ldab,const i_type *ipiv,d_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(dgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),
                        _rc(&ldab),_rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrs<s_type>(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,
                   const s_type *ab, i_type ldab,const i_type *ipiv,s_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(sgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),
                        _rc(&ldab),_rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrs<c_type>(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,
                   const c_type *ab, i_type ldab,const i_type *ipiv,c_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(cgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),
                        _rc(&ldab),_rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gbtrs<z_type>(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,
                   const z_type *ab, i_type ldab,const i_type *ipiv,z_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(zgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),
                        _rc(&ldab),_rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void cgbtrs(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const c_type *ab,
            i_type ldab,const i_type *ipiv,c_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(cgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
        _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void sgbtrs(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const s_type *ab,
            i_type ldab,const i_type *ipiv,s_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(sgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
        _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void dgbtrs(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const d_type *ab,
            i_type ldab,const i_type *ipiv,d_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(dgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
        _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE
void zgbtrs(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const z_type *ab,
            i_type ldab,const i_type *ipiv,z_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(zgbtrs)(_rc(trans),_rc(&n),_rc(&kl),_rc(&ku),_rc(&nrhs),_rc(ab),_rc(&ldab),
        _rc(ipiv),_rc(b),_rc(&ldb),_rc(info));
};

//-----------------------------------------------------------------------
//                          TBTRS
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void tbtrs<c_type>(const char *uplo,const char *trans,const char *diag,i_type n,
                   i_type kd,i_type nrhs, const c_type *ab,i_type ldab,c_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(ctbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void tbtrs<d_type>(const char *uplo,const char *trans,const char *diag,i_type n,
                   i_type kd,i_type nrhs, const d_type *ab,i_type ldab,d_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(dtbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void tbtrs<s_type>(const char *uplo,const char *trans,const char *diag,i_type n,
                   i_type kd,i_type nrhs, const s_type *ab,i_type ldab,s_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(stbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void tbtrs<z_type>(const char *uplo,const char *trans,const char *diag,i_type n,
                   i_type kd,i_type nrhs, const z_type *ab,i_type ldab,z_type *b,
                   i_type ldb,i_type *info)
{
    LAPACK_NAME(ztbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void ctbtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,
            i_type nrhs, const c_type *ab,i_type ldab,c_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(ctbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void stbtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,
            i_type nrhs, const s_type *ab,i_type ldab,s_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(stbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dtbtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,
            i_type nrhs, const d_type *ab,i_type ldab,d_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(dtbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void ztbtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,
            i_type nrhs, const z_type *ab,i_type ldab,z_type *b,i_type ldb,i_type *info)
{
    LAPACK_NAME(ztbtrs)(_rc(uplo),_rc(trans),_rc(diag),_rc(&n),_rc(&kd),_rc(&nrhs),_rc(ab),
                    _rc(&ldab),_rc(b),_rc(&ldb),_rc(info));
};

//-----------------------------------------------------------------------
//                          GEES
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gees<s_type>(const char *jobvs,const char *sort,sel_fun selctg,
                  i_type n, s_type *a,i_type lda,i_type *sdim, 
                  c_type *eig, s_type *vs,i_type ldvs,
                  s_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(sgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),nullptr,nullptr,_rc(vs),_rc(&ldvs),
                          _rc(work),_rc(&lwork),_rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;

        workspace eigr(n+1);
        workspace eigi(n+1);

        LAPACK_NAME(sgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),_rc(&eigr[0]),_rc(&eigi[0]),_rc(vs),_rc(&ldvs),
                          _rc(work),_rc(&lwork),_rcl(bwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            eig[i] = c_type(eigr[i],eigi[i]);

        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gees<d_type>(const char *jobvs,const char *sort,sel_fun selctg,
                  i_type n, d_type *a,i_type lda,i_type *sdim, 
                  z_type *eig,d_type *vs,i_type ldvs,
                  d_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(dgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),nullptr,nullptr,_rc(vs),_rc(&ldvs), 
                          _rc(work),_rc(&lwork),_rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        workspace eigr(n+1);
        workspace eigi(n+1);

        LAPACK_NAME(dgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),_rc(&eigr[0]),_rc(&eigi[0]),_rc(vs),_rc(&ldvs),
                          _rc(work),_rc(&lwork),_rcl(bwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            eig[i] = z_type(eigr[i],eigi[i]);

        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gees<c_type>(const char *jobvs,const char *sort,sel_fun selctg,
                  i_type n, c_type *a,i_type lda,i_type *sdim, 
                  c_type *eig, c_type *vs,i_type ldvs,
                  c_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(cgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                      _rc(sdim),_rc(eig),_rc(vs),_rc(&ldvs),
                      _rc(work),_rc(&lwork),nullptr,_rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;

        workspace rwork(n+1);

        LAPACK_NAME(cgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                      _rc(sdim),_rc(eig),_rc(vs),_rc(&ldvs),
                      _rc(work),_rc(&lwork),_rc(&rwork[0]),_rcl(bwork),_rc(info));
        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gees<z_type>(const char *jobvs,const char *sort,sel_fun selctg,
                  i_type n, z_type *a,i_type lda,i_type *sdim, 
                  z_type *eig, z_type *vs,i_type ldvs,
                  z_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(zgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                      _rc(sdim),_rc(eig),_rc(vs),_rc(&ldvs),
                      _rc(work),_rc(&lwork),nullptr,_rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        workspace rwork(n+1);
        LAPACK_NAME(zgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                      _rc(sdim),_rc(eig),_rc(vs),_rc(&ldvs),
                      _rc(work),_rc(&lwork),_rc(&rwork[0]),_rcl(bwork),_rc(info));
    };
};

BLAS_EXPORT INLINE_TYPE 
void cgees(const char *jobvs,const char *sort,sel_fun selctg,
           i_type n, c_type *a,i_type lda,i_type *sdim, 
           c_type *eig, c_type *vs,i_type ldvs, 
           c_type *work,i_type lwork,s_type *rwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(cgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),_rc(eig), _rc(vs),_rc(&ldvs), 
                          _rc(work),_rc(&lwork),_rc(rwork),_rcl(bwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sgees(const char *jobvs,const char *sort,sel_fun selctg,
           i_type n,s_type *a,i_type lda,i_type *sdim,
           s_type *eigr,s_type *eigi,s_type *vs,i_type ldvs,
           s_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(sgees)(_rc(jobvs),_rc(sort),_rc(selctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),_rc(eigr),_rc(eigi),_rc(vs),_rc(&ldvs), 
                          _rc(work),_rc(&lwork),_rcl(bwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dgees(const char *jobvs,const char *sort,sel_fun delctg,
           i_type n,d_type *a,i_type lda,i_type *sdim,
           d_type *eigr,d_type *eigi,d_type *vs,i_type ldvs,
           d_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(dgees)(_rc(jobvs),_rc(sort),_rc(delctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),_rc(eigr),_rc(eigi),_rc(vs),_rc(&ldvs), 
                          _rc(work),_rc(&lwork),_rcl(bwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zgees(const char *jobvs,const char *sort,sel_fun delctg,
           i_type n, z_type *a,i_type lda,i_type *sdim, 
           z_type *eig,  z_type *vs,i_type ldvs,  
           z_type *work,i_type lwork,d_type *rwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(zgees)(_rc(jobvs),_rc(sort),_rc(delctg),_rc(&n),_rc(a),_rc(&lda),
                          _rc(sdim),_rc(eig), _rc(vs),_rc(&ldvs), 
                          _rc(work),_rc(&lwork),_rc(rwork),_rcl(bwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          GGES
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gges<s_type>(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun selctg,
                  i_type n, s_type *a,i_type lda, s_type *b,i_type ldb,i_type *sdim, 
                  c_type *alpha, s_type *beta, s_type *vsl,i_type ldvsl, s_type *vsr,
                  i_type ldvsr, s_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(sgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),nullptr,nullptr,_rc(beta),
                           _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),_rc(&lwork),
                           _rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;

        workspace alphar(n+1), alphai(n+1);
        LAPACK_NAME(sgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(&alphar[0]),_rc(&alphai[0]),
                           _rc(beta),_rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),
                           _rc(&lwork),_rcl(bwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            alpha[i] = c_type(alphar[i],alphai[i]);

        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gges<d_type>(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun selctg,
                  i_type n, d_type *a,i_type lda, d_type *b,i_type ldb,i_type *sdim, 
                  z_type *alpha, d_type *beta, d_type *vsl,i_type ldvsl, d_type *vsr,i_type ldvsr, 
                  d_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(dgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),nullptr,nullptr,_rc(beta),
                           _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),_rc(&lwork),
                           _rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        workspace alphar(n+1), alphai(n+1);

        LAPACK_NAME(dgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(alphar.ptr()),
                           _rc(alphai.ptr()),_rc(beta),_rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), 
                          _rc(work),_rc(&lwork),_rcl(bwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            alpha[i] = z_type(alphar[i],alphai[i]);

        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gges<c_type>(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun selctg,
                  i_type n, c_type *a,i_type lda, c_type *b,i_type ldb,i_type *sdim, 
                  c_type *alpha, c_type *beta, c_type *vsl,i_type ldvsl, c_type *vsr,i_type ldvsr, 
                  c_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(cgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(alpha),_rc(beta),
                           _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),
                           _rc(&lwork),nullptr,_rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;

        workspace rwork(8*n+1);
        LAPACK_NAME(cgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(alpha),_rc(beta),
                           _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),_rc(&lwork),
                           _rc(&rwork[0]),_rcl(bwork),_rc(info));
        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gges<z_type>(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun selctg,
                  i_type n, z_type *a,i_type lda, z_type *b,i_type ldb,i_type *sdim, 
                  z_type *alpha, z_type *beta, z_type *vsl,i_type ldvsl, z_type *vsr,i_type ldvsr, 
                  z_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(zgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(alpha),_rc(beta),
                           _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),_rc(&lwork),
                           nullptr,_rcl(bwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        workspace rwork(8*n+1);

        LAPACK_NAME(zgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                           _rc(&lda),_rc(b),_rc(&ldb),_rc(sdim),_rc(alpha),_rc(beta),
                           _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),_rc(&lwork),
                           _rc(rwork.ptr()),_rcl(bwork),_rc(info));
    };
};

BLAS_EXPORT INLINE_TYPE 
void cgges(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun selctg,
           i_type n, c_type *a,i_type lda, c_type *b,i_type ldb,i_type *sdim, 
           c_type *alpha, c_type *beta, c_type *vsl,i_type ldvsl, c_type *vsr,i_type ldvsr, 
           c_type *work,i_type lwork,s_type *rwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(cgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                       _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(alpha), _rc(beta),
                       _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),_rc(&lwork),
                       _rc(rwork),_rcl(bwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sgges(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun selctg,
           i_type n,s_type *a,i_type lda,s_type *b,i_type ldb,i_type *sdim,
           s_type *alphar,s_type *alphai,s_type *beta,s_type *vsl,i_type ldvsl,
           s_type *vsr,i_type ldvsr,s_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(sgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(selctg),_rc(&n),_rc(a),
                       _rc(&lda),_rc(b),_rc(&ldb),_rc(sdim),_rc(alphar),_rc(alphai),
                       _rc(beta),_rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),
                       _rc(&lwork),_rcl(bwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dgges(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun delctg,
           i_type n,d_type *a,i_type lda,d_type *b,i_type ldb,i_type *sdim,
           d_type *alphar,d_type *alphai,d_type *beta,d_type *vsl,i_type ldvsl,
           d_type *vsr,i_type ldvsr,d_type *work,i_type lwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(dgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(delctg),_rc(&n),_rc(a),
                       _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(alphar),_rc(alphai),
                       _rc(beta),_rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),
                       _rc(&lwork),_rcl(bwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zgges(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun delctg,
           i_type n, z_type *a,i_type lda, z_type *b,i_type ldb,i_type *sdim, 
           z_type *alpha, z_type *beta, z_type *vsl,i_type ldvsl, z_type *vsr,i_type ldvsr, 
           z_type *work,i_type lwork,d_type *rwork,l_type *bwork,i_type *info)
{
    LAPACK_NAME(zgges)(_rc(jobvsl),_rc(jobvsr),_rc(sort),_rc(delctg),_rc(&n),_rc(a),
                       _rc(&lda),_rc(b),_rc(&ldb), _rc(sdim),_rc(alpha), _rc(beta),
                       _rc(vsl),_rc(&ldvsl),_rc(vsr),_rc(&ldvsr), _rc(work),_rc(&lwork),
                       _rc(rwork),_rcl(bwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          HEEVR / SYEVR
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void heevr<c_type>(const char *jobv,const char *range,const char *uplo, i_type n, c_type *a,
                   i_type lda, s_type vl, s_type vu, i_type il, i_type iu, s_type abstol,
                   i_type *m, s_type* w, c_type* z, i_type ldz, i_type *isuppz,
                   c_type *work,i_type lwork, s_type *rwork, i_type lrwork,
                   i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(cheevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(rwork), _rc(&lrwork), _rc(iwork),
                        _rc(&liwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void heevr<z_type>(const char *jobv,const char *range,const char *uplo, i_type n, z_type *a,
                   i_type lda, d_type vl, d_type vu,i_type il, i_type iu, d_type abstol,
                   i_type *m, d_type* w, z_type* z, i_type ldz, i_type *isuppz,
                   z_type *work,i_type lwork, d_type *rwork, i_type lrwork,
                   i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(zheevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(rwork), _rc(&lrwork), _rc(iwork),
                        _rc(&liwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void heevr<s_type>(const char *jobv,const char *range,const char *uplo, i_type n, s_type *a,
                   i_type lda, s_type vl, s_type vu, i_type il, i_type iu, s_type abstol,
                   i_type *m, s_type* w, s_type* z, i_type ldz, i_type *isuppz,
                   s_type *work,i_type lwork, s_type *rwork, i_type lrwork,
                   i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ssyevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(iwork),_rc(&liwork),_rc(info));

    if (lrwork == -1 || lwork == -1 || liwork == -1)
        *rwork = 0;
}

template<> BLAS_EXPORT INLINE_TYPE
void heevr<d_type>(const char *jobv,const char *range,const char *uplo, i_type n, d_type *a,
                   i_type lda, d_type vl, d_type vu, i_type il, i_type iu, d_type abstol,
                   i_type *m, d_type* w, d_type* z, i_type ldz, i_type *isuppz,
                   d_type *work,i_type lwork, d_type *rwork, i_type lrwork,
                   i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(dsyevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(iwork),_rc(&liwork),_rc(info));

    if (lrwork == -1 || lwork == -1 || liwork == -1)
        *rwork = 0;
}

BLAS_EXPORT INLINE_TYPE 
void cheevr(const char *jobv,const char *range,const char *uplo,i_type n, c_type *a,
            i_type lda, s_type vl, s_type vu, i_type il, i_type iu, s_type abstol,
            i_type *m, s_type* w, c_type* z, i_type ldz, i_type *isuppz,
            c_type *work,i_type lwork, s_type *rwork, i_type lrwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(cheevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(rwork), _rc(&lrwork),
                          _rc(iwork),_rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ssyevr(const char *jobv,const char *range,const char *uplo,i_type n, s_type *a,
            i_type lda, s_type vl, s_type vu, i_type il, i_type iu, s_type abstol,
            i_type *m, s_type* w, s_type* z, i_type ldz, i_type *isuppz,
            s_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ssyevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(iwork),_rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dsyevr(const char *jobv,const char *range,const char *uplo,i_type n, d_type *a,
            i_type lda, d_type vl, d_type vu, i_type il, i_type iu, d_type abstol,
            i_type *m, d_type* w, d_type* z, i_type ldz, i_type *isuppz,
            d_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(dsyevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(iwork),_rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void zheevr(const char *jobv,const char *range,const char *uplo,i_type n, z_type *a,
            i_type lda, d_type vl, d_type vu, i_type il, i_type iu, d_type abstol,
            i_type *m, d_type* w, z_type* z, i_type ldz, i_type *isuppz,
            z_type *work,i_type lwork, d_type *rwork, i_type lrwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(zheevr)(_rc(jobv),_rc(range),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(&vl), _rc(&vu), _rc(&il), _rc(&iu), _rc(&abstol),
                          _rc(m), _rc(w), _rc(z), _rc(&ldz), _rc(isuppz),
                          _rc(work),_rc(&lwork), _rc(rwork), _rc(&lrwork),
                          _rc(iwork),_rc(&liwork),_rc(info));
}

//-----------------------------------------------------------------------
//                          HEEV / SYEV
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void heev<c_type>(const char *jobv,const char *uplo, i_type n, c_type *a,i_type lda, 
                  s_type* w, c_type *work,i_type lwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(cheev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                              _rc(w),_rc(work),_rc(&lwork), nullptr,_rc(info));
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;

        workspace rwork(lapack::maximum(1,3*n-2));
        LAPACK_NAME(cheev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                              _rc(w),_rc(work),_rc(&lwork), _rc(&rwork[0]),_rc(info));
    }
}

template<> BLAS_EXPORT INLINE_TYPE
void heev<z_type>(const char *jobv,const char *uplo, i_type n, z_type *a,i_type lda, 
                  d_type* w, z_type *work,i_type lwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(zheev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                              _rc(w),_rc(work),_rc(&lwork), nullptr,_rc(info));
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        workspace rwork(lapack::maximum(1,3*n-2));
        LAPACK_NAME(zheev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                              _rc(w),_rc(work),_rc(&lwork), _rc(&rwork[0]),_rc(info));
    }
}

template<> BLAS_EXPORT INLINE_TYPE
void heev<s_type>(const char *jobv,const char *uplo, i_type n, s_type *a,i_type lda, 
                  s_type* w, s_type *work,i_type lwork, i_type *info)
{
    LAPACK_NAME(ssyev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(w), _rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void heev<d_type>(const char *jobv,const char *uplo, i_type n, d_type *a,i_type lda, 
                  d_type* w, d_type *work,i_type lwork, i_type *info)
{
    LAPACK_NAME(dsyev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(w),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void cheev(const char *jobv,const char *uplo,i_type n, c_type *a,i_type lda, 
           s_type* w, c_type *work,i_type lwork, s_type *rwork, i_type *info)
{
    LAPACK_NAME(cheev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(w), _rc(work),_rc(&lwork), _rc(rwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ssyev(const char *jobv,const char *uplo,i_type n, s_type *a,i_type lda, 
            s_type* w, s_type *work,i_type lwork, i_type *info)
{
    LAPACK_NAME(ssyev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(w), _rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dsyev(const char *jobv,const char *uplo,i_type n, d_type *a,i_type lda, 
            d_type* w, d_type *work,i_type lwork, i_type *info)
{
    LAPACK_NAME(dsyev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(w),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void zheev(const char *jobv,const char *uplo,i_type n, z_type *a,i_type lda, 
            d_type* w, z_type *work,i_type lwork, d_type *rwork, i_type *info)
{
    LAPACK_NAME(zheev)(_rc(jobv),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                          _rc(w), _rc(work),_rc(&lwork), _rc(rwork),_rc(info));
}

//-----------------------------------------------------------------------
//                          SYGVD/HEGVD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void hegvd<c_type>(i_type itype, const char* jobz,const char *uplo, i_type n, 
                    c_type *a,i_type lda, c_type *b,i_type ldb, s_type * w,
                    c_type *work,i_type lwork, s_type* rwork, i_type lrwork,
                    i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(chegvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda),  
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(rwork), 
                        _rc(&lrwork), _rc(iwork),_rc(&liwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hegvd<s_type>(i_type itype, const char* jobz,const char *uplo, i_type n, 
                    s_type *a,i_type lda, s_type *b,i_type ldb, s_type * w,
                    s_type *work,i_type lwork, s_type* rwork, i_type lrwork,
                    i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ssygvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda),
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(iwork),
                        _rc(&liwork),_rc(info));

    if (lrwork == -1 || lwork == -1 || liwork == -1)
        *rwork = 0;
}

template<> BLAS_EXPORT INLINE_TYPE
void hegvd<d_type>(i_type itype, const char* jobz,const char *uplo, i_type n, 
                    d_type *a,i_type lda, d_type *b,i_type ldb, d_type *w,
                    d_type *work,i_type lwork, d_type* rwork, i_type lrwork,
                    i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(dsygvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(iwork),
                        _rc(&liwork),_rc(info));

    if (lrwork == -1 || lwork == -1 || liwork == -1)
        *rwork = 0;
}

template<> BLAS_EXPORT INLINE_TYPE
void hegvd<z_type>(i_type itype, const char* jobz,const char *uplo, i_type n, 
                    z_type *a,i_type lda, z_type *b,i_type ldb, d_type *w,
                    z_type *work,i_type lwork, d_type* rwork, i_type lrwork,
                    i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(zhegvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(rwork),
                        _rc(&lrwork), _rc(iwork),_rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void chegvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
            c_type *a,i_type lda, c_type *b,i_type ldb, s_type* w,
            c_type *work,i_type lwork, s_type* rwork, i_type lrwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(chegvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda),
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(rwork),
                        _rc(&lrwork), _rc(iwork),_rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ssygvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
            s_type *a,i_type lda, s_type *b,i_type ldb, s_type *w,
            s_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ssygvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(iwork),
                        _rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dsygvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
            d_type *a,i_type lda, d_type *b,i_type ldb, d_type * w,
            d_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(dsygvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(iwork),
                        _rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void zhegvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
            z_type *a,i_type lda, z_type *b,i_type ldb, d_type * w,
            z_type *work,i_type lwork, d_type* rwork, i_type lrwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(zhegvd)(_rc(&itype),_rc(jobz),_rc(uplo),_rc(&n), _rc(a),_rc(&lda), 
                        _rc(b),_rc(&ldb), _rc(w), _rc(work),_rc(&lwork), _rc(rwork),
                        _rc(&lrwork), _rc(iwork),_rc(&liwork),_rc(info));
}

//-----------------------------------------------------------------------
//                          GESVD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gesvd<s_type>(const char *jobu,const char *jobvt,i_type m,i_type n,s_type *a,
                   i_type lda, s_type *s,s_type *u,i_type ldu,s_type *vt,i_type ldvt,
                    s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),
                        _rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),
                        _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesvd<d_type>(const char *jobu,const char *jobvt,i_type m,i_type n,d_type *a,
                   i_type lda, d_type *s,d_type *u,i_type ldu,d_type *vt,i_type ldvt,
                   d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),
                        _rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),
                        _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesvd<c_type>(const char *jobu,const char *jobvt,i_type m,i_type n,c_type *a,
                   i_type lda, s_type *s,c_type *u,i_type ldu,c_type *vt,i_type ldvt,
                   c_type *work,i_type lwork,i_type *info)
{    
    if (lwork == -1)
    {
        LAPACK_NAME(cgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),
                            _rc(s),_rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),
                            _rc(&lwork),nullptr,_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;

        i_type K = (m<n)?m:n;
        workspace rwork(5*K+1);

        LAPACK_NAME(cgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),
                            _rc(s),_rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),
                            _rc(&lwork),_rc(&rwork[0]),_rc(info));
        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesvd<z_type>(const char *jobu,const char *jobvt,i_type m,i_type n,z_type *a,
                   i_type lda,d_type *s,z_type *u,i_type ldu,z_type *vt,i_type ldvt,
                   z_type *work,i_type lwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(zgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),
                            _rc(s),_rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),
                            _rc(&lwork),nullptr,_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        i_type K = (m<n)?m:n;
        workspace rwork(5*K+1);

        LAPACK_NAME(zgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),
                            _rc(s),_rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),
                            _rc(&lwork),_rc(&rwork[0]),_rc(info));
        return;
    };
};

BLAS_EXPORT INLINE_TYPE 
void cgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,c_type *a,i_type lda,
            s_type *s,c_type *u,i_type ldu,c_type *vt,i_type ldvt,
            c_type *work,i_type lwork,s_type *rwork, i_type *info)
{
    LAPACK_NAME(cgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),
                        _rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),
                        _rc(rwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,s_type *a,i_type lda,
            s_type *s,s_type *u,i_type ldu,s_type *vt,i_type ldvt,
            s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),
                        _rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,d_type *a,i_type lda,
            d_type *s,d_type *u,i_type ldu,d_type *vt,i_type ldvt,
            d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),
                        _rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,z_type *a,i_type lda,
            d_type *s, z_type *u,i_type ldu,z_type *vt,i_type ldvt,
            z_type *work,i_type lwork,d_type *rwork,i_type *info)
{
    LAPACK_NAME(zgesvd)(_rc(jobu),_rc(jobvt),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),
                        _rc(u),_rc(&ldu), _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),
                        _rc(rwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          GESDD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void gesdd<s_type>(const char *jobu,i_type m,i_type n,s_type *a,i_type lda,
                    s_type *s,s_type *u,i_type ldu,s_type *vt,i_type ldvt,
                    s_type *work,i_type lwork, i_type *iwork,i_type *info)
{
    LAPACK_NAME(sgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesdd<d_type>(const char *jobu,i_type m,i_type n,d_type *a,i_type lda,
                    d_type *s,d_type *u,i_type ldu,d_type *vt,i_type ldvt,
                    d_type *work,i_type lwork, i_type *iwork,i_type *info)
{
    LAPACK_NAME(dgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesdd<c_type>(const char *jobu,i_type m,i_type n,c_type *a,i_type lda,
                    s_type *s,c_type *u,i_type ldu,c_type *vt,i_type ldvt,
                    c_type *work,i_type lwork, i_type *iwork,i_type *info)
{    
    if (lwork == -1)
    {
        LAPACK_NAME(cgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),nullptr,_rc(iwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;

        i_type MX   = (m>n)?m:n;
        i_type MN   = (m<n)?m:n;
        i_type LW;

        if (jobu == nullptr || jobu[0] == 'N' || jobu[0] == 'n')
            LW      = 5*MN+1;
        else
            LW      = maximum(5*MN*MN + 5*MN, 2*MX*MN + 2*MN*MN + MN) + 1;

        workspace rwork(LW);

        LAPACK_NAME(cgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(&rwork[0]),_rc(iwork),_rc(info));
        return;
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void gesdd<z_type>(const char *jobu,i_type m,i_type n,z_type *a,i_type lda,
                    d_type *s,z_type *u,i_type ldu,z_type *vt,i_type ldvt,
                    z_type *work,i_type lwork, i_type *iwork,i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(zgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),nullptr,_rc(iwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        i_type MX   = (m>n)?m:n;
        i_type MN   = (m<n)?m:n;
        i_type LW;

        if (jobu == nullptr || jobu[0] == 'N' || jobu[0] == 'n')
            LW      = 5*MN+1;
        else
            LW      = maximum(5*MN*MN + 5*MN, 2*MX*MN + 2*MN*MN + MN) + 1;

        workspace rwork(LW);

        LAPACK_NAME(zgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(&rwork[0]),_rc(iwork),_rc(info));
        return;
    };
};

BLAS_EXPORT INLINE_TYPE 
void lapack::cgesdd(const char *jobu,i_type m,i_type n,c_type *a,i_type lda,
                    s_type *s,c_type *u,i_type ldu,c_type *vt,i_type ldvt,
                    c_type *work,i_type lwork,s_type *rwork, i_type *iwork, i_type *info)
{
    LAPACK_NAME(cgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(rwork),_rc(iwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void lapack::sgesdd(const char *jobu,i_type m,i_type n,s_type *a,i_type lda,
                    s_type *s,s_type *u,i_type ldu,s_type *vt,i_type ldvt,
                    s_type *work,i_type lwork, i_type *iwork,i_type *info)
{
    LAPACK_NAME(sgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void lapack::dgesdd(const char *jobu,i_type m,i_type n,d_type *a,i_type lda,
                    d_type *s,d_type *u,i_type ldu,d_type *vt,i_type ldvt,
                    d_type *work,i_type lwork, i_type *iwork,i_type *info)
{
    LAPACK_NAME(dgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void lapack::zgesdd(const char *jobu,i_type m,i_type n,z_type *a,i_type lda,
                    d_type *s, z_type *u,i_type ldu,z_type *vt,i_type ldvt,
                    z_type *work,i_type lwork,d_type *rwork, i_type *iwork,i_type *info)
{
    LAPACK_NAME(zgesdd)(_rc(jobu),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(s),_rc(u),_rc(&ldu),
                _rc(vt),_rc(&ldvt),_rc(work),_rc(&lwork),_rc(rwork),_rc(iwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          TGSEN
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void tgsen<s_type>(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
                    s_type *a,i_type lda,s_type *b,i_type ldb, 
                    c_type *alpha, s_type *beta, s_type *q,i_type ldq, s_type *z,
                    i_type ldz,i_type *m,s_type *pl,s_type *pr,s_type *dif, 
                    s_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    if (liwork == -1 || lwork == -1)
    {
        LAPACK_NAME(stgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),
                            _rc(a),_rc(&lda),_rc(b), _rc(&ldb),nullptr,nullptr,_rc(beta),
                            _rc(q),_rc(&ldq),_rc(z),_rc(&ldz),_rc(m), _rc(pl),_rc(pr),
                            _rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;
        workspace alphar(n+1), alphai(n+1);

        LAPACK_NAME(stgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                            _rc(&lda),_rc(b), _rc(&ldb),&alphar[0],&alphai[0],_rc(beta),
                            _rc(q),_rc(&ldq),_rc(z),_rc(&ldz),_rc(m), _rc(pl),_rc(pr),
                            _rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            alpha[i] = c_type(alphar[i],alphai[i]);

    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void tgsen<d_type>(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
                    d_type *a,i_type lda,d_type *b,i_type ldb, 
                    z_type *alpha, d_type *beta, d_type *q,i_type ldq, d_type *z,
                    i_type ldz,i_type *m,d_type *pl,d_type *pr,d_type *dif, 
                    d_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    if (liwork == -1 || lwork == -1)
    {
        LAPACK_NAME(dtgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),
                            _rc(a),_rc(&lda),_rc(b), _rc(&ldb),nullptr,nullptr,_rc(beta),
                            _rc(q),_rc(&ldq),_rc(z),_rc(&ldz),_rc(m), _rc(pl),_rc(pr),
                            _rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        workspace alphar(n+1), alphai(n+1);

        LAPACK_NAME(dtgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                            _rc(&lda),_rc(b), _rc(&ldb),&alphar[0],&alphai[0],_rc(beta),
                            _rc(q),_rc(&ldq),_rc(z),_rc(&ldz),_rc(m), _rc(pl),_rc(pr),
                            _rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            alpha[i] = z_type(alphar[i],alphai[i]);
    };
};

template<> BLAS_EXPORT INLINE_TYPE 
void tgsen<c_type>(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
                    c_type *a,i_type lda,c_type *b,i_type ldb, 
                    c_type *alpha, c_type *beta, c_type *q,i_type ldq, c_type *z,
                    i_type ldz,i_type *m,s_type *pl,s_type *pr,s_type *dif, 
                    c_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ctgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                        _rc(&lda),_rc(b),_rc(&ldb),_rc(alpha),_rc(beta),_rc(q),_rc(&ldq),
                        _rc(z),_rc(&ldz),_rc(m),_rc(pl),_rc(pr),_rc(dif),_rc(work),
                        _rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void tgsen<z_type>(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
                    z_type *a,i_type lda,z_type *b,i_type ldb, 
                    z_type *alpha, z_type *beta, z_type *q,i_type ldq, z_type *z,
                    i_type ldz,i_type *m,d_type *pl,d_type *pr,d_type *dif, 
                    z_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ztgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                        _rc(&lda),_rc(b), _rc(&ldb),_rc(alpha),_rc(beta),_rc(q),_rc(&ldq),
                        _rc(z),_rc(&ldz),_rc(m), _rc(pl),_rc(pr),_rc(dif),_rc(work),
                        _rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void ctgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
            c_type *a,i_type lda,c_type *b,i_type ldb, c_type *alpha, 
            c_type *beta, c_type *q,i_type ldq, c_type *z,i_type ldz,i_type *m,
            s_type *pl,s_type *pr,s_type *dif, c_type *work,i_type lwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ctgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                        _rc(&lda),_rc(b), _rc(&ldb),_rc(alpha),_rc(beta),_rc(q),_rc(&ldq),
                        _rc(z),_rc(&ldz),_rc(m), _rc(pl),_rc(pr),_rc(dif),_rc(work),
                        _rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void stgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
            s_type *a,i_type lda,s_type *b,i_type ldb,s_type *alphar,s_type *alphai,
            s_type *beta,s_type *q,i_type ldq,s_type *z,i_type ldz,i_type *m,
            s_type *pl,s_type *pr,s_type *dif,s_type *work,i_type lwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(stgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                        _rc(&lda),_rc(b), _rc(&ldb),_rc(alphar),_rc(alphai),_rc(beta),_rc(q),
                        _rc(&ldq),_rc(z),_rc(&ldz),_rc(m), _rc(pl),_rc(pr),_rc(dif),_rc(work),
                        _rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dtgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
            d_type *a,i_type lda,d_type *b,i_type ldb,d_type *alphar,d_type *alphai,
            d_type *beta,d_type *q,i_type ldq,d_type *z,i_type ldz,i_type *m,
            d_type *pl,d_type *pr,d_type *dif,d_type *work,i_type lwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(dtgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                        _rc(&lda),_rc(b), _rc(&ldb),_rc(alphar),_rc(alphai),_rc(beta),
                        _rc(q),_rc(&ldq),_rc(z),_rc(&ldz),_rc(m),_rc(pl),_rc(pr),_rc(dif),
                        _rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void ztgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n, 
            z_type *a,i_type lda, z_type *b,i_type ldb, z_type *alpha, 
            z_type *beta, z_type *q,i_type ldq, z_type *z,i_type ldz,i_type *m,
            d_type *pl,d_type *pr,d_type *dif, z_type *work,i_type lwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(ztgsen)(_rc(&ijob),_rc(&wantq),_rc(&wantz),_rc(select),_rc(&n),_rc(a),
                        _rc(&lda),_rc(b),_rc(&ldb),_rc(alpha),_rc(beta),_rc(q),_rc(&ldq),
                        _rc(z),_rc(&ldz),_rc(m),_rc(pl),_rc(pr),_rc(dif),_rc(work),
                        _rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          TRSEN
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void trsen<c_type>(const char * job,const char * compq,const i_type *select,i_type n,
                    c_type *t,i_type ldt,c_type *q,i_type ldq, c_type *w, 
                    i_type *m, s_type *s,
                    s_type *sep,c_type *work,i_type lwork,
                    i_type *,i_type ,i_type *info)
{
    LAPACK_NAME(ctrsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                           _rc(&ldq),_rc(w),_rc(m),
                           _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void trsen<z_type>(const char * job,const char * compq,const i_type *select,i_type n,
                    z_type *t,i_type ldt,z_type *q,i_type ldq, z_type *w, 
                    i_type *m, d_type *s,
                    d_type *sep,z_type *work,i_type lwork,
                    i_type *,i_type ,i_type *info)
{
    LAPACK_NAME(ztrsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                           _rc(&ldq),_rc(w),_rc(m),
                           _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void trsen<d_type>(const char * job,const char * compq,const i_type *select,i_type n,
                    d_type *t,i_type ldt,d_type *q,i_type ldq, z_type *w, 
                    i_type *m, d_type *s,
                    d_type *sep,d_type *work,i_type lwork,
                    i_type *iwork,i_type liwork,i_type *info)
{
    if (liwork == -1 || lwork == -1)
    {
        LAPACK_NAME(dtrsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                           _rc(&ldq),nullptr,nullptr,_rc(m),
                           _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),
                            _rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<d_type>;

        workspace wr(n+1), wi(n+1);
        
        LAPACK_NAME(dtrsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                           _rc(&ldq),_rc(&wr[0]),_rc(&wi[0]),_rc(m),
                           _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            w[i] = z_type(wr[i],wi[i]);
    };
};

template<> BLAS_EXPORT INLINE_TYPE
void trsen<s_type>(const char * job,const char * compq,const i_type *select,i_type n,
                    s_type *t,i_type ldt,s_type *q,i_type ldq,c_type *w, 
                    i_type *m, s_type *s,
                    s_type *sep,s_type *work,i_type lwork,
                    i_type *iwork,i_type liwork,i_type *info)
{
    if (liwork == -1 || lwork == -1)
    {
        LAPACK_NAME(strsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                           _rc(&ldq),nullptr,nullptr,_rc(m),
                           _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
        return;
    }
    else
    {
        using workspace     = matcl::pod_workspace<s_type>;
        workspace wr(n+1), wi(n+1);
        
        LAPACK_NAME(strsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                           _rc(&ldq),_rc(&wr[0]),_rc(&wi[0]),_rc(m),
                           _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));

        for (i_type i = 0; i < n; ++i)
            w[i] = z_type(wr[i],wi[i]);
    };
};

BLAS_EXPORT INLINE_TYPE 
void ctrsen(const char * job,const char * compq,const i_type *select,i_type n,
            c_type *t,i_type ldt,c_type *q,i_type ldq, c_type *w,
            i_type *m, s_type *s,s_type *sep,c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(ctrsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                           _rc(&ldq),_rc(w),_rc(m),
                           _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void strsen(const char * job,const char * compq,const i_type *select,i_type n,
            s_type *t,i_type ldt,s_type *q,i_type ldq,s_type *wr,s_type *wi,
            i_type *m, s_type *s,s_type *sep,s_type *work,i_type lwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(strsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                        _rc(&ldq),_rc(wr),_rc(wi),_rc(m),
                        _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),
                        _rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dtrsen(const char * job,const char * compq,const i_type *select,i_type n,
            d_type *t,i_type ldt,d_type *q,i_type ldq,d_type *wr,d_type *wi,
            i_type *m, d_type *s,d_type *sep,d_type *work,i_type lwork,
            i_type *iwork,i_type liwork,i_type *info)
{
    LAPACK_NAME(dtrsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                        _rc(&ldq),_rc(wr),_rc(wi),_rc(m),
                        _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(iwork),_rc(&liwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ztrsen(const char * job,const char * compq,const i_type *select,i_type n, 
            z_type *t,i_type ldt, z_type *q,i_type ldq, z_type *w, 
            i_type *m, d_type *s,d_type *sep, z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(ztrsen)(_rc(job),_rc(compq),_rc(select),_rc(&n),_rc(t),_rc(&ldt),_rc(q),
                        _rc(&ldq),_rc(w),_rc(m),
                        _rc(s),_rc(sep),_rc(work),_rc(&lwork),_rc(info));
}

//-----------------------------------------------------------------------
//                          TRSYL
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void trsyl<c_type>(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                const c_type *a,i_type lda, const c_type *b,i_type ldb, c_type *c,i_type ldc,
                s_type *scale, i_type *info)
{
    LAPACK_NAME(ctrsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void trsyl<s_type>(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                const s_type *a,i_type lda, const s_type *b,i_type ldb, s_type *c,i_type ldc,
                s_type *scale, i_type *info)
{
    LAPACK_NAME(strsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void trsyl<d_type>(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                const d_type *a,i_type lda, const d_type *b,i_type ldb, d_type *c,i_type ldc,
                d_type *scale, i_type *info)
{
    LAPACK_NAME(dtrsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void trsyl<z_type>(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                    const z_type *a,i_type lda, const z_type *b,i_type ldb, z_type *c,i_type ldc,
                    d_type *scale, i_type *info)
{
    LAPACK_NAME(ztrsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void lapack::ctrsyl(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                    const c_type *a,i_type lda, const c_type *b,i_type ldb, c_type *c,i_type ldc,
                    s_type *scale, i_type *info)
{
    LAPACK_NAME(ctrsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void lapack::strsyl(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                    const s_type *a,i_type lda, const s_type *b,i_type ldb, s_type *c,i_type ldc,
                    s_type *scale, i_type *info)
{
    LAPACK_NAME(strsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void lapack::dtrsyl(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                    const d_type *a,i_type lda, const d_type *b,i_type ldb, d_type *c,i_type ldc,
                    d_type *scale, i_type *info)
{
    LAPACK_NAME(dtrsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void lapack::ztrsyl(const char * trana,const char * tranb,const i_type isgn,i_type m,i_type n,
                    const z_type *a,i_type lda, const z_type *b,i_type ldb, z_type *c,i_type ldc,
                    d_type *scale, i_type *info)
{
    LAPACK_NAME(ztrsyl)(_rc(trana),_rc(tranb),_rc(&isgn),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(scale),_rc(info));
}

//-----------------------------------------------------------------------
//                          TGSYL
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void tgsyl<c_type>(const char* trans, const i_type ijob, i_type m, i_type n, 
                   const c_type *a, i_type lda, const c_type *b, i_type ldb, 
                   c_type *c, i_type ldc, const c_type *d, i_type ldd, const c_type *e,
                   i_type lde, c_type *f, i_type ldf, s_type *scale, s_type *dif,
                   c_type *work, i_type lwork, i_type* i_work, i_type *info)
{
    LAPACK_NAME(ctgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n), _rc(a),_rc(&lda),_rc(b),
                        _rc(&ldb),_rc(c),_rc(&ldc),_rc(d),_rc(&ldd),_rc(e),_rc(&lde),
                        _rc(f),_rc(&ldf), _rc(scale),_rc(dif),_rc(work),_rc(&lwork),
                        i_work,_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void tgsyl<s_type>(const char* trans, i_type ijob, i_type m, i_type n, const s_type *a, 
                   i_type lda, const s_type* b, i_type ldb, s_type* c, i_type ldc, 
                   const s_type* d, i_type ldd, const s_type* e, i_type lde, s_type* f,i_type ldf,
                   s_type* scale, s_type* dif, s_type* work, i_type lwork, i_type* i_work,
                   i_type *info)
{
    LAPACK_NAME(stgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb),
                        _rc(c),_rc(&ldc), _rc(d),_rc(&ldd),_rc(e),_rc(&lde),_rc(f),_rc(&ldf), 
                        _rc(scale),_rc(dif),_rc(work),_rc(&lwork), i_work,_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void tgsyl<d_type>(const char* trans, i_type ijob, i_type m, i_type n, const d_type* a,
                   i_type lda, const d_type* b, i_type ldb, d_type *c, i_type ldc, const d_type* d,
                   i_type ldd, const d_type* e, i_type lde, d_type* f, i_type ldf, d_type *scale,
                   d_type *dif, d_type *work, i_type lwork, i_type* i_work, i_type *info)
{
    LAPACK_NAME(dtgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb),
                        _rc(c),_rc(&ldc), _rc(d),_rc(&ldd),_rc(e),_rc(&lde),_rc(f),_rc(&ldf),
                        _rc(scale),_rc(dif),_rc(work),_rc(&lwork), i_work,_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void tgsyl<z_type>(const char* trans, i_type ijob, i_type m, i_type n, const z_type* a,
                   i_type lda, const z_type* b, i_type ldb, z_type* c, i_type ldc, const z_type* d,
                   i_type ldd, const z_type* e, i_type lde, z_type* f, i_type ldf, d_type *scale, 
                   d_type *dif, z_type *work, i_type lwork, i_type* i_work, i_type* info)
{
    LAPACK_NAME(ztgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb),
                        _rc(c),_rc(&ldc), _rc(d),_rc(&ldd),_rc(e),_rc(&lde),_rc(f),_rc(&ldf),
                        _rc(scale),_rc(dif),_rc(work),_rc(&lwork),i_work, _rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ctgsyl(const char * trans,const i_type ijob,i_type m,i_type n,
            c_type *a,i_type lda, c_type *b,i_type ldb, c_type *c,i_type ldc,
            c_type *d,i_type ldd, c_type *e,i_type lde, c_type *f,i_type ldf,
            s_type *scale, s_type *dif,  
            c_type *work, i_type lwork, i_type *iwork, i_type *info)
{
    LAPACK_NAME(ctgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(d),_rc(&ldd),_rc(e),_rc(&lde),_rc(f),_rc(&ldf),
                           _rc(scale),_rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void stgsyl(const char * trans,const i_type ijob,i_type m,i_type n,
            s_type *a,i_type lda, s_type *b,i_type ldb, s_type *c,i_type ldc,
            s_type *d,i_type ldd, s_type *e,i_type lde, s_type *f,i_type ldf,
            s_type *scale, s_type *dif,  
            s_type *work, i_type lwork, i_type *iwork, i_type *info)
{
    LAPACK_NAME(stgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n),   
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(d),_rc(&ldd),_rc(e),_rc(&lde),_rc(f),_rc(&ldf),
                           _rc(scale),_rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dtgsyl(const char * trans,const i_type ijob,i_type m,i_type n,
            d_type *a,i_type lda, d_type *b,i_type ldb, d_type *c,i_type ldc,
            d_type *d,i_type ldd, d_type *e,i_type lde, d_type *f,i_type ldf,
            d_type *scale, d_type *dif,  
            d_type *work, i_type lwork, i_type *iwork, i_type *info)
{
    LAPACK_NAME(dtgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(d),_rc(&ldd),_rc(e),_rc(&lde),_rc(f),_rc(&ldf),
                           _rc(scale),_rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ztgsyl(const char * trans,const i_type ijob,i_type m,i_type n,
            z_type *a,i_type lda, z_type *b,i_type ldb, z_type *c,i_type ldc,
            z_type *d,i_type ldd, z_type *e,i_type lde, z_type *f,i_type ldf,
            d_type *scale, d_type *dif,  
            z_type *work, i_type lwork, i_type *iwork, i_type *info)
{
    LAPACK_NAME(ztgsyl)(_rc(trans),_rc(&ijob),_rc(&m),_rc(&n),
                           _rc(a),_rc(&lda),_rc(b),_rc(&ldb),_rc(c),_rc(&ldc),
                           _rc(d),_rc(&ldd),_rc(e),_rc(&lde),_rc(f),_rc(&ldf),
                           _rc(scale),_rc(dif),_rc(work),_rc(&lwork),_rc(iwork),_rc(info));
}

//-----------------------------------------------------------------------
//                          GEQRF
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void geqrf<c_type>(i_type m,i_type n,c_type *a,i_type lda,c_type *tau,c_type *work,
                    i_type lwork,i_type *info)
{
    LAPACK_NAME(cgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void geqrf<d_type>(i_type m,i_type n,d_type *a,i_type lda,d_type *tau,d_type *work,
                    i_type lwork,i_type *info)
{
    LAPACK_NAME(dgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void geqrf<s_type>(i_type m,i_type n,s_type *a,i_type lda,s_type *tau,s_type *work,
                    i_type lwork,i_type *info)
{
    LAPACK_NAME(sgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void geqrf<z_type>(i_type m,i_type n,z_type *a,i_type lda,z_type *tau,z_type *work,
                    i_type lwork,i_type *info)
{
    LAPACK_NAME(zgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void cgeqrf(i_type m,i_type n,c_type *a,i_type lda,c_type *tau, c_type *work,i_type lwork,
            i_type *info)
{
    LAPACK_NAME(cgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sgeqrf(i_type m,i_type n,s_type *a,i_type lda,s_type *tau, s_type *work,i_type lwork,
            i_type *info)
{
    LAPACK_NAME(sgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dgeqrf(i_type m,i_type n,d_type *a,i_type lda,d_type *tau, d_type *work,i_type lwork,
            i_type *info)
{
    LAPACK_NAME(dgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zgeqrf(i_type m,i_type n,z_type *a,i_type lda,z_type *tau, z_type *work,i_type lwork,
            i_type *info)
{
    LAPACK_NAME(zgeqrf)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(tau),_rc(work),_rc(&lwork),
                        _rc(info));
};

//-----------------------------------------------------------------------
//                          GEQP3
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
geqp3<d_type>(i_type m, i_type n, d_type *a, i_type lda, i_type* jpvt, d_type *tau, d_type *work, 
              i_type lwork, i_type *info)
{
    LAPACK_NAME(dgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
        _rc(tau),_rc(work),_rc(&lwork),_rc(info));
};

template<> BLAS_EXPORT void
geqp3<s_type>(i_type m, i_type n, s_type *a, i_type lda, i_type* jpvt, s_type *tau, s_type *work, 
              i_type lwork, i_type *info)
{
    LAPACK_NAME(sgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
        _rc(tau),_rc(work),_rc(&lwork),_rc(info));
};

template<> BLAS_EXPORT void
geqp3<c_type>(i_type m, i_type n, c_type *a, i_type lda, i_type* jpvt, c_type *tau, c_type *work, 
              i_type lwork, i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(cgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
            _rc(tau),_rc(work),_rc(&lwork), (s_type*)nullptr, _rc(info));
        work[0] = work[0] + c_type(s_type(1 * n));
    }
    else
    {
        s_type* rwork   = reinterpret_cast<s_type*>(work);
        work            = work + n;
        lwork           = lwork - n;

        LAPACK_NAME(cgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
            _rc(tau),_rc(work),_rc(&lwork), _rc(rwork), _rc(info));
        work[0] = work[0] + c_type(s_type(1 * n));
    };
};

template<> BLAS_EXPORT void
geqp3<z_type>(i_type m, i_type n, z_type *a, i_type lda, i_type* jpvt, z_type *tau, z_type *work, 
              i_type lwork, i_type *info)
{
    if (lwork == -1)
    {
        LAPACK_NAME(zgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
            _rc(tau),_rc(work),_rc(&lwork), (d_type*)nullptr, _rc(info));
        work[0] = work[0] + z_type(d_type(1 * n));
    }
    else
    {
        d_type* rwork   = reinterpret_cast<d_type*>(work);
        work            = work + n;
        lwork           = lwork - n;

        LAPACK_NAME(zgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
            _rc(tau),_rc(work),_rc(&lwork), _rc(rwork), _rc(info));
        work[0] = work[0] + z_type(d_type(1 * n));
    };
};

BLAS_EXPORT INLINE_TYPE 
void cgeqp3(i_type m,i_type n,c_type *a,i_type lda,i_type* jpvt,c_type *tau,
                           c_type *work,i_type lwork,s_type* rwork, i_type *info)
{
    LAPACK_NAME(cgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
        _rc(tau),_rc(work),_rc(&lwork),_rc(rwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sgeqp3(i_type m,i_type n,s_type *a,i_type lda,i_type* jpvt,s_type *tau,
                           s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
        _rc(tau),_rc(work),_rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dgeqp3(i_type m,i_type n,d_type *a,i_type lda,i_type* jpvt,d_type *tau,
                           d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
        _rc(tau),_rc(work),_rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zgeqp3(i_type m,i_type n,z_type *a,i_type lda,i_type* jpvt,z_type *tau,
                           z_type *work,i_type lwork,d_type* rwork, i_type *info)
{
    LAPACK_NAME(zgeqp3)(_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(jpvt),
        _rc(tau),_rc(work),_rc(&lwork),_rc(rwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          ORGQR
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void orgqr<d_type>(i_type m, i_type n, i_type k, d_type *a, i_type lda, const d_type *tau, 
                   d_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(dorgqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void orgqr<s_type>(i_type m, i_type n, i_type k, s_type *a, i_type lda, const s_type *tau,
                   s_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(sorgqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void orgqr<c_type>(i_type m, i_type n, i_type k, c_type *a, i_type lda, const c_type *tau,
                   c_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(cungqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void orgqr<z_type>(i_type m, i_type n, i_type k, z_type *a, i_type lda, const z_type *tau, 
                   z_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(zungqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dorgqr(i_type m, i_type n, i_type k, d_type *a, i_type lda, const d_type *tau, d_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(dorgqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sorgqr(i_type m, i_type n, i_type k, s_type *a, i_type lda, const s_type *tau, s_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(sorgqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          UNGQR
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void ungqr<z_type>(i_type m, i_type n, i_type k, z_type *a, i_type lda, const z_type *tau,
                   z_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(zungqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zungqr(i_type m, i_type n, i_type k, z_type *a, i_type lda, const z_type *tau, z_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(zungqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void ungqr<c_type>(i_type m, i_type n, i_type k, c_type *a, i_type lda, const c_type *tau, 
                   c_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(cungqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void cungqr(i_type m, i_type n, i_type k, c_type *a, i_type lda, const c_type *tau, 
            c_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(cungqr)(_rc(&m),_rc(&n),_rc(&k), _rc(a),_rc(&lda),_rc(tau),_rc(work),
                        _rc(&lwork),_rc(info));
};

//-----------------------------------------------------------------------
//                          GEHRD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void gehrd<c_type>(i_type n,i_type ilo, i_type ihi, c_type *a, 
                    i_type lda, c_type *tau,c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(cgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void gehrd<s_type>(i_type n,i_type ilo, i_type ihi, s_type *a, 
                    i_type lda, s_type *tau,s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void gehrd<d_type>(i_type n,i_type ilo, i_type ihi, d_type *a, 
                    i_type lda, d_type *tau,d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void gehrd<z_type>(i_type n,i_type ilo, i_type ihi, z_type *a, 
                    i_type lda, z_type *tau,z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(zgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void cgehrd(i_type n,i_type ilo, i_type ihi, c_type *a,i_type lda,c_type *tau,
            c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(cgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void sgehrd(i_type n,i_type ilo, i_type ihi, s_type *a,i_type lda,s_type *tau,
            s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void dgehrd(i_type n,i_type ilo, i_type ihi, d_type *a,i_type lda,d_type *tau,
            d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void zgehrd(i_type n,i_type ilo, i_type ihi, z_type *a,i_type lda,z_type *tau,
            z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(zgehrd)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

//-----------------------------------------------------------------------
//                          SYTRD/HETRD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void sytrd<s_type>(const char* uplo, i_type n,s_type *a, i_type lda, s_type *d, s_type *e,
        s_type *tau,s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(ssytrd)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(d), _rc(e),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void sytrd<d_type>(const char* uplo, i_type n,d_type *a, i_type lda, d_type *d, d_type *e,
        d_type *tau,d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dsytrd)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(d), _rc(e),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void ssytrd(const char* uplo, i_type n, s_type *a,i_type lda,
                           s_type *d, s_type *e, s_type *tau,
                           s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(ssytrd)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(d), _rc(e),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void dsytrd(const char* uplo, i_type n, d_type *a,i_type lda,
                           d_type *d, d_type *e, d_type *tau,
                           d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dsytrd)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(d), _rc(e),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void chetrd(const char* uplo, i_type n, c_type *a,i_type lda,
                           s_type *d, s_type *e, c_type *tau,
                           c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(chetrd)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(d), _rc(e),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void zhetrd(const char* uplo, i_type n, z_type *a,i_type lda,
                           d_type *d, d_type *e, z_type *tau,
                           z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(zhetrd)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(d), _rc(e),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

//-----------------------------------------------------------------------
//                          UNGTR/ORGTR
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void orgtr<s_type>(const char* uplo, i_type n,s_type *a, i_type lda,
        s_type *tau,s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sorgtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), 
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void orgtr<d_type>(const char* uplo, i_type n,d_type *a, i_type lda, 
        d_type *tau,d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dorgtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void sorgtr(const char* uplo, i_type n, s_type *a,i_type lda,
                           s_type *tau, s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sorgtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void dorgtr(const char* uplo, i_type n, d_type *a,i_type lda,
                           d_type *tau, d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dorgtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), 
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void ungtr<c_type>(const char* uplo, i_type n,c_type *a, i_type lda, 
        c_type *tau,c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(cungtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void ungtr<z_type>(const char* uplo, i_type n,z_type *a, i_type lda, 
        z_type *tau,z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(zungtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void cungtr(const char* uplo, i_type n, c_type *a,i_type lda,
                           c_type *tau,
                           c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(cungtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE void zungtr(const char* uplo, i_type n, z_type *a,i_type lda,
                           z_type *tau,
                           z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(zungtr)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda),
        _rc(tau), _rc(work), _rc(&lwork), _rc(info));
};

//-----------------------------------------------------------------------
//                          ORGHR
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void orghr<d_type>(i_type n,i_type ilo, i_type ihi, d_type *a, 
                    i_type lda, d_type *tau,d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dorghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void dorghr(i_type n,i_type ilo, i_type ihi, d_type *a, 
            i_type lda, d_type *tau,d_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(dorghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void orghr<s_type>(i_type n,i_type ilo, i_type ihi, s_type *a,i_type lda,s_type *tau,
                    s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sorghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void sorghr(i_type n,i_type ilo, i_type ihi, s_type *a,i_type lda,s_type *tau,
            s_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(sorghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

//-----------------------------------------------------------------------
//                          UNGHR
//-----------------------------------------------------------------------
template<> BLAS_EXPORT INLINE_TYPE 
void unghr<z_type>(i_type n,i_type ilo, i_type ihi, z_type *a,i_type lda,z_type *tau,
                    z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(zunghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void zunghr(i_type n,i_type ilo, i_type ihi, z_type *a,i_type lda,z_type *tau,
            z_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(zunghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE 
void unghr<c_type>( i_type n,i_type ilo, i_type ihi, c_type *a, 
                    i_type lda, c_type *tau,c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(cunghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

BLAS_EXPORT INLINE_TYPE 
void cunghr(i_type n,i_type ilo, i_type ihi, c_type *a, 
            i_type lda, c_type *tau,c_type *work,i_type lwork,i_type *info)
{
    LAPACK_NAME(cunghr)(_rc(&n), _rc(&ilo), _rc(&ihi), _rc(a), _rc(&lda), _rc (tau),
        _rc(work), _rc(&lwork), _rc(info));
};

//-----------------------------------------------------------------------
//                          POTRF
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void potrf<c_type>(const char *uplo, i_type n, c_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(cpotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void potrf<s_type>(const char *uplo, i_type n, s_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(spotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void potrf<d_type>(const char *uplo, i_type n, d_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(dpotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void potrf<z_type>(const char *uplo, i_type n, z_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(zpotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void cpotrf(const char *uplo, i_type n, c_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(cpotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void spotrf(const char *uplo, i_type n, s_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(spotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void dpotrf(const char *uplo, i_type n, d_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(dpotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void zpotrf(const char *uplo, i_type n, z_type *a, i_type lda, i_type *info)
{
	LAPACK_NAME(zpotrf)(_rc(uplo), _rc(&n), _rc(a), _rc(&lda), _rc(info));
};

//-----------------------------------------------------------------------
//                          PBTRF
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void pbtrf<c_type>(const char *uplo, i_type n, i_type kd, c_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(cpbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void pbtrf<s_type>(const char *uplo, i_type n, i_type kd, s_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(spbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void pbtrf<d_type>(const char *uplo, i_type n, i_type kd, d_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(dpbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

template<> BLAS_EXPORT INLINE_TYPE
void pbtrf<z_type>(const char *uplo, i_type n, i_type kd, z_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(zpbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void cpbtrf(const char *uplo, i_type n, i_type kd, c_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(cpbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void spbtrf(const char *uplo, i_type n, i_type kd, s_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(spbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void dpbtrf(const char *uplo, i_type n, i_type kd, d_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(dpbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

BLAS_EXPORT INLINE_TYPE
void zpbtrf(const char *uplo, i_type n, i_type kd, z_type *ab, i_type ldab, i_type *info)
{
	LAPACK_NAME(zpbtrf)(_rc(uplo), _rc(&n), _rc(&kd), _rc(ab), _rc(&ldab), _rc(info));
};

//-----------------------------------------------------------------------
//                          SYTRF/HETRF
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void sytrf<s_type>(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void sytrf<d_type>(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void sytrf<c_type>(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(csytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void sytrf<z_type>(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(zsytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ssytrf(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dsytrf(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void csytrf(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(csytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void zsytrf(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(zsytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf<s_type>(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
                    i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf<d_type>(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
                    i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf<c_type>(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
                i_type lwork, i_type *info)
{
    LAPACK_NAME(chetrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void shetrf(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dhetrf(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf<z_type>(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
                i_type lwork, i_type *info)
{
    LAPACK_NAME(zhetrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void chetrf(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(chetrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void zhetrf(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(zhetrf)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

//-----------------------------------------------------------------------
//                          SYTRF_ROOK/HETRF_ROOK
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE
void sytrf_rook<s_type>(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void sytrf_rook<d_type>(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void sytrf_rook<c_type>(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(csytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void sytrf_rook<z_type>(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
                      i_type lwork, i_type *info)
{
    LAPACK_NAME(zsytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void ssytrf_rook(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dsytrf_rook(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void csytrf_rook(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(csytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void zsytrf_rook(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(zsytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf_rook<s_type>(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
                    i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf_rook<d_type>(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
                    i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf_rook<c_type>(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
                i_type lwork, i_type *info)
{
    LAPACK_NAME(chetrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void shetrf_rook(const char *uplo, i_type n, s_type *a, i_type lda, i_type *ipiv, s_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(ssytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dhetrf_rook(const char *uplo, i_type n, d_type *a, i_type lda, i_type *ipiv, d_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(dsytrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

template<> BLAS_EXPORT INLINE_TYPE
void hetrf_rook<z_type>(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
                i_type lwork, i_type *info)
{
    LAPACK_NAME(zhetrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void chetrf_rook(const char *uplo, i_type n, c_type *a, i_type lda, i_type *ipiv, c_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(chetrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void zhetrf_rook(const char *uplo, i_type n, z_type *a, i_type lda, i_type *ipiv, z_type *work, 
            i_type lwork, i_type *info)
{
    LAPACK_NAME(zhetrf_rook)(_rc(uplo),_rc(&n),_rc(a),_rc(&lda),_rc(ipiv),_rc(work),_rc(&lwork),_rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dgeqr2(i_type m, i_type n, d_type *a, i_type lda, d_type *tau, d_type *work, i_type *info)
{
    LAPACK_NAME(dgeqr2)(_rc(&m), _rc(&n), _rc(a), _rc(&lda), _rc(tau), _rc(work), _rc(info));
}

BLAS_EXPORT INLINE_TYPE 
void dlaset( const char *uplo, i_type m, i_type n, d_type *alpha, d_type *beta, d_type *a, i_type lda)
{
    LAPACK_NAME(dlaset)(_rc(uplo), _rc(&m), _rc(&n), _rc(alpha), _rc(beta), _rc(a),
        _rc(&lda));
}

BLAS_EXPORT INLINE_TYPE 
void dormqr( const char *side, const char *trans, i_type m, i_type n, 
            i_type k, d_type *a, i_type lda, d_type *tau, d_type *c,
            i_type ldc, d_type *work, i_type lwork, i_type *info)
{
    LAPACK_NAME(dormqr)(_rc(side), _rc(trans), _rc(&m), _rc(&n), _rc(&k), _rc(a), _rc(&lda),
        _rc(tau), _rc(c), _rc(&ldc), _rc(work), _rc(&lwork), _rc(info));
}

//-----------------------------------------------------------------------
//                          laln2
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
laln2<s_type>(bool LTRANS, i_type NA, i_type NW, s_type SMIN, s_type CA, const s_type* A, 
              i_type LDA, s_type D1, s_type D2, const s_type* B, i_type LDB, s_type WR,
              s_type WI, s_type* X, i_type LDX, s_type* SCALE, s_type* XNORM, i_type* INFO )
{
    i_type ILTRANS = LTRANS;
    LAPACK_NAME(slaln2)(&ILTRANS, &NA, &NW, &SMIN, &CA, const_cast<s_type*>(A), &LDA, 
                        &D1, &D2, const_cast<s_type*>(B), &LDB, &WR, &WI, X, &LDX, SCALE, XNORM, INFO);
};

template<> BLAS_EXPORT void
laln2<d_type>(bool LTRANS, i_type NA, i_type NW, d_type SMIN, d_type CA, const d_type* A, 
              i_type LDA, d_type D1, d_type D2, const d_type* B, i_type LDB, d_type WR, d_type WI,
              d_type* X, i_type LDX, d_type* SCALE, d_type* XNORM, i_type* INFO )
{
    i_type ILTRANS = LTRANS;
    LAPACK_NAME(dlaln2)(&ILTRANS, &NA, &NW, &SMIN, &CA, const_cast<d_type*>(A), &LDA, 
                        &D1, &D2, const_cast<d_type*>(B), &LDB, &WR, &WI, X, &LDX, SCALE, XNORM, INFO);
};

BLAS_EXPORT void slaln2(bool LTRANS, i_type NA, i_type NW, s_type SMIN, s_type CA, const s_type* A, 
                        i_type LDA, s_type D1, s_type D2, const s_type* B, i_type LDB, s_type WR, 
                        s_type WI, s_type* X, i_type LDX, s_type* SCALE, s_type* XNORM, i_type* INFO )
{
    i_type ILTRANS = LTRANS;
    LAPACK_NAME(slaln2)(&ILTRANS, &NA, &NW, &SMIN, &CA, const_cast<s_type*>(A), &LDA, 
                        &D1, &D2, const_cast<s_type*>(B), &LDB, &WR, &WI, X, &LDX, SCALE, XNORM, INFO);
};

BLAS_EXPORT void dlaln2(bool LTRANS, i_type NA, i_type NW, d_type SMIN, d_type CA, const d_type* A, 
                        i_type LDA, d_type D1, d_type D2, const d_type* B, i_type LDB, d_type WR, 
                        d_type WI, d_type* X, i_type LDX, d_type* SCALE, d_type* XNORM, i_type* INFO )
{
    i_type ILTRANS = LTRANS;
    LAPACK_NAME(dlaln2)(&ILTRANS, &NA, &NW, &SMIN, &CA, const_cast<d_type*>(A), &LDA, 
                        &D1, &D2, const_cast<d_type*>(B), &LDB, &WR, &WI, X, &LDX, SCALE, XNORM, INFO);
};

//-----------------------------------------------------------------------
//                          DLAGV2
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lagv2<s_type>(s_type* A, i_type LDA, s_type* B, i_type LDB, s_type* ALPHAR, s_type* ALPHAI, 
                s_type* BETA, s_type& CSL, s_type& SNL, s_type& CSR, s_type& SNR )
{
    LAPACK_NAME(slagv2)(A,&LDA, B, &LDB, ALPHAR, ALPHAI, BETA, &CSL, &SNL, &CSR, &SNR);
};

template<> BLAS_EXPORT void
lagv2<d_type>(d_type* A, i_type LDA, d_type* B, i_type LDB, d_type* ALPHAR, d_type* ALPHAI, 
                d_type* BETA, d_type& CSL, d_type& SNL, d_type& CSR, d_type& SNR )
{
    LAPACK_NAME(dlagv2)(A,&LDA, B, &LDB, ALPHAR, ALPHAI, BETA, &CSL, &SNL, &CSR, &SNR);
};

//-----------------------------------------------------------------------
//                          DSTEVD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT
void
stevd<d_type>(const char * JOBZ, i_type N, d_type* D, d_type* E, d_type* Z, i_type LDZ, d_type* WORK, 
      i_type LWORK, i_type* IWORK, i_type LIWORK, i_type* INFO )
{
    LAPACK_NAME(dstevd)(_rc(JOBZ), &N, D, E, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, INFO);
};

template<> BLAS_EXPORT
void
stevd<s_type>(const char * JOBZ, i_type N, s_type* D, s_type* E, s_type* Z, i_type LDZ, s_type* WORK, 
      i_type LWORK, i_type* IWORK, i_type LIWORK, i_type* INFO )
{
    LAPACK_NAME(sstevd)(_rc(JOBZ), &N, D, E, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, INFO);
};

BLAS_EXPORT void sstevd(const char * JOBZ, i_type N, s_type* D, s_type* E, s_type* Z, 
                        i_type LDZ, s_type* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK,
                        i_type* INFO )
{
    LAPACK_NAME(sstevd)(_rc(JOBZ), &N, D, E, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, INFO);
};

BLAS_EXPORT void dstevd(const char * JOBZ, i_type N, d_type* D, d_type* E, d_type* Z, 
                        i_type LDZ, d_type* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK,
                        i_type* INFO )
{
    LAPACK_NAME(dstevd)(_rc(JOBZ), &N, D, E, Z, &LDZ, WORK, &LWORK, IWORK, &LIWORK, INFO);
};

//-----------------------------------------------------------------------
//                          DLARFT
//-----------------------------------------------------------------------

template<> BLAS_EXPORT
void
larft<d_type>(const char *direct, const char *storev, i_type n, i_type k, 
        const d_type *v, i_type ldv, const d_type* tau, d_type* t, i_type ldt)
{
    LAPACK_NAME(dlarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), t, &ldt);
};

template<> BLAS_EXPORT
void
larft<s_type>(const char *direct, const char *storev, i_type n, i_type k, 
        const s_type *v, i_type ldv, const s_type* tau, s_type* t, i_type ldt)
{
    LAPACK_NAME(slarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), t, &ldt);
};

template<> BLAS_EXPORT
void
larft<c_type>(const char *direct, const char *storev, i_type n, i_type k, 
        const c_type *v, i_type ldv, const c_type* tau, c_type* t, i_type ldt)
{
    LAPACK_NAME(clarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), _rc(t), &ldt);
};

template<> BLAS_EXPORT
void
larft<z_type>(const char *direct, const char *storev, i_type n, i_type k, 
        const z_type *v, i_type ldv, const z_type* tau, z_type* t, i_type ldt)
{
    LAPACK_NAME(zlarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), _rc(t), &ldt);
};

BLAS_EXPORT void
lapack::dlarft(const char *direct, const char *storev, i_type n, i_type k, 
        const d_type *v, i_type ldv, const d_type* tau, d_type* t, i_type ldt)
{
    LAPACK_NAME(dlarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), t, &ldt);
};

BLAS_EXPORT void
lapack::slarft(const char *direct, const char *storev, i_type n, i_type k, 
        const s_type *v, i_type ldv, const s_type* tau, s_type* t, i_type ldt)
{
    LAPACK_NAME(slarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), t, &ldt);
};

BLAS_EXPORT void
lapack::clarft(const char *direct, const char *storev, i_type n, i_type k, 
        const c_type *v, i_type ldv, const c_type* tau, c_type* t, i_type ldt)
{
    LAPACK_NAME(clarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), _rc(t), &ldt);
};

BLAS_EXPORT void
lapack::zlarft(const char *direct, const char *storev, i_type n, i_type k, 
        const z_type *v, i_type ldv, const z_type* tau, z_type* t, i_type ldt)
{
    LAPACK_NAME(zlarft)(_rc(direct), _rc(storev), &n, &k, _rc(v), &ldv, _rc(tau), _rc(t), &ldt);
};

//-----------------------------------------------------------------------
//                          DLARFB
//-----------------------------------------------------------------------
template<> BLAS_EXPORT
void
larfb<d_type>(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, const d_type* v, i_type ldv,  const d_type* t, i_type ldt, 
    d_type* c__, i_type ldc, d_type* work, i_type ldwork)
{
    LAPACK_NAME(dlarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, _rc(v), &ldv, 
                        _rc(t), &ldt, c__, &ldc, work, &ldwork);
};

template<> BLAS_EXPORT
void
larfb<s_type>(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, const s_type* v, i_type ldv,  const s_type* t, i_type ldt, 
    s_type* c__, i_type ldc, s_type* work, i_type ldwork)
{
    LAPACK_NAME(slarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, _rc(v), &ldv, 
                        _rc(t), &ldt, c__, &ldc, work, &ldwork);
};

template<> BLAS_EXPORT
void
larfb<c_type>(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, const c_type* v, i_type ldv,  const c_type* t, i_type ldt, 
    c_type* c__, i_type ldc, c_type* work, i_type ldwork)
{
    LAPACK_NAME(clarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, _rc(v), &ldv, 
                        _rc(t), &ldt, _rc(c__), &ldc, _rc(work), &ldwork);
};

template<> BLAS_EXPORT
void
larfb<z_type>(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, const z_type* v, i_type ldv,  const z_type* t, i_type ldt, 
    z_type* c__, i_type ldc, z_type* work, i_type ldwork)
{
    LAPACK_NAME(zlarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, _rc(v), &ldv, 
                        _rc(t), &ldt, _rc(c__), &ldc, _rc(work), &ldwork);
};

BLAS_EXPORT void
dlarfb(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, d_type* v, i_type ldv,  d_type* t, i_type ldt, 
    d_type* c__, i_type ldc, d_type* work, i_type ldwork)
{
    LAPACK_NAME(dlarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, v, &ldv, 
                        t, &ldt, c__, &ldc, work, &ldwork);
};

BLAS_EXPORT void
slarfb(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, s_type* v, i_type ldv,  s_type* t, i_type ldt, 
    s_type* c__, i_type ldc, s_type* work, i_type ldwork)
{
    LAPACK_NAME(slarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, v, &ldv, 
                        t, &ldt, c__, &ldc, work, &ldwork);
};

BLAS_EXPORT void
clarfb(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, c_type* v, i_type ldv,  c_type* t, i_type ldt, 
    c_type* c__, i_type ldc, c_type* work, i_type ldwork)
{
    LAPACK_NAME(clarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, _rc(v), &ldv, 
                        _rc(t), &ldt, _rc(c__), &ldc, _rc(work), &ldwork);
};

BLAS_EXPORT void
zlarfb(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, z_type* v, i_type ldv,  z_type* t, i_type ldt, 
    z_type* c__, i_type ldc, z_type* work, i_type ldwork)
{
    LAPACK_NAME(zlarfb)(_rc(side), _rc(trans), _rc(direct), _rc(storev), &m, &n, &k, _rc(v), &ldv, 
                        _rc(t), &ldt, _rc(c__), &ldc, _rc(work), &ldwork);
};

//-----------------------------------------------------------------------
//                          DGEQR2
//-----------------------------------------------------------------------

template<> BLAS_EXPORT
void
geqr2<d_type>(i_type m, i_type n, d_type* a, i_type lda, d_type* tau, d_type* work, i_type& info)
{
    LAPACK_NAME(dgeqr2)(&m, &n, a, &lda, tau, work, &info);
};

template<> BLAS_EXPORT
void
geqr2<s_type>(i_type m, i_type n, s_type* a, i_type lda, s_type* tau, s_type* work, i_type& info)
{
    LAPACK_NAME(sgeqr2)(&m, &n, a, &lda, tau, work, &info);
};

template<> BLAS_EXPORT
void
geqr2<c_type>(i_type m, i_type n, c_type* a, i_type lda, c_type* tau, c_type* work, i_type& info)
{
    LAPACK_NAME(cgeqr2)(&m, &n, _rc(a), &lda, _rc(tau), _rc(work), &info);
};

template<> BLAS_EXPORT
void
geqr2<z_type>(i_type m, i_type n, z_type* a, i_type lda, z_type* tau, z_type* work, i_type& info)
{
    LAPACK_NAME(zgeqr2)(&m, &n, _rc(a), &lda, _rc(tau), _rc(work), &info);
};

BLAS_EXPORT void
dgeqr2(i_type m, i_type n, d_type* a, i_type lda, d_type* tau, d_type* work, i_type& info)
{
    LAPACK_NAME(dgeqr2)(&m, &n, a, &lda, tau, work, &info);
};

BLAS_EXPORT void
sgeqr2(i_type m, i_type n, s_type* a, i_type lda, s_type* tau, s_type* work, i_type& info)
{
    LAPACK_NAME(sgeqr2)(&m, &n, a, &lda, tau, work, &info);
};

BLAS_EXPORT void
cgeqr2(i_type m, i_type n, c_type* a, i_type lda, c_type* tau, c_type* work, i_type& info)
{
    LAPACK_NAME(cgeqr2)(&m, &n, _rc(a), &lda, _rc(tau), _rc(work), &info);
};

BLAS_EXPORT void
zgeqr2(i_type m, i_type n, z_type* a, i_type lda, z_type* tau, z_type* work, i_type& info)
{
    LAPACK_NAME(zgeqr2)(&m, &n, _rc(a), &lda, _rc(tau), _rc(work), &info);
};

//-----------------------------------------------------------------------
//                          DLARFG
//-----------------------------------------------------------------------

template<> BLAS_EXPORT
void
larfg<d_type>(i_type n, d_type* alpha, d_type* x, i_type incx, d_type* tau)
{
    LAPACK_NAME(dlarfg)(&n, alpha, x, &incx, tau);
};

template<> BLAS_EXPORT
void
larfg<s_type>(i_type n, s_type* alpha, s_type* x, i_type incx, s_type* tau)
{
    LAPACK_NAME(slarfg)(&n, alpha, x, &incx, tau);
};

template<> BLAS_EXPORT
void
larfg<c_type>(i_type n, c_type* alpha, c_type* x, i_type incx, c_type* tau)
{
    LAPACK_NAME(clarfg)(&n, _rc(alpha), _rc(x), &incx, _rc(tau));
};

template<> BLAS_EXPORT
void
larfg<z_type>(i_type n, z_type* alpha, z_type* x, i_type incx, z_type* tau)
{
    LAPACK_NAME(zlarfg)(&n, _rc(alpha), _rc(x), &incx, _rc(tau));
};

BLAS_EXPORT void
dlarfg(i_type n, d_type* alpha, d_type* x, i_type incx, d_type* tau)
{
    LAPACK_NAME(dlarfg)(&n, alpha, x, &incx, tau);
};

BLAS_EXPORT void
slarfg(i_type n, s_type* alpha, s_type* x, i_type incx, s_type* tau)
{
    LAPACK_NAME(slarfg)(&n, alpha, x, &incx, tau);
};

BLAS_EXPORT void
clarfg(i_type n, c_type* alpha, c_type* x, i_type incx, c_type* tau)
{
    LAPACK_NAME(clarfg)(&n, _rc(alpha), _rc(x), &incx, _rc(tau));
};

BLAS_EXPORT void
zlarfg(i_type n, z_type* alpha, z_type* x, i_type incx, z_type* tau)
{
    LAPACK_NAME(zlarfg)(&n, _rc(alpha), _rc(x), &incx, _rc(tau));
};

//-----------------------------------------------------------------------
//                          DLARF
//-----------------------------------------------------------------------

template<> BLAS_EXPORT
void
larf<d_type>(const char *side, i_type m, i_type n, d_type* v, i_type incv, const d_type* tau, 
             d_type* c__, i_type ldc, d_type* work)
{
    LAPACK_NAME(dlarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

template<> BLAS_EXPORT
void
larf<s_type>(const char *side, i_type m, i_type n, s_type* v, i_type incv, const s_type* tau, 
             s_type* c__, i_type ldc, s_type* work)
{
    LAPACK_NAME(slarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

template<> BLAS_EXPORT
void
larf<c_type>(const char *side, i_type m, i_type n, c_type* v, i_type incv, const c_type* tau, 
             c_type* c__, i_type ldc, c_type* work)
{
    LAPACK_NAME(clarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

template<> BLAS_EXPORT
void
larf<z_type>(const char *side, i_type m, i_type n, z_type* v, i_type incv, const z_type* tau, 
             z_type* c__, i_type ldc, z_type* work)
{
    LAPACK_NAME(zlarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

BLAS_EXPORT void
dlarf(const char *side, i_type m, i_type n, d_type* v, i_type incv, const d_type* tau, 
             d_type* c__, i_type ldc, d_type* work)
{
    LAPACK_NAME(dlarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

BLAS_EXPORT void
slarf(const char *side, i_type m, i_type n, s_type* v, i_type incv, const s_type* tau, 
             s_type* c__, i_type ldc, s_type* work)
{
    LAPACK_NAME(slarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

BLAS_EXPORT void
clarf(const char *side, i_type m, i_type n, c_type* v, i_type incv, const c_type* tau, 
             c_type* c__, i_type ldc, c_type* work)
{
    LAPACK_NAME(clarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

BLAS_EXPORT void
zlarf(const char *side, i_type m, i_type n, z_type* v, i_type incv, const z_type* tau, 
             z_type* c__, i_type ldc, z_type* work)
{
    LAPACK_NAME(zlarf)(_rc(side), &m, &n, _rc(v), &incv, _rc(tau), _rc(c__), &ldc, _rc(work));
};

//-----------------------------------------------------------------------
//                          DGEBRD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
gebrd<d_type>(i_type M, i_type N, d_type* A, i_type LDA, d_type* D, d_type* E, d_type* tauq, 
              d_type* taup, d_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(dgebrd)(&M, &N, A, &LDA, D, E, tauq, taup, WORK, &lwork, &info);
};

template<> BLAS_EXPORT void
gebrd<s_type>(i_type M, i_type N, s_type* A, i_type LDA, s_type* D, s_type* E, s_type* tauq, 
              s_type* taup, s_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(sgebrd)(&M, &N, A, &LDA, D, E, tauq, taup, WORK, &lwork, &info);
};

template<> BLAS_EXPORT void
gebrd<c_type>(i_type M, i_type N, c_type* A, i_type LDA, s_type* D, s_type* E, c_type* tauq, 
              c_type* taup, c_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(cgebrd)(&M, &N, _rc(A), &LDA, _rc(D), _rc(E), _rc(tauq), _rc(taup), _rc(WORK),
                        &lwork, &info);
};

template<> BLAS_EXPORT void
gebrd<z_type>(i_type M, i_type N, z_type* A, i_type LDA, d_type* D, d_type* E, z_type* tauq, 
              z_type* taup, z_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(zgebrd)(&M, &N, _rc(A), &LDA, _rc(D), _rc(E), _rc(tauq), _rc(taup),
                        _rc(WORK), &lwork, &info);
};

BLAS_EXPORT void 
lapack::dgebrd(i_type M, i_type N, d_type* A, i_type LDA, d_type* D, d_type* E, d_type* tauq, 
       d_type* taup, d_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(dgebrd)(&M, &N, _rc(A), &LDA, _rc(D), _rc(E), _rc(tauq), _rc(taup), 
                        _rc(WORK), &lwork, &info);
};

BLAS_EXPORT void 
lapack::sgebrd(i_type M, i_type N, s_type* A, i_type LDA, s_type* D, s_type* E, s_type* tauq, 
              s_type* taup, s_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(sgebrd)(&M, &N, A, &LDA, D, E, tauq, taup, WORK, &lwork, &info);
};

BLAS_EXPORT void 
lapack::cgebrd(i_type M, i_type N, c_type* A, i_type LDA, s_type* D, s_type* E, c_type* tauq, 
              c_type* taup, c_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(cgebrd)(&M, &N, _rc(A), &LDA, _rc(D), _rc(E), _rc(tauq), _rc(taup), _rc(WORK), 
                        &lwork, &info);
};

BLAS_EXPORT void 
lapack::zgebrd(i_type M, i_type N, z_type* A, i_type LDA, d_type* D, d_type* E, z_type* tauq, 
              z_type* taup, z_type* WORK, i_type lwork, i_type& info)
{
    LAPACK_NAME(zgebrd)(&M, &N, _rc(A), &LDA, _rc(D), _rc(E), _rc(tauq), _rc(taup), _rc(WORK),
                        &lwork, &info);
};

//-----------------------------------------------------------------------
//                          DGGHRD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::gghrd(const char* COMPQ, const char* COMPZ, i_type N, i_type ILO, i_type IHI, d_type* A, 
      i_type LDA, d_type* B, i_type LDB, d_type* Q, i_type LDQ, d_type* Z, i_type LDZ, i_type& INFO )
{
    LAPACK_NAME(dgghrd)(_rc(COMPQ), _rc(COMPZ), &N, &ILO, &IHI, _rc(A), &LDA, _rc(B), &LDB,
                        _rc(Q), &LDQ, _rc(Z), &LDZ, &INFO );
};

template<> BLAS_EXPORT void
lapack::gghrd(const char* COMPQ, const char* COMPZ, i_type N, i_type ILO, i_type IHI, s_type* A, 
      i_type LDA, s_type* B, i_type LDB, s_type* Q, i_type LDQ, s_type* Z, i_type LDZ, i_type& INFO )
{
    LAPACK_NAME(sgghrd)(_rc(COMPQ), _rc(COMPZ), &N, &ILO, &IHI, _rc(A), &LDA, _rc(B), &LDB, 
                        _rc(Q), &LDQ, _rc(Z), &LDZ, &INFO );
};

template<> BLAS_EXPORT void
lapack::gghrd(const char* COMPQ, const char* COMPZ, i_type N, i_type ILO, i_type IHI, c_type* A, 
      i_type LDA, c_type* B, i_type LDB, c_type* Q, i_type LDQ, c_type* Z, i_type LDZ, i_type& INFO )
{
    LAPACK_NAME(cgghrd)(_rc(COMPQ), _rc(COMPZ), &N, &ILO, &IHI, _rc(A), &LDA, _rc(B), &LDB, _rc(Q), 
                        &LDQ, _rc(Z), &LDZ, &INFO );
};

template<> BLAS_EXPORT void
lapack::gghrd(const char* COMPQ, const char* COMPZ, i_type N, i_type ILO, i_type IHI, z_type* A, 
      i_type LDA, z_type* B, i_type LDB, z_type* Q, i_type LDQ, z_type* Z, i_type LDZ, i_type& INFO )
{
    LAPACK_NAME(zgghrd)(_rc(COMPQ), _rc(COMPZ), &N, &ILO, &IHI, _rc(A), &LDA, _rc(B), &LDB, _rc(Q),
                        &LDQ, _rc(Z), &LDZ, &INFO );
};

//-----------------------------------------------------------------------
//                          sbtrd
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::sbtrd<d_type>(const char *VECT, const char *UPLO, i_type N, i_type KD, d_type *AB, i_type LDAB, 
      d_type* D, d_type* E, d_type* Q, i_type LDQ, d_type* work, i_type& info)
{
    LAPACK_NAME(dsbtrd)(_rc(VECT), _rc(UPLO), &N, &KD, _rc(AB), &LDAB, _rc(D), _rc(E), _rc(Q), &LDQ, 
                        _rc(work), &info);
};

template<> BLAS_EXPORT void
lapack::sbtrd<s_type>(const char *VECT, const char *UPLO, i_type N, i_type KD, s_type *AB, i_type LDAB, 
      s_type* D, s_type* E, s_type* Q, i_type LDQ, s_type* work, i_type& info)
{
    LAPACK_NAME(ssbtrd)(_rc(VECT), _rc(UPLO), &N, &KD, _rc(AB), &LDAB, _rc(D), _rc(E), _rc(Q), &LDQ, 
                        _rc(work), &info);
};

template<> BLAS_EXPORT void
lapack::sbtrd<c_type>(const char *VECT, const char *UPLO, i_type N, i_type KD, c_type *AB, i_type LDAB, 
      s_type* D, s_type* E, c_type* Q, i_type LDQ, c_type* work, i_type& info)
{
    LAPACK_NAME(chbtrd)(_rc(VECT), _rc(UPLO), &N, &KD, _rc(AB), &LDAB, _rc(D), _rc(E), _rc(Q), &LDQ, 
                        _rc(work), &info);
};

template<> BLAS_EXPORT void
lapack::sbtrd<z_type>(const char *VECT, const char *UPLO, i_type N, i_type KD, z_type *AB, i_type LDAB, 
      d_type* D, d_type* E, z_type* Q, i_type LDQ, z_type* work, i_type& info)
{
    LAPACK_NAME(zhbtrd)(_rc(VECT), _rc(UPLO), &N, &KD, _rc(AB), &LDAB, _rc(D), _rc(E), _rc(Q), &LDQ, 
                        _rc(work), &info);
};

//-----------------------------------------------------------------------
//                          DSYEVD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::syevd<d_type>(const char *jobv,const char *uplo, i_type n, d_type *a,i_type lda, d_type* w, 
            d_type *work, i_type lwork, i_type* iwork, i_type liwork, i_type *info)
{
    LAPACK_NAME(dsyevd)(_rc(jobv), _rc(uplo), &n, _rc(a), &lda, _rc(w), _rc(work), &lwork, 
                        iwork, &liwork, info);
};

template<> BLAS_EXPORT void
lapack::syevd<s_type>(const char *jobv,const char *uplo, i_type n, s_type *a,i_type lda, s_type* w, 
            s_type *work, i_type lwork, i_type* iwork, i_type liwork, i_type *info)
{
    LAPACK_NAME(ssyevd)(_rc(jobv), _rc(uplo), &n, _rc(a), &lda, _rc(w), _rc(work), &lwork, 
                        iwork, &liwork, info);
};

template<> BLAS_EXPORT void
lapack::heevd<c_type>(const char *jobv,const char *uplo, i_type n, c_type *a,i_type lda, s_type* w, 
                      c_type *work,i_type lwork, s_type* rwork, i_type lrwork, i_type* iwork, 
                      i_type liwork, i_type *info)
{
    LAPACK_NAME(cheevd)(_rc(jobv), _rc(uplo), &n, _rc(a), &lda, _rc(w), _rc(work), &lwork, _rc(rwork),
                        &lrwork, iwork, &liwork, info);
};

template<> BLAS_EXPORT void
lapack::heevd<z_type>(const char *jobv,const char *uplo, i_type n, z_type *a,i_type lda, d_type* w, 
                     z_type *work,i_type lwork, d_type* rwork, i_type lrwork, i_type* iwork, 
                      i_type liwork, i_type *info)
{
    LAPACK_NAME(zheevd)(_rc(jobv), _rc(uplo), &n, _rc(a), &lda, _rc(w), _rc(work), &lwork, 
                        _rc(rwork), &lrwork, iwork, &liwork, info);
};

//-----------------------------------------------------------------------
//                          hbev
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::sbev<d_type>(const char *jobv,const char *uplo, i_type n, i_type LD, d_type *a, i_type lda, 
     d_type* w, d_type* Z, i_type LDZ, d_type *work, i_type *info)
{
    LAPACK_NAME(dsbev)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, 
                       _rc(work), info);
};

template<> BLAS_EXPORT void
lapack::sbev<s_type>(const char *jobv,const char *uplo, i_type n, i_type LD, s_type *a, i_type lda, 
     s_type* w, s_type* Z, i_type LDZ, s_type *work, i_type *info)
{
    LAPACK_NAME(ssbev)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, 
                       _rc(work), info);
};

template<> BLAS_EXPORT void
lapack::hbev<c_type>(const char *jobv,const char *uplo, i_type n, i_type LD, c_type *a, i_type lda, 
     s_type* w, c_type* Z, i_type LDZ, c_type *work, s_type *rwork, i_type *info)
{
    LAPACK_NAME(chbev)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, 
                       _rc(work), _rc(rwork), info);
};

template<> BLAS_EXPORT void
lapack::hbev<z_type>(const char *jobv,const char *uplo, i_type n, i_type LD, z_type *a, i_type lda, 
     d_type* w, z_type* Z, i_type LDZ, z_type *work, d_type *rwork, i_type *info)
{
    LAPACK_NAME(zhbev)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, 
                       _rc(work), _rc(rwork), info);
};

//-----------------------------------------------------------------------
//                          DSBEVD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::sbevd<d_type>(const char *jobv,const char *uplo, i_type n, i_type LD, d_type *a,
                      i_type lda, d_type* w, d_type* Z, i_type LDZ, d_type *work, i_type lwork, 
                      i_type* iwork, i_type liwork, i_type *info)
{
    LAPACK_NAME(dsbevd)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, _rc(work), 
                        &lwork, iwork, &liwork, info);
};

template<> BLAS_EXPORT void
lapack::sbevd<s_type>(const char *jobv,const char *uplo, i_type n, i_type LD, s_type *a,
                      i_type lda, s_type* w, s_type* Z, i_type LDZ, s_type *work, i_type lwork, 
                      i_type* iwork, i_type liwork, i_type *info)
{
    LAPACK_NAME(ssbevd)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, _rc(work), 
                        &lwork, iwork, &liwork, info);
};

template<> BLAS_EXPORT void
lapack::hbevd<c_type>(const char *jobv,const char *uplo, i_type n, i_type LD, c_type *a,
                      i_type lda, s_type* w, c_type* Z, i_type LDZ, c_type *work, i_type lwork,
                      s_type *rwork, i_type lrwork, i_type* iwork, i_type liwork, i_type *info)
{
    LAPACK_NAME(chbevd)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, _rc(work), 
                        &lwork, _rc(rwork), &lrwork, iwork, &liwork, info);
};

template<> BLAS_EXPORT void
lapack::hbevd<z_type>(const char *jobv,const char *uplo, i_type n, i_type LD, z_type *a,
                      i_type lda, d_type* w, z_type* Z, i_type LDZ, z_type *work, i_type lwork, 
                      d_type *rwork, i_type lrwork, i_type* iwork, i_type liwork, i_type *info)
{
    LAPACK_NAME(zhbevd)(_rc(jobv), _rc(uplo), &n, &LD, _rc(a), &lda, _rc(w), _rc(Z), &LDZ, _rc(work), 
                        &lwork, _rc(rwork), &lrwork, iwork, &liwork, info);
};

//-----------------------------------------------------------------------
//                          DTREVC
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::trevc<d_type>(const char* side, const char* howmny, i_type* select, i_type n, d_type * t, 
    i_type ldt, d_type *vl, i_type ldvl, d_type* vr, i_type ldvr, i_type mm, i_type m, d_type* work, 
    i_type* info)
{
    LAPACK_NAME(dtrevc)(_rc(side), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl,
                        _rc(vr), &ldvr, &mm, &m, _rc(work), info);
};

template<> BLAS_EXPORT void
lapack::trevc<s_type>(const char* side, const char* howmny, i_type* select, i_type n, s_type * t, 
    i_type ldt, s_type *vl, i_type ldvl, s_type* vr, i_type ldvr, i_type mm, i_type m, s_type* work, 
    i_type* info)
{
    LAPACK_NAME(strevc)(_rc(side), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl, 
                        _rc(vr), &ldvr, &mm, &m, _rc(work), info);
};

template<> BLAS_EXPORT void
lapack::trevc<c_type>(const char* side, const char* howmny, i_type* select, i_type n, c_type * t, 
    i_type ldt, c_type *vl, i_type ldvl, c_type* vr, i_type ldvr, i_type mm, i_type m, c_type* work, 
    i_type* info)
{
    // work must contain 3N elements, we split this array into complex work of size 2n and real work
    // of size n
    using VR    = s_type;
    VR* rwork   = reinterpret_cast<VR*>(work + 2*n);

    LAPACK_NAME(ctrevc)(_rc(side), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl, _rc(vr), &ldvr,
                        &mm, &m, _rc(work), _rc(rwork), info);
};

template<> BLAS_EXPORT void
lapack::trevc<z_type>(const char* side, const char* howmny, i_type* select, i_type n, z_type * t, 
    i_type ldt, z_type *vl, i_type ldvl, z_type* vr, i_type ldvr, i_type mm, i_type m, z_type* work, 
    i_type* info)
{
    // work must contain 3N elements, we split this array into complex work of size 2n and real work
    // of size n
    using VR    = d_type;
    VR* rwork   = reinterpret_cast<VR*>(work + 2*n);

    LAPACK_NAME(ztrevc)(_rc(side), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl, _rc(vr),
                        &ldvr, &mm, &m, _rc(work), _rc(rwork), info);
};

//-----------------------------------------------------------------------
//                          DTRSNA
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::trsna<d_type>(const char *job, const char *howmny, const i_type *select, i_type n, const d_type* t, 
        i_type ldt, const d_type* vl, i_type ldvl, const d_type* vr, i_type ldvr, d_type* s, d_type * sep, 
        i_type mm, i_type m, d_type* work, i_type ldwork, i_type* iwork, d_type* rwork, i_type* info)
{
    (void)rwork;
    LAPACK_NAME(dtrsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl, _rc(vr), 
                        &ldvr, _rc(s), _rc(sep), &mm, &m, _rc(work), &ldwork, iwork, info);
};

template<> BLAS_EXPORT void
lapack::trsna<s_type>(const char *job, const char *howmny, const i_type *select, i_type n, const s_type* t, 
        i_type ldt, const s_type* vl, i_type ldvl, const s_type* vr, i_type ldvr, s_type* s, s_type * sep, 
        i_type mm, i_type m, s_type* work, i_type ldwork, i_type* iwork, s_type* rwork, i_type* info)
{
    (void)rwork;
    LAPACK_NAME(strsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl, _rc(vr), 
                        &ldvr, _rc(s), _rc(sep), &mm, &m, _rc(work), &ldwork, iwork, info);
};

template<> BLAS_EXPORT void
lapack::trsna<c_type>(const char *job, const char *howmny, const i_type *select, i_type n, const c_type* t, 
        i_type ldt, const c_type* vl, i_type ldvl, const c_type* vr, i_type ldvr, s_type* s, s_type * sep, 
        i_type mm, i_type m, c_type* work, i_type ldwork, i_type* iwork, s_type* rwork, i_type* info)
{
    (void)iwork;
    LAPACK_NAME(ctrsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl, _rc(vr), 
                        &ldvr, _rc(s), _rc(sep), &mm, &m, _rc(work), &ldwork, rwork, info);
};

template<> BLAS_EXPORT void
lapack::trsna<z_type>(const char *job, const char *howmny, const i_type *select, i_type n, const z_type* t, 
        i_type ldt, const z_type* vl, i_type ldvl, const z_type* vr, i_type ldvr, d_type* s, d_type * sep, 
        i_type mm, i_type m, z_type* work, i_type ldwork, i_type* iwork, d_type* rwork, i_type* info)
{
    (void)iwork;
    LAPACK_NAME(ztrsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(t), &ldt, _rc(vl), &ldvl, _rc(vr), 
                        &ldvr, _rc(s), _rc(sep), &mm, &m, _rc(work), &ldwork, rwork, info);
};

//-----------------------------------------------------------------------
//                          DLACN2
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::lacn2<d_type>(i_type n, d_type* v, d_type* x, i_type* isgn, d_type* est, i_type* kase, i_type* isave)
{
    LAPACK_NAME(dlacn2)(&n, _rc(v), _rc(x), _rc(isgn), _rc(est), kase, isave);
};

template<> BLAS_EXPORT void
lapack::lacn2<s_type>(i_type n, s_type* v, s_type* x, i_type* isgn, s_type* est, i_type* kase, i_type* isave)
{
    LAPACK_NAME(slacn2)(&n, _rc(v), _rc(x), _rc(isgn), _rc(est), kase, isave);
};

template<> BLAS_EXPORT void
lapack::lacn2<c_type>(i_type n, c_type* v, c_type* x, i_type* isgn, s_type* est, i_type* kase, i_type* isave)
{
    (void)isgn;
    LAPACK_NAME(clacn2)(&n, _rc(v), _rc(x), _rc(est), kase, isave);
};

template<> BLAS_EXPORT void
lapack::lacn2<z_type>(i_type n, z_type* v, z_type* x, i_type* isgn, d_type* est, i_type* kase, i_type* isave)
{
    (void)isgn;
    LAPACK_NAME(zlacn2)(&n, _rc(v), _rc(x), _rc(est), kase, isave);
};

//-----------------------------------------------------------------------
//                          HSEIN
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::hsein<d_type>(const char *side, const char *eigsrc, const char *initv, i_type* select, i_type n, 
      const d_type* h, i_type ldh, d_type* wr, d_type* wi, d_type* vl, i_type ldvl, d_type* vr, i_type ldvr,
      i_type mm, i_type m, d_type* work, i_type* ifaill, i_type* ifailr, i_type* info)
{
    LAPACK_NAME(dhsein)(_rc(side), _rc(eigsrc), _rc(initv), select, &n, _rc(h), &ldh, _rc(wr), _rc(wi), 
                        _rc(vl), &ldvl, _rc(vr), &ldvr, &mm, &m, _rc(work), ifaill, ifailr, info);
};

template<> BLAS_EXPORT void
lapack::hsein<s_type>(const char *side, const char *eigsrc, const char *initv, i_type* select, i_type n, 
      const s_type* h, i_type ldh, s_type* wr, s_type* wi, s_type* vl, i_type ldvl, s_type* vr, i_type ldvr,
      i_type mm, i_type m, s_type* work, i_type* ifaill, i_type* ifailr, i_type* info)
{
    LAPACK_NAME(shsein)(_rc(side), _rc(eigsrc), _rc(initv), select, &n, _rc(h), &ldh, _rc(wr), _rc(wi), 
                        _rc(vl), &ldvl, _rc(vr), &ldvr, &mm, &m, _rc(work), ifaill, ifailr, info);
};

template<> BLAS_EXPORT void
lapack::hsein<c_type>(const char *side, const char *eigsrc, const char *initv, i_type* select, i_type n, 
      const c_type* h, i_type ldh, c_type* w, c_type* vl, i_type ldvl, c_type* vr, i_type ldvr, i_type mm, 
      i_type m, c_type* work, s_type* rwork, i_type* ifaill, i_type* ifailr, i_type* info)
{
    LAPACK_NAME(chsein)(_rc(side), _rc(eigsrc), _rc(initv), select, &n, _rc(h), &ldh, _rc(w),
                        _rc(vl), &ldvl, _rc(vr), &ldvr, &mm, &m, _rc(work), _rc(rwork), ifaill, ifailr, info);
};

template<> BLAS_EXPORT void
lapack::hsein<z_type>(const char *side, const char *eigsrc, const char *initv, i_type* select, i_type n, 
      const z_type* h, i_type ldh, z_type* w, z_type* vl, i_type ldvl, z_type* vr, i_type ldvr, i_type mm, 
      i_type m, z_type* work, d_type* rwork, i_type* ifaill, i_type* ifailr, i_type* info)
{
    LAPACK_NAME(zhsein)(_rc(side), _rc(eigsrc), _rc(initv), select, &n, _rc(h), &ldh, _rc(w),
                        _rc(vl), &ldvl, _rc(vr), &ldvr, &mm, &m, _rc(work), _rc(rwork), ifaill, ifailr, info);
};

//-----------------------------------------------------------------------
//                          DSTEIN
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lapack::stein<d_type>(i_type n, const d_type* d, const d_type* e, i_type m, const d_type* w, 
                      const i_type* iblock, const i_type* isplit, d_type* z, i_type ldz, 
                      d_type * work, i_type* iwork, i_type* ifail, i_type* info)
{
    LAPACK_NAME(dstein)(&n, _rc(d), _rc(e), &m, _rc(w), _rc(iblock), _rc(isplit), _rc(z), &ldz, _rc(work), 
                        iwork, ifail, info);
};

template<> BLAS_EXPORT void
lapack::stein<s_type>(i_type n, const s_type* d, const s_type* e, i_type m, const s_type* w,
                      const i_type* iblock, const i_type* isplit, s_type* z, i_type ldz, 
                      s_type * work, i_type* iwork, i_type* ifail, i_type* info)
{
    LAPACK_NAME(sstein)(&n, _rc(d), _rc(e), &m, _rc(w), _rc(iblock), _rc(isplit), _rc(z), &ldz,
                        _rc(work), iwork, ifail, info);
};

template<> BLAS_EXPORT void
lapack::stein<c_type>(i_type n, const s_type* d, const s_type* e, i_type m, const s_type* w, 
                      const i_type* iblock, const i_type* isplit, c_type* z, i_type ldz,
                      s_type * work, i_type* iwork, i_type* ifail, i_type* info)
{
    LAPACK_NAME(cstein)(&n, _rc(d), _rc(e), &m, _rc(w), _rc(iblock), _rc(isplit), _rc(z), &ldz,
                        _rc(work), iwork, ifail, info);
};

template<> BLAS_EXPORT void
lapack::stein<z_type>(i_type n, const d_type* d, const d_type* e, i_type m, const d_type* w, 
                      const i_type* iblock, const i_type* isplit, z_type* z, i_type ldz, 
                      d_type * work, i_type* iwork, i_type* ifail, i_type* info)
{
    LAPACK_NAME(zstein)(&n, _rc(d), _rc(e), &m, _rc(w), _rc(iblock), _rc(isplit), _rc(z), &ldz,
                        _rc(work), iwork, ifail, info);
};

//-----------------------------------------------------------------------
//                          DSTEIN
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
tgevc<d_type>(const char *side, const char *howmny, i_type *select, i_type n, const d_type* s, 
              i_type lds, const d_type* p, i_type ldp, d_type* vl, i_type ldvl, d_type* vr, 
              i_type ldvr, i_type mm, i_type m, d_type* work, i_type* info)
{
    LAPACK_NAME(dtgevc)(_rc(side), _rc(howmny), select, &n, _rc(s), &lds, _rc(p), &ldp, _rc(vl), &ldvl,
            _rc(vr), &ldvr, &mm, &m, _rc(work), info);
};

template<> BLAS_EXPORT void
tgevc<s_type>(const char *side, const char *howmny, i_type *select, i_type n, const s_type* s,
              i_type lds, const s_type* p, i_type ldp, s_type* vl, i_type ldvl, s_type* vr, 
              i_type ldvr, i_type mm, i_type m, s_type* work, i_type* info)
{
    LAPACK_NAME(stgevc)(_rc(side), _rc(howmny), select, &n, _rc(s), &lds, _rc(p), &ldp, _rc(vl), &ldvl,
            _rc(vr), &ldvr, &mm, &m, _rc(work), info);
};

template<> BLAS_EXPORT void
tgevc<c_type>(const char *side, const char *howmny, i_type *select, i_type n, const c_type* s, 
              i_type lds, const c_type* p, i_type ldp, c_type* vl, i_type ldvl, c_type* vr, 
              i_type ldvr, i_type mm, i_type m, c_type* work, i_type* info)
{
    using VR = s_type;

    //complex version requires less storage, we separate WORK into WORK and RWORK
    VR* RWORK   = reinterpret_cast<VR*>(work + 2*n);

    LAPACK_NAME(ctgevc)(_rc(side), _rc(howmny), select, &n, _rc(s), &lds, _rc(p), &ldp, _rc(vl), &ldvl,
            _rc(vr), &ldvr, &mm, &m, _rc(work), _rc(RWORK), info);
};

template<> BLAS_EXPORT void
tgevc<z_type>(const char *side, const char *howmny, i_type *select, i_type n, const z_type* s, 
              i_type lds, const z_type* p, i_type ldp, z_type* vl, i_type ldvl, z_type* vr,
              i_type ldvr, i_type mm, i_type m, z_type* work, i_type* info)
{
    using VR = d_type;

    //complex version requires less storage, we separate WORK into WORK and RWORK
    VR* RWORK   = reinterpret_cast<VR*>(work + 2*n);

    LAPACK_NAME(ztgevc)(_rc(side), _rc(howmny), select, &n, _rc(s), &lds, _rc(p), &ldp, 
                        _rc(vl), &ldvl, _rc(vr), &ldvr, &mm, &m, _rc(work), _rc(RWORK), info);
};

//-----------------------------------------------------------------------
//                          DTGSNA
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
tgsna<d_type>(const char *job, const char *howmny, const i_type *select, i_type n, 
              const d_type* a, i_type lda, const d_type* b, i_type ldb, d_type* vl, 
              i_type ldvl, d_type* vr, i_type ldvr, d_type* s, d_type* dif, i_type mm, 
              i_type m, d_type* work, i_type lwork, i_type* iwork, i_type *info)
{
    LAPACK_NAME(dtgsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(a), &lda, _rc(b), &ldb,
                        _rc(vl), &ldvl, _rc(vr), &ldvr, _rc(s), dif, &mm, &m, _rc(work), &lwork,
                        iwork, info);
};

template<> BLAS_EXPORT void
tgsna<s_type>(const char *job, const char *howmny, const i_type *select, i_type n, const s_type* a,
              i_type lda, const s_type* b, i_type ldb, s_type* vl, i_type ldvl, s_type* vr, i_type ldvr, 
      s_type* s, s_type* dif, i_type mm, i_type m, s_type* work, i_type lwork, i_type* iwork, i_type *info)
{
    LAPACK_NAME(stgsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(a), &lda, _rc(b), &ldb, _rc(vl), &ldvl,
            _rc(vr), &ldvr, _rc(s), dif, &mm, &m, _rc(work), &lwork, iwork, info);
};

template<> BLAS_EXPORT void
tgsna<c_type>(const char *job, const char *howmny, const i_type *select, i_type n, const c_type* a, 
              i_type lda, const c_type* b, i_type ldb, c_type* vl, i_type ldvl, c_type* vr, 
              i_type ldvr, s_type* s, s_type* dif, i_type mm, i_type m, c_type* work, i_type lwork, 
              i_type* iwork, i_type *info)
{
    LAPACK_NAME(ctgsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(a), &lda, _rc(b), &ldb, _rc(vl), &ldvl,
            _rc(vr), &ldvr, _rc(s), dif, &mm, &m, _rc(work), &lwork, iwork, info);
};

template<> BLAS_EXPORT void
tgsna<z_type>(const char *job, const char *howmny, const i_type *select, i_type n, const z_type* a,
              i_type lda, const z_type* b, i_type ldb, z_type* vl, i_type ldvl, z_type* vr, i_type ldvr, 
      d_type* s, d_type* dif, i_type mm, i_type m, z_type* work, i_type lwork, i_type* iwork, i_type *info)
{
    LAPACK_NAME(ztgsna)(_rc(job), _rc(howmny), _rc(select), &n, _rc(a), &lda, _rc(b), &ldb, _rc(vl), &ldvl,
            _rc(vr), &ldvr, _rc(s), dif, &mm, &m, _rc(work), &lwork, iwork, info);
};

//-----------------------------------------------------------------------
//                          DLASSQ
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lassq<d_type>(i_type n, const d_type* x, i_type incx, d_type& scale, d_type& sumsq)
{
    LAPACK_NAME(dlassq)(&n, _rc(x), &incx, _rc(&scale), _rc(&sumsq));
};

template<> BLAS_EXPORT void
lassq<s_type>(i_type n, const s_type* x, i_type incx, s_type& scale, s_type& sumsq)
{
    LAPACK_NAME(slassq)(&n, _rc(x), &incx, _rc(&scale), _rc(&sumsq));
};

template<> BLAS_EXPORT void
lassq<c_type>(i_type n, const c_type* x, i_type incx, s_type& scale, s_type& sumsq)
{
    LAPACK_NAME(classq)(&n, _rc(x), &incx, _rc(&scale), _rc(&sumsq));
};

template<> BLAS_EXPORT void
lassq<z_type>(i_type n, const z_type* x, i_type incx, d_type& scale, d_type& sumsq)
{
    LAPACK_NAME(zlassq)(&n, _rc(x), &incx, _rc(&scale), _rc(&sumsq));
};

//-----------------------------------------------------------------------
//                          DLASCL
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lascl<d_type>(const char *type, i_type kl, i_type ku, d_type cfrom, d_type cto, i_type m, i_type n, 
    d_type* a, i_type lda, i_type& info)
{
    LAPACK_NAME(dlascl)(_rc(type), &kl, &ku, &cfrom, &cto, &m, &n, a, &lda, &info);
};

template<> BLAS_EXPORT void
lascl<s_type>(const char *type, i_type kl, i_type ku, s_type cfrom, s_type cto, i_type m, i_type n, 
    s_type* a, i_type lda, i_type& info)
{
    LAPACK_NAME(slascl)(_rc(type), &kl, &ku, &cfrom, &cto, &m, &n, a, &lda, &info);
};

template<> BLAS_EXPORT void
lascl<c_type>(const char *type, i_type kl, i_type ku, s_type cfrom, s_type cto, i_type m, i_type n, 
    c_type* a, i_type lda, i_type& info)
{
    LAPACK_NAME(clascl)(_rc(type), &kl, &ku, &cfrom, &cto, &m, &n, _rc(a), &lda, &info);
};

template<> BLAS_EXPORT void
lascl<z_type>(const char *type, i_type kl, i_type ku, d_type cfrom, d_type cto, i_type m, i_type n, 
    z_type* a, i_type lda, i_type& info)
{
    LAPACK_NAME(zlascl)(_rc(type), &kl, &ku, &cfrom, &cto, &m, &n, _rc(a), &lda, &info);
};

//-----------------------------------------------------------------------
//                          DGGBAL
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
ggbal<d_type>(const char *job, i_type n, d_type* a, i_type lda, d_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, d_type* lscale, d_type* rscale, d_type* work, i_type& info)
{
    LAPACK_NAME(dggbal)(_rc(job), &n, _rc(a), &lda, _rc(b), &ldb, &ilo, &ihi, _rc(lscale), _rc(rscale), 
            _rc(work), &info);
};

template<> BLAS_EXPORT void
ggbal<s_type>(const char *job, i_type n, s_type* a, i_type lda, s_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, s_type* lscale, s_type* rscale, s_type* work, i_type& info)
{
    LAPACK_NAME(sggbal)(_rc(job), &n, _rc(a), &lda, _rc(b), &ldb, &ilo, &ihi, _rc(lscale), _rc(rscale), 
            _rc(work), &info);
};

template<> BLAS_EXPORT void
ggbal<c_type>(const char *job, i_type n, c_type* a, i_type lda, c_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, s_type* lscale, s_type* rscale, s_type* work, i_type& info)
{
    LAPACK_NAME(cggbal)(_rc(job), &n, _rc(a), &lda, _rc(b), &ldb, &ilo, &ihi, _rc(lscale), _rc(rscale), 
            _rc(work), &info);
};

template<> BLAS_EXPORT void
ggbal<z_type>(const char *job, i_type n, z_type* a, i_type lda, z_type* b, i_type ldb, i_type& ilo,
        i_type& ihi, d_type* lscale, d_type* rscale, d_type* work, i_type& info)
{
    LAPACK_NAME(zggbal)(_rc(job), &n, _rc(a), &lda, _rc(b), &ldb, &ilo, &ihi, _rc(lscale), _rc(rscale), 
            _rc(work), &info);
};

//-----------------------------------------------------------------------
//                          DGGBAK
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
ggbak<d_type>(const char* job, const char* side, i_type n, i_type ilo, i_type ihi, const d_type* lscale, 
    const d_type* rscale, i_type m, d_type* v, i_type ldv, i_type& info)
{
    LAPACK_NAME(dggbak)(_rc(job), _rc(side), &n, &ilo, &ihi, _rc(lscale), _rc(rscale), &m, _rc(v), &ldv,
            &info);
};

template<> BLAS_EXPORT void
ggbak<s_type>(const char* job, const char* side, i_type n, i_type ilo, i_type ihi, const s_type* lscale, 
    const s_type* rscale, i_type m, s_type* v, i_type ldv, i_type& info)
{
    LAPACK_NAME(sggbak)(_rc(job), _rc(side), &n, &ilo, &ihi, _rc(lscale), _rc(rscale), &m, _rc(v), &ldv,
            &info);
};

template<> BLAS_EXPORT void
ggbak<c_type>(const char* job, const char* side, i_type n, i_type ilo, i_type ihi, const s_type* lscale, 
    const s_type* rscale, i_type m, c_type* v, i_type ldv, i_type& info)
{
    LAPACK_NAME(cggbak)(_rc(job), _rc(side), &n, &ilo, &ihi, _rc(lscale), _rc(rscale), &m, _rc(v), &ldv,
            &info);
};

template<> BLAS_EXPORT void
ggbak<z_type>(const char* job, const char* side, i_type n, i_type ilo, i_type ihi, const d_type* lscale, 
    const d_type* rscale, i_type m, z_type* v, i_type ldv, i_type& info)
{
    LAPACK_NAME(zggbak)(_rc(job), _rc(side), &n, &ilo, &ihi, _rc(lscale), _rc(rscale), &m, _rc(v), &ldv,
            &info);
};

//-----------------------------------------------------------------------
//                          DORMQR
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
ormqr<d_type>(const char* side, const char* trans, i_type m, i_type n, i_type k, d_type* a, i_type lda, 
        const d_type* tau, d_type* c, i_type ldc, d_type* work, i_type lwork, i_type& info)
{
    LAPACK_NAME(dormqr)(_rc(side), _rc(trans), &m, &n, &k, _rc(a), &lda, _rc(tau), _rc(c), &ldc, _rc(work),
            &lwork, &info);
};

template<> BLAS_EXPORT void
ormqr<s_type>(const char* side, const char* trans, i_type m, i_type n, i_type k, s_type* a, i_type lda, 
        const s_type* tau, s_type* c, i_type ldc, s_type* work, i_type lwork, i_type& info)
{
    LAPACK_NAME(sormqr)(_rc(side), _rc(trans), &m, &n, &k, _rc(a), &lda, _rc(tau), _rc(c), &ldc, _rc(work),
            &lwork, &info);
};

template<> BLAS_EXPORT void
ormqr<c_type>(const char* side, const char* trans, i_type m, i_type n, i_type k, c_type* a, i_type lda, 
        const c_type* tau, c_type* c, i_type ldc, c_type* work, i_type lwork, i_type& info)
{
    LAPACK_NAME(cunmqr)(_rc(side), _rc(trans), &m, &n, &k, _rc(a), &lda, _rc(tau), _rc(c), &ldc, _rc(work),
            &lwork, &info);
};

template<> BLAS_EXPORT void
ormqr<z_type>(const char* side, const char* trans, i_type m, i_type n, i_type k, z_type* a, i_type lda, 
        const z_type* tau, z_type* c, i_type ldc, z_type* work, i_type lwork, i_type& info)
{
    LAPACK_NAME(zunmqr)(_rc(side), _rc(trans), &m, &n, &k, _rc(a), &lda, _rc(tau), _rc(c), &ldc, _rc(work),
            &lwork, &info);
};

//-----------------------------------------------------------------------
//                          DHGEQZ
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
hgeqz<d_type>(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
    d_type* h, i_type ldh, d_type* t, i_type ldt, d_type* alphar, d_type* alphai, d_type* beta, d_type* q, 
    i_type ldq, d_type* z, i_type ldz, d_type* work, i_type lwork, i_type& info)
{
    LAPACK_NAME(dhgeqz)(_rc(job), _rc(compq), _rc(compz), &n, &ilo, &ihi, _rc(h), &ldh, _rc(t), &ldt, 
    _rc(alphar), _rc(alphai), _rc(beta), _rc(q), &ldq, _rc(z), &ldz, _rc(work), &lwork, &info);
};

template<> BLAS_EXPORT void
hgeqz<s_type>(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
    s_type* h, i_type ldh, s_type* t, i_type ldt, s_type* alphar, s_type* alphai, s_type* beta, s_type* q, 
    i_type ldq, s_type* z, i_type ldz, s_type* work, i_type lwork, i_type& info)
{
    LAPACK_NAME(shgeqz)(_rc(job), _rc(compq), _rc(compz), &n, &ilo, &ihi, _rc(h), &ldh, _rc(t), &ldt, 
    _rc(alphar), _rc(alphai), _rc(beta), _rc(q), &ldq, _rc(z), &ldz, _rc(work), &lwork, &info);
};

template<> BLAS_EXPORT void
hgeqz<c_type>(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
    c_type* h, i_type ldh, c_type* t, i_type ldt, c_type* alpha, c_type* beta, c_type* q, 
    i_type ldq, c_type* z, i_type ldz, c_type* work, i_type lwork, s_type* rwork, i_type& info)
{
    LAPACK_NAME(chgeqz)(_rc(job), _rc(compq), _rc(compz), &n, &ilo, &ihi, _rc(h), &ldh, _rc(t), &ldt, 
    _rc(alpha), _rc(beta), _rc(q), &ldq, _rc(z), &ldz, _rc(work), &lwork, _rc(rwork), &info);
};

template<> BLAS_EXPORT void
hgeqz<z_type>(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
    z_type* h, i_type ldh, z_type* t, i_type ldt, z_type* alpha, z_type* beta, z_type* q, 
    i_type ldq, z_type* z, i_type ldz, z_type* work, i_type lwork, d_type* rwork, i_type& info)
{
    LAPACK_NAME(zhgeqz)(_rc(job), _rc(compq), _rc(compz), &n, &ilo, &ihi, _rc(h), &ldh, _rc(t), &ldt, 
    _rc(alpha), _rc(beta), _rc(q), &ldq, _rc(z), &ldz, _rc(work), &lwork, _rc(rwork), &info);
};

//-----------------------------------------------------------------------
//                          DLABAD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
labad<d_type>(d_type& small, d_type& large)
{
    LAPACK_NAME(dlabad)(&small, &large);
};

template<> BLAS_EXPORT void
labad<s_type>(s_type& small, s_type& large)
{
    LAPACK_NAME(slabad)(&small, &large);
};

//-----------------------------------------------------------------------
//                          DLANHS
//-----------------------------------------------------------------------

template<> BLAS_EXPORT d_type
lanhs<d_type>(const char* NORM, i_type N, const d_type* A, i_type LDA, d_type* WORK )
{
    return LAPACK_NAME(dlanhs)(_rc(NORM), &N, _rc(A), &LDA, WORK );
};

template<> BLAS_EXPORT s_type
lanhs<s_type>(const char* NORM, i_type N, const s_type* A, i_type LDA, s_type* WORK )
{
    return LAPACK_NAME(slanhs)(_rc(NORM), &N, _rc(A), &LDA, WORK );
};

//-----------------------------------------------------------------------
//                          DLAG2
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
lag2<d_type>(const d_type* A, i_type LDA, const d_type* B, i_type LDB, d_type SAFMIN, d_type& SCALE1, 
             d_type& SCALE2, d_type& WR1, d_type& WR2, d_type& WI )
{
    LAPACK_NAME(dlag2)(_rc(A), &LDA, _rc(B), &LDB, &SAFMIN, &SCALE1, &SCALE2, &WR1, &WR2, &WI );
};

template<> BLAS_EXPORT void
lag2<s_type>(const s_type* A, i_type LDA, const s_type* B, i_type LDB, s_type SAFMIN, s_type& SCALE1, 
             s_type& SCALE2, s_type& WR1, s_type& WR2, s_type& WI )
{
    LAPACK_NAME(slag2)(_rc(A), &LDA, _rc(B), &LDB, &SAFMIN, &SCALE1, &SCALE2, &WR1, &WR2, &WI );
};

//-----------------------------------------------------------------------
//                          DSBGVD
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
hbgvd<d_type>(const char* JOBZ, const char* UPLO, i_type N, i_type KA, i_type KB,
              d_type* AB, i_type LDAB, d_type* BB, i_type LDBB, d_type* W, d_type* Z, 
              i_type LDZ, d_type* WORK, i_type LWORK, d_type* RWORK, i_type LRWORK, 
              i_type* IWORK, i_type LIWORK, i_type& INFO )
{
    LAPACK_NAME(dsbgvd)(_rc(JOBZ), _rc(UPLO), &N, &KA, &KB, _rc(AB), &LDAB, _rc(BB), &LDBB, _rc(W), 
                _rc(Z), &LDZ, _rc(WORK), &LWORK, IWORK, &LIWORK, &INFO);

    if (LRWORK == -1 || LWORK == -1 || LIWORK == -1)
    {
        *RWORK = 0;
    }
};

template<> BLAS_EXPORT void
hbgvd<s_type>(const char* JOBZ, const char* UPLO, i_type N, i_type KA, i_type KB, 
              s_type* AB, i_type LDAB, s_type* BB, i_type LDBB, s_type* W, s_type* Z,
              i_type LDZ, s_type* WORK, i_type LWORK, s_type* RWORK, i_type LRWORK, 
              i_type* IWORK, i_type LIWORK, i_type& INFO )
{
    LAPACK_NAME(ssbgvd)(_rc(JOBZ), _rc(UPLO), &N, &KA, &KB, _rc(AB), &LDAB, _rc(BB), &LDBB, _rc(W), 
                _rc(Z), &LDZ, _rc(WORK), &LWORK, IWORK, &LIWORK, &INFO);

    if (LRWORK == -1 || LWORK == -1 || LIWORK == -1)
    {
        *RWORK = 0;
    }
};

template<> BLAS_EXPORT void
hbgvd<c_type>(const char* JOBZ, const char* UPLO, i_type N, i_type KA, i_type KB, 
              c_type* AB, i_type LDAB, c_type* BB, i_type LDBB, s_type* W, c_type* Z,
              i_type LDZ, c_type* WORK, i_type LWORK, 
              s_type* RWORK, i_type LRWORK, i_type* IWORK, i_type LIWORK, i_type& INFO )
{
    LAPACK_NAME(chbgvd)(_rc(JOBZ), _rc(UPLO), &N, &KA, &KB, _rc(AB), &LDAB, _rc(BB), &LDBB, _rc(W), 
                _rc(Z), &LDZ, _rc(WORK), &LWORK, _rc(RWORK), &LRWORK, IWORK, &LIWORK, &INFO);
};

template<> BLAS_EXPORT void
hbgvd<z_type>(const char* JOBZ, const char* UPLO, i_type N, i_type KA, i_type KB, 
              z_type* AB, i_type LDAB, z_type* BB, i_type LDBB, d_type* W, z_type* Z, 
              i_type LDZ, z_type* WORK, i_type LWORK, d_type* RWORK, i_type LRWORK, 
              i_type* IWORK, i_type LIWORK, i_type& INFO )
{
    LAPACK_NAME(zhbgvd)(_rc(JOBZ), _rc(UPLO), &N, &KA, &KB, _rc(AB), &LDAB, _rc(BB), &LDBB, _rc(W), 
                _rc(Z), &LDZ, _rc(WORK), &LWORK, _rc(RWORK), &LRWORK, IWORK, &LIWORK, &INFO);
};

//-----------------------------------------------------------------------
//                          DNRM2
//-----------------------------------------------------------------------

template<> BLAS_EXPORT d_type
nrm2(i_type N, const d_type* X, i_type INCX)
{
    return LAPACK_NAME(dnrm2)(&N, _rc(X), &INCX);
};

template<> BLAS_EXPORT s_type
nrm2(i_type N, const s_type* X, i_type INCX)
{
    return (s_type)LAPACK_NAME(snrm2)(&N, _rc(X), &INCX);
};

template<> BLAS_EXPORT d_type
nrm2(i_type N, const z_type* X, i_type INCX)
{
    return LAPACK_NAME(dznrm2)(&N, _rc(X), &INCX);
};

template<> BLAS_EXPORT s_type
nrm2(i_type N, const c_type* X, i_type INCX)
{
    return (s_type)LAPACK_NAME(scnrm2)(&N, _rc(X), &INCX);
};

//-----------------------------------------------------------------------
//                          pttrf
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
pttrf<d_type>(i_type N, d_type* D, d_type* E, i_type& info)
{
    LAPACK_NAME(dpttrf)(&N, _rc(D), _rc(E), &info);
};

template<> BLAS_EXPORT void
pttrf<s_type>(i_type N, s_type* D, s_type* E, i_type& info)
{
    LAPACK_NAME(spttrf)(&N, _rc(D), _rc(E), &info);
};

template<> BLAS_EXPORT void
pttrf<c_type>(i_type N, s_type* D, c_type* E, i_type& info)
{
    LAPACK_NAME(cpttrf)(&N, _rc(D), _rc(E), &info);
};

template<> BLAS_EXPORT void
pttrf<z_type>(i_type N, d_type* D, z_type* E, i_type& info)
{
    LAPACK_NAME(zpttrf)(&N, _rc(D), _rc(E), &info);
};

//-----------------------------------------------------------------------
//                          PTTRS
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
pttrs<d_type>(const char*, i_type N, i_type NRHS, const d_type* D, const d_type* E, 
              d_type* B, i_type LDB, i_type& INFO )
{
    LAPACK_NAME(dpttrs)(&N, &NRHS, _rc(D), _rc(E), _rc(B), &LDB, &INFO );
};

template<> BLAS_EXPORT void
pttrs<s_type>(const char*, i_type N, i_type NRHS, const s_type* D, const s_type* E, 
              s_type* B, i_type LDB, i_type& INFO )
{
    LAPACK_NAME(spttrs)(&N, &NRHS, _rc(D), _rc(E), _rc(B), &LDB, &INFO );
};

template<> BLAS_EXPORT void
pttrs<z_type>(const char* UPLO, i_type N, i_type NRHS, const d_type* D, const z_type* E, z_type* B, 
              i_type LDB, i_type& INFO )
{
    LAPACK_NAME(zpttrs)(_rc(UPLO), &N, &NRHS, _rc(D), _rc(E), _rc(B), &LDB, &INFO );
};

template<> BLAS_EXPORT void
pttrs<c_type>(const char* UPLO, i_type N, i_type NRHS, const s_type* D, const c_type* E, c_type* B, 
              i_type LDB, i_type& INFO )
{
    LAPACK_NAME(cpttrs)(_rc(UPLO), &N, &NRHS, _rc(D), _rc(E), _rc(B), &LDB, &INFO );
};

//-----------------------------------------------------------------------
//                          DTRTRI
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
trtri<d_type>(const char* UPLO, const char* DIAG, i_type N, d_type* A, i_type LDA, i_type& INFO)
{
    LAPACK_NAME(dtrtri)(_rc(UPLO), _rc(DIAG), &N, _rc(A), &LDA, &INFO);
};

template<> BLAS_EXPORT void
trtri<s_type>(const char* UPLO, const char* DIAG, i_type N, s_type* A, i_type LDA, i_type& INFO)
{
    LAPACK_NAME(strtri)(_rc(UPLO), _rc(DIAG), &N, _rc(A), &LDA, &INFO);
};

template<> BLAS_EXPORT void
trtri<z_type>(const char* UPLO, const char* DIAG, i_type N, z_type* A, i_type LDA, i_type& INFO)
{
    LAPACK_NAME(ztrtri)(_rc(UPLO), _rc(DIAG), &N, _rc(A), &LDA, &INFO);
};

template<> BLAS_EXPORT void
trtri<c_type>(const char* UPLO, const char* DIAG, i_type N, c_type* A, i_type LDA, i_type& INFO)
{
    LAPACK_NAME(ctrtri)(_rc(UPLO), _rc(DIAG), &N, _rc(A), &LDA, &INFO);
};

//-----------------------------------------------------------------------
//                          DGETRI
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
getri<d_type>(i_type N, d_type* A, i_type LDA, const i_type* IPIV, d_type* WORK, i_type LWORK, 
              i_type& INFO )
{
    LAPACK_NAME(dgetri)(&N, _rc(A), &LDA, _rc(IPIV), _rc(WORK), &LWORK, &INFO);
};

template<> BLAS_EXPORT void
getri<s_type>(i_type N, s_type* A, i_type LDA, const i_type* IPIV, s_type* WORK, i_type LWORK, 
              i_type& INFO )
{
    LAPACK_NAME(sgetri)(&N, _rc(A), &LDA, _rc(IPIV), _rc(WORK), &LWORK, &INFO);
};

template<> BLAS_EXPORT void
getri<z_type>(i_type N, z_type* A, i_type LDA, const i_type* IPIV, z_type* WORK, i_type LWORK, 
              i_type& INFO )
{
    LAPACK_NAME(zgetri)(&N, _rc(A), &LDA, _rc(IPIV), _rc(WORK), &LWORK, &INFO);
};

template<> BLAS_EXPORT void
getri<c_type>(i_type N, c_type* A, i_type LDA, const i_type* IPIV, c_type* WORK, i_type LWORK, 
              i_type& INFO )
{
    LAPACK_NAME(cgetri)(&N, _rc(A), &LDA, _rc(IPIV), _rc(WORK), &LWORK, &INFO);
};

//-----------------------------------------------------------------------
//                          gttrf
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
gttrf<d_type>(i_type N, d_type* DL, d_type* D, d_type* DU, d_type* DU2, i_type* IPIV, i_type& INFO)
{
    LAPACK_NAME(dgttrf)(&N, _rc(DL), _rc(D), _rc(DU), _rc(DU2), IPIV, &INFO);
};

template<> BLAS_EXPORT void
gttrf<s_type>(i_type N, s_type* DL, s_type* D, s_type* DU, s_type* DU2, i_type* IPIV, i_type& INFO)
{
    LAPACK_NAME(sgttrf)(&N, _rc(DL), _rc(D), _rc(DU), _rc(DU2), IPIV, &INFO);
};

template<> BLAS_EXPORT void
gttrf<z_type>(i_type N, z_type* DL, z_type* D, z_type* DU, z_type* DU2, i_type* IPIV, i_type& INFO)
{
    LAPACK_NAME(zgttrf)(&N, _rc(DL), _rc(D), _rc(DU), _rc(DU2), IPIV, &INFO);
};

template<> BLAS_EXPORT void
gttrf<c_type>(i_type N, c_type* DL, c_type* D, c_type* DU, c_type* DU2, i_type* IPIV, i_type& INFO)
{
    LAPACK_NAME(cgttrf)(&N, _rc(DL), _rc(D), _rc(DU), _rc(DU2), IPIV, &INFO);
};

//-----------------------------------------------------------------------
//                          gttrs
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
gttrs<d_type>(const char* TRANS, i_type N, i_type NRHS, const d_type* DL, const d_type* D, 
              const d_type* DU, const d_type* DU2, const i_type* IPIV, d_type* B, i_type LDB,
              i_type& INFO)
{
    LAPACK_NAME(dgttrs)(_rc(TRANS), &N, &NRHS, _rc(DL), _rc(D), _rc(DU), _rc(DU2), _rc(IPIV), _rc(B), 
                        &LDB, &INFO);
};

template<> BLAS_EXPORT void
gttrs<s_type>(const char* TRANS, i_type N, i_type NRHS, const s_type* DL, const s_type* D, 
              const s_type* DU, const s_type* DU2, const i_type* IPIV, s_type* B, i_type LDB, 
              i_type& INFO)
{
    LAPACK_NAME(sgttrs)(_rc(TRANS), &N, &NRHS, _rc(DL), _rc(D), _rc(DU), _rc(DU2), _rc(IPIV), _rc(B), 
                        &LDB, &INFO);
};

template<> BLAS_EXPORT void
gttrs<z_type>(const char* TRANS, i_type N, i_type NRHS, const z_type* DL, const z_type* D, 
              const z_type* DU, const z_type* DU2, const i_type* IPIV, z_type* B, i_type LDB, 
              i_type& INFO)
{
    LAPACK_NAME(zgttrs)(_rc(TRANS), &N, &NRHS, _rc(DL), _rc(D), _rc(DU), _rc(DU2), _rc(IPIV), _rc(B), 
                        &LDB, &INFO);
};

template<> BLAS_EXPORT void
gttrs<c_type>(const char* TRANS, i_type N, i_type NRHS, const c_type* DL, const c_type* D, 
              const c_type* DU, const c_type* DU2, const i_type* IPIV, c_type* B, i_type LDB,
              i_type& INFO)
{
    LAPACK_NAME(cgttrs)(_rc(TRANS), &N, &NRHS, _rc(DL), _rc(D), _rc(DU), _rc(DU2), _rc(IPIV), _rc(B), 
                        &LDB, &INFO);
};

//-----------------------------------------------------------------------
//                          gttrs
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
stebz<d_type>(const char* RANGE, const char* ORDER, i_type N, d_type VL, d_type VU, i_type IL,
              i_type IU, d_type ABSTOL, const d_type* D, const d_type* E, i_type& M, 
              i_type& NSPLIT, d_type* W, i_type* IBLOCK, i_type* ISPLIT, d_type* WORK, 
              i_type* IWORK, i_type& INFO )
{
    LAPACK_NAME(dstebz)(_rc(RANGE), _rc(ORDER), &N, &VL, &VU, &IL, &IU, &ABSTOL, _rc(D), _rc(E), &M, 
                        &NSPLIT, _rc(W), _rc(IBLOCK), _rc(ISPLIT), _rc(WORK), IWORK, &INFO );
};

template<> BLAS_EXPORT void
stebz<s_type>(const char* RANGE, const char* ORDER, i_type N, s_type VL, s_type VU, i_type IL, 
              i_type IU, s_type ABSTOL, const s_type* D, const s_type* E, i_type& M, 
              i_type& NSPLIT, s_type* W, i_type* IBLOCK, i_type* ISPLIT, s_type* WORK, 
              i_type* IWORK, i_type& INFO )
{
    LAPACK_NAME(sstebz)(_rc(RANGE), _rc(ORDER), &N, &VL, &VU, &IL, &IU, &ABSTOL, _rc(D), _rc(E), &M, 
                        &NSPLIT, _rc(W), _rc(IBLOCK), _rc(ISPLIT), _rc(WORK), IWORK, &INFO );
};

//-----------------------------------------------------------------------
//                          ggsvd3
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
ggsvd3<d_type>(const char* JOBU, const char* JOBV, const char* JOBQ, i_type M, i_type N, 
               i_type P, i_type& K, i_type& L, d_type* ptr_A, i_type LDA, d_type* ptr_B, 
               i_type LDB, d_type* ptr_alpha, d_type* ptr_beta, d_type* ptr_U, i_type LDU,
               d_type* ptr_V, i_type LDV, d_type* ptr_Q, i_type LDQ, d_type* work, i_type lwork,
               d_type* rwork, i_type* iwork, i_type& info)
{
    (void)rwork;
    LAPACK_NAME(dggsvd3)(_rc(JOBU), _rc(JOBV), _rc(JOBQ), &M, &N, &P, &K, &L, 
        _rc(ptr_A), &LDA, _rc(ptr_B), &LDB, _rc(ptr_alpha), _rc(ptr_beta), _rc(ptr_U), 
        &LDU, _rc(ptr_V), &LDV, _rc(ptr_Q), &LDQ, _rc(work), &lwork, _rc(iwork), &info);
};

template<> BLAS_EXPORT void
ggsvd3<s_type>(const char* JOBU, const char* JOBV, const char* JOBQ, i_type M, i_type N, i_type P,
               i_type& K, i_type& L, s_type* ptr_A, i_type LDA, s_type* ptr_B, i_type LDB, 
               s_type* ptr_alpha, s_type* ptr_beta, s_type* ptr_U, i_type LDU, s_type* ptr_V,
               i_type LDV, s_type* ptr_Q, i_type LDQ, s_type* work, i_type lwork, s_type* rwork, 
               i_type* iwork, i_type& info)
{
    (void)rwork;
    LAPACK_NAME(sggsvd3)(_rc(JOBU), _rc(JOBV), _rc(JOBQ), &M, &N, &P, &K, &L, 
        _rc(ptr_A), &LDA, _rc(ptr_B), &LDB, _rc(ptr_alpha), _rc(ptr_beta), _rc(ptr_U), 
        &LDU, _rc(ptr_V), &LDV, _rc(ptr_Q), &LDQ, _rc(work), &lwork, _rc(iwork), &info);
};

template<> BLAS_EXPORT void
ggsvd3<c_type>(const char* JOBU, const char* JOBV, const char* JOBQ, i_type M, i_type N, i_type P, 
               i_type& K, i_type& L, c_type* ptr_A, i_type LDA, c_type* ptr_B, i_type LDB, 
               s_type* ptr_alpha, s_type* ptr_beta, c_type* ptr_U, i_type LDU, c_type* ptr_V,
               i_type LDV, c_type* ptr_Q, i_type LDQ, c_type* work, i_type lwork, s_type* rwork, 
               i_type* iwork, i_type& info)
{
    LAPACK_NAME(cggsvd3)(_rc(JOBU), _rc(JOBV), _rc(JOBQ), &M, &N, &P, &K, &L, 
        _rc(ptr_A), &LDA, _rc(ptr_B), &LDB, _rc(ptr_alpha), _rc(ptr_beta), _rc(ptr_U), 
        &LDU, _rc(ptr_V), &LDV, _rc(ptr_Q), &LDQ, _rc(work), &lwork, _rc(rwork), _rc(iwork), &info);
};

template<> BLAS_EXPORT void
ggsvd3<z_type>(const char* JOBU, const char* JOBV, const char* JOBQ, i_type M, i_type N, i_type P,
               i_type& K, i_type& L, z_type* ptr_A, i_type LDA, z_type* ptr_B, i_type LDB,
               d_type* ptr_alpha, d_type* ptr_beta, z_type* ptr_U, i_type LDU, z_type* ptr_V, 
               i_type LDV, z_type* ptr_Q, i_type LDQ, z_type* work, i_type lwork, d_type* rwork, 
               i_type* iwork, i_type& info)
{
    LAPACK_NAME(zggsvd3)(_rc(JOBU), _rc(JOBV), _rc(JOBQ), &M, &N, &P, &K, &L, 
        _rc(ptr_A), &LDA, _rc(ptr_B), &LDB, _rc(ptr_alpha), _rc(ptr_beta), _rc(ptr_U), 
        &LDU, _rc(ptr_V), &LDV, _rc(ptr_Q), &LDQ, _rc(work), &lwork, _rc(rwork), _rc(iwork), &info);
};

//-----------------------------------------------------------------------
//                          DGESVDX
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
gesvdx<d_type>(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N,
               d_type* A, i_type LDA, d_type VL, d_type VU, i_type IL, i_type IU, i_type& NS,
               d_type* S, d_type* U, i_type LDU, d_type* VT, i_type LDVT, d_type* WORK, 
               i_type LWORK, d_type* RWORK, i_type* IWORK, i_type& INFO)
{
    (void)RWORK;
    LAPACK_NAME(dgesvdx)(_rc(JOBU), _rc(JOBVT), _rc(RANGE), &M, &N, _rc(A), &LDA, &VL,
                         &VU, &IL, &IU, &NS, _rc(S), _rc(U), &LDU, _rc(VT), &LDVT, 
                         _rc(WORK), &LWORK, IWORK, &INFO);
};

template<> BLAS_EXPORT void
gesvdx<s_type>(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N,
               s_type* A, i_type LDA, s_type VL, s_type VU, i_type IL, i_type IU, i_type& NS,
               s_type* S, s_type* U, i_type LDU, s_type* VT, i_type LDVT, s_type* WORK, 
               i_type LWORK, s_type* RWORK, i_type* IWORK, i_type& INFO)
{
    (void)RWORK;
    LAPACK_NAME(sgesvdx)(_rc(JOBU), _rc(JOBVT), _rc(RANGE), &M, &N, _rc(A), &LDA, &VL,
                         &VU, &IL, &IU, &NS, _rc(S), _rc(U), &LDU, _rc(VT), &LDVT,
                         _rc(WORK), &LWORK, IWORK, &INFO);
};

template<> BLAS_EXPORT void
gesvdx<c_type>(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N, 
               c_type* A, i_type LDA, s_type VL, s_type VU, i_type IL, i_type IU, i_type& NS,
               s_type* S, c_type* U, i_type LDU, c_type* VT, i_type LDVT, c_type* WORK,
               i_type LWORK, s_type* RWORK, i_type* IWORK, i_type& INFO)
{
    LAPACK_NAME(cgesvdx)(_rc(JOBU), _rc(JOBVT), _rc(RANGE), &M, &N, _rc(A), &LDA, &VL,
                         &VU, &IL, &IU, &NS, _rc(S), _rc(U), &LDU, _rc(VT), &LDVT,
                         _rc(WORK), &LWORK, _rc(RWORK), IWORK, &INFO);
};

template<> BLAS_EXPORT void
gesvdx<z_type>(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, i_type N,
               z_type* A, i_type LDA, d_type VL, d_type VU, i_type IL, i_type IU, i_type& NS,
               d_type* S, z_type* U, i_type LDU, z_type* VT, i_type LDVT, z_type* WORK,
               i_type LWORK, d_type* RWORK, i_type* IWORK, i_type& INFO)
{
    LAPACK_NAME(zgesvdx)(_rc(JOBU), _rc(JOBVT), _rc(RANGE), &M, &N, _rc(A), &LDA, &VL,
                         &VU, &IL, &IU, &NS, _rc(S), _rc(U), &LDU, _rc(VT), &LDVT, 
                         _rc(WORK), &LWORK, _rc(RWORK), IWORK, &INFO);
};


//-----------------------------------------------------------------------
//                          DBDSVDX
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
bdsvdx<d_type>(const char* UPLO, const char* JOBZ, const char* RANGE, i_type N, 
               const d_type* D, const d_type* E, d_type VL, d_type VU, i_type IL, i_type IU,
               i_type& NS, d_type* S, d_type* Z, i_type LDZ, d_type* WORK, i_type* IWORK, 
               i_type& INFO )
{
    LAPACK_NAME(dbdsvdx)(_rc(UPLO), _rc(JOBZ), _rc(RANGE), &N, _rc(D), _rc(E), _rc(&VL),
                         _rc(&VU), &IL, &IU, &NS, _rc(S), _rc(Z), &LDZ, _rc(WORK), 
                         IWORK, &INFO );
};

template<> BLAS_EXPORT void
bdsvdx<s_type>(const char* UPLO, const char* JOBZ, const char* RANGE, i_type N, const s_type* D, 
               const s_type* E, s_type VL, s_type VU, i_type IL, i_type IU, i_type& NS, s_type* S, 
               s_type* Z, i_type LDZ, s_type* WORK, i_type* IWORK, i_type& INFO )
{
    LAPACK_NAME(sbdsvdx)(_rc(UPLO), _rc(JOBZ), _rc(RANGE), &N, _rc(D), _rc(E), _rc(&VL), 
                         _rc(&VU), &IL, &IU, &NS, _rc(S), _rc(Z), &LDZ, _rc(WORK), IWORK, 
                         &INFO );
};

//-----------------------------------------------------------------------
//                          DORMBR
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
ormbr<d_type>(const char* VECT, const char* SIDE, const char* TRANS, i_type M, i_type N, 
              i_type K, const d_type* A, i_type LDA, const d_type* TAU, d_type* C, i_type LDC, 
              d_type* WORK, i_type LWORK, i_type& INFO )
{
    LAPACK_NAME(dormbr)(_rc(VECT), _rc(SIDE), _rc(TRANS), &M, &N, &K, _rc(A), &LDA, _rc(TAU), 
                        _rc(C), &LDC, _rc(WORK), &LWORK, &INFO );
};

template<> BLAS_EXPORT void
ormbr<s_type>(const char* VECT, const char* SIDE, const char* TRANS, i_type M, i_type N, 
              i_type K, const s_type* A, i_type LDA, const s_type* TAU, s_type* C, i_type LDC,
              s_type* WORK, i_type LWORK, i_type& INFO )
{
    LAPACK_NAME(sormbr)(_rc(VECT), _rc(SIDE), _rc(TRANS), &M, &N, &K, _rc(A), &LDA, _rc(TAU), 
                        _rc(C), &LDC, _rc(WORK), &LWORK, &INFO );
};

//-----------------------------------------------------------------------
//                          DGELQF
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
gelqf<d_type>(i_type M, i_type N, d_type* A, i_type LDA, d_type* TAU, d_type* WORK, 
              i_type LWORK, i_type& INFO )
{
    LAPACK_NAME(dgelqf)(&M, &N, _rc(A), &LDA, _rc(TAU), _rc(WORK), &LWORK, &INFO);
};

template<> BLAS_EXPORT void
gelqf<s_type>(i_type M, i_type N, s_type* A, i_type LDA, s_type* TAU, s_type* WORK,
              i_type LWORK, i_type& INFO )
{
    LAPACK_NAME(sgelqf)(&M, &N, _rc(A), &LDA, _rc(TAU), _rc(WORK), &LWORK, &INFO);
};

//-----------------------------------------------------------------------
//                          DORMLQ
//-----------------------------------------------------------------------

template<> BLAS_EXPORT void
ormlq<d_type>(const char* SIDE, const char* TRANS, i_type M, i_type N, i_type K, d_type* A, 
              i_type LDA, const d_type* TAU, d_type* C, i_type LDC, d_type* WORK, i_type LWORK,
              i_type& INFO )
{
    LAPACK_NAME(dormlq)(_rc(SIDE), _rc(TRANS), &M, &N, &K, _rc(A), &LDA, _rc(TAU), _rc(C), &LDC, 
                        _rc(WORK), &LWORK, &INFO );
};

template<> BLAS_EXPORT void
ormlq<s_type>(const char* SIDE, const char* TRANS, i_type M, i_type N, i_type K, s_type* A,
              i_type LDA, const s_type* TAU, s_type* C, i_type LDC, s_type* WORK, i_type LWORK, 
              i_type& INFO )
{
    LAPACK_NAME(sormlq)(_rc(SIDE), _rc(TRANS), &M, &N, &K, _rc(A), &LDA, _rc(TAU), _rc(C), &LDC, 
                        _rc(WORK), &LWORK, &INFO );
};

//-----------------------------------------------------------------------
//                          lacpy
//-----------------------------------------------------------------------

template<> BLAS_EXPORT INLINE_TYPE 
void lacpy<s_type>(const char *uplo,i_type m,i_type n, const s_type *a,i_type lda, s_type *b,
                   i_type ldb)
{
    lacpy2<s_type>::eval(uplo,m,n,a,lda,b,ldb);
    //LAPACK_NAME(slacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void lacpy<d_type>(const char *uplo,i_type m,i_type n, const d_type *a,i_type lda, d_type *b,
                   i_type ldb)
{
    lacpy2<d_type>::eval(uplo,m,n,a,lda,b,ldb);
    //LAPACK_NAME(dlacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void lacpy<c_type>(const char *uplo,i_type m,i_type n, const c_type *a,i_type lda, c_type *b,
                   i_type ldb)
{
    lacpy2<c_type>::eval(uplo,m,n,a,lda,b,ldb);
    //LAPACK_NAME(clacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

template<> BLAS_EXPORT INLINE_TYPE 
void lacpy<z_type>(const char *uplo,i_type m,i_type n, const z_type *a,i_type lda, z_type *b,
                   i_type ldb)
{
    //lacpy2<z_type>::eval(uplo,m,n,a,lda,b,ldb);
    LAPACK_NAME(zlacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void clacpy(const char *uplo,i_type m,i_type n, const c_type *a,i_type lda, c_type *b,
            i_type ldb)
{
    lacpy2<c_type>::eval(uplo,m,n,a,lda,b,ldb);
    //LAPACK_NAME(clacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void slacpy(const char *uplo,i_type m,i_type n, const s_type *a,i_type lda, s_type *b,
            i_type ldb)
{
    lacpy2<s_type>::eval(uplo,m,n,a,lda,b,ldb);
    //LAPACK_NAME(slacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void dlacpy(const char *uplo,i_type m,i_type n, const d_type *a,i_type lda, d_type *b,
            i_type ldb)
{
    lacpy2<d_type>::eval(uplo,m,n,a,lda,b,ldb);
    //LAPACK_NAME(dlacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

BLAS_EXPORT INLINE_TYPE
void zlacpy(const char *uplo,i_type m,i_type n, const z_type *a,i_type lda, z_type *b,
            i_type ldb)
{
    lacpy2<z_type>::eval(uplo,m,n,a,lda,b,ldb);
    //LAPACK_NAME(zlacpy)(_rc(uplo),_rc(&m),_rc(&n),_rc(a),_rc(&lda),_rc(b),_rc(&ldb));
};

};};

#undef INLINE_TYPE
