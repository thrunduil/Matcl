/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-blas-lapack/blas/details/blas_utils.h"

namespace matcl { namespace lapack
{

// BLAS 1, 2, 3 functions

// BLAS Level 0
//-----------------------------------------------------------------
//                          utils
//-----------------------------------------------------------------
// print error on std::cerr and continue
BLAS_EXPORT void xerbla(const char* srname, int info );

// maximum of two values
template<class T> inline T maximum(T a, T b)
{
    return (a > b)? a : b;
};

// minimum of two values
template<class T> inline T minimum(T a, T b)
{
    return (a < b)? a : b;
};

// bitwise integer and
inline i_type iand(i_type a, i_type b)
{
    return a & b;
};

// compute x*x
template<class V>
inline V pow_2(const V& arg)    { return arg * arg;};

//-----------------------------------------------------------------
//                          REAL
//-----------------------------------------------------------------
// return real part of complex value, do nothing for real values
template<class T>
inline typename details::real_type<T>::type	real(const T& val);

template<> inline s_type real(const s_type& val)    { return val; };
template<> inline d_type real(const d_type& val)    { return val; };
template<> inline s_type real(const c_type& val)    { return std::real(val); };
template<> inline d_type real(const z_type& val)    { return std::real(val); };

//-----------------------------------------------------------------
//                          IMAG
//-----------------------------------------------------------------
// return imaginary part of complex value, return zero for real values
template<class T>
inline typename details::real_type<T>::type	imag(const T& val);

template<> inline s_type imag(const s_type& )       { return 0; };
template<> inline d_type imag(const d_type& )       { return 0; };
template<> inline s_type imag(const c_type& val)    { return std::imag(val); };
template<> inline d_type imag(const z_type& val)    { return std::imag(val); };

//-----------------------------------------------------------------
//                          ABS
//-----------------------------------------------------------------
// return true if value represent NaN
template<class T> bool isnan(T a);

// return absolute value
template<class T> 
typename details::real_type<T>::type abs(T a);

template<> inline s_type abs(s_type a) { return std::abs(a); };
template<> inline d_type abs(d_type a) { return std::abs(a); };
template<> inline s_type abs(c_type a) { return std::abs(a); };
template<> inline d_type abs(z_type a) { return std::abs(a); };

//-----------------------------------------------------------------
//                          ABS2
//-----------------------------------------------------------------
// return square of absolute value
template<class T> 
inline typename details::real_type<T>::type	abs2(const T& val);

template<> inline s_type abs2(const s_type& val)    { return val*val; };
template<> inline d_type abs2(const d_type& val)    { return val*val; };
template<> inline s_type abs2(const c_type& val)    { return abs2(std::real(val))+abs2(std::imag(val)); };
template<> inline d_type abs2(const z_type& val)    { return abs2(std::real(val))+abs2(std::imag(val)); };

//-----------------------------------------------------------------
//                          CONJ
//-----------------------------------------------------------------
// return complex conjugation, do nothing for real values
template<class T> 
inline T conj(const T& val);

template<> inline s_type conj(const s_type& val)    { return val; };
template<> inline d_type conj(const d_type& val)    { return val; };
template<> inline c_type conj(const c_type& val)    { return c_type(real(val), -imag(val)); };
template<> inline z_type conj(const z_type& val)    { return z_type(real(val), -imag(val)); };

// BLAS Level1
//-----------------------------------------------------------------
//                          SCAL
//-----------------------------------------------------------------
// scales a vector by a constant
template<class V1, class V2> BLAS_EXPORT
typename details::enable_if_equal_or_real1<void,V1,V2>::type
scal(i_type n, V1 a, V2 *x, i_type incx);

BLAS_EXPORT void sscal (i_type n, s_type a, s_type *x, i_type incx);
BLAS_EXPORT void cscal (i_type n, c_type a, c_type *x, i_type incx); 
BLAS_EXPORT void csscal(i_type n, s_type a, c_type *x, i_type incx); 
BLAS_EXPORT void dscal (i_type n, d_type a, d_type *x, i_type incx);
BLAS_EXPORT void zdscal(i_type n, d_type a, z_type *x, i_type incx); 
BLAS_EXPORT void zscal (i_type n, z_type a, z_type *x, i_type incx); 

//-----------------------------------------------------------------
//                          AXPY
//-----------------------------------------------------------------
// constant times a vector plus a vector.
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
axpy(i_type n, V alpha, const V *x, i_type incx, V *y, i_type incy);

BLAS_EXPORT void saxpy(i_type n, s_type alpha, const s_type *x, i_type incx, 
                       s_type *y, i_type incy);
BLAS_EXPORT void caxpy(i_type n, c_type alpha, const c_type *x, i_type incx, 
                       c_type *y, i_type incy); 
BLAS_EXPORT void daxpy(i_type n, d_type alpha, const d_type *x, i_type incx, 
                       d_type *y, i_type incy);
BLAS_EXPORT void zaxpy(i_type n, z_type alpha, const z_type *x, i_type incx, 
                       z_type *y, i_type incy); 

//-----------------------------------------------------------------
//                          AMAX
//-----------------------------------------------------------------
// finds the index of element having max. absolute value.
template<class V> BLAS_EXPORT
typename details::enable_if_valid<i_type,V>::type
amax(i_type n, const V *x, i_type incx);

BLAS_EXPORT i_type  idamax(i_type n, const d_type *x, i_type incx);
BLAS_EXPORT i_type  isamax(i_type n, const s_type *x, i_type incx);
BLAS_EXPORT i_type  icamax(i_type n, const c_type *x, i_type incx); 
BLAS_EXPORT i_type  izamax(i_type n, const z_type *x, i_type incx); 

//-----------------------------------------------------------------
//                          SWAP
//-----------------------------------------------------------------
// interchanges two vectors.
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
swap(i_type n, V *x, i_type incx, V *y, i_type incy);

BLAS_EXPORT void   sswap(i_type n, s_type *x, i_type incx, s_type *y, i_type incy);
BLAS_EXPORT void   cswap(i_type n, c_type *x, i_type incx, c_type *y, i_type incy); 
BLAS_EXPORT void   dswap(i_type n, d_type *x, i_type incx, d_type *y, i_type incy);
BLAS_EXPORT void   zswap(i_type n, z_type *x, i_type incx, z_type *y, i_type incy); 

//-----------------------------------------------------------------
//                          COPY
//-----------------------------------------------------------------
// copies a vector, x, to a vector, y.
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
copy(i_type n, const V *x, i_type incx, V *y, i_type incy);

BLAS_EXPORT void   scopy(i_type n, const s_type *x, i_type incx, s_type *y, 
                         i_type incy);
BLAS_EXPORT void   ccopy(i_type n, const c_type *x, i_type incx, c_type *y, 
                         i_type incy); 
BLAS_EXPORT void   dcopy(i_type n, const d_type *x, i_type incx, d_type *y, 
                         i_type incy);
BLAS_EXPORT void   zcopy(i_type n, const z_type *x, i_type incx, z_type *y, 
                         i_type incy); 

//-----------------------------------------------------------------
//                          DOT
//-----------------------------------------------------------------
// dot product of two vectors.
template<class V> BLAS_EXPORT
typename details::enable_if_valid<V,V>::type
dot(i_type n, const V *x, i_type incx, const V *y, i_type incy);

BLAS_EXPORT s_type sdot(i_type n, const s_type *x, i_type incx, const s_type *y, 
                        i_type incy);
BLAS_EXPORT d_type ddot(i_type n, const d_type *x, i_type incx, const d_type *y, 
                        i_type incy);
BLAS_EXPORT c_type cdotu(i_type n, const c_type *x, i_type incx, const c_type *y,
                         i_type incy);
BLAS_EXPORT z_type zdotu(i_type n, const z_type *x, i_type incx, const z_type *y,
                         i_type incy);

//-----------------------------------------------------------------
//                          DOT
//-----------------------------------------------------------------
// form dot product of two vectors conjugating the first
template<class V> BLAS_EXPORT
typename details::enable_if_valid<V,V>::type
dotc(i_type n, const V *x, i_type incx, const V *y, i_type incy);

BLAS_EXPORT c_type cdotc(i_type n, const c_type *x, i_type incx, 
                         const c_type *y, i_type incy);
BLAS_EXPORT z_type zdotc(i_type n, const z_type *x, i_type incx, 
                         const z_type *y, i_type incy);

//-----------------------------------------------------------------
//                          ROT
//-----------------------------------------------------------------
// apply givens rotation
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
rot(i_type n, V *x, i_type incx, V *y, i_type incy, 
    typename details::real_type<V>::type C, V S);

BLAS_EXPORT void srot(i_type n, s_type *x, i_type incx, s_type *y, 
                      i_type incy, s_type C, s_type S);
BLAS_EXPORT void drot(i_type n, d_type *x, i_type incx, d_type *y, 
                      i_type incy, d_type C, d_type S);
BLAS_EXPORT void crot(i_type n, c_type *x, i_type incx, c_type *y,
                      i_type incy, s_type C, c_type S);
BLAS_EXPORT void zrot(i_type n, z_type *x, i_type incx, z_type *y, 
                      i_type incy, d_type C, z_type S);

// BLAS Level2
//-----------------------------------------------------------------
//                          GER
//-----------------------------------------------------------------
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
geru(i_type m, i_type n, V alpha, const V *x, i_type incx, const V *y, 
     i_type incy, V *a, i_type lda);

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gerc(i_type m, i_type n, V alpha, const V *x, i_type incx, const V *y, 
     i_type incy, V *a, i_type lda);

BLAS_EXPORT void sger(i_type m, i_type n, s_type alpha, const s_type *x, 
                      i_type incx, const s_type *y, i_type incy, s_type *a,
                      i_type lda);
BLAS_EXPORT void dger(i_type m, i_type n, d_type alpha, const d_type *x, 
                      i_type incx, const d_type *y, i_type incy, d_type *a,
                      i_type lda);
BLAS_EXPORT void cgerc(i_type m, i_type n,c_type alpha, const c_type *x, 
                       i_type incx, const c_type *y, i_type incy, c_type *a,
                       i_type lda); 
BLAS_EXPORT void cgeru(i_type m, i_type n,c_type alpha, const c_type *x, 
                       i_type incx, const c_type *y, i_type incy, c_type *a, 
                       i_type lda); 
BLAS_EXPORT void zgerc(i_type m, i_type n,z_type alpha, const z_type *x, 
                       i_type incx, const z_type *y, i_type incy, z_type *a, 
                       i_type lda); 
BLAS_EXPORT void zgeru(i_type m, i_type n,z_type alpha, const z_type *x, 
                       i_type incx, const z_type *y, i_type incy, z_type *a,
                       i_type lda); 

//-----------------------------------------------------------------
//                          GEMV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGEMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading m by n part of the array A must
*           contain the matrix of coefficients.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, m ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry with BETA non-zero, the incremented array Y
*           must contain the vector y. On exit, Y is overwritten by the
*           updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gemv(const char *trans, i_type m, i_type n, V alpha, const V *a, i_type lda,
     const V *x, i_type incx, V beta, V *y, i_type incy);

BLAS_EXPORT void sgemv(const char *trans, i_type m, i_type n, s_type alpha, 
                       const s_type *a, i_type lda, const s_type *x, i_type incx,
                       s_type beta, s_type *y, i_type incy);
BLAS_EXPORT void cgemv(const char *trans, i_type m, i_type n, c_type alpha,
                       const c_type *a, i_type lda, const c_type *x, i_type incx,
                       c_type beta, c_type *y, i_type incy); 
BLAS_EXPORT void dgemv(const char *trans, i_type m, i_type n, d_type alpha,
                       const d_type *a, i_type lda, const d_type *x, i_type incx, 
                       d_type beta, d_type *y, i_type incy);
BLAS_EXPORT void zgemv(const char *trans, i_type m, i_type n, z_type alpha,
                       const z_type *a, i_type lda, const z_type *x, i_type incx, 
                       z_type beta, z_type *y, i_type incy); 

//-----------------------------------------------------------------
//                          TRMV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTRMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular matrix and the strictly lower triangular part of
*           A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular matrix and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u', the diagonal elements of
*           A are not referenced either, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trmv(const char *uplo, const char *transa, const char *diag, i_type n, 
     const V *a, i_type lda, V *b, i_type incx);

BLAS_EXPORT void strmv(const char *uplo, const char *transa, const char *diag,
                       i_type n, const s_type *a, i_type lda, s_type *b, 
                       i_type incx);
BLAS_EXPORT void ctrmv(const char *uplo, const char *transa, const char *diag,
                       i_type n, const c_type *a,  i_type lda, c_type *b, 
                       i_type incx); 
BLAS_EXPORT void dtrmv(const char *uplo, const char *transa, const char *diag, 
                       i_type n, const d_type *a, i_type lda, d_type *b, 
                       i_type incx);
BLAS_EXPORT void ztrmv(const char *uplo, const char *transa, const char *diag, 
                       i_type n, const z_type *a, i_type lda, z_type *b, 
                       i_type incx); 

//-----------------------------------------------------------------
//                          SYMV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSYMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the symmetric matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the symmetric matrix and the strictly
*           upper triangular part of A is not referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
symv(const char *uplo, i_type n, V alpha, const V *a, i_type lda, 
     const V *x, i_type incx, V beta, V *y, i_type incy);

BLAS_EXPORT void ssymv(const char *uplo, i_type n, s_type alpha, 
                       const s_type *a, i_type lda, const s_type *x, 
                       i_type incx, s_type beta, s_type *y, i_type incy);
BLAS_EXPORT void dsymv(const char *uplo, i_type n, d_type alpha,
                       const d_type *a, i_type lda, const d_type *x, 
                       i_type incx, d_type beta, d_type *y, i_type incy);
BLAS_EXPORT void csymv(const char *uplo, i_type n, c_type alpha, 
                       const c_type *a, i_type lda, const c_type *x, 
                       i_type incx, c_type beta, c_type *y, i_type incy);
BLAS_EXPORT void zsymv(const char *uplo, i_type n, z_type alpha, 
                       const z_type *a, i_type lda, const z_type *x, 
                       i_type incx, z_type beta, z_type *y, i_type incy);

//-----------------------------------------------------------------
//                          HEMV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHEMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n hermitian matrix.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the array A is to be referenced as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of A
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of A
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, n ).
*           Before entry with  UPLO = 'U' or 'u', the leading n by n
*           upper triangular part of the array A must contain the upper
*           triangular part of the hermitian matrix and the strictly
*           lower triangular part of A is not referenced.
*           Before entry with UPLO = 'L' or 'l', the leading n by n
*           lower triangular part of the array A must contain the lower
*           triangular part of the hermitian matrix and the strictly
*           upper triangular part of A is not referenced.
*           Note that the imaginary parts of the diagonal elements need
*           not be set and are assumed to be zero.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           max( 1, n ).
*           Unchanged on exit.
*
*  X      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - COMPLEX*16       array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the n
*           element vector y. On exit, Y is overwritten by the updated
*           vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
hemv(const char *uplo, i_type n, V alpha, const V *a, i_type lda,
     const V *x, i_type incx,  V beta, V *y, i_type incy); 

BLAS_EXPORT void chemv(const char *uplo, i_type n, c_type alpha, 
                       const c_type *a, i_type lda, const c_type *x, 
                       i_type incx, c_type beta, c_type *y, i_type incy); 
BLAS_EXPORT void zhemv(const char *uplo, i_type n, z_type alpha, 
                       const z_type *a, i_type lda, const z_type *x, 
                       i_type incx, z_type beta, z_type *y, i_type incy); 

//-----------------------------------------------------------------
//                          GMBV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGBMV  performs one of the matrix-vector operations
*
*     y := alpha*A*x + beta*y,   or   y := alpha*A'*x + beta*y,
*
*  where alpha and beta are scalars, x and y are vectors and A is an
*  m by n band matrix, with kl sub-diagonals and ku super-diagonals.
*
*  Parameters
*  ==========
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   y := alpha*A*x + beta*y.
*
*              TRANS = 'T' or 't'   y := alpha*A'*x + beta*y.
*
*              TRANS = 'C' or 'c'   y := alpha*A'*x + beta*y.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of the matrix A.
*           M must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  KL     - INTEGER.
*           On entry, KL specifies the number of sub-diagonals of the
*           matrix A. KL must satisfy  0 .le. KL.
*           Unchanged on exit.
*
*  KU     - INTEGER.
*           On entry, KU specifies the number of super-diagonals of the
*           matrix A. KU must satisfy  0 .le. KU.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry, the leading ( kl + ku + 1 ) by n part of the
*           array A must contain the matrix of coefficients, supplied
*           column by column, with the leading diagonal of the matrix in
*           row ( ku + 1 ) of the array, the first super-diagonal
*           starting at position 2 in row ku, the first sub-diagonal
*           starting at position 1 in row ( ku + 2 ), and so on.
*           Elements in the array A that do not correspond to elements
*           in the band matrix (such as the top left ku by ku triangle)
*           are not referenced.
*           The following program segment will transfer a band matrix
*           from conventional full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    K = KU + 1 - J
*                    DO 10, I = MAX( 1, J - KU ), MIN( M, J + KL )
*                       A( K + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( kl + ku + 1 ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( m - 1 )*abs( INCX ) ) otherwise.
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta. When BETA is
*           supplied as zero then Y need not be set on input.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( m - 1 )*abs( INCY ) ) when TRANS = 'N' or 'n'
*           and at least
*           ( 1 + ( n - 1 )*abs( INCY ) ) otherwise.
*           Before entry, the incremented array Y must contain the
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
gbmv(const char *trans, i_type m, i_type n, i_type kl, i_type ku, V alpha,
     const V *a, i_type lda, const V *x, i_type incx, V beta, V *y, 
     i_type incy);

BLAS_EXPORT void sgbmv(const char *trans, i_type m, i_type n, i_type kl,
                       i_type ku, s_type alpha,  const s_type *a, i_type lda, 
                       const s_type *x, i_type incx, s_type beta, s_type *y, 
                       i_type incy);
BLAS_EXPORT void cgbmv(const char *trans, i_type m, i_type n, i_type kl,
                       i_type ku, c_type alpha, const c_type *a, i_type lda, 
                       const c_type *x, i_type incx, c_type beta, c_type *y, 
                       i_type incy); 
BLAS_EXPORT void dgbmv(const char *trans, i_type m, i_type n, i_type kl, 
                       i_type ku, d_type alpha, const d_type *a, i_type lda, 
                       const d_type *x, i_type incx, d_type beta, d_type *y, 
                       i_type incy);
BLAS_EXPORT void zgbmv(const char *trans, i_type m, i_type n, i_type kl, 
                       i_type ku, z_type alpha, const z_type *a, i_type lda, 
                       const z_type *x, i_type incx, z_type beta, z_type *y, 
                       i_type incy); 

//-----------------------------------------------------------------
//                          TBMV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTBMV  performs one of the matrix-vector operations
*
*     x := A*x,   or   x := A'*x,
*
*  where x is an n element vector and  A is an n by n unit, or non-unit,
*  upper or lower triangular band matrix, with ( k + 1 ) diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   x := A*x.
*
*              TRANS = 'T' or 't'   x := A'*x.
*
*              TRANS = 'C' or 'c'   x := A'*x.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U' or 'u', K specifies the number of
*           super-diagonals of the matrix A.
*           On entry with UPLO = 'L' or 'l', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer an upper
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer a lower
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Note that when DIAG = 'U' or 'u' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element vector x. On exit, X is overwritten with the
*           tranformed vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
tbmv(const char *uplo, const char *trans, const char *diag, i_type n, 
     i_type k, const V *a, i_type lda, V *x, i_type incx);

BLAS_EXPORT void stbmv(const char *uplo, const char *trans, const char *diag,
                       i_type n, i_type k, const s_type *a, i_type lda, 
                       s_type *x, i_type incx);
BLAS_EXPORT void ctbmv(const char *uplo, const char *trans, const char *diag,
                       i_type n, i_type k, const c_type *a, i_type lda, 
                       c_type *x, i_type incx); 
BLAS_EXPORT void dtbmv(const char *uplo, const char *trans, const char *diag,
                       i_type n, i_type k, const d_type *a, i_type lda, 
                       d_type *x, i_type incx);
BLAS_EXPORT void ztbmv(const char *uplo, const char *trans, const char *diag,
                       i_type n, i_type k, const z_type *a, i_type lda, 
                       z_type *x, i_type incx); 

//-----------------------------------------------------------------
//                          SBMV
//-----------------------------------------------------------------
/*
*  DSBMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n symmetric band matrix, with k super-diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the band matrix A is being supplied as
*           follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  being supplied.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  being supplied.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry, K specifies the number of super-diagonals of the
*           matrix A. K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the symmetric matrix, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer the upper
*           triangular part of a symmetric band matrix from conventional
*           full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the symmetric matrix, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer the lower
*           triangular part of a symmetric band matrix from conventional
*           full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  Y      - DOUBLE PRECISION array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
sbmv(const char *uplo, i_type n, i_type k, V alpha, const V *a, i_type lda,
     const V *x, i_type incx, V beta, V *y, i_type incy); 

BLAS_EXPORT void ssbmv(const char *uplo, i_type n, i_type k, s_type alpha, 
                       const s_type *a, i_type lda, const s_type *x, i_type incx,
                       s_type beta, s_type *y, i_type incy);
BLAS_EXPORT void dsbmv(const char *uplo, i_type n, i_type k, d_type alpha, 
                       const d_type *a, i_type lda, const d_type *x, i_type incx,
                       d_type beta, d_type *y, i_type incy);
BLAS_EXPORT void csbmv(const char *uplo, i_type n, i_type k, c_type alpha, 
                       const c_type *a, i_type lda, const c_type *x, i_type incx,
                       c_type beta, c_type *y, i_type incy);
BLAS_EXPORT void zsbmv(const char *uplo, i_type n, i_type k, z_type alpha, 
                       const z_type *a, i_type lda, const z_type *x, i_type incx, 
                       z_type beta, z_type *y, i_type incy);

//-----------------------------------------------------------------
//                          HBMV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  CHBMV  performs the matrix-vector  operation
*
*     y := alpha*A*x + beta*y,
*
*  where alpha and beta are scalars, x and y are n element vectors and
*  A is an n by n hermitian band matrix, with k super-diagonals.
*
*  Parameters
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the upper or lower
*           triangular part of the band matrix A is being supplied as
*           follows:
*
*              UPLO = 'U' or 'u'   The upper triangular part of A is
*                                  being supplied.
*
*              UPLO = 'L' or 'l'   The lower triangular part of A is
*                                  being supplied.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry, K specifies the number of super-diagonals of the
*           matrix A. K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX         .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX          array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the hermitian matrix, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer the upper
*           triangular part of a hermitian band matrix from conventional
*           full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the hermitian matrix, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer the lower
*           triangular part of a hermitian band matrix from conventional
*           full matrix storage to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Note that the imaginary parts of the diagonal elements need
*           not be set and are assumed to be zero.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - COMPLEX          array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the
*           vector x.
*           Unchanged on exit.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*
*  BETA   - COMPLEX         .
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  Y      - COMPLEX          array of DIMENSION at least
*           ( 1 + ( n - 1 )*abs( INCY ) ).
*           Before entry, the incremented array Y must contain the
*           vector y. On exit, Y is overwritten by the updated vector y.
*
*  INCY   - INTEGER.
*           On entry, INCY specifies the increment for the elements of
*           Y. INCY must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
hbmv(const char *uplo, i_type n, i_type k, V alpha, const V *a, i_type lda, 
     const V *x, i_type incx, V beta, V *y, i_type incy); 

BLAS_EXPORT void chbmv(const char *uplo, i_type n, i_type k, c_type alpha,
                       const c_type *a, i_type lda, const c_type *x, i_type incx,
                       c_type beta, c_type *y, i_type incy); 
BLAS_EXPORT void zhbmv(const char *uplo, i_type n, i_type k, z_type alpha, 
                       const z_type *a, i_type lda, const z_type *x, i_type incx, 
                       z_type beta, z_type *y, i_type incy); 

//-----------------------------------------------------------------
//                          TBSV
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTBSV  solves one of the systems of equations
*
*     A*x = b,   or   A**T*x = b,
*
*  where b and x are n element vectors and A is an n by n unit, or
*  non-unit, upper or lower triangular band matrix, with ( k + 1 )
*  diagonals.
*
*  No test for singularity or near-singularity is included in this
*  routine. Such tests must be performed before calling this routine.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry, TRANS specifies the equations to be solved as
*           follows:
*
*              TRANS = 'N' or 'n'   A*x = b.
*
*              TRANS = 'T' or 't'   A**T*x = b.
*
*              TRANS = 'C' or 'c'   A**T*x = b.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit
*           triangular as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the order of the matrix A.
*           N must be at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with UPLO = 'U' or 'u', K specifies the number of
*           super-diagonals of the matrix A.
*           On entry with UPLO = 'L' or 'l', K specifies the number of
*           sub-diagonals of the matrix A.
*           K must satisfy  0 .le. K.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, n ).
*           Before entry with UPLO = 'U' or 'u', the leading ( k + 1 )
*           by n part of the array A must contain the upper triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row
*           ( k + 1 ) of the array, the first super-diagonal starting at
*           position 2 in row k, and so on. The top left k by k triangle
*           of the array A is not referenced.
*           The following program segment will transfer an upper
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = K + 1 - J
*                    DO 10, I = MAX( 1, J - K ), J
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Before entry with UPLO = 'L' or 'l', the leading ( k + 1 )
*           by n part of the array A must contain the lower triangular
*           band part of the matrix of coefficients, supplied column by
*           column, with the leading diagonal of the matrix in row 1 of
*           the array, the first sub-diagonal starting at position 1 in
*           row 2, and so on. The bottom right k by k triangle of the
*           array A is not referenced.
*           The following program segment will transfer a lower
*           triangular band matrix from conventional full matrix storage
*           to band storage:
*
*                 DO 20, J = 1, N
*                    M = 1 - J
*                    DO 10, I = J, MIN( N, J + K )
*                       A( M + I, J ) = matrix( I, J )
*              10    CONTINUE
*              20 CONTINUE
*
*           Note that when DIAG = 'U' or 'u' the elements of the array A
*           corresponding to the diagonal elements of the matrix are not
*           referenced, but are assumed to be unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. LDA must be at least
*           ( k + 1 ).
*           Unchanged on exit.
*
*  X      - DOUBLE PRECISION array of dimension at least
*           ( 1 + ( n - 1 )*abs( INCX ) ).
*           Before entry, the incremented array X must contain the n
*           element right-hand side vector b. On exit, X is overwritten
*           with the solution vector x.
*
*  INCX   - INTEGER.
*           On entry, INCX specifies the increment for the elements of
*           X. INCX must not be zero.
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
tbsv(const char *uplo, const char *trans, const char *diag, i_type n, i_type k, 
     const V *a, i_type lda, V *x, i_type incx);

BLAS_EXPORT void stbsv(const char *uplo, const char *trans, const char *diag, 
                       i_type n, i_type k, const s_type *a, i_type lda, s_type *x,
                       i_type incx);
BLAS_EXPORT void dtbsv(const char *uplo, const char *trans, const char *diag,
                       i_type n, i_type k, const d_type *a, i_type lda, d_type *x,
                       i_type incx);
BLAS_EXPORT void ctbsv(const char *uplo, const char *trans, const char *diag, 
                       i_type n, i_type k, const c_type *a, i_type lda, c_type *x, 
                       i_type incx);
BLAS_EXPORT void ztbsv(const char *uplo, const char *trans, const char *diag, 
                       i_type n, i_type k, const z_type *a, i_type lda, z_type *x, 
                       i_type incx);

// BLAS Level3
//-----------------------------------------------------------------
//                          TRSM
//-----------------------------------------------------------------
/*  Purpose
*  =======
*
*  DTRSM  solves one of the matrix equations
*
*     op( A )*X = alpha*B,   or   X*op( A ) = alpha*B,
*
*  where alpha is a scalar, X and B are m by n matrices, A is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  The matrix X is overwritten on B.
*
*  Arguments
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry, SIDE specifies whether op( A ) appears on the left
*           or right of X as follows:
*
*              SIDE = 'L' or 'l'   op( A )*X = alpha*B.
*
*              SIDE = 'R' or 'r'   X*op( A ) = alpha*B.
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain  the  right-hand  side  matrix  B,  and  on exit  is
*           overwritten by the solution matrix  X.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trsm(const char *side, const char *uplo, const char *transa, const char *diag,
    i_type m, i_type n, V alpha, const V *a, i_type lda, V *b, i_type ldb);

BLAS_EXPORT void strsm(const char *side, const char *uplo, const char *transa, 
                       const char *diag, i_type m, i_type n, s_type alpha, 
                       const s_type *a, i_type lda, s_type *b, i_type ldb);
BLAS_EXPORT void ctrsm(const char *side, const char *uplo, const char *transa,
                       const char *diag, i_type m, i_type n, c_type alpha, 
                       const c_type *a, i_type lda, c_type *b, i_type ldb); 
BLAS_EXPORT void dtrsm(const char *side, const char *uplo, const char *transa,
                       const char *diag, i_type m, i_type n, d_type alpha, 
                       const d_type *a, i_type lda, d_type *b, i_type ldb);
BLAS_EXPORT void ztrsm(const char *side, const char *uplo, const char *transa,
                       const char *diag, i_type m, i_type n, z_type alpha, 
                       const z_type *a, i_type lda, z_type *b, i_type ldb); 

//-----------------------------------------------------------------
//                          GEMM
//-----------------------------------------------------------------
/*  Purpose
*  =======
*
*  DGEMM  performs one of the matrix-matrix operations
*
*     C := alpha*op( A )*op( B ) + beta*C,
*
*  where  op( X ) is one of
*
*     op( X ) = X   or   op( X ) = X',
*
*  alpha and beta are scalars, and A, B and C are matrices, with op( A )
*  an m by k matrix,  op( B )  a  k by n matrix and  C an m by n matrix.
*
*  Arguments
*  ==========
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n',  op( A ) = A.
*
*              TRANSA = 'T' or 't',  op( A ) = A'.
*
*              TRANSA = 'C' or 'c',  op( A ) = A'.
*
*           Unchanged on exit.
*
*  TRANSB - CHARACTER*1.
*           On entry, TRANSB specifies the form of op( B ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSB = 'N' or 'n',  op( B ) = B.
*
*              TRANSB = 'T' or 't',  op( B ) = B'.
*
*              TRANSB = 'C' or 'c',  op( B ) = B'.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies  the number  of rows  of the  matrix
*           op( A )  and of the  matrix  C.  M  must  be at least  zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N  specifies the number  of columns of the matrix
*           op( B ) and the number of columns of the matrix C. N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry,  K  specifies  the number of columns of the matrix
*           op( A ) and the number of rows of the matrix op( B ). K must
*           be at least  zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANSA = 'N' or 'n',  and is  m  otherwise.
*           Before entry with  TRANSA = 'N' or 'n',  the leading  m by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by m  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program. When  TRANSA = 'N' or 'n' then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, k ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, kb ), where kb is
*           n  when  TRANSB = 'N' or 'n',  and is  k  otherwise.
*           Before entry with  TRANSB = 'N' or 'n',  the leading  k by n
*           part of the array  B  must contain the matrix  B,  otherwise
*           the leading  n by k  part of the array  B  must contain  the
*           matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in the calling (sub) program. When  TRANSB = 'N' or 'n' then
*           LDB must be at least  max( 1, k ), otherwise  LDB must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n  matrix
*           ( alpha*op( A )*op( B ) + beta*C ).
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gemm(const char *transa, const char *transb, i_type m, i_type n, i_type k,
     V alpha, const V *a, i_type lda, const V *b, i_type ldb, V beta, V *c, 
     i_type ldc);

BLAS_EXPORT void sgemm(const char *transa, const char *transb, i_type m, 
                       i_type n, i_type k, s_type alpha, const s_type *a, 
                       i_type lda, const s_type *b, i_type ldb, s_type beta,
                       s_type *c, i_type ldc);
BLAS_EXPORT void cgemm(const char *transa, const char *transb, i_type m, 
                       i_type n, i_type k, c_type alpha, const c_type *a, 
                       i_type lda, const c_type *b, i_type ldb, c_type beta,
                       c_type *c, i_type ldc); 
BLAS_EXPORT void dgemm(const char *transa, const char *transb, i_type m, 
                       i_type n, i_type k, d_type alpha, const d_type *a, 
                       i_type lda, const d_type *b, i_type ldb, d_type beta, 
                       d_type *c, i_type ldc);
BLAS_EXPORT void zgemm(const char *transa, const char *transb, i_type m, 
                       i_type n, i_type k, z_type alpha, const z_type *a, 
                       i_type lda, const z_type *b, i_type ldb, z_type beta, 
                       z_type *c, i_type ldc); 


//-----------------------------------------------------------------
//                          HERK
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHERK  performs one of the hermitian rank k operations
*
*     C := alpha*A*conjg( A' ) + beta*C,
*
*  or
*
*     C := alpha*conjg( A' )*A + beta*C,
*
*  where  alpha and beta  are  real scalars,  C is an  n by n  hermitian
*  matrix and  A  is an  n by k  matrix in the  first case and a  k by n
*  matrix in the second case.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*conjg( A' ) + beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*conjg( A' )*A + beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns   of  the   matrix   A,   and  on   entry   with
*           TRANS = 'C' or 'c',  K  specifies  the number of rows of the
*           matrix A.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION            .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - COMPLEX*16          array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  hermitian matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  hermitian matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*           Note that the imaginary parts of the diagonal elements need
*           not be set,  they are assumed to be zero,  and on exit they
*           are set to zero.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
herk(const char *uplo, const char *trans, i_type n, i_type k, 
     typename details::real_type<V>::type alpha, const V *a, i_type lda, 
     typename details::real_type<V>::type beta, V* c, i_type ldc);

BLAS_EXPORT void cherk(const char *uplo, const char *trans, i_type n,
                       i_type k, s_type alpha, const c_type *a, i_type lda,
                       s_type beta, c_type* c, i_type ldc);
BLAS_EXPORT void zherk(const char *uplo, const char *trans, i_type n, 
                       i_type k, d_type alpha,  const z_type *a, i_type lda, 
                       d_type beta, z_type* c, i_type ldc);

//-----------------------------------------------------------------
//                          SYRK
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSYRK  performs one of the symmetric rank k operations
*
*     C := alpha*A*A**T + beta*C,
*
*  or
*
*     C := alpha*A**T*A + beta*C,
*
*  where  alpha and beta  are scalars, C is an  n by n  symmetric matrix
*  and  A  is an  n by k  matrix in the first case and a  k by n  matrix
*  in the second case.
*
*  Arguments
*  ==========
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of the  array  C  is to be  referenced  as
*           follows:
*
*              UPLO = 'U' or 'u'   Only the  upper triangular part of  C
*                                  is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the  lower triangular part of  C
*                                  is to be referenced.
*
*           Unchanged on exit.
*
*  TRANS  - CHARACTER*1.
*           On entry,  TRANS  specifies the operation to be performed as
*           follows:
*
*              TRANS = 'N' or 'n'   C := alpha*A*A**T + beta*C.
*
*              TRANS = 'T' or 't'   C := alpha*A**T*A + beta*C.
*
*              TRANS = 'C' or 'c'   C := alpha*A**T*A + beta*C.
*
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry,  N specifies the order of the matrix C.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  K      - INTEGER.
*           On entry with  TRANS = 'N' or 'n',  K  specifies  the number
*           of  columns   of  the   matrix   A,   and  on   entry   with
*           TRANS = 'T' or 't' or 'C' or 'c',  K  specifies  the  number
*           of rows of the matrix  A.  K must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           k  when  TRANS = 'N' or 'n',  and is  n  otherwise.
*           Before entry with  TRANS = 'N' or 'n',  the  leading  n by k
*           part of the array  A  must contain the matrix  A,  otherwise
*           the leading  k by n  part of the array  A  must contain  the
*           matrix A.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in  the  calling  (sub)  program.   When  TRANS = 'N' or 'n'
*           then  LDA must be at least  max( 1, n ), otherwise  LDA must
*           be at least  max( 1, k ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry, BETA specifies the scalar beta.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry  with  UPLO = 'U' or 'u',  the leading  n by n
*           upper triangular part of the array C must contain the upper
*           triangular part  of the  symmetric matrix  and the strictly
*           lower triangular part of C is not referenced.  On exit, the
*           upper triangular part of the array  C is overwritten by the
*           upper triangular part of the updated matrix.
*           Before entry  with  UPLO = 'L' or 'l',  the leading  n by n
*           lower triangular part of the array C must contain the lower
*           triangular part  of the  symmetric matrix  and the strictly
*           upper triangular part of C is not referenced.  On exit, the
*           lower triangular part of the array  C is overwritten by the
*           lower triangular part of the updated matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, n ).
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
syrk(const char *uplo, const char *trans, i_type n, i_type k, 
     const V& alpha, const V *a, i_type lda, const V& beta, V* c, i_type ldc);

BLAS_EXPORT void ssyrk(const char *uplo, const char *trans, i_type n, i_type k,
                       const s_type& alpha, const s_type *a, i_type lda, 
                       const s_type& beta, s_type* c, i_type ldc);
BLAS_EXPORT void dsyrk(const char *uplo, const char *trans, i_type n, i_type k, 
                       const d_type& alpha, const d_type *a, i_type lda, 
                       const d_type& beta, d_type* c, i_type ldc);
BLAS_EXPORT void csyrk(const char *uplo, const char *trans, i_type n, i_type k, 
                       const c_type& alpha, const c_type *a, i_type lda, 
                       const c_type& beta, c_type* c, i_type ldc);
BLAS_EXPORT void zsyrk(const char *uplo, const char *trans, i_type n, i_type k, 
                       const z_type& alpha, const z_type *a, i_type lda, 
                       const z_type& beta, z_type* c, i_type ldc);

//-----------------------------------------------------------------
//                          TRMM
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTRMM  performs one of the matrix-matrix operations
*
*     B := alpha*op( A )*B,   or   B := alpha*B*op( A ),
*
*  where  alpha  is a scalar,  B  is an m by n matrix,  A  is a unit, or
*  non-unit,  upper or lower triangular matrix  and  op( A )  is one  of
*
*     op( A ) = A   or   op( A ) = A'.
*
*  Arguments
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE specifies whether  op( A ) multiplies B from
*           the left or right as follows:
*
*              SIDE = 'L' or 'l'   B := alpha*op( A )*B.
*
*              SIDE = 'R' or 'r'   B := alpha*B*op( A ).
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On entry, UPLO specifies whether the matrix A is an upper or
*           lower triangular matrix as follows:
*
*              UPLO = 'U' or 'u'   A is an upper triangular matrix.
*
*              UPLO = 'L' or 'l'   A is a lower triangular matrix.
*
*           Unchanged on exit.
*
*  TRANSA - CHARACTER*1.
*           On entry, TRANSA specifies the form of op( A ) to be used in
*           the matrix multiplication as follows:
*
*              TRANSA = 'N' or 'n'   op( A ) = A.
*
*              TRANSA = 'T' or 't'   op( A ) = A'.
*
*              TRANSA = 'C' or 'c'   op( A ) = A'.
*
*           Unchanged on exit.
*
*  DIAG   - CHARACTER*1.
*           On entry, DIAG specifies whether or not A is unit triangular
*           as follows:
*
*              DIAG = 'U' or 'u'   A is assumed to be unit triangular.
*
*              DIAG = 'N' or 'n'   A is not assumed to be unit
*                                  triangular.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry, M specifies the number of rows of B. M must be at
*           least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of B.  N must be
*           at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry,  ALPHA specifies the scalar  alpha. When  alpha is
*           zero then  A is not referenced and  B need not be set before
*           entry.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, k ), where k is m
*           when  SIDE = 'L' or 'l'  and is  n  when  SIDE = 'R' or 'r'.
*           Before entry  with  UPLO = 'U' or 'u',  the  leading  k by k
*           upper triangular part of the array  A must contain the upper
*           triangular matrix  and the strictly lower triangular part of
*           A is not referenced.
*           Before entry  with  UPLO = 'L' or 'l',  the  leading  k by k
*           lower triangular part of the array  A must contain the lower
*           triangular matrix  and the strictly upper triangular part of
*           A is not referenced.
*           Note that when  DIAG = 'U' or 'u',  the diagonal elements of
*           A  are not referenced either,  but are assumed to be  unity.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA  must be at least  max( 1, m ),  when  SIDE = 'R' or 'r'
*           then LDA must be at least max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry,  the leading  m by n part of the array  B must
*           contain the matrix  B,  and  on exit  is overwritten  by the
*           transformed matrix.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
trmm(const char *side, const char *uplo, const char *transa, const char *diag, 
     i_type m, i_type n, V alpha, const V *a, i_type lda, V *b, i_type ldb);

BLAS_EXPORT void strmm(const char *side, const char *uplo, const char *transa,
                       const char *diag, i_type m, i_type n, s_type alpha, 
                       const s_type *a, i_type lda, s_type *b, i_type ldb);
BLAS_EXPORT void ctrmm(const char *side, const char *uplo, const char *transa, 
                       const char *diag, i_type m, i_type n, c_type alpha, 
                       const c_type *a, i_type lda, c_type *b, i_type ldb); 
BLAS_EXPORT void dtrmm(const char *side, const char *uplo, const char *transa, 
                       const char *diag, i_type m, i_type n, d_type alpha, 
                       const d_type *a, i_type lda, d_type *b, i_type ldb);
BLAS_EXPORT void ztrmm(const char *side, const char *uplo, const char *transa,
                       const char *diag, i_type m, i_type n, z_type alpha, 
                       const z_type *a, i_type lda, z_type *b, i_type ldb); 

//-----------------------------------------------------------------
//                          SYMM
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSYMM  performs one of the matrix-matrix operations
*
*     C := alpha*A*B + beta*C,
*
*  or
*
*     C := alpha*B*A + beta*C,
*
*  where alpha and beta are scalars,  A is a symmetric matrix and  B and
*  C are  m by n matrices.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE  specifies whether  the  symmetric matrix  A
*           appears on the  left or right  in the  operation as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
*
*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of  the  symmetric  matrix   A  is  to  be
*           referenced as follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of the
*                                  symmetric matrix is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of the
*                                  symmetric matrix is to be referenced.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies the number of rows of the matrix  C.
*           M  must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix C.
*           N  must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - DOUBLE PRECISION.
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - DOUBLE PRECISION array of DIMENSION ( LDA, ka ), where ka is
*           m  when  SIDE = 'L' or 'l'  and is  n otherwise.
*           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading m by m upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  m by m  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
*           the array  A  must contain the  symmetric matrix,  such that
*           when  UPLO = 'U' or 'u', the leading n by n upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  symmetric matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  n by n  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  symmetric
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the calling (sub) program.  When  SIDE = 'L' or 'l'  then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least  max( 1, n ).
*           Unchanged on exit.
*
*  B      - DOUBLE PRECISION array of DIMENSION ( LDB, n ).
*           Before entry, the leading  m by n part of the array  B  must
*           contain the matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*  BETA   - DOUBLE PRECISION.
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - DOUBLE PRECISION array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/
template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
symm(const char *side, const char *uplo, i_type m, i_type n, V alpha, 
     const V *a, i_type lda, const V *b, i_type ldb, V beta, V *c, 
     i_type ldc);

BLAS_EXPORT void ssymm(const char *side, const char *uplo, i_type m,
                       i_type n, s_type alpha, const s_type *a, i_type lda, 
                       const s_type *b, i_type ldb, s_type beta, s_type *c, 
                       i_type ldc);
BLAS_EXPORT void csymm(const char *side, const char *uplo, i_type m, i_type n,
                       c_type alpha, const c_type *a, i_type lda,const c_type *b, 
                       i_type ldb, c_type beta, c_type *c, i_type ldc); 
BLAS_EXPORT void dsymm(const char *side, const char *uplo, i_type m, i_type n,
                       d_type alpha, const d_type *a, i_type lda, const d_type *b,
                       i_type ldb, d_type beta, d_type *c, i_type ldc);
BLAS_EXPORT void zsymm(const char *side, const char *uplo, i_type m, i_type n, 
                       z_type alpha, const z_type *a, i_type lda, const z_type *b, 
                       i_type ldb, z_type beta, z_type *c, i_type ldc); 

//-----------------------------------------------------------------
//                          HEMM
//-----------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHEMM  performs one of the matrix-matrix operations
*
*     C := alpha*A*B + beta*C,
*
*  or
*
*     C := alpha*B*A + beta*C,
*
*  where alpha and beta are scalars, A is an hermitian matrix and  B and
*  C are m by n matrices.
*
*  Parameters
*  ==========
*
*  SIDE   - CHARACTER*1.
*           On entry,  SIDE  specifies whether  the  hermitian matrix  A
*           appears on the  left or right  in the  operation as follows:
*
*              SIDE = 'L' or 'l'   C := alpha*A*B + beta*C,
*
*              SIDE = 'R' or 'r'   C := alpha*B*A + beta*C,
*
*           Unchanged on exit.
*
*  UPLO   - CHARACTER*1.
*           On  entry,   UPLO  specifies  whether  the  upper  or  lower
*           triangular  part  of  the  hermitian  matrix   A  is  to  be
*           referenced as follows:
*
*              UPLO = 'U' or 'u'   Only the upper triangular part of the
*                                  hermitian matrix is to be referenced.
*
*              UPLO = 'L' or 'l'   Only the lower triangular part of the
*                                  hermitian matrix is to be referenced.
*
*           Unchanged on exit.
*
*  M      - INTEGER.
*           On entry,  M  specifies the number of rows of the matrix  C.
*           M  must be at least zero.
*           Unchanged on exit.
*
*  N      - INTEGER.
*           On entry, N specifies the number of columns of the matrix C.
*           N  must be at least zero.
*           Unchanged on exit.
*
*  ALPHA  - COMPLEX*16      .
*           On entry, ALPHA specifies the scalar alpha.
*           Unchanged on exit.
*
*  A      - COMPLEX*16       array of DIMENSION ( LDA, ka ), where ka is
*           m  when  SIDE = 'L' or 'l'  and is n  otherwise.
*           Before entry  with  SIDE = 'L' or 'l',  the  m by m  part of
*           the array  A  must contain the  hermitian matrix,  such that
*           when  UPLO = 'U' or 'u', the leading m by m upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  hermitian matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  m by m  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  hermitian
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Before entry  with  SIDE = 'R' or 'r',  the  n by n  part of
*           the array  A  must contain the  hermitian matrix,  such that
*           when  UPLO = 'U' or 'u', the leading n by n upper triangular
*           part of the array  A  must contain the upper triangular part
*           of the  hermitian matrix and the  strictly  lower triangular
*           part of  A  is not referenced,  and when  UPLO = 'L' or 'l',
*           the leading  n by n  lower triangular part  of the  array  A
*           must  contain  the  lower triangular part  of the  hermitian
*           matrix and the  strictly upper triangular part of  A  is not
*           referenced.
*           Note that the imaginary parts  of the diagonal elements need
*           not be set, they are assumed to be zero.
*           Unchanged on exit.
*
*  LDA    - INTEGER.
*           On entry, LDA specifies the first dimension of A as declared
*           in the  calling (sub) program. When  SIDE = 'L' or 'l'  then
*           LDA must be at least  max( 1, m ), otherwise  LDA must be at
*           least max( 1, n ).
*           Unchanged on exit.
*
*  B      - COMPLEX*16       array of DIMENSION ( LDB, n ).
*           Before entry, the leading  m by n part of the array  B  must
*           contain the matrix B.
*           Unchanged on exit.
*
*  LDB    - INTEGER.
*           On entry, LDB specifies the first dimension of B as declared
*           in  the  calling  (sub)  program.   LDB  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*
*  BETA   - COMPLEX*16      .
*           On entry,  BETA  specifies the scalar  beta.  When  BETA  is
*           supplied as zero then C need not be set on input.
*           Unchanged on exit.
*
*  C      - COMPLEX*16       array of DIMENSION ( LDC, n ).
*           Before entry, the leading  m by n  part of the array  C must
*           contain the matrix  C,  except when  beta  is zero, in which
*           case C need not be set on entry.
*           On exit, the array  C  is overwritten by the  m by n updated
*           matrix.
*
*  LDC    - INTEGER.
*           On entry, LDC specifies the first dimension of C as declared
*           in  the  calling  (sub)  program.   LDC  must  be  at  least
*           max( 1, m ).
*           Unchanged on exit.
*/

template<class V> BLAS_EXPORT 
typename details::enable_if_valid<void,V>::type
hemm(const char *side, const char *uplo, i_type m, i_type n, V alpha, 
     const V *a, i_type lda, const V *b, i_type ldb, V beta, V *c, i_type ldc); 

BLAS_EXPORT void chemm(const char *side, const char *uplo, i_type m, i_type n, 
                       c_type alpha, const c_type *a, i_type lda, const c_type *b, 
                       i_type ldb, c_type beta, c_type *c, i_type ldc); 
BLAS_EXPORT void zhemm(const char *side, const char *uplo, i_type m, i_type n,
                       z_type alpha, const z_type *a, i_type lda, const z_type *b,
                       i_type ldb, z_type beta, z_type *c, i_type ldc); 

};};
