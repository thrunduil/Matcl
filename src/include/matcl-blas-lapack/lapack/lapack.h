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

#pragma once

#include "matcl-blas-lapack/blas/details/blas_utils.h"

namespace matcl { namespace lapack
{

//-----------------------------------------------------------------------
//                          ILAENV
//-----------------------------------------------------------------------
// get optimal value of given parameter
BLAS_EXPORT i_type    ilaenv(i_type ispec, const char *name, const char *opts, 
                           i_type n1, i_type n2, i_type n3, i_type n4);

//-----------------------------------------------------------------------
//                          laswp
//-----------------------------------------------------------------------
/* Purpose
*  =======
*
*  DLASWP performs a series of row interchanges on the matrix A.
*  One row interchange is initiated for each of rows K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix of column dimension N to which the row
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a row interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a row interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies rows K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
laswp(i_type n, V *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,
      i_type incx);

BLAS_EXPORT void    claswp(i_type n, c_type *a,i_type lda,i_type k1,i_type k2,
                           const i_type *ipiv,i_type incx);
BLAS_EXPORT void    slaswp(i_type n, s_type *a,i_type lda,i_type k1,i_type k2,
                           const i_type *ipiv,i_type incx);
BLAS_EXPORT void    dlaswp(i_type n, d_type *a,i_type lda,i_type k1,i_type k2,
                           const i_type *ipiv,i_type incx);
BLAS_EXPORT void    zlaswp(i_type n, z_type *a,i_type lda,i_type k1,i_type k2,
                           const i_type *ipiv,i_type incx);

//-----------------------------------------------------------------------
//                          LACPY
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLACPY copies all or part of a two-dimensional matrix A to another
*  matrix B.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be copied to B.
*          = 'U':      Upper triangular part
*          = 'L':      Lower triangular part
*          Otherwise:  All of the matrix A
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.  If UPLO = 'U', only the upper triangle
*          or trapezoid is accessed; if UPLO = 'L', only the lower
*          triangle or trapezoid is accessed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  B       (output) DOUBLE PRECISION array, dimension (LDB,N)
*          On exit, B = A in the locations specified by UPLO.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,M).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
lacpy(const char *uplo,i_type m,i_type n, const V *a,i_type lda, V *b,
      i_type ldb);

BLAS_EXPORT void    clacpy(const char *uplo,i_type m,i_type n, const c_type *a,
                           i_type lda, c_type *b,i_type ldb);
BLAS_EXPORT void    slacpy(const char *uplo,i_type m,i_type n, const s_type *a,
                           i_type lda, s_type *b,i_type ldb);
BLAS_EXPORT void    dlacpy(const char *uplo,i_type m,i_type n, const d_type *a,
                           i_type lda, d_type *b,i_type ldb);
BLAS_EXPORT void    zlacpy(const char *uplo,i_type m,i_type n, const z_type *a,
                           i_type lda, z_type *b,i_type ldb);

//-----------------------------------------------------------------------
//                          LASET
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLASET initializes an m-by-n matrix A to BETA on the diagonal and
*  ALPHA on the offdiagonals.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the part of the matrix A to be set.
*          = 'U':      Upper triangular part is set; the strictly lower
*                      triangular part of A is not changed.
*          = 'L':      Lower triangular part is set; the strictly upper
*                      triangular part of A is not changed.
*          Otherwise:  All of the matrix A is set.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  ALPHA   (input) DOUBLE PRECISION
*          The constant to which the offdiagonal elements are to be set.
*
*  BETA    (input) DOUBLE PRECISION
*          The constant to which the diagonal elements are to be set.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On exit, the leading m-by-n submatrix of A is set as follows:
*
*          if UPLO = 'U', A(i,j) = ALPHA, 1<=i<=j-1, 1<=j<=n,
*          if UPLO = 'L', A(i,j) = ALPHA, j+1<=i<=m, 1<=j<=n,
*          otherwise,     A(i,j) = ALPHA, 1<=i<=m, 1<=j<=n, i.ne.j,
*
*          and, for all UPLO, A(i,i) = BETA, 1<=i<=min(m,n).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
laset(const char *uplo, i_type m, i_type n, V alpha, V beta, V* a, i_type lda);

BLAS_EXPORT void    claset(const char *uplo, i_type m, i_type n, c_type alpha, 
                           c_type beta, c_type* a, i_type lda);
BLAS_EXPORT void    slaset(const char *uplo, i_type m, i_type n, s_type alpha, 
                           s_type beta, s_type* a, i_type lda);
BLAS_EXPORT void    dlaset(const char *uplo, i_type m, i_type n, d_type alpha, 
                           d_type beta, d_type* a, i_type lda);
BLAS_EXPORT void    zlaset(const char *uplo, i_type m, i_type n, z_type alpha, 
                           z_type beta, z_type* a, i_type lda);

//-----------------------------------------------------------------------
//                          LARTG
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLARTG generate a plane rotation so that
*
*     [  CS  SN  ]  .  [ F ]  =  [ R ]   where CS**2 + SN**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a slower, more accurate version of the BLAS1 routine DROTG,
*  with the following other differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations (saves work in DBDSQR when
*        there are zeros on the diagonal).
*
*  If F exceeds G in magnitude, CS will be positive.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The first component of vector to be rotated.
*
*  G       (input) DOUBLE PRECISION
*          The second component of vector to be rotated.
*
*  CS      (output) DOUBLE PRECISION
*          The cosine of the rotation.
*
*  SN      (output) DOUBLE PRECISION
*          The sine of the rotation.
*
*  R       (output) DOUBLE PRECISION
*          The nonzero component of the rotated vector.
*
*  This version has a few statements commented out for thread safety
*  (machine parameters are computed on each entry). 10 feb 03, SJH.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
lartg(V f, V g, typename details::real_type<V>::type *cs, V *sn,  V *r);

BLAS_EXPORT void    clartg(c_type f, c_type g, s_type *cs, c_type *sn,  c_type *r);
BLAS_EXPORT void    dlartg(d_type f, d_type g, d_type *cs, d_type *sn, d_type *r);
BLAS_EXPORT void    slartg(s_type f, s_type g, s_type *cs, s_type *sn, s_type *r);
BLAS_EXPORT void    zlartg(z_type f, z_type g, d_type * cs, z_type *sn, z_type *r);

//-----------------------------------------------------------------------
//                          LANGE
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLANGE  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  real matrix A.
*
*  Description
*  ===========
*
*  DLANGE returns the value
*
*     DLANGE = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANGE as described
*          above.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.  When M = 0,
*          DLANGE is set to zero.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.  When N = 0,
*          DLANGE is set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The m by n matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(M,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= M when NORM = 'I'; otherwise, WORK is not
*          referenced.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<typename details::real_type<V>::type,V>::type
lange(const char *norm,i_type m,i_type n,const V *a,i_type lda,
      typename details::real_type<V>::type *work);

BLAS_EXPORT s_type clange(const char *norm,i_type m,i_type n,const c_type *a,
                          i_type lda,s_type *work);
BLAS_EXPORT s_type slange(const char *norm,i_type m,i_type n,const s_type *a,
                          i_type lda,s_type *work);
BLAS_EXPORT d_type dlange(const char *norm,i_type m,i_type n,const d_type *a,
                          i_type lda,d_type *work);
BLAS_EXPORT d_type zlange(const char *norm,i_type m,i_type n,const z_type *a,
                          i_type lda,d_type *work);

//-----------------------------------------------------------------------
//                          LAMCH
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by DLAMCH:
*          = 'E' or 'e',   DLAMCH := eps
*          = 'S' or 's ,   DLAMCH := sfmin
*          = 'B' or 'b',   DLAMCH := base
*          = 'P' or 'p',   DLAMCH := eps*base
*          = 'N' or 'n',   DLAMCH := t
*          = 'R' or 'r',   DLAMCH := rnd
*          = 'M' or 'm',   DLAMCH := emin
*          = 'U' or 'u',   DLAMCH := rmin
*          = 'L' or 'l',   DLAMCH := emax
*          = 'O' or 'o',   DLAMCH := rmax
*
*          where
*
*          eps   = relative machine precision
*          sfmin = safe minimum, such that 1/sfmin does not overflow
*          base  = base of the machine
*          prec  = eps*base
*          t     = number of (base) digits in the mantissa
*          rnd   = 1.0 when rounding occurs in addition, 0.0 otherwise
*          emin  = minimum exponent before (gradual) underflow
*          rmin  = underflow threshold - base**(emin-1)
*          emax  = largest exponent before overflow
*          rmax  = overflow threshold  - (base**emax)*(1-eps)
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<typename details::real_type<V>::type,V>::type
lamch(const char *cmach);

BLAS_EXPORT s_type  slamch(const char *cmach);
BLAS_EXPORT d_type  dlamch(const char *cmach);

//-----------------------------------------------------------------------
//                          LAE2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLAE2  computes the eigenvalues of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, and RT2
*  is the eigenvalue of smaller absolute value.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) and (2,1) elements of the 2-by-2 matrix.
*
*  C       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  RT1     (output) DOUBLE PRECISION
*          The eigenvalue of larger absolute value.
*
*  RT2     (output) DOUBLE PRECISION
*          The eigenvalue of smaller absolute value.
*
*  Further Details
*  ===============
*
*  RT1 is accurate to a few ulps barring over/underflow.
*
*  RT2 may be inaccurate if there is massive cancellation in the
*  determinant A*C-B*B; higher precision or correctly rounded or
*  correctly truncated arithmetic would be needed to compute RT2
*  accurately in all cases.
*
*  Overflow is possible only if RT1 is within a factor of 5 of overflow.
*  Underflow is harmless if the input data is 0 or exceeds
*     underflow_threshold / macheps.
*/
template<class V>
typename details::enable_if_valid<void,V>::type
lae2(typename details::real_type<V>::type a, V b, typename details::real_type<V>::type c, 
     typename details::real_type<V>::type& rt1, typename details::real_type<V>::type& rt2);

//-----------------------------------------------------------------------
//                          LACGV
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZLACGV conjugates a complex vector of length N.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The length of the vector X.  N >= 0.
*
*  X       (input/output) COMPLEX*16 array, dimension
*                         (1+(N-1)*abs(INCX))
*          On entry, the vector of length N to be conjugated.
*          On exit, X is overwritten with conjg(X).
*
*  INCX    (input) INTEGER
*          The spacing between successive elements of X.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
lacgv(i_type N, V* X, i_type INCX);

BLAS_EXPORT void clacgv(i_type N, c_type* X, i_type INCX);
BLAS_EXPORT void dlacgv(i_type N, d_type* X, i_type INCX);
BLAS_EXPORT void slacgv(i_type N, s_type* X, i_type INCX);
BLAS_EXPORT void zlacgv(i_type N, z_type* X, i_type INCX);

//-----------------------------------------------------------------------
//                          GBTRF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGBTRF computes an LU factorization of a real m-by-n band matrix A
*  using partial pivoting with row interchanges.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(kl+ku+1+i-j,j) = A(i,j) for max(1,j-ku)<=i<=min(m,j+kl)
*
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (min(M,N))
*          The pivot indices; for 1 <= i <= min(M,N), row i of the
*          matrix was interchanged with row IPIV(i).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = +i, U(i,i) is exactly zero. The factorization
*               has been completed, but the factor U is exactly
*               singular, and division by zero will occur if it is used
*               to solve a system of equations.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gbtrf(i_type m,i_type n,i_type kl,i_type ku, V *ab,i_type ldab,i_type *ipiv,
      i_type *info);

BLAS_EXPORT void    cgbtrf(i_type m,i_type n,i_type kl,i_type ku, c_type *ab,
                           i_type ldab,i_type *ipiv,i_type *info);
BLAS_EXPORT void    sgbtrf(i_type m,i_type n,i_type kl,i_type ku, s_type *ab,
                           i_type ldab,i_type *ipiv,i_type *info);
BLAS_EXPORT void    dgbtrf(i_type m,i_type n,i_type kl,i_type ku, d_type *ab,
                           i_type ldab,i_type *ipiv,i_type *info);
BLAS_EXPORT void    zgbtrf(i_type m,i_type n,i_type kl,i_type ku, z_type *ab,
                           i_type ldab,i_type *ipiv,i_type *info);

//-----------------------------------------------------------------------
//                          GESV
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGESV computes the solution to a real system of linear equations
*     A * X = B,
*  where A is an N-by-N matrix and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as
*     A = P * L * U,
*  where P is a permutation matrix, L is unit lower triangular, and U is
*  upper triangular.  The factored form of A is then used to solve the
*  system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N coefficient matrix A.
*          On exit, the factors L and U from the factorization
*          A = P*L*U; the unit diagonal elements of L are not stored.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS matrix of right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, so the solution could not be computed.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gesv(i_type n,i_type nrhs, V *a,i_type lda,i_type *ipiv,V* b,i_type ldb,
     i_type *info);

BLAS_EXPORT void    dgesv(i_type n,i_type nrhs, d_type *a,i_type lda,i_type *ipiv,
                          d_type* b,i_type ldb,i_type *info);
BLAS_EXPORT void    cgesv(i_type n,i_type nrhs, c_type *a,i_type lda,i_type *ipiv,
                          c_type* b,i_type ldb,i_type *info);
BLAS_EXPORT void    sgesv(i_type n,i_type nrhs, s_type *a,i_type lda,i_type *ipiv,
                          s_type* b,i_type ldb,i_type *info);
BLAS_EXPORT void    zgesv(i_type n,i_type nrhs, z_type *a,i_type lda,i_type *ipiv,
                          z_type* b,i_type ldb,i_type *info);

//-----------------------------------------------------------------------
//                          GETRS
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The factors L and U from the factorization A = P*L*U
*          as computed by DGETRF.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
getrs(const char *trans,i_type n,i_type nrhs,const V *a,i_type lda,
      const i_type *ipiv,V *b,i_type ldb, i_type *info);

BLAS_EXPORT void    cgetrs(const char *trans,i_type n,i_type nrhs,
                           const c_type *a,i_type lda, const i_type *ipiv,
                           c_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    sgetrs(const char *trans,i_type n,i_type nrhs,
                           const s_type *a,i_type lda, const i_type *ipiv,
                           s_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    dgetrs(const char *trans,i_type n,i_type nrhs,
                           const d_type *a,i_type lda, const i_type *ipiv,
                           d_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    zgetrs(const char *trans,i_type n,i_type nrhs,
                           const z_type *a,i_type lda, const i_type *ipiv,
                           z_type *b,i_type ldb,i_type *info);

//-----------------------------------------------------------------------
//                          TRTRS
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTRTRS solves a triangular system of the form
*
*     A * X = B  or  A**T * X = B,
*
*  where A is a triangular matrix of order N, and B is an N-by-NRHS
*  matrix.  A check is made to verify that A is nonsingular.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The triangular matrix A.  If UPLO = 'U', the leading N-by-N
*          upper triangular part of the array A contains the upper
*          triangular matrix, and the strictly lower triangular part of
*          A is not referenced.  If UPLO = 'L', the leading N-by-N lower
*          triangular part of the array A contains the lower triangular
*          matrix, and the strictly upper triangular part of A is not
*          referenced.  If DIAG = 'U', the diagonal elements of A are
*          also not referenced and are assumed to be 1.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, if INFO = 0, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, the i-th diagonal element of A is zero,
*               indicating that the matrix is singular and the solutions
*               X have not been computed.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trtrs(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
      const V *a,i_type lda,V *b,i_type ldb,i_type *info);

BLAS_EXPORT void    dtrtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type nrhs, const d_type *a,i_type lda,
                           d_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    ctrtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type nrhs, const c_type *a,i_type lda,
                           c_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    strtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type nrhs, const s_type *a,i_type lda,
                           s_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    ztrtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type nrhs, const z_type *a,i_type lda,
                           z_type *b,i_type ldb,i_type *info);

//-----------------------------------------------------------------------
//                          GBSV
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGBSV computes the solution to a real system of linear equations
*  A * X = B, where A is a band matrix of order N with KL subdiagonals
*  and KU superdiagonals, and X and B are N-by-NRHS matrices.
*
*  The LU decomposition with partial pivoting and row interchanges is
*  used to factor A as A = L * U, where L is a product of permutation
*  and unit lower triangular matrices with KL subdiagonals, and U is
*  upper triangular with KL+KU superdiagonals.  The factored form of A
*  is then used to solve the system of equations A * X = B.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of linear equations, i.e., the order of the
*          matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the matrix A in band storage, in rows KL+1 to
*          2*KL+KU+1; rows 1 to KL of the array need not be set.
*          The j-th column of A is stored in the j-th column of the
*          array AB as follows:
*          AB(KL+KU+1+i-j,j) = A(i,j) for max(1,j-KU)<=i<=min(N,j+KL)
*          On exit, details of the factorization: U is stored as an
*          upper triangular band matrix with KL+KU superdiagonals in
*          rows 1 to KL+KU+1, and the multipliers used during the
*          factorization are stored in rows KL+KU+2 to 2*KL+KU+1.
*          See below for further details.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices that define the permutation matrix P;
*          row i of the matrix was interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the N-by-NRHS right hand side matrix B.
*          On exit, if INFO = 0, the N-by-NRHS solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero.  The factorization
*                has been completed, but the factor U is exactly
*                singular, and the solution has not been computed.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gbsv(i_type n,i_type kl,i_type ku,i_type nrhs,V *ab,i_type ldab, i_type *ipiv,
     V *b,i_type ldb,i_type *info);

BLAS_EXPORT void    cgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,c_type *ab,
                          i_type ldab, i_type *ipiv,c_type *b,i_type ldb,
                          i_type *info);
BLAS_EXPORT void    sgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,s_type *ab,
                          i_type ldab, i_type *ipiv,s_type *b,i_type ldb,
                          i_type *info);
BLAS_EXPORT void    dgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,d_type *ab,
                          i_type ldab, i_type *ipiv,d_type *b,i_type ldb,
                          i_type *info);
BLAS_EXPORT void    zgbsv(i_type n,i_type kl,i_type ku,i_type nrhs,z_type *ab,
                          i_type ldab, i_type *ipiv,z_type *b,i_type ldb,
                          i_type *info);

//-----------------------------------------------------------------------
//                          GBTRS
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGBTRS solves a system of linear equations
*     A * X = B  or  A' * X = B
*  with a general band matrix A using the LU factorization computed
*  by DGBTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KL      (input) INTEGER
*          The number of subdiagonals within the band of A.  KL >= 0.
*
*  KU      (input) INTEGER
*          The number of superdiagonals within the band of A.  KU >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          Details of the LU factorization of the band matrix A, as
*          computed by DGBTRF.  U is stored as an upper triangular band
*          matrix with KL+KU superdiagonals in rows 1 to KL+KU+1, and
*          the multipliers used during the factorization are stored in
*          rows KL+KU+2 to 2*KL+KU+1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= 2*KL+KU+1.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= N, row i of the matrix was
*          interchanged with row IPIV(i).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gbtrs(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const V *ab,
      i_type ldab,const i_type *ipiv,V *b,i_type ldb,i_type *info);

BLAS_EXPORT void    cgbtrs(const char *trans,i_type n,i_type kl,i_type ku,
                           i_type nrhs,const c_type *ab, i_type ldab,
                           const i_type *ipiv,c_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    sgbtrs(const char *trans,i_type n,i_type kl,i_type ku,
                           i_type nrhs,const s_type *ab, i_type ldab,
                           const i_type *ipiv,s_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    dgbtrs(const char *trans,i_type n,i_type kl,i_type ku,
                           i_type nrhs,const d_type *ab, i_type ldab,
                           const i_type *ipiv,d_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    zgbtrs(const char *trans,i_type n,i_type kl,i_type ku,
                           i_type nrhs,const z_type *ab, i_type ldab,
                           const i_type *ipiv,z_type *b,i_type ldb,i_type *info);

//-----------------------------------------------------------------------
//                          TBTRS
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTBTRS solves a triangular system of the form
*
*     A * X = B  or  A**T * X = B,
*
*  where A is a triangular band matrix of order N, and B is an
*  N-by NRHS matrix.  A check is made to verify that A is nonsingular.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form the system of equations:
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A**T * X = B  (Transpose)
*          = 'C':  A**H * X = B  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals or subdiagonals of the
*          triangular band matrix A.  KD >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  AB      (input) DOUBLE PRECISION array, dimension (LDAB,N)
*          The upper or lower triangular band matrix A, stored in the
*          first kd+1 rows of AB.  The j-th column of A is stored
*          in the j-th column of the array AB as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*          If DIAG = 'U', the diagonal elements of A are not referenced
*          and are assumed to be 1.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side matrix B.
*          On exit, if INFO = 0, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the i-th diagonal element of A is zero,
*                indicating that the matrix is singular and the
*                solutions X have not been computed.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
tbtrs(const char *uplo,const char *trans,const char *diag,i_type n,
      i_type kd,i_type nrhs, const V *ab,i_type ldab,V *b,i_type ldb,
      i_type *info);

BLAS_EXPORT void    ctbtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type kd,i_type nrhs, const c_type *ab,
                           i_type ldab,c_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    stbtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type kd,i_type nrhs, const s_type *ab,
                           i_type ldab,s_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    dtbtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type kd,i_type nrhs, const d_type *ab,
                           i_type ldab,d_type *b,i_type ldb,i_type *info);
BLAS_EXPORT void    ztbtrs(const char *uplo,const char *trans,const char *diag,
                           i_type n,i_type kd,i_type nrhs, const z_type *ab,
                           i_type ldab,z_type *b,i_type ldb,i_type *info);

//-----------------------------------------------------------------------
//                          GEES
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
* GEES computes for an N-by-N real nonsymmetric matrix A, the
* eigenvalues, the real Schur form T, and, optionally, the matrix of
* Schur vectors Z.  This gives the Schur factorization A = Z*T*(Z**T).
*
* Optionally, it also orders the eigenvalues on the diagonal of the
* real Schur form so that selected eigenvalues are at the top left.
* The leading columns of Z then form an orthonormal basis for the
* invariant subspace corresponding to the selected eigenvalues.
*
* A matrix is in real Schur form if it is upper quasi-triangular with
* 1-by-1 and 2-by-2 blocks. 2-by-2 blocks will be standardized in the
* form
*         [  a  b  ]
*         [  c  a  ]
*
* where b*c < 0. The eigenvalues of such a block are a +- sqrt(bc).
*
*  Arguments:
*  ==========
*
* [in] JOBVS
*          JOBVS is CHARACTER*1
*          = 'N': Schur vectors are not computed;
*          = 'V': Schur vectors are computed.
*
* [in] SORT
*          SORT is CHARACTER*1
*          Specifies whether or not to order the eigenvalues on the
*          diagonal of the Schur form.
*          = 'N': Eigenvalues are not ordered;
*          = 'S': Eigenvalues are ordered (see SELECT).
*
* [in] SELECT
*          SELECT is a LOGICAL FUNCTION of two DOUBLE PRECISION arguments
*          SELECT must be declared EXTERNAL in the calling subroutine.
*          If SORT = 'S', SELECT is used to select eigenvalues to sort
*          to the top left of the Schur form.
*          If SORT = 'N', SELECT is not referenced.
*          An eigenvalue WR(j)+sqrt(-1)*WI(j) is selected if
*          SELECT(WR(j),WI(j)) is true; i.e., if either one of a complex
*          conjugate pair of eigenvalues is selected, then both complex
*          eigenvalues are selected.
*          Note that a selected complex eigenvalue may no longer
*          satisfy SELECT(WR(j),WI(j)) = .TRUE. after ordering, since
*          ordering may change the value of complex eigenvalues
*          (especially if the eigenvalue is ill-conditioned); in this
*          case INFO is set to N+2 (see INFO below).
*
* [in] N
*          N is INTEGER
*          The order of the matrix A. N >= 0.
*
* [in,out] A
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N matrix A.
*          On exit, A has been overwritten by its real Schur form T.
*
* [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
* [out] SDIM
*          SDIM is INTEGER
*          If SORT = 'N', SDIM = 0.
*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*                         for which SELECT is true. (Complex conjugate
*                         pairs for which SELECT is true for either
*                         eigenvalue count as 2.)
*
* [out] WR
*          WR is DOUBLE PRECISION array, dimension (N)
*
* [out] WI
*          WI is DOUBLE PRECISION array, dimension (N)
*          WR and WI contain the real and imaginary parts,
*          respectively, of the computed eigenvalues in the same order
*          that they appear on the diagonal of the output Schur form T.
*          Complex conjugate pairs of eigenvalues will appear
*          consecutively with the eigenvalue having the positive
*          imaginary part first.
*
* [out] VS
*          VS is DOUBLE PRECISION array, dimension (LDVS,N)
*          If JOBVS = 'V', VS contains the orthogonal matrix Z of Schur
*          vectors.
*          If JOBVS = 'N', VS is not referenced.
*
* [in] LDVS
*          LDVS is INTEGER
*          The leading dimension of the array VS.  LDVS >= 1; if
*          JOBVS = 'V', LDVS >= N.
*
* [out] WORK
*          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) contains the optimal LWORK.
*
* [in] LWORK
*          LWORK is INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,3*N).
*          For good performance, LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
* [out] BWORK
*          BWORK is LOGICAL array, dimension (N)
*          Not referenced if SORT = 'N'.
*
* [out] INFO
*          INFO is INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*          > 0: if INFO = i, and i is
*             <= N: the QR algorithm failed to compute all the
*                   eigenvalues; elements 1:ILO-1 and i+1:N of WR and WI
*                   contain those eigenvalues which have converged; if
*                   JOBVS = 'V', VS contains the matrix which reduces A
*                   to its partially converged Schur form.
*             = N+1: the eigenvalues could not be reordered because some
*                   eigenvalues were too close to separate (the problem
*                   is very ill-conditioned);
*             = N+2: after reordering, roundoff changed values of some
*                   complex eigenvalues so that leading eigenvalues in
*                   the Schur form no longer satisfy SELECT=.TRUE.  This
*                   could also be caused by underflow due to scaling.
*
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gees(const char *jobvs,const char *sort,sel_fun selctg, i_type n, V *a,
     i_type lda,i_type *sdim, typename details::complex_type<V>::type *eig, 
     V *vs,i_type ldvs, V *work,i_type lwork,l_type *bwork,i_type *info);

BLAS_EXPORT void    cgees(const char *jobvs,const char *sort,sel_fun selctg,
                          i_type n, c_type *a,i_type lda,i_type *sdim, 
                          c_type *eig, c_type *vs,i_type ldvs, 
                          c_type *work,i_type lwork,s_type *rwork,l_type *bwork,
                          i_type *info);
BLAS_EXPORT void    sgees(const char *jobvs,const char *sort,sel_fun selctg,
                          i_type n,s_type *a,i_type lda,i_type *sdim,
                          s_type *eigr,s_type *eigi,s_type *vs,i_type ldvs,
                          s_type *work,i_type lwork,l_type *bwork,i_type *info);
BLAS_EXPORT void    dgees(const char *jobvs,const char *sort,sel_fun delctg,
                          i_type n,d_type *a,i_type lda,i_type *sdim,
                          d_type *eigr,d_type *eigi,d_type *vs,i_type ldvs,
                          d_type *work,i_type lwork,l_type *bwork,i_type *info);
BLAS_EXPORT void    zgees(const char *jobvs,const char *sort,sel_fun delctg,
                          i_type n, z_type *a,i_type lda,i_type *sdim, 
                          z_type *eig, z_type *vs,i_type ldvs, z_type *work,
                          i_type lwork,d_type *rwork,l_type *bwork,i_type *info);

//-----------------------------------------------------------------------
//                          GGES
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGGES computes for a pair of N-by-N real nonsymmetric matrices (A,B),
*  the generalized eigenvalues, the generalized real Schur form (S,T),
*  optionally, the left and/or right matrices of Schur vectors (VSL and
*  VSR). This gives the generalized Schur factorization
*
*           (A,B) = ( (VSL)*S*(VSR)**T, (VSL)*T*(VSR)**T )
*
*  Optionally, it also orders the eigenvalues so that a selected cluster
*  of eigenvalues appears in the leading diagonal blocks of the upper
*  quasi-triangular matrix S and the upper triangular matrix T.The
*  leading columns of VSL and VSR then form an orthonormal basis for the
*  corresponding left and right eigenspaces (deflating subspaces).
*
*  (If only the generalized eigenvalues are needed, use the driver
*  DGGEV instead, which is faster.)
*
*  A generalized eigenvalue for a pair of matrices (A,B) is a scalar w
*  or a ratio alpha/beta = w, such that  A - w*B is singular.  It is
*  usually represented as the pair (alpha,beta), as there is a
*  reasonable interpretation for beta=0 or both being zero.
*
*  A pair of matrices (S,T) is in generalized real Schur form if T is
*  upper triangular with non-negative diagonal and S is block upper
*  triangular with 1-by-1 and 2-by-2 blocks.  1-by-1 blocks correspond
*  to real generalized eigenvalues, while 2-by-2 blocks of S will be
*  "standardized" by making the corresponding elements of T have the
*  form:
*          [  a  0  ]
*          [  0  b  ]
*
*  and the pair of corresponding 2-by-2 blocks in S and T will have a
*  complex conjugate pair of generalized eigenvalues.
*
*
*  Arguments
*  =========
*
*  JOBVSL  (input) CHARACTER*1
*          = 'N':  do not compute the left Schur vectors;
*          = 'V':  compute the left Schur vectors.
*
*  JOBVSR  (input) CHARACTER*1
*          = 'N':  do not compute the right Schur vectors;
*          = 'V':  compute the right Schur vectors.
*
*  SORT    (input) CHARACTER*1
*          Specifies whether or not to order the eigenvalues on the
*          diagonal of the generalized Schur form.
*          = 'N':  Eigenvalues are not ordered;
*          = 'S':  Eigenvalues are ordered (see SELCTG);
*
*  SELCTG  (external procedure) LOGICAL FUNCTION of three DOUBLE PRECISION arguments
*          SELCTG must be declared EXTERNAL in the calling subroutine.
*          If SORT = 'N', SELCTG is not referenced.
*          If SORT = 'S', SELCTG is used to select eigenvalues to sort
*          to the top left of the Schur form.
*          An eigenvalue (ALPHAR(j)+ALPHAI(j))/BETA(j) is selected if
*          SELCTG(ALPHAR(j),ALPHAI(j),BETA(j)) is true; i.e. if either
*          one of a complex conjugate pair of eigenvalues is selected,
*          then both complex eigenvalues are selected.
*
*          Note that in the ill-conditioned case, a selected complex
*          eigenvalue may no longer satisfy SELCTG(ALPHAR(j),ALPHAI(j),
*          BETA(j)) = .TRUE. after ordering. INFO is to be set to N+2
*          in this case.
*
*  N       (input) INTEGER
*          The order of the matrices A, B, VSL, and VSR.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the first of the pair of matrices.
*          On exit, A has been overwritten by its generalized Schur
*          form S.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the second of the pair of matrices.
*          On exit, B has been overwritten by its generalized Schur
*          form T.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  LDB >= max(1,N).
*
*  SDIM    (output) INTEGER
*          If SORT = 'N', SDIM = 0.
*          If SORT = 'S', SDIM = number of eigenvalues (after sorting)
*          for which SELCTG is true.  (Complex conjugate pairs for which
*          SELCTG is true for either eigenvalue count as 2.)
*
*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
*          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i,
*          and  BETA(j),j=1,...,N are the diagonals of the complex Schur
*          form (S,T) that would result if the 2-by-2 diagonal blocks of
*          the real Schur form of (A,B) were further reduced to
*          triangular form using 2-by-2 complex unitary transformations.
*          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
*          positive, then the j-th and (j+1)-st eigenvalues are a
*          complex conjugate pair, with ALPHAI(j+1) negative.
*
*          Note: the quotients ALPHAR(j)/BETA(j) and ALPHAI(j)/BETA(j)
*          may easily over- or underflow, and BETA(j) may even be zero.
*          Thus, the user should avoid naively computing the ratio.
*          However, ALPHAR and ALPHAI will be always less than and
*          usually comparable with norm(A) in magnitude, and BETA always
*          less than and usually comparable with norm(B).
*
*  VSL     (output) DOUBLE PRECISION array, dimension (LDVSL,N)
*          If JOBVSL = 'V', VSL will contain the left Schur vectors.
*          Not referenced if JOBVSL = 'N'.
*
*  LDVSL   (input) INTEGER
*          The leading dimension of the matrix VSL. LDVSL >=1, and
*          if JOBVSL = 'V', LDVSL >= N.
*
*  VSR     (output) DOUBLE PRECISION array, dimension (LDVSR,N)
*          If JOBVSR = 'V', VSR will contain the right Schur vectors.
*          Not referenced if JOBVSR = 'N'.
*
*  LDVSR   (input) INTEGER
*          The leading dimension of the matrix VSR. LDVSR >= 1, and
*          if JOBVSR = 'V', LDVSR >= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N = 0, LWORK >= 1, else LWORK >= 8*N+16.
*          For good performance , LWORK must generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  BWORK   (workspace) LOGICAL array, dimension (N)
*          Not referenced if SORT = 'N'.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1,...,N:
*                The QZ iteration failed.  (A,B) are not in Schur
*                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
*                be correct for j=INFO+1,...,N.
*          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
*                =N+2: after reordering, roundoff changed values of
*                      some complex eigenvalues so that leading
*                      eigenvalues in the Generalized Schur form no
*                      longer satisfy SELCTG=.TRUE.  This could also
*                      be caused due to scaling.
*                =N+3: reordering failed in DTGSEN.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gges(const char *jobvsl,const char *jobvsr,const char *sort,sel_fun selctg,
                          i_type n, V *a,i_type lda, V *b,i_type ldb,i_type *sdim, 
                          typename details::complex_type<V>::type *alpha, V *beta, 
                          V *vsl,i_type ldvsl, V *vsr,i_type ldvsr, 
                          V *work,i_type lwork,l_type *bwork,i_type *info);

BLAS_EXPORT void    cgges(const char *jobvsl,const char *jobvsr,const char *sort,
                          sel_fun selctg, i_type n, c_type *a,i_type lda, c_type *b,
                          i_type ldb,i_type *sdim, c_type *alpha, c_type *beta, 
                          c_type *vsl,i_type ldvsl, c_type *vsr,i_type ldvsr, 
                          c_type *work,i_type lwork,s_type *rwork,l_type *bwork,
                          i_type *info);
BLAS_EXPORT void    sgges(const char *jobvsl,const char *jobvsr,const char *sort,
                          sel_fun selctg,i_type n,s_type *a,i_type lda,s_type *b,
                          i_type ldb,i_type *sdim,s_type *alphar,s_type *alphai,
                          s_type *beta,s_type *vsl,i_type ldvsl,s_type *vsr,
                          i_type ldvsr,s_type *work,i_type lwork,l_type *bwork,
                          i_type *info);
BLAS_EXPORT void    dgges(const char *jobvsl,const char *jobvsr,const char *sort,
                          sel_fun delctg,i_type n,d_type *a,i_type lda,d_type *b,
                          i_type ldb,i_type *sdim,d_type *alphar,d_type *alphai,
                          d_type *beta,d_type *vsl,i_type ldvsl,d_type *vsr,
                          i_type ldvsr,d_type *work,i_type lwork,l_type *bwork,
                          i_type *info);
BLAS_EXPORT void    zgges(const char *jobvsl,const char *jobvsr,const char *sort,
                          sel_fun delctg,i_type n, z_type *a,i_type lda, z_type *b,
                          i_type ldb,i_type *sdim, z_type *alpha, z_type *beta, 
                          z_type *vsl,i_type ldvsl, z_type *vsr,i_type ldvsr, 
                          z_type *work,i_type lwork,d_type *rwork,l_type *bwork,
                          i_type *info);

//-----------------------------------------------------------------------
//                          HEEVR / SYEVR
//-----------------------------------------------------------------------
/*
* Purpose:
*
* ZHEEVR computes selected eigenvalues and, optionally, eigenvectors
* of a complex Hermitian matrix A.  Eigenvalues and eigenvectors can
* be selected by specifying either a range of values or a range of
* indices for the desired eigenvalues.
*
* ZHEEVR first reduces the matrix A to tridiagonal form T with a call
* to ZHETRD.  Then, whenever possible, ZHEEVR calls ZSTEMR to compute
* eigenspectrum using Relatively Robust Representations.  ZSTEMR
* computes eigenvalues by the dqds algorithm, while orthogonal
* eigenvectors are computed from various "good" L D L^T representations
* (also known as Relatively Robust Representations). Gram-Schmidt
* orthogonalization is avoided as far as possible. More specifically,
* the various steps of the algorithm are as follows.
*
* For each unreduced block (submatrix) of T,
*    (a) Compute T - sigma I  = L D L^T, so that L and D
*        define all the wanted eigenvalues to high relative accuracy.
*        This means that small relative changes in the entries of D and L
*        cause only small relative changes in the eigenvalues and
*        eigenvectors. The standard (unfactored) representation of the
*        tridiagonal matrix T does not have this property in general.
*    (b) Compute the eigenvalues to suitable accuracy.
*        If the eigenvectors are desired, the algorithm attains full
*        accuracy of the computed eigenvalues only right before
*        the corresponding vectors have to be computed, see steps c) and d).
*    (c) For each cluster of close eigenvalues, select a new
*        shift close to the cluster, find a new factorization, and refine
*        the shifted eigenvalues to suitable accuracy.
*    (d) For each eigenvalue with a large enough relative separation compute
*        the corresponding eigenvector by forming a rank revealing twisted
*        factorization. Go back to (c) for any clusters that remain.
*
* The desired accuracy of the output can be specified by the input
* parameter ABSTOL.
*
* For more details, see DSTEMR's documentation and:
* - Inderjit S. Dhillon and Beresford N. Parlett: "Multiple representations
*   to compute orthogonal eigenvectors of symmetric tridiagonal matrices,"
*   Linear Algebra and its Applications, 387(1), pp. 1-28, August 2004.
* - Inderjit Dhillon and Beresford Parlett: "Orthogonal Eigenvectors and
*   Relative Gaps," SIAM Journal on Matrix Analysis and Applications, Vol. 25,
*   2004.  Also LAPACK Working Note 154.
* - Inderjit Dhillon: "A new O(n^2) algorithm for the symmetric
*   tridiagonal eigenvalue/eigenvector problem",
*   Computer Science Division Technical Report No. UCB/CSD-97-971,
*   UC Berkeley, May 1997.
*
*
* Note 1 : ZHEEVR calls ZSTEMR when the full spectrum is requested
* on machines which conform to the ieee-754 floating point standard.
* ZHEEVR calls DSTEBZ and ZSTEIN on non-ieee machines and
* when partial spectrum requests are made.
*
* Normal execution of ZSTEMR may create NaNs and infinities and
* hence may abort due to a floating point exception in environments
* which do not handle NaNs and infinities in the ieee standard default
* manner.
*
*  Arguments:
*  ==========
*
* [in] JOBZ
*          JOBZ is CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
* [in] RANGE
*          RANGE is CHARACTER*1
*          = 'A': all eigenvalues will be found.
*          = 'V': all eigenvalues in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th eigenvalues will be found.
*          For RANGE = 'V' or 'I' and IU - IL < N - 1, DSTEBZ and
*          ZSTEIN are called
*
* [in] UPLO
*          UPLO is CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
* [in] N
*          N is INTEGER
*          The order of the matrix A.  N >= 0.
*
* [in,out] A
*          A is COMPLEX*16 array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, the lower triangle (if UPLO='L') or the upper
*          triangle (if UPLO='U') of A, including the diagonal, is
*          destroyed.
*
* [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
* [in] VL
*          VL is DOUBLE PRECISION
*
* [in] VU
*          VU is DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues. VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
* [in] IL
*          IL is INTEGER
*
* [in] IU
*          IU is INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.

*
* [in] ABSTOL
*          ABSTOL is DOUBLE PRECISION
*          The absolute error tolerance for the eigenvalues.
*          An approximate eigenvalue is accepted as converged
*          when it is determined to lie in an interval [a,b]
*          of width less than or equal to
*
*                  ABSTOL + EPS *  max( |a|,|b| ) ,
*
*          where EPS is the machine precision.  If ABSTOL is less than
*          or equal to zero, then  EPS*|T|  will be used in its place,
*          where |T| is the 1-norm of the tridiagonal matrix obtained
*          by reducing A to tridiagonal form.
*
*          See "Computing Small Singular Values of Bidiagonal Matrices
*          with Guaranteed High Relative Accuracy," by Demmel and
*          Kahan, LAPACK Working Note #3.
*
*          If high relative accuracy is important, set ABSTOL to
*          DLAMCH( 'Safe minimum' ).  Doing so will guarantee that
*          eigenvalues are computed to high relative accuracy when
*          possible in future releases.  The current code does not
*          make any guarantees about high relative accuracy, but
*          furutre releases will. See J. Barlow and J. Demmel,
*          "Computing Accurate Eigensystems of Scaled Diagonally
*          Dominant Matrices", LAPACK Working Note #7, for a discussion
*          of which matrices define their eigenvalues to high relative
*          accuracy.
*
* [out] M
*          M is INTEGER
*          The total number of eigenvalues found.  0 <= M <= N.
*          If RANGE = 'A', M = N, and if RANGE = 'I', M = IU-IL+1.
*
* [out] W
*          W is DOUBLE PRECISION array, dimension (N)
*          The first M elements contain the selected eigenvalues in
*          ascending order.
*
* [out] Z
*          Z is COMPLEX*16 array, dimension (LDZ, max(1,M))
*          If JOBZ = 'V', then if INFO = 0, the first M columns of Z
*          contain the orthonormal eigenvectors of the matrix A
*          corresponding to the selected eigenvalues, with the i-th
*          column of Z holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*          Note: the user must ensure that at least max(1,M) columns are
*          supplied in the array Z; if RANGE = 'V', the exact value of M
*          is not known in advance and an upper bound must be used.
*
* [in] LDZ
*          LDZ is INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
* [out] ISUPPZ
*          ISUPPZ is INTEGER array, dimension ( 2*max(1,M) )
*          The support of the eigenvectors in Z, i.e., the indices
*          indicating the nonzero elements in Z. The i-th eigenvector
*          is nonzero only in elements ISUPPZ( 2*i-1 ) through
*          ISUPPZ( 2*i ).
*          Implemented only for RANGE = 'A' or 'I' and IU - IL = N - 1
*
* [out] WORK
*          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
* [in] LWORK
*          LWORK is INTEGER
*          The length of the array WORK.  LWORK >= max(1,2*N).
*          For optimal efficiency, LWORK >= (NB+1)*N,
*          where NB is the max of the blocksize for ZHETRD and for
*          ZUNMTR as returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK, RWORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
* [out] RWORK
*          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
*          On exit, if INFO = 0, RWORK(1) returns the optimal
*          (and minimal) LRWORK.
*
* [in] LRWORK
*          LRWORK is INTEGER
*          The length of the array RWORK.  LRWORK >= max(1,24*N).
*
*          If LRWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
* [out] IWORK
*          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal
*          (and minimal) LIWORK.
*
* [in] LIWORK
*          LIWORK is INTEGER
*          The dimension of the array IWORK.  LIWORK >= max(1,10*N).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
* [out] INFO
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  Internal error
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
heevr(const char *jobv,const char *range,const char *uplo, i_type n, V *a,i_type lda,
      typename details::real_type<V>::type vl,typename details::real_type<V>::type vu,
      i_type il, i_type iu, typename details::real_type<V>::type abstol, i_type *m, 
      typename details::real_type<V>::type* w, V* z, i_type ldz, i_type *isuppz, V *work,
      i_type lwork, typename details::real_type<V>::type* rwork, i_type lrwork, 
      i_type *iwork,i_type liwork,i_type *info);

BLAS_EXPORT void    cheevr(const char *jobv,const char *range,const char *uplo,i_type n, 
                           c_type *a,i_type lda, s_type vl, s_type vu, i_type il, i_type iu,
                           s_type abstol, i_type *m, s_type* w, c_type* z, i_type ldz, 
                           i_type *isuppz, c_type *work,i_type lwork, s_type *rwork, 
                           i_type lrwork, i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    ssyevr(const char *jobv,const char *range,const char *uplo,i_type n,
                           s_type *a,i_type lda, s_type vl, s_type vu, i_type il, i_type iu, 
                           s_type abstol,i_type *m, s_type* w, s_type* z, i_type ldz, 
                           i_type *isuppz,s_type *work,i_type lwork, i_type *iwork,
                           i_type liwork,i_type *info);
BLAS_EXPORT void    dsyevr(const char *jobv,const char *range,const char *uplo,i_type n,
                           d_type *a,i_type lda, d_type vl, d_type vu, i_type il, i_type iu, 
                           d_type, i_type *m, d_type* w, d_type* z, i_type ldz, i_type *isuppz,
                          d_type *work,i_type lwork, i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    zheevr(const char *jobv,const char *range,const char *uplo,i_type n, 
                           z_type *a,i_type lda, d_type vl, d_type vu, i_type il, i_type iu, 
                           d_type abstol, i_type *m, d_type* w, z_type* z, i_type ldz, 
                           i_type *isuppz, z_type *work,i_type lwork, d_type *rwork, 
                           i_type lrwork, i_type *iwork,i_type liwork,i_type *info);

//-----------------------------------------------------------------------
//                          SYEV/HEEV
//-----------------------------------------------------------------------

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
heev(const char *jobv,const char *uplo, i_type n, V *a,i_type lda, 
                          typename details::real_type<V>::type* w, 
                          V *work,i_type lwork,i_type *info);

BLAS_EXPORT void    cheev(const char *jobv,const char *uplo,i_type n, c_type *a,
                          i_type lda, s_type* w, c_type *work,i_type lwork, 
                          s_type *rwork, i_type *info);
BLAS_EXPORT void    ssyev(const char *jobv,const char *uplo,i_type n, s_type *a,
                          i_type lda, s_type* w, s_type *work,i_type lwork, 
                          i_type *info);
BLAS_EXPORT void    dsyev(const char *jobv,const char *uplo,i_type n, d_type *a,
                          i_type lda, d_type* w, d_type *work,i_type lwork, 
                          i_type *info);
BLAS_EXPORT void    zheev(const char *jobv,const char *uplo,i_type n, z_type *a,
                          i_type lda, d_type* w, z_type *work,i_type lwork, 
                          d_type *rwork, i_type *info);

//-----------------------------------------------------------------------
//                          SYGVD/HEGVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
* DSYGVD computes all the eigenvalues, and optionally, the eigenvectors
* of a real generalized symmetric-definite eigenproblem, of the form
* A*x=(lambda)*B*x,  A*Bx=(lambda)*x,  or B*A*x=(lambda)*x.  Here A and
* B are assumed to be symmetric and B is also positive definite.
* If eigenvectors are desired, it uses a divide and conquer algorithm.
*
* The divide and conquer algorithm makes very mild assumptions about
* floating point arithmetic. It will work on machines with a guard
* digit in add/subtract, or on those binary machines without guard
* digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
* Cray-2. It could conceivably fail on hexadecimal or decimal machines
* without guard digits, but we know of none.
*
*  Arguments:
*  ==========
*
*              [in] ITYPE
*          ITYPE is INTEGER
*          Specifies the problem type to be solved:
*          = 1:  A*x = (lambda)*B*x
*          = 2:  A*B*x = (lambda)*x
*          = 3:  B*A*x = (lambda)*x
*
*              [in] JOBZ
*          JOBZ is CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*              [in] UPLO
*          UPLO is CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*              [in] N
*          N is INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*              [in,out] A
*          A is DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          matrix Z of eigenvectors.  The eigenvectors are normalized
*          as follows:
*          if ITYPE = 1 or 2, Z**T*B*Z = I;
*          if ITYPE = 3, Z**T*inv(B)*Z = I.
*          If JOBZ = 'N', then on exit the upper triangle (if UPLO='U')
*          or the lower triangle (if UPLO='L') of A, including the
*          diagonal, is destroyed.
*
*              [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*              [in,out] B
*          B is DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the symmetric matrix B.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of B contains the
*          upper triangular part of the matrix B.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of B contains
*          the lower triangular part of the matrix B.
*
*          On exit, if INFO <= N, the part of B containing the matrix is
*          overwritten by the triangular factor U or L from the Cholesky
*          factorization B = U**T*U or B = L*L**T.
*
*              [in] LDB
*          LDB is INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*              [out] W
*          W is DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*              [out] WORK
*          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*              [in] LWORK
*          LWORK is INTEGER
*          The dimension of the array WORK.
*          If N <= 1,               LWORK >= 1.
*          If JOBZ = 'N' and N > 1, LWORK >= 2*N+1.
*          If JOBZ = 'V' and N > 1, LWORK >= 1 + 6*N + 2*N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK and IWORK
*          arrays, returns these values as the first entries of the WORK
*          and IWORK arrays, and no error message related to LWORK or
*          LIWORK is issued by XERBLA.
*
*              [out] IWORK
*          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*              [in] LIWORK
*          LIWORK is INTEGER
*          The dimension of the array IWORK.
*          If N <= 1,                LIWORK >= 1.
*          If JOBZ  = 'N' and N > 1, LIWORK >= 1.
*          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK and IWORK arrays, and no error message related to
*          LWORK or LIWORK is issued by XERBLA.
*
*              [out] INFO
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  DPOTRF or DSYEVD returned an error code:
*             <= N:  if INFO = i and JOBZ = 'N', then the algorithm
*                    failed to converge; i off-diagonal elements of an
*                    intermediate tridiagonal form did not converge to
*                    zero;
*                    if INFO = i and JOBZ = 'V', then the algorithm
*                    failed to compute an eigenvalue while working on
*                    the submatrix lying in rows and columns INFO/(N+1)
*                    through mod(INFO,N+1);
*             > N:   if INFO = N + i, for 1 <= i <= N, then the leading
*                    minor of order i of B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
*
* Further Details:
*  =====================
*
*  Modified so that no backsubstitution is performed if DSYEVD fails to
*  converge (NEIG in old code could be greater than N causing out of
*  bounds reference to A - reported by Ralf Meyer).  Also corrected the
*  description of INFO and the test on ITYPE. Sven, 16 Feb 05.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
hegvd(i_type itype, const char* jobz,const char *uplo, i_type n, V *a,i_type lda,
      V *b,i_type ldb, typename details::real_type<V>::type *w, V *work,i_type lwork,
      typename details::real_type<V>::type* rwork, i_type lrwork, i_type *iwork,
      i_type liwork,i_type *info);

BLAS_EXPORT void    chegvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
                          c_type *a,i_type lda, c_type *b,i_type ldb, s_type* w,
                          c_type *work,i_type lwork, s_type* rwork, i_type lrwork,
                          i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    ssygvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
                          s_type *a,i_type lda, s_type *b,i_type ldb, s_type w,
                          s_type *work,i_type lwork, 
                          i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    dsygvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
                          d_type *a,i_type lda, d_type *b,i_type ldb, d_type * w,
                          d_type *work,i_type lwork, 
                          i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    zhegvd(i_type itype, const char* jobz,const char *uplo, i_type n, 
                          z_type *a,i_type lda, z_type *b,i_type ldb, d_type * w,
                          z_type *work,i_type lwork, d_type* rwork, i_type lrwork,
                          i_type *iwork,i_type liwork,i_type *info);

//-----------------------------------------------------------------------
//                          GESVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGESVD computes the singular value decomposition (SVD) of a real
*  M-by-N matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written
*
*       A = U * SIGMA * transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
*
*  Note that the routine returns V**T, not V.
*
*  Arguments
*  =========
*
*  JOBU    (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix U:
*          = 'A':  all M columns of U are returned in array U:
*          = 'S':  the first min(m,n) columns of U (the left singular
*                  vectors) are returned in the array U;
*          = 'O':  the first min(m,n) columns of U (the left singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no columns of U (no left singular vectors) are
*                  computed.
*
*  JOBVT   (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix
*          V**T:
*          = 'A':  all N rows of V**T are returned in the array VT;
*          = 'S':  the first min(m,n) rows of V**T (the right singular
*                  vectors) are returned in the array VT;
*          = 'O':  the first min(m,n) rows of V**T (the right singular
*                  vectors) are overwritten on the array A;
*          = 'N':  no rows of V**T (no right singular vectors) are
*                  computed.
*
*          JOBVT and JOBU cannot both be 'O'.
*
*  M       (input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the input matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*          if JOBU = 'O',  A is overwritten with the first min(m,n)
*                          columns of U (the left singular vectors,
*                          stored columnwise);
*          if JOBVT = 'O', A is overwritten with the first min(m,n)
*                          rows of V**T (the right singular vectors,
*                          stored rowwise);
*          if JOBU .ne. 'O' and JOBVT .ne. 'O', the contents of A
*                          are destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The singular values of A, sorted so that S(i) >= S(i+1).
*
*  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
*          (LDU,M) if JOBU = 'A' or (LDU,min(M,N)) if JOBU = 'S'.
*          If JOBU = 'A', U contains the M-by-M orthogonal matrix U;
*          if JOBU = 'S', U contains the first min(m,n) columns of U
*          (the left singular vectors, stored columnwise);
*          if JOBU = 'N' or 'O', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= 1; if
*          JOBU = 'S' or 'A', LDU >= M.
*
*  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
*          If JOBVT = 'A', VT contains the N-by-N orthogonal matrix
*          V**T;
*          if JOBVT = 'S', VT contains the first min(m,n) rows of
*          V**T (the right singular vectors, stored rowwise);
*          if JOBVT = 'N' or 'O', VT is not referenced.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.  LDVT >= 1; if
*          JOBVT = 'A', LDVT >= N; if JOBVT = 'S', LDVT >= min(M,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
*          if INFO > 0, WORK(2:MIN(M,N)) contains the unconverged
*          superdiagonal elements of an upper bidiagonal matrix B
*          whose diagonal is in S (not necessarily sorted). B
*          satisfies A = U * B * VT, so it has the same singular values
*          as A, and singular vectors related by U and VT.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          LWORK >= MAX(1,3*MIN(M,N)+MAX(M,N),5*MIN(M,N)).
*          For good performance, LWORK should generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if DBDSQR did not converge, INFO specifies how many
*                superdiagonals of an intermediate bidiagonal form B
*                did not converge to zero. See the description of WORK
*                above for details.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gesvd(const char *jobu,const char *jobvt,i_type m,i_type n,V *a,i_type lda,
      typename details::real_type<V>::type *s,V *u,i_type ldu,V *vt,i_type ldvt,
      V *work,i_type lwork,i_type *info);

BLAS_EXPORT void    cgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,
                           c_type *a,i_type lda, s_type *s,c_type *u,i_type ldu,
                           c_type *vt,i_type ldvt, c_type *work,i_type lwork,
                           s_type *rwork, i_type *info);
BLAS_EXPORT void    sgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,
                           s_type *a,i_type lda, s_type *s,s_type *u,i_type ldu,
                           s_type *vt,i_type ldvt, s_type *work,i_type lwork,
                           i_type *info);
BLAS_EXPORT void    dgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,
                           d_type *a,i_type lda, d_type *s,d_type *u,i_type ldu,
                           d_type *vt,i_type ldvt, d_type *work,i_type lwork,
                           i_type *info);
BLAS_EXPORT void    zgesvd(const char *jobu,const char *jobvt,i_type m,i_type n,
                           z_type *a,i_type lda, d_type *s, z_type *u,i_type ldu,
                           z_type *vt,i_type ldvt, z_type *work,i_type lwork,
                           d_type *rwork,i_type *info);

//-----------------------------------------------------------------------
//                          DGESDD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGESDD computes the singular value decomposition (SVD) of a real
*  M-by-N matrix A, optionally computing the left and right singular
*  vectors.  If singular vectors are desired, it uses a
*  divide-and-conquer algorithm.
*
*  The SVD is written
*
*       A = U * SIGMA * transpose(V)
*
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
*
*  Note that the routine returns VT = V**T, not V.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          Specifies options for computing all or part of the matrix U:
*          = 'A':  all M columns of U and all N rows of V**T are
*                  returned in the arrays U and VT;
*          = 'S':  the first min(M,N) columns of U and the first
*                  min(M,N) rows of V**T are returned in the arrays U
*                  and VT;
*          = 'O':  If M >= N, the first N columns of U are overwritten
*                  on the array A and all rows of V**T are returned in
*                  the array VT;
*                  otherwise, all columns of U are returned in the
*                  array U and the first M rows of V**T are overwritten
*                  in the array A;
*          = 'N':  no columns of U or rows of V**T are computed.
*
*  M       (input) INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the input matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit,
*          if JOBZ = 'O',  A is overwritten with the first N columns
*                          of U (the left singular vectors, stored
*                          columnwise) if M >= N;
*                          A is overwritten with the first M rows
*                          of V**T (the right singular vectors, stored
*                          rowwise) otherwise.
*          if JOBZ .ne. 'O', the contents of A are destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  S       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The singular values of A, sorted so that S(i) >= S(i+1).
*
*  U       (output) DOUBLE PRECISION array, dimension (LDU,UCOL)
*          UCOL = M if JOBZ = 'A' or JOBZ = 'O' and M < N;
*          UCOL = min(M,N) if JOBZ = 'S'.
*          If JOBZ = 'A' or JOBZ = 'O' and M < N, U contains the M-by-M
*          orthogonal matrix U;
*          if JOBZ = 'S', U contains the first min(M,N) columns of U
*          (the left singular vectors, stored columnwise);
*          if JOBZ = 'O' and M >= N, or JOBZ = 'N', U is not referenced.
*
*  LDU     (input) INTEGER
*          The leading dimension of the array U.  LDU >= 1; if
*          JOBZ = 'S' or 'A' or JOBZ = 'O' and M < N, LDU >= M.
*
*  VT      (output) DOUBLE PRECISION array, dimension (LDVT,N)
*          If JOBZ = 'A' or JOBZ = 'O' and M >= N, VT contains the
*          N-by-N orthogonal matrix V**T;
*          if JOBZ = 'S', VT contains the first min(M,N) rows of
*          V**T (the right singular vectors, stored rowwise);
*          if JOBZ = 'O' and M < N, or JOBZ = 'N', VT is not referenced.
*
*  LDVT    (input) INTEGER
*          The leading dimension of the array VT.  LDVT >= 1; if
*          JOBZ = 'A' or JOBZ = 'O' and M >= N, LDVT >= N;
*          if JOBZ = 'S', LDVT >= min(M,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= 1.
*          If JOBZ = 'N',
*            LWORK >= 3*min(M,N) + max(max(M,N),7*min(M,N)).
*          If JOBZ = 'O',
*            LWORK >= 3*min(M,N)*min(M,N) + 
*                     max(max(M,N),5*min(M,N)*min(M,N)+4*min(M,N)).
*          If JOBZ = 'S' or 'A'
*            LWORK >= 3*min(M,N)*min(M,N) +
*                     max(max(M,N),4*min(M,N)*min(M,N)+4*min(M,N)).
*          For good performance, LWORK should generally be larger.
*          If LWORK = -1 but other input arguments are legal, WORK(1)
*          returns the optimal LWORK.
*
*  IWORK   (workspace) INTEGER array, dimension (8*min(M,N))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  DBDSDC did not converge, updating process failed.
*
*  Further Details
*  ===============
*
*  Based on contributions by
*     Ming Gu and Huan Ren, Computer Science Division, University of
*     California at Berkeley, USA
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gesdd(const char *jobu, i_type m,i_type n,V *a,i_type lda,
      typename details::real_type<V>::type *s,V *u,i_type ldu,V *vt,i_type ldvt,
      V *work,i_type lwork, i_type *iwork, i_type *info);

BLAS_EXPORT void    cgesdd(const char *jobu,i_type m,i_type n,c_type *a,i_type lda,
                           s_type *s,c_type *u,i_type ldu,c_type *vt,i_type ldvt,
                           c_type *work,i_type lwork,s_type *rwork, i_type *iwork, 
                           i_type *info);
BLAS_EXPORT void    sgesdd(const char *jobu,i_type m,i_type n,s_type *a,i_type lda,
                           s_type *s,s_type *u,i_type ldu,s_type *vt,i_type ldvt,
                           s_type *work,i_type lwork, i_type *iwork,i_type *info);
BLAS_EXPORT void    dgesdd(const char *jobu,i_type m,i_type n,d_type *a,i_type lda,
                           d_type *s,d_type *u,i_type ldu,d_type *vt,i_type ldvt,
                           d_type *work,i_type lwork, i_type *iwork,i_type *info);
BLAS_EXPORT void    zgesdd(const char *jobu,i_type m,i_type n,z_type *a,i_type lda,
                           d_type *s, z_type *u,i_type ldu,z_type *vt,i_type ldvt,
                           z_type *work,i_type lwork,d_type *rwork, i_type *iwork,
                           i_type *info);

//-----------------------------------------------------------------------
//                          TGSEN
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTGSEN reorders the generalized real Schur decomposition of a real
*  matrix pair (A, B) (in terms of an orthonormal equivalence trans-
*  formation Q' * (A, B) * Z), so that a selected cluster of eigenvalues
*  appears in the leading diagonal blocks of the upper quasi-triangular
*  matrix A and the upper triangular B. The leading columns of Q and
*  Z form orthonormal bases of the corresponding left and right eigen-
*  spaces (deflating subspaces). (A, B) must be in generalized real
*  Schur canonical form (as returned by DGGES), i.e. A is block upper
*  triangular with 1-by-1 and 2-by-2 diagonal blocks. B is upper
*  triangular.
*
*  DTGSEN also computes the generalized eigenvalues
*
*              w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
*
*  of the reordered matrix pair (A, B).
*
*  Optionally, DTGSEN computes the estimates of reciprocal condition
*  numbers for eigenvalues and eigenspaces. These are Difu[(A11,B11),
*  (A22,B22)] and Difl[(A11,B11), (A22,B22)], i.e. the separation(s)
*  between the matrix pairs (A11, B11) and (A22,B22) that correspond to
*  the selected cluster and the eigenvalues outside the cluster, resp.,
*  and norms of "projections" onto left and right eigenspaces w.r.t.
*  the selected cluster in the (1,1)-block.
*
*  Arguments
*  =========
*
*  IJOB    (input) INTEGER
*          Specifies whether condition numbers are required for the
*          cluster of eigenvalues (PL and PR) or the deflating subspaces
*          (Difu and Difl):
*           =0: Only reorder w.r.t. SELECT. No extras.
*           =1: Reciprocal of norms of "projections" onto left and right
*               eigenspaces w.r.t. the selected cluster (PL and PR).
*           =2: Upper bounds on Difu and Difl. F-norm-based estimate
*               (DIF(1:2)).
*           =3: Estimate of Difu and Difl. 1-norm-based estimate
*               (DIF(1:2)).
*               About 5 times as expensive as IJOB = 2.
*           =4: Compute PL, PR and DIF (i.e. 0, 1 and 2 above): Economic
*               version to get it all.
*           =5: Compute PL, PR and DIF (i.e. 0, 1 and 3 above)
*
*  WANTQ   (input) LOGICAL
*          .TRUE. : update the left transformation matrix Q;
*          .FALSE.: do not update Q.
*
*  WANTZ   (input) LOGICAL
*          .TRUE. : update the right transformation matrix Z;
*          .FALSE.: do not update Z.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          SELECT specifies the eigenvalues in the selected cluster.
*          To select a real eigenvalue w(j), SELECT(j) must be set to
*          .TRUE.. To select a complex conjugate pair of eigenvalues
*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*          either SELECT(j) or SELECT(j+1) or both must be set to
*          .TRUE.; a complex conjugate pair of eigenvalues must be
*          either both included in the cluster or both excluded.
*
*  N       (input) INTEGER
*          The order of the matrices A and B. N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension(LDA,N)
*          On entry, the upper quasi-triangular matrix A, with (A, B) in
*          generalized real Schur canonical form.
*          On exit, A is overwritten by the reordered matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension(LDB,N)
*          On entry, the upper triangular matrix B, with (A, B) in
*          generalized real Schur canonical form.
*          On exit, B is overwritten by the reordered matrix B.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          On exit, (ALPHAR(j) + ALPHAI(j)*i)/BETA(j), j=1,...,N, will
*          be the generalized eigenvalues.  ALPHAR(j) + ALPHAI(j)*i
*          and BETA(j),j=1,...,N  are the diagonals of the complex Schur
*          form (S,T) that would result if the 2-by-2 diagonal blocks of
*          the real generalized Schur form of (A,B) were further reduced
*          to triangular form using complex unitary transformations.
*          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
*          positive, then the j-th and (j+1)-st eigenvalues are a
*          complex conjugate pair, with ALPHAI(j+1) negative.
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if WANTQ = .TRUE., Q is an N-by-N matrix.
*          On exit, Q has been postmultiplied by the left orthogonal
*          transformation matrix which reorder (A, B); The leading M
*          columns of Q form orthonormal bases for the specified pair of
*          left eigenspaces (deflating subspaces).
*          If WANTQ = .FALSE., Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= 1;
*          and if WANTQ = .TRUE., LDQ >= N.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          On entry, if WANTZ = .TRUE., Z is an N-by-N matrix.
*          On exit, Z has been postmultiplied by the left orthogonal
*          transformation matrix which reorder (A, B); The leading M
*          columns of Z form orthonormal bases for the specified pair of
*          left eigenspaces (deflating subspaces).
*          If WANTZ = .FALSE., Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z. LDZ >= 1;
*          If WANTZ = .TRUE., LDZ >= N.
*
*  M       (output) INTEGER
*          The dimension of the specified pair of left and right eigen-
*          spaces (deflating subspaces). 0 <= M <= N.
*
*  PL      (output) DOUBLE PRECISION
*  PR      (output) DOUBLE PRECISION
*          If IJOB = 1, 4 or 5, PL, PR are lower bounds on the
*          reciprocal of the norm of "projections" onto left and right
*          eigenspaces with respect to the selected cluster.
*          0 < PL, PR <= 1.
*          If M = 0 or M = N, PL = PR  = 1.
*          If IJOB = 0, 2 or 3, PL and PR are not referenced.
*
*  DIF     (output) DOUBLE PRECISION array, dimension (2).
*          If IJOB >= 2, DIF(1:2) store the estimates of Difu and Difl.
*          If IJOB = 2 or 4, DIF(1:2) are F-norm-based upper bounds on
*          Difu and Difl. If IJOB = 3 or 5, DIF(1:2) are 1-norm-based
*          estimates of Difu and Difl.
*          If M = 0 or N, DIF(1:2) = F-norm([A, B]).
*          If IJOB = 0 or 1, DIF is not referenced.
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*          dimension (MAX(1,LWORK)) 
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >=  4*N+16.
*          If IJOB = 1, 2 or 4, LWORK >= MAX(4*N+16, 2*M*(N-M)).
*          If IJOB = 3 or 5, LWORK >= MAX(4*N+16, 4*M*(N-M)).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          IF IJOB = 0, IWORK is not referenced.  Otherwise,
*          on exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK. LIWORK >= 1.
*          If IJOB = 1, 2 or 4, LIWORK >=  N+6.
*          If IJOB = 3 or 5, LIWORK >= MAX(2*M*(N-M), N+6).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*            =0: Successful exit.
*            <0: If INFO = -i, the i-th argument had an illegal value.
*            =1: Reordering of (A, B) failed because the transformed
*                matrix pair (A, B) would be too far from generalized
*                Schur form; the problem is very ill-conditioned.
*                (A, B) may have been partially reordered.
*                If requested, 0 is returned in DIF(*), PL and PR.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
tgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,i_type n,
      V *a,i_type lda,V *b,i_type ldb, typename details::complex_type<V>::type *alpha, 
      V *beta, V *q,i_type ldq, V *z,i_type ldz,i_type *m, 
      typename details::real_type<V>::type *pl, typename details::real_type<V>::type *pr,
      typename details::real_type<V>::type *dif, V *work,i_type lwork, i_type *iwork,
      i_type liwork,i_type *info);

BLAS_EXPORT void    ctgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,
                           i_type n, c_type *a,i_type lda,c_type *b,i_type ldb, 
                           c_type *alpha, c_type *beta, c_type *q,i_type ldq, c_type *z,
                           i_type ldz,i_type *m, s_type *pl,s_type *pr,s_type *dif, 
                           c_type *work,i_type lwork, i_type *iwork,i_type liwork,
                           i_type *info);
BLAS_EXPORT void    stgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,
                           i_type n, s_type *a,i_type lda,s_type *b,i_type ldb,
                           s_type *alphar,s_type *alphai, s_type *beta,s_type *q,
                           i_type ldq,s_type *z,i_type ldz,i_type *m, s_type *pl,
                           s_type *pr,s_type *dif,s_type *work,i_type lwork, 
                           i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    dtgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,
                           i_type n, d_type *a,i_type lda,d_type *b,i_type ldb,
                           d_type *alphar,d_type *alphai, d_type *beta,d_type *q,
                           i_type ldq,d_type *z,i_type ldz,i_type *m, d_type *pl,
                           d_type *pr,d_type *dif,d_type *work,i_type lwork, 
                           i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    ztgsen(i_type ijob,i_type wantq,i_type wantz,const i_type *select,
                           i_type n, z_type *a,i_type lda, z_type *b,i_type ldb, 
                           z_type *alpha, z_type *beta, z_type *q,i_type ldq, z_type *z,
                           i_type ldz,i_type *m, d_type *pl,d_type *pr,d_type *dif, 
                           z_type *work,i_type lwork, i_type *iwork,i_type liwork,
                           i_type *info);

//-----------------------------------------------------------------------
//                          TRSEN
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
* DTRSEN reorders the real Schur factorization of a real matrix
* A = Q*T*Q**T, so that a selected cluster of eigenvalues appears in
* the leading diagonal blocks of the upper quasi-triangular matrix T,
* and the leading columns of Q form an orthonormal basis of the
* corresponding right invariant subspace.
*
* Optionally the routine computes the reciprocal condition numbers of
* the cluster of eigenvalues and/or the invariant subspace.
*
* T must be in Schur canonical form (as returned by DHSEQR), that is,
* block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
* 2-by-2 diagonal block has its diagonal elements equal and its
* off-diagonal elements of opposite sign.
*
*  Arguments:
*  ==========
*
* [in] JOB
*          JOB is CHARACTER*1
*          Specifies whether condition numbers are required for the
*          cluster of eigenvalues (S) or the invariant subspace (SEP):
*          = 'N': none;
*          = 'E': for eigenvalues only (S);
*          = 'V': for invariant subspace only (SEP);
*          = 'B': for both eigenvalues and invariant subspace (S and
*                 SEP).
*
* [in] COMPQ
*          COMPQ is CHARACTER*1
*          = 'V': update the matrix Q of Schur vectors;
*          = 'N': do not update Q.
*
* [in] SELECT
*          SELECT is LOGICAL array, dimension (N)
*          SELECT specifies the eigenvalues in the selected cluster. To
*          select a real eigenvalue w(j), SELECT(j) must be set to
*          .TRUE.. To select a complex conjugate pair of eigenvalues
*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*          either SELECT(j) or SELECT(j+1) or both must be set to
*          .TRUE.; a complex conjugate pair of eigenvalues must be
*          either both included in the cluster or both excluded.
*
* [in] N
*          N is INTEGER
*          The order of the matrix T. N >= 0.
*
* [in,out] T
*          T is DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          canonical form.
*          On exit, T is overwritten by the reordered matrix T, again in
*          Schur canonical form, with the selected eigenvalues in the
*          leading diagonal blocks.
*
* [in] LDT
*          LDT is INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
* [in,out] Q
*          Q is DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*          orthogonal transformation matrix which reorders T; the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If COMPQ = 'N', Q is not referenced.
*
* [in] LDQ
*          LDQ is INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1; and if COMPQ = 'V', LDQ >= N.
*
* [out] WR
*          WR is DOUBLE PRECISION array, dimension (N)
* [out] WI
*          WI is DOUBLE PRECISION array, dimension (N)
*
*          The real and imaginary parts, respectively, of the reordered
*          eigenvalues of T. The eigenvalues are stored in the same
*          order as on the diagonal of T, with WR(i) = T(i,i) and, if
*          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
*          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
*          sufficiently ill-conditioned, then its value may differ
*          significantly from its value before reordering.
*
* [out] M
*          M is INTEGER
*          The dimension of the specified invariant subspace.
*          0 < = M <= N.
*
* [out] S
*          S is DOUBLE PRECISION
*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
*          condition number for the selected cluster of eigenvalues.
*          S cannot underestimate the true reciprocal condition number
*          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
*          If JOB = 'N' or 'V', S is not referenced.
*
* [out] SEP
*          SEP is DOUBLE PRECISION
*          If JOB = 'V' or 'B', SEP is the estimated reciprocal
*          condition number of the specified invariant subspace. If
*          M = 0 or N, SEP = norm(T).
*          If JOB = 'N' or 'E', SEP is not referenced.
*
* [out] WORK
*          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
* [in] LWORK
*          LWORK is INTEGER
*          The dimension of the array WORK.
*          If JOB = 'N', LWORK >= max(1,N);
*          if JOB = 'E', LWORK >= max(1,M*(N-M));
*          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
* [out] IWORK
*          IWORK is INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
* [in] LIWORK
*          LIWORK is INTEGER
*          The dimension of the array IWORK.
*          If JOB = 'N' or 'E', LIWORK >= 1;
*          if JOB = 'V' or 'B', LIWORK >= max(1,M*(N-M)).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
* [out] INFO
*          INFO is INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: reordering of T failed because some eigenvalues are too
*               close to separate (the problem is very ill-conditioned);
*               T may have been partially reordered, and WR and WI
*               contain the eigenvalues in the same order as in T; S and
*               SEP (if requested) are set to zero.
*
* Further Details:
*  =====================
*
*
*  DTRSEN first collects the selected eigenvalues by computing an
*  orthogonal transformation Z to move them to the top left corner of T.
*  In other words, the selected eigenvalues are the eigenvalues of T11
*  in:
*
*          Z**T * T * Z = ( T11 T12 ) n1
*                         (  0  T22 ) n2
*                            n1  n2
*
*  where N = n1+n2 and Z**T means the transpose of Z. The first n1 columns
*  of Z span the specified invariant subspace of T.
*
*  If T has been obtained from the real Schur factorization of a matrix
*  A = Q*T*Q**T, then the reordered real Schur factorization of A is given
*  by A = (Q*Z)*(Z**T*T*Z)*(Q*Z)**T, and the first n1 columns of Q*Z span
*  the corresponding invariant subspace of A.
*
*  The reciprocal condition number of the average of the eigenvalues of
*  T11 may be returned in S. S lies between 0 (very badly conditioned)
*  and 1 (very well conditioned). It is computed as follows. First we
*  compute R so that
*
*                         P = ( I  R ) n1
*                             ( 0  0 ) n2
*                               n1 n2
*
*  is the projector on the invariant subspace associated with T11.
*  R is the solution of the Sylvester equation:
*
*                        T11*R - R*T22 = T12.
*
*  Let F-norm(M) denote the Frobenius-norm of M and 2-norm(M) denote
*  the two-norm of M. Then S is computed as the lower bound
*
*                      (1 + F-norm(R)**2)**(-1/2)
*
*  on the reciprocal of 2-norm(P), the true reciprocal condition number.
*  S cannot underestimate 1 / 2-norm(P) by more than a factor of
*  sqrt(N).
*
*  An approximate error bound for the computed average of the
*  eigenvalues of T11 is
*
*                         EPS * norm(T) / S
*
*  where EPS is the machine precision.
*
*  The reciprocal condition number of the right invariant subspace
*  spanned by the first n1 columns of Z (or of Q*Z) is returned in SEP.
*  SEP is defined as the separation of T11 and T22:
*
*                     sep( T11, T22 ) = sigma-min( C )
*
*  where sigma-min(C) is the smallest singular value of the
*  n1*n2-by-n1*n2 matrix
*
*     C  = kprod( I(n2), T11 ) - kprod( transpose(T22), I(n1) )
*
*  I(m) is an m by m identity matrix, and kprod denotes the Kronecker
*  product. We estimate sigma-min(C) by the reciprocal of an estimate of
*  the 1-norm of inverse(C). The true reciprocal 1-norm of inverse(C)
*  cannot differ from sigma-min(C) by more than a factor of sqrt(n1*n2).
*
*  When SEP is small, small changes in T can cause large changes in
*  the invariant subspace. An approximate bound on the maximum angular
*  error in the computed right invariant subspace is
*
*                      EPS * norm(T) / SEP
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trsen(const char * job,const char * compq,const i_type *select,i_type n,
      V *t,i_type ldt,V *q,i_type ldq, typename details::complex_type<V>::type *w, 
      i_type *m, typename details::real_type<V>::type *s, 
      typename details::real_type<V>::type *sep,V *work,i_type lwork, i_type *iwork,
      i_type liwork,i_type *info);

BLAS_EXPORT void    ctrsen(const char * job,const char * compq,const i_type *select,
                           i_type n, c_type *t,i_type ldt,c_type *q,i_type ldq, c_type *w,
                           i_type *m, s_type *s,s_type *sep,c_type *work,i_type lwork,
                           i_type *info);
BLAS_EXPORT void    strsen(const char * job,const char * compq,const i_type *select,
                           i_type n, s_type *t,i_type ldt,s_type *q,i_type ldq,s_type *wr,
                           s_type *wi, i_type *m, s_type *s,s_type *sep,s_type *work,
                           i_type lwork, i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    dtrsen(const char * job,const char * compq,const i_type *select,
                           i_type n, d_type *t,i_type ldt,d_type *q,i_type ldq,d_type *wr,
                           d_type *wi, i_type *m, d_type *s,d_type *sep,d_type *work,
                           i_type lwork, i_type *iwork,i_type liwork,i_type *info);
BLAS_EXPORT void    ztrsen(const char * job,const char * compq,const i_type *select,
                           i_type n, z_type *t,i_type ldt, z_type *q,i_type ldq, z_type *w, 
                           i_type *m, d_type *s,d_type *sep, z_type *work,i_type lwork,
                           i_type *info);

//-----------------------------------------------------------------------
//                          TRSYL
//-----------------------------------------------------------------------
/*
* DTRSYL solves the real Sylvester matrix equation:
*
*    op(A)*X + X*op(B) = scale*C or
*    op(A)*X - X*op(B) = scale*C,
*
* where op(A) = A or A**T, and  A and B are both upper quasi-
* triangular. A is M-by-M and B is N-by-N; the right hand side C and
* the solution X are M-by-N; and scale is an output scale factor, set
* <= 1 to avoid overflow in X.
*
* A and B must be in Schur canonical form (as returned by DHSEQR), that
* is, block upper triangular with 1-by-1 and 2-by-2 diagonal blocks;
* each 2-by-2 diagonal block has its diagonal elements equal and its
* off-diagonal elements of opposite sign.
*
*  Arguments:
*  ==========
*
* [in] TRANA
*          TRANA is CHARACTER*1
*          Specifies the option op(A):
*          = 'N': op(A) = A    (No transpose)
*          = 'T': op(A) = A**T (Transpose)
*          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
*
* [in] TRANB
*          TRANB is CHARACTER*1
*          Specifies the option op(B):
*          = 'N': op(B) = B    (No transpose)
*          = 'T': op(B) = B**T (Transpose)
*          = 'C': op(B) = B**H (Conjugate transpose = Transpose)
*
* [in] ISGN
*          ISGN is INTEGER
*          Specifies the sign in the equation:
*          = +1: solve op(A)*X + X*op(B) = scale*C
*          = -1: solve op(A)*X - X*op(B) = scale*C
*
* [in] M
*          M is INTEGER
*          The order of the matrix A, and the number of rows in the
*          matrices X and C. M >= 0.
*
* [in] N
*          N is INTEGER
*          The order of the matrix B, and the number of columns in the
*          matrices X and C. N >= 0.
*
* [in] A
*          A is DOUBLE PRECISION array, dimension (LDA,M)
*          The upper quasi-triangular matrix A, in Schur canonical form.
*
* [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
* [in] B
*          B is DOUBLE PRECISION array, dimension (LDB,N)
*          The upper quasi-triangular matrix B, in Schur canonical form.
*
* [in] LDB
*          LDB is INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
* [in,out] C
*          C is DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N right hand side matrix C.
*          On exit, C is overwritten by the solution matrix X.
*
* [in] LDC
*          LDC is INTEGER
*          The leading dimension of the array C. LDC >= max(1,M)
*
* [out] SCALE
*          SCALE is DOUBLE PRECISION
*          The scale factor, scale, set <= 1 to avoid overflow in X.
*
* [out] INFO
*          INFO is INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: A and B have common or very close eigenvalues; perturbed
*               values were used to solve the equation (but the matrices
*               A and B are unchanged).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trsyl(const char * trana, const char * tranb,const i_type isgn,i_type m,i_type n,
      const V *a,i_type lda, const V *b,i_type ldb, V *c,i_type ldc,
      typename details::real_type<V>::type *scale, i_type *info);

BLAS_EXPORT void    ctrsyl(const char * trana,const char * tranb,
                           const i_type isgn,i_type m,i_type n, const c_type *a,
                           i_type lda, const c_type *b,i_type ldb, c_type *c,
                           i_type ldc, s_type *scale, i_type *info);
BLAS_EXPORT void    strsyl(const char * trana,const char * tranb,
                           const i_type isgn,i_type m,i_type n, const s_type *a,
                           i_type lda, const s_type *b,i_type ldb, s_type *c,
                           i_type ldc, s_type *scale, i_type *info);
BLAS_EXPORT void    dtrsyl(const char * trana,const char * tranb,const i_type isgn,
                           i_type m,i_type n, const d_type *a,i_type lda,
                           const d_type *b,i_type ldb, d_type *c,i_type ldc,
                           d_type *scale, i_type *info);
BLAS_EXPORT void    ztrsyl(const char * trana,const char * tranb,const i_type isgn,
                           i_type m,i_type n,const z_type *a,i_type lda, 
                           const z_type *b,i_type ldb, z_type *c,i_type ldc,
                           d_type *scale, i_type *info);

//-----------------------------------------------------------------------
//                          TGSYL
//-----------------------------------------------------------------------
/*
*
* DTGSYL solves the generalized Sylvester equation:
*
*             A * R - L * B = scale * C                 (1)
*             D * R - L * E = scale * F
*
* where R and L are unknown m-by-n matrices, (A, D), (B, E) and
* (C, F) are given matrix pairs of size m-by-m, n-by-n and m-by-n,
* respectively, with real entries. (A, D) and (B, E) must be in
* generalized (real) Schur canonical form, i.e. A, B are upper quasi
* triangular and D, E are upper triangular.
*
* The solution (R, L) overwrites (C, F). 0 <= SCALE <= 1 is an output
* scaling factor chosen to avoid overflow.
*
* In matrix notation (1) is equivalent to solve  Zx = scale b, where
* Z is defined as
*
*            Z = [ kron(In, A)  -kron(B**T, Im) ]         (2)
*                [ kron(In, D)  -kron(E**T, Im) ].
*
* Here Ik is the identity matrix of size k and X**T is the transpose of
* X. kron(X, Y) is the Kronecker product between the matrices X and Y.
*
* If TRANS = 'T', DTGSYL solves the transposed system Z**T*y = scale*b,
* which is equivalent to solve for R and L in
*
*             A**T * R + D**T * L = scale * C           (3)
*             R * B**T + L * E**T = scale * -F
*
* This case (TRANS = 'T') is used to compute an one-norm-based estimate
* of Dif[(A,D), (B,E)], the separation between the matrix pairs (A,D)
* and (B,E), using DLACON.
*
* If IJOB >= 1, DTGSYL computes a Frobenius norm-based estimate
* of Dif[(A,D),(B,E)]. That is, the reciprocal of a lower bound on the
* reciprocal of the smallest singular value of Z. See [1-2] for more
* information.
*
* This is a level 3 BLAS algorithm.
*
*  Arguments:
*  ==========
*
* [in] TRANS
*          TRANS is CHARACTER*1
*          = 'N', solve the generalized Sylvester equation (1).
*          = 'T', solve the 'transposed' system (3).
*
* [in] IJOB
*          IJOB is INTEGER
*          Specifies what kind of functionality to be performed.
*           =0: solve (1) only.
*           =1: The functionality of 0 and 3.
*           =2: The functionality of 0 and 4.
*           =3: Only an estimate of Dif[(A,D), (B,E)] is computed.
*               (look ahead strategy IJOB  = 1 is used).
*           =4: Only an estimate of Dif[(A,D), (B,E)] is computed.
*               ( DGECON on sub-systems is used ).
*          Not referenced if TRANS = 'T'.
*
* [in] M
*          M is INTEGER
*          The order of the matrices A and D, and the row dimension of
*          the matrices C, F, R and L.
*
* [in] N
*          N is INTEGER
*          The order of the matrices B and E, and the column dimension
*          of the matrices C, F, R and L.
*
* [in] A
*          A is DOUBLE PRECISION array, dimension (LDA, M)
*          The upper quasi triangular matrix A.
*
* [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A. LDA >= max(1, M).
*
* [in] B
*          B is DOUBLE PRECISION array, dimension (LDB, N)
*          The upper quasi triangular matrix B.
*
* [in] LDB
*          LDB is INTEGER
*          The leading dimension of the array B. LDB >= max(1, N).
*
* [in,out] C
*          C is DOUBLE PRECISION array, dimension (LDC, N)
*          On entry, C contains the right-hand-side of the first matrix
*          equation in (1) or (3).
*          On exit, if IJOB = 0, 1 or 2, C has been overwritten by
*          the solution R. If IJOB = 3 or 4 and TRANS = 'N', C holds R,
*          the solution achieved during the computation of the
*          Dif-estimate.
*
* [in] LDC
*          LDC is INTEGER
*          The leading dimension of the array C. LDC >= max(1, M).
*
* [in] D
*          D is DOUBLE PRECISION array, dimension (LDD, M)
*          The upper triangular matrix D.
*
* [in] LDD
*          LDD is INTEGER
*          The leading dimension of the array D. LDD >= max(1, M).
*
* [in] E
*          E is DOUBLE PRECISION array, dimension (LDE, N)
*          The upper triangular matrix E.
*
* [in] LDE
*          LDE is INTEGER
*          The leading dimension of the array E. LDE >= max(1, N).
*
* [in,out] F
*          F is DOUBLE PRECISION array, dimension (LDF, N)
*          On entry, F contains the right-hand-side of the second matrix
*          equation in (1) or (3).
*          On exit, if IJOB = 0, 1 or 2, F has been overwritten by
*          the solution L. If IJOB = 3 or 4 and TRANS = 'N', F holds L,
*          the solution achieved during the computation of the
*          Dif-estimate.
*
* [in] LDF
*          LDF is INTEGER
*          The leading dimension of the array F. LDF >= max(1, M).
*
* [out] DIF
*          DIF is DOUBLE PRECISION
*          On exit DIF is the reciprocal of a lower bound of the
*          reciprocal of the Dif-function, i.e. DIF is an upper bound of
*          Dif[(A,D), (B,E)] = sigma_min(Z), where Z as in (2).
*          IF IJOB = 0 or TRANS = 'T', DIF is not touched.
*
* [out] SCALE
*          SCALE is DOUBLE PRECISION
*          On exit SCALE is the scaling factor in (1) or (3).
*          If 0 < SCALE < 1, C and F hold the solutions R and L, resp.,
*          to a slightly perturbed system but the input matrices A, B, D
*          and E have not been changed. If SCALE = 0, C and F hold the
*          solutions R and L, respectively, to the homogeneous system
*          with C = F = 0. Normally, SCALE = 1.
*
* [out] WORK
*          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
* [in] LWORK
*          LWORK is INTEGER
*          The dimension of the array WORK. LWORK > = 1.
*          If IJOB = 1 or 2 and TRANS = 'N', LWORK >= max(1,2*M*N).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
* [out] IWORK
*          IWORK is INTEGER array, dimension (M+N+6)
*
* [out] INFO
*          INFO is INTEGER
*            =0: successful exit
*            <0: If INFO = -i, the i-th argument had an illegal value.
*            >0: (A, D) and (B, E) have common or close eigenvalues.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
tgsyl(const char * trans,const i_type ijob,i_type m,i_type n, const V *a,
      i_type lda, const V *b,i_type ldb, V *c, i_type ldc, const V *d,i_type ldd, 
      const V *e,i_type lde, V *f,i_type ldf, 
      typename details::real_type<V>::type *scale, 
      typename details::real_type<V>::type *dif, V *work, i_type lwork, 
      i_type* iwork, i_type *info);

BLAS_EXPORT void    ctgsyl(const char * trans,const i_type ijob,i_type m,
                           i_type n, c_type *a,i_type lda, c_type *b,i_type ldb, 
                           c_type *c,i_type ldc, c_type *d,i_type ldd, c_type *e,
                           i_type lde, c_type *f,i_type ldf, s_type *scale, s_type *dif,  
                           c_type *work, i_type lwork, i_type *iwork, i_type *info);
BLAS_EXPORT void    stgsyl(const char * trans,const i_type ijob,i_type m,i_type n,
                           s_type *a,i_type lda, s_type *b,i_type ldb, s_type *c,
                           i_type ldc, s_type *d,i_type ldd, s_type *e,i_type lde, 
                           s_type *f,i_type ldf, s_type *scale, s_type *dif,  
                           s_type *work, i_type lwork, i_type *iwork, i_type *info);
BLAS_EXPORT void    dtgsyl(const char * trans,const i_type ijob,i_type m,i_type n,
                           d_type *a,i_type lda, d_type *b,i_type ldb, d_type *c,
                           i_type ldc, d_type *d,i_type ldd, d_type *e,i_type lde, 
                           d_type *f,i_type ldf, d_type *scale, d_type *dif,  
                           d_type *work, i_type lwork, i_type *iwork, i_type *info);
BLAS_EXPORT void    ztgsyl(const char * trans,const i_type ijob,i_type m,i_type n,
                           z_type *a,i_type lda, z_type *b,i_type ldb, z_type *c,
                           i_type ldc, z_type *d,i_type ldd, z_type *e,i_type lde, 
                           z_type *f,i_type ldf, d_type *scale, d_type *dif,  
                           z_type *work, i_type lwork, i_type *iwork, i_type *info);

//-----------------------------------------------------------------------
//                          GEQRF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGEQRF computes a QR factorization of a real M-by-N matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(M,N)-by-N upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of min(m,n) elementary reflectors (see Further
*          Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is
*          the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
geqrf(i_type m,i_type n,V *a,i_type lda,V *tau,V *work,i_type lwork,i_type *info);

BLAS_EXPORT void    cgeqrf(i_type m,i_type n,c_type *a,i_type lda,c_type *tau,
                           c_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    sgeqrf(i_type m,i_type n,s_type *a,i_type lda,s_type *tau,
                           s_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    dgeqrf(i_type m,i_type n,d_type *a,i_type lda,d_type *tau,
                           d_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    zgeqrf(i_type m,i_type n,z_type *a,i_type lda,z_type *tau,
                           z_type *work,i_type lwork,i_type *info);

//-----------------------------------------------------------------------
//                          GEQP3
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  Computes the QR factorization of a general m-by-n matrix
*  with column pivoting using level 3 BLAS.
*  A * P = Q * R.
*
*  Arguments
*  =========
*
* * [in] M
* 
*          M is INTEGER
*          The number of rows of the matrix A. M >= 0.
* 
*
* [in] N
* 
*          N is INTEGER
*          The number of columns of the matrix A.  N >= 0.
* 
*
* [in,out] A
* 
*          A is COMPLEX array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the upper triangle of the array contains the
*          min(M,N)-by-N upper trapezoidal matrix R; the elements below
*          the diagonal, together with the array TAU, represent the
*          unitary matrix Q as a product of min(M,N) elementary
*          reflectors.
* 
*
* [in] LDA
* 
*          LDA is INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
* 
*
* [in,out] JPVT
* 
*          JPVT is INTEGER array, dimension (N)
*          On entry, if JPVT(J).ne.0, the J-th column of A is permuted
*          to the front of A*P (a leading column); if JPVT(J)=0,
*          the J-th column of A is a free column.
*          On exit, if JPVT(J)=K, then the J-th column of A*P was the
*          the K-th column of A.
* 
*
* [out] TAU
* 
*          TAU is COMPLEX array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors.
* 
*
* [out] WORK
* 
*          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
* 
*
* [in] LWORK
* 
*          LWORK is INTEGER
*          The dimension of the array WORK. LWORK >= N+1.
*          For optimal performance LWORK >= 2*N + ( N+1 )*NB, where NB
*          is the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
* 
* [out] INFO
* 
*          INFO is INTEGER
*          = 0: successful exit.
*          < 0: if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v**H
*
*  where tau is a real/complex scalar, and v is a real/complex vector
*  with v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in
*  A(i+1:m,i), and tau in TAU(i).
* 
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
geqp3(i_type m, i_type n, V *a, i_type lda, i_type* jpvt, V *tau, V *work, 
      i_type lwork, i_type *info);

BLAS_EXPORT void    cgeqp3(i_type m,i_type n,c_type *a,i_type lda,i_type* jpvt,
                           c_type *tau,c_type *work,i_type lwork,s_type* rwork,
                           i_type *info);
BLAS_EXPORT void    sgeqp3(i_type m,i_type n,s_type *a,i_type lda,i_type* jpvt,
                           s_type *tau, s_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    dgeqp3(i_type m,i_type n,d_type *a,i_type lda,i_type* jpvt,
                           d_type *tau, d_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    zgeqp3(i_type m,i_type n,z_type *a,i_type lda,i_type* jpvt,
                           z_type *tau, z_type *work,i_type lwork,d_type* rwork, 
                           i_type *info);

//-----------------------------------------------------------------------
//                          ORGQR/UNGQR
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*   DORGQR generates an M-by-N real matrix Q with orthonormal columns,
* which is defined as the first N columns of a product of K elementary
* reflectors of order M
*
*       Q  =  H(1) H(2) . . . H(k)
*
* as returned by DGEQRF.
**  Arguments:
*  ==========
*
* [in] M
* 
*          M is INTEGER
*          The number of rows of the matrix Q. M >= 0.
* 
*
* [in] N
* 
*          N is INTEGER
*          The number of columns of the matrix Q. M >= N >= 0.
* 
*
* [in] K
* 
*          K is INTEGER
*          The number of elementary reflectors whose product defines the
*          matrix Q. N >= K >= 0.
* 
*
* [in,out] A
* 
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the i-th column must contain the vector which
*          defines the elementary reflector H(i), for i = 1,2,...,k, as
*          returned by DGEQRF in the first k columns of its array
*          argument A.
*          On exit, the M-by-N matrix Q.
* 
*
* [in] LDA
* 
*          LDA is INTEGER
*          The first dimension of the array A. LDA >= max(1,M).
* 
*
* [in] TAU
* 
*          TAU is DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
* 
*
* [out] WORK
* 
*          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
* 
*
* [in] LWORK
* 
*          LWORK is INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
* 
*
* [out] INFO
* 
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument has an illegal value
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
orgqr(i_type m, i_type n, i_type k, V *a, i_type lda, const V *tau, V *work, 
      i_type lwork, i_type *info);

BLAS_EXPORT void    sorgqr(i_type m, i_type n, i_type k, s_type *a, i_type lda,
                           const s_type *tau, s_type *work, i_type lwork, 
                           i_type *info);
BLAS_EXPORT void    dorgqr(i_type m, i_type n, i_type k, d_type *a, i_type lda, 
                           const d_type *tau, d_type *work, i_type lwork, 
                           i_type *info);

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ungqr(i_type m, i_type n, i_type k, V *a, i_type lda, const V *tau, V *work, 
      i_type lwork, i_type *info);

BLAS_EXPORT void    cungqr(i_type m, i_type n, i_type k, c_type *a, i_type lda,
                           const c_type *tau, c_type *work, i_type lwork, 
                           i_type *info);
BLAS_EXPORT void    zungqr(i_type m, i_type n, i_type k, z_type *a, i_type lda, 
                           const z_type *tau, z_type *work, i_type lwork,
                           i_type *info);

//-----------------------------------------------------------------------
//                          DGEHRD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGEHRD reduces a real general matrix A to upper Hessenberg form H by
*  an orthogonal similarity transformation:  Q**T * A * Q = H .
*
*  Arguments
*  =========
*
* [in] N
* 
*          N is INTEGER
*          The order of the matrix A.  N >= 0.
*
* [in] ILO
* 
*          ILO is INTEGER
*
* [in] IHI
* 
*          IHI is INTEGER
*
*          It is assumed that A is already upper triangular in rows
*          and columns 1:ILO-1 and IHI+1:N. ILO and IHI are normally
*          set by a previous call to DGEBAL; otherwise they should be
*          set to 1 and N respectively. See Further Details.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
* [in,out] A
* 
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the N-by-N general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the upper Hessenberg matrix H, and the
*          elements below the first subdiagonal, with the array TAU,
*          represent the orthogonal matrix Q as a product of elementary
*          reflectors. See Further Details.
*
* [in] LDA
* 
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
* [out] TAU
* 
*          TAU is DOUBLE PRECISION array, dimension (N-1)
*          The scalar factors of the elementary reflectors (see Further
*          Details). Elements 1:ILO-1 and IHI:N-1 of TAU are set to
*          zero.
*
* [out] WORK
* 
*          WORK is DOUBLE PRECISION array, dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
* [in] LWORK
* 
*          LWORK is INTEGER
*          The length of the array WORK.  LWORK >= max(1,N).
*          For optimum performance LWORK >= N*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
* [out] INFO
* 
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
* 
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of (ihi-ilo) elementary
*  reflectors
*
*     Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v**T
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i) = 0, v(i+1) = 1 and v(ihi+1:n) = 0; v(i+2:ihi) is stored on
*  exit in A(i+2:ihi,i), and tau in TAU(i).
*
*  The contents of A are illustrated by the following example, with
*  n = 7, ilo = 2 and ihi = 6:
*
*  on entry,                        on exit,
*
*  ( a   a   a   a   a   a   a )    (  a   a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      a   h   h   h   h   a )
*  (     a   a   a   a   a   a )    (      h   h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  h   h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  h   h   h   h )
*  (     a   a   a   a   a   a )    (      v2  v3  v4  h   h   h )
*  (                         a )    (                          a )
*
*  where a denotes an element of the original matrix A, h denotes a
*  modified element of the upper Hessenberg matrix H, and vi denotes an
*  element of the vector defining H(i).
*
*  This file is a slight modification of LAPACK-3.0's DGEHRD
*  subroutine incorporating improvements proposed by Quintana-Orti and
*  Van de Geijn (2006). (See DLAHR2.)
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gehrd(i_type n,i_type ilo, i_type ihi, V *a, i_type lda, V *tau,V *work,
      i_type lwork,i_type *info);

BLAS_EXPORT void    cgehrd(i_type n,i_type ilo, i_type ihi, c_type *a,
                           i_type lda,c_type *tau, c_type *work,i_type lwork,
                           i_type *info);
BLAS_EXPORT void    sgehrd(i_type n,i_type ilo, i_type ihi, s_type *a,i_type lda,
                           s_type *tau, s_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    dgehrd(i_type n,i_type ilo, i_type ihi, d_type *a,i_type lda,
                           d_type *tau, d_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    zgehrd(i_type n,i_type ilo, i_type ihi, z_type *a,i_type lda,
                           z_type *tau, z_type *work,i_type lwork,i_type *info);

//-----------------------------------------------------------------------
//                          DSYTRD/ZHETRD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*   
*   DSYTRD reduces a real symmetric matrix A to real symmetric
*   tridiagonal form T by an orthogonal similarity transformation:
*   Q**T * A * Q = T.
*   
*  Arguments
*  =========
*   [in] UPLO
*   
*            UPLO is CHARACTER*1
*            = 'U':  Upper triangle of A is stored;
*            = 'L':  Lower triangle of A is stored.
*   
*  
*   [in] N
*   
*            N is INTEGER
*            The order of the matrix A.  N >= 0.
*   
*  
*   [in,out] A
*   
*            A is DOUBLE PRECISION array, dimension (LDA,N)
*            On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*            N-by-N upper triangular part of A contains the upper
*            triangular part of the matrix A, and the strictly lower
*            triangular part of A is not referenced.  If UPLO = 'L', the
*            leading N-by-N lower triangular part of A contains the lower
*            triangular part of the matrix A, and the strictly upper
*            triangular part of A is not referenced.
*            On exit, if UPLO = 'U', the diagonal and first superdiagonal
*            of A are overwritten by the corresponding elements of the
*            tridiagonal matrix T, and the elements above the first
*            superdiagonal, with the array TAU, represent the orthogonal
*            matrix Q as a product of elementary reflectors; if UPLO
*            = 'L', the diagonal and first subdiagonal of A are over-
*            written by the corresponding elements of the tridiagonal
*            matrix T, and the elements below the first subdiagonal, with
*            the array TAU, represent the orthogonal matrix Q as a product
*            of elementary reflectors. See Further Details.
*   
*  
*   [in] LDA
*   
*            LDA is INTEGER
*            The leading dimension of the array A.  LDA >= max(1,N).
*   
*  
*   [out] D
*   
*            D is DOUBLE PRECISION array, dimension (N)
*            The diagonal elements of the tridiagonal matrix T:
*            D(i) = A(i,i).
*   
*  
*   [out] E
*   
*            E is DOUBLE PRECISION array, dimension (N-1)
*            The off-diagonal elements of the tridiagonal matrix T:
*            E(i) = A(i,i+1) if UPLO = 'U', E(i) = A(i+1,i) if UPLO = 'L'.
*   
*  
*   [out] TAU
*   
*            TAU is DOUBLE PRECISION array, dimension (N-1)
*            The scalar factors of the elementary reflectors (see Further
*            Details).
*   
*  
*   [out] WORK
*   
*            WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*   
*  
*   [in] LWORK
*   
*            LWORK is INTEGER
*            The dimension of the array WORK.  LWORK >= 1.
*            For optimum performance LWORK >= N*NB, where NB is the
*            optimal blocksize.
*  
*            If LWORK = -1, then a workspace query is assumed; the routine
*            only calculates the optimal size of the WORK array, returns
*            this value as the first entry of the WORK array, and no error
*            message related to LWORK is issued by XERBLA.
*   
*  
*   [out] INFO
*   
*            INFO is INTEGER
*            = 0:  successful exit
*            < 0:  if INFO = -i, the i-th argument had an illegal value
*   
*   Further Details:
* 
*    If UPLO = 'U', the matrix Q is represented as a product of elementary
*    reflectors
*  
*       Q = H(n-1) . . . H(2) H(1).
*  
*    Each H(i) has the form
*  
*       H(i) = I - tau * v * v**T
*  
*    where tau is a real scalar, and v is a real vector with
*    v(i+1:n) = 0 and v(i) = 1; v(1:i-1) is stored on exit in
*    A(1:i-1,i+1), and tau in TAU(i).
*  
*    If UPLO = 'L', the matrix Q is represented as a product of elementary
*    reflectors
*  
*       Q = H(1) H(2) . . . H(n-1).
*  
*    Each H(i) has the form
*  
*       H(i) = I - tau * v * v**T
*  
*    where tau is a real scalar, and v is a real vector with
*    v(1:i) = 0 and v(i+1) = 1; v(i+2:n) is stored on exit in A(i+2:n,i),
*    and tau in TAU(i).
*  
*    The contents of A on exit are illustrated by the following examples
*    with n = 5:
*  
*    if UPLO = 'U':                       if UPLO = 'L':
*  
*      (  d   e   v2  v3  v4 )              (  d                  )
*      (      d   e   v3  v4 )              (  e   d              )
*      (          d   e   v4 )              (  v1  e   d          )
*      (              d   e  )              (  v1  v2  e   d      )
*      (                  d  )              (  v1  v2  v3  e   d  )
*  
*    where d and e denote diagonal and off-diagonal elements of T, and vi
*    denotes an element of the vector defining H(i).
*   
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
sytrd(const char* uplo, i_type n,V *a, i_type lda, V *d, V *e, V *tau,V *work,
      i_type lwork,i_type *info);

BLAS_EXPORT void    ssytrd(const char* uplo, i_type n, s_type *a,i_type lda,
                           s_type *d, s_type *e, s_type *tau, s_type *work,
                           i_type lwork,i_type *info);
BLAS_EXPORT void    dsytrd(const char* uplo, i_type n, d_type *a,i_type lda,
                           d_type *d, d_type *e, d_type *tau, d_type *work,
                           i_type lwork,i_type *info);

BLAS_EXPORT void    chetrd(const char* uplo, i_type n, c_type *a,i_type lda,
                           s_type *d, s_type *e, c_type *tau, c_type *work,
                           i_type lwork,i_type *info);
BLAS_EXPORT void    zhetrd(const char* uplo, i_type n, z_type *a,i_type lda,
                           d_type *d, d_type *e, z_type *tau, z_type *work,
                           i_type lwork,i_type *info);

//-----------------------------------------------------------------------
//                          ORGTR/UNGTR
//-----------------------------------------------------------------------
/*
*   ZUNGTR generates a complex unitary matrix Q which is defined as the
*   product of n-1 elementary reflectors of order N, as returned by
*   ZHETRD:
*  
*   if UPLO = 'U', Q = H(n-1) . . . H(2) H(1),
*  
*   if UPLO = 'L', Q = H(1) H(2) . . . H(n-1).
*
*   Arguments
* 
*   [in] UPLO
*   
*            UPLO is CHARACTER*1
*            = 'U': Upper triangle of A contains elementary reflectors
*                   from ZHETRD;
*            = 'L': Lower triangle of A contains elementary reflectors
*                   from ZHETRD.
*   
*  
*   [in] N
*   
*            N is INTEGER
*            The order of the matrix Q. N >= 0.
*   
*  
*   [in,out] A
*   
*            A is COMPLEX*16 array, dimension (LDA,N)
*            On entry, the vectors which define the elementary reflectors,
*            as returned by ZHETRD.
*            On exit, the N-by-N unitary matrix Q.
*   
*  
*   [in] LDA
*   
*            LDA is INTEGER
*            The leading dimension of the array A. LDA >= N.
*   
*  
*   [in] TAU
*   
*            TAU is COMPLEX*16 array, dimension (N-1)
*            TAU(i) must contain the scalar factor of the elementary
*            reflector H(i), as returned by ZHETRD.
*   
*  
*   [out] WORK
*   
*            WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
*            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*   
*  
*   [in] LWORK
*   
*            LWORK is INTEGER
*            The dimension of the array WORK. LWORK >= N-1.
*            For optimum performance LWORK >= (N-1)*NB, where NB is
*            the optimal blocksize.
*  
*            If LWORK = -1, then a workspace query is assumed; the routine
*            only calculates the optimal size of the WORK array, returns
*            this value as the first entry of the WORK array, and no error
*            message related to LWORK is issued by XERBLA.
*   
*  
*   [out] INFO
*   
*            INFO is INTEGER
*            = 0:  successful exit
*            < 0:  if INFO = -i, the i-th argument had an illegal value
*   
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
orgtr(const char* uplo, i_type n,V *a, i_type lda, V *tau,V *work,
      i_type lwork,i_type *info);

BLAS_EXPORT void    sorgtr(const char* uplo, i_type n, s_type *a,i_type lda,
                           s_type *tau, s_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    dorgtr(const char* uplo, i_type n, d_type *a,i_type lda,
                           d_type *tau, d_type *work,i_type lwork,i_type *info);

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ungtr(const char* uplo, i_type n,V *a, i_type lda, V *tau,V *work,i_type lwork,
      i_type *info);

BLAS_EXPORT void    cungtr(const char* uplo, i_type n, c_type *a,i_type lda,
                           c_type *tau, c_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    zungtr(const char* uplo, i_type n, z_type *a,i_type lda,
                           z_type *tau, z_type *work,i_type lwork,i_type *info);

//-----------------------------------------------------------------------
//                          ORGHR/UNGHR
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*   CUNGHR generates a complex unitary matrix Q which is defined as the
*   product of IHI-ILO elementary reflectors of order N, as returned by
*   CGEHRD:
*  
*   Q = H(ilo) H(ilo+1) . . . H(ihi-1).
*   
**  Arguments:
*  ==========
*
*   [in] N
*   
*            N is INTEGER
*            The order of the matrix Q. N >= 0.
*   
*  
*   [in] ILO
*   
*            ILO is INTEGER
*   
*  
*   [in] IHI
*   
*            IHI is INTEGER
*  
*            ILO and IHI must have the same values as in the previous call
*            of CGEHRD. Q is equal to the unit matrix except in the
*            submatrix Q(ilo+1:ihi,ilo+1:ihi).
*            1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*   
*  
*   [in,out] A
*   
*            A is COMPLEX array, dimension (LDA,N)
*            On entry, the vectors which define the elementary reflectors,
*            as returned by CGEHRD.
*            On exit, the N-by-N unitary matrix Q.
*   
*  
*   [in] LDA
*   
*            LDA is INTEGER
*            The leading dimension of the array A. LDA >= max(1,N).
*   
*  
*   [in] TAU
*   
*            TAU is COMPLEX array, dimension (N-1)
*            TAU(i) must contain the scalar factor of the elementary
*            reflector H(i), as returned by CGEHRD.
*   
*  
*   [out] WORK
*   
*            WORK is COMPLEX array, dimension (MAX(1,LWORK))
*            On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*   
*  
*   [in] LWORK
*   
*            LWORK is INTEGER
*            The dimension of the array WORK. LWORK >= IHI-ILO.
*            For optimum performance LWORK >= (IHI-ILO)*NB, where NB is
*            the optimal blocksize.
*  
*            If LWORK = -1, then a workspace query is assumed; the routine
*            only calculates the optimal size of the WORK array, returns
*            this value as the first entry of the WORK array, and no error
*            message related to LWORK is issued by XERBLA.
*   
*  
*   [out] INFO
*   
*            INFO is INTEGER
*            = 0:  successful exit
*            < 0:  if INFO = -i, the i-th argument had an illegal value
*/   

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
orghr(i_type n,i_type ilo, i_type ihi, V *a, i_type lda, V *tau,V *work,i_type lwork,
      i_type *info);

BLAS_EXPORT void    sorghr(i_type n,i_type ilo, i_type ihi, s_type *a,i_type lda,
                           s_type *tau, s_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    dorghr(i_type n,i_type ilo, i_type ihi, d_type *a,i_type lda,
                           d_type *tau, d_type *work,i_type lwork,i_type *info);

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
unghr(i_type m, i_type n, i_type k, V *a, i_type lda, V *tau, V *work, i_type lwork,
      i_type *info);

BLAS_EXPORT void    cunghr(i_type n,i_type ilo, i_type ihi, c_type *a,i_type lda,
                           c_type *tau, c_type *work,i_type lwork,i_type *info);
BLAS_EXPORT void    zunghr(i_type n,i_type ilo, i_type ihi, z_type *a,i_type lda,
                           z_type *tau, z_type *work,i_type lwork,i_type *info);

//-----------------------------------------------------------------------
//                          POTRF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*  DPOTRF computes the Cholesky factorization of a real symmetric
*  positive definite matrix A.
*
*  The factorization has the form
*     A = U**T * U,  if UPLO = 'U', or
*     A = L  * L**T,  if UPLO = 'L',
*  where U is an upper triangular matrix and L is lower triangular.
*
*  This is the block version of the algorithm, calling Level 3 BLAS.
*
*  Arguments:
*  ==========
*
*  [in] UPLO
*          UPLO is CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  [in] N
*          N is INTEGER
*          The order of the matrix A.  N >= 0.
*
*  [in,out] A
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, if INFO = 0, the factor U or L from the Cholesky
*          factorization A = U**T*U or A = L*L**T.
*
*  [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  [out] INFO
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*                positive definite, and the factorization could not be
*                completed.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
potrf(const char *uplo, i_type n, V *a, i_type lda, i_type *info);

BLAS_EXPORT void    cpotrf(const char *uplo, i_type n, c_type *a, i_type lda,
                           i_type *info);
BLAS_EXPORT void    spotrf(const char *uplo, i_type n, s_type *a, i_type lda,
                           i_type *info);
BLAS_EXPORT void    dpotrf(const char *uplo, i_type n, d_type *a, i_type lda,
                           i_type *info);
BLAS_EXPORT void    zpotrf(const char *uplo, i_type n, z_type *a, i_type lda,
                           i_type *info);

//-----------------------------------------------------------------------
//                          PBTRF
//-----------------------------------------------------------------------
/*
*	Purpose
*	=======
*
*	DPBTRF computes the Cholesky factorization of a real symmetric
*	positive definite band matrix A.
*
*	The factorization has the form
*		A = U**T * U,  if UPLO = 'U', or
*		A = L  * L**T,  if UPLO = 'L',
*	where U is an upper triangular matrix and L is lower triangular.
*
*
*	Arguments:
*	==========
*
*	[in] UPLO
*
*          UPLO is CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*	[in] N
*          N is INTEGER
*          The order of the matrix A.  N >= 0.
*
*	[in] KD
*          KD is INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*	[in,out] AB
*          AB is DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, if INFO = 0, the triangular factor U or L from the
*          Cholesky factorization A = U**T*U or A = L*L**T of the band
*          matrix A, in the same storage format as A.
*
*	[in] LDAB
*          LDAB is INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*	[out] INFO
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the leading minor of order i is not
*	              positive definite, and the factorization could not be
*                completed.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
pbtrf(const char *uplo, i_type n, i_type kd, V *ab, i_type ldab, i_type *info);

BLAS_EXPORT void    cpbtrf(const char *uplo, i_type n, i_type kd, c_type *ab, 
                           i_type ldab, i_type *info);
BLAS_EXPORT void    spbtrf(const char *uplo, i_type n, i_type kd, s_type *ab,
                           i_type ldab, i_type *info);
BLAS_EXPORT void    dpbtrf(const char *uplo, i_type n, i_type kd, d_type *ab, 
                           i_type ldab, i_type *info);
BLAS_EXPORT void    zpbtrf(const char *uplo, i_type n, i_type kd, z_type *ab, 
                           i_type ldab, i_type *info);

//-----------------------------------------------------------------------
//                          SYTRF_ROOK
//-----------------------------------------------------------------------
/*> \par Purpose:
 *  =============
 *>
 *> \verbatim
 *>
 *> DSYTRF_ROOK computes the factorization of a real symmetric matrix A
 *> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
 *> The form of the factorization is
 *>
 *>    A = U*D*U**T  or  A = L*D*L**T
 *>
 *> where U (or L) is a product of permutation and unit upper (lower)
 *> triangular matrices, and D is symmetric and block diagonal with
 *> 1-by-1 and 2-by-2 diagonal blocks.
 *>
 *> This is the blocked version of the algorithm, calling Level 3 BLAS.
 *> \endverbatim
 *
 *  Arguments:
 *  ==========
 *
 *> \param[in] UPLO
 *> \verbatim
 *>          UPLO is CHARACTER*1
 *>          = 'U':  Upper triangle of A is stored;
 *>          = 'L':  Lower triangle of A is stored.
 *> \endverbatim
 *>
 *> \param[in] N
 *> \verbatim
 *>          N is INTEGER
 *>          The order of the matrix A.  N >= 0.
 *> \endverbatim
 *>
 *> \param[in,out] A
 *> \verbatim
 *>          A is DOUBLE PRECISION array, dimension (LDA,N)
 *>          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
 *>          N-by-N upper triangular part of A contains the upper
 *>          triangular part of the matrix A, and the strictly lower
 *>          triangular part of A is not referenced.  If UPLO = 'L', the
 *>          leading N-by-N lower triangular part of A contains the lower
 *>          triangular part of the matrix A, and the strictly upper
 *>          triangular part of A is not referenced.
 *>
 *>          On exit, the block diagonal matrix D and the multipliers used
 *>          to obtain the factor U or L (see below for further details).
 *> \endverbatim
 *>
 *> \param[in] LDA
 *> \verbatim
 *>          LDA is INTEGER
 *>          The leading dimension of the array A.  LDA >= max(1,N).
 *> \endverbatim
 *>
 *> \param[out] IPIV
 *> \verbatim
 *>          IPIV is INTEGER array, dimension (N)
 *>          Details of the interchanges and the block structure of D.
 *>
 *>          If UPLO = 'U':
 *>               If IPIV(k) > 0, then rows and columns k and IPIV(k)
 *>               were interchanged and D(k,k) is a 1-by-1 diagonal block.
 *>
 *>               If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
 *>               columns k and -IPIV(k) were interchanged and rows and
 *>               columns k-1 and -IPIV(k-1) were inerchaged,
 *>               D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
 *>
 *>          If UPLO = 'L':
 *>               If IPIV(k) > 0, then rows and columns k and IPIV(k)
 *>               were interchanged and D(k,k) is a 1-by-1 diagonal block.
 *>
 *>               If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
 *>               columns k and -IPIV(k) were interchanged and rows and
 *>               columns k+1 and -IPIV(k+1) were inerchaged,
 *>               D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
 *> \endverbatim
 *>
 *> \param[out] WORK
 *> \verbatim
 *>          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK)).
 *>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *> \endverbatim
 *>
 *> \param[in] LWORK
 *> \verbatim
 *>          LWORK is INTEGER
 *>          The length of WORK.  LWORK >=1.  For best performance
 *>          LWORK >= N*NB, where NB is the block size returned by ILAENV.
 *>
 *>          If LWORK = -1, then a workspace query is assumed; the routine
 *>          only calculates the optimal size of the WORK array, returns
 *>          this value as the first entry of the WORK array, and no error
 *>          message related to LWORK is issued by XERBLA.
 *> \endverbatim
 *>
 *> \param[out] INFO
 *> \verbatim
 *>          INFO is INTEGER
 *>          = 0:  successful exit
 *>          < 0:  if INFO = -i, the i-th argument had an illegal value
 *>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
 *>                has been completed, but the block diagonal matrix D is
 *>                exactly singular, and division by zero will occur if it
 *>                is used to solve a system of equations.
 *> \endverbatim
 *
 *  Authors:
 *  ========
 *
 *> \author Univ. of Tennessee
 *> \author Univ. of California Berkeley
 *> \author Univ. of Colorado Denver
 *> \author NAG Ltd.
 *
 *> \date April 2012
 *
 *> \ingroup doubleSYcomputational
 *
 *> \par Further Details:
 *  =====================
 *>
 *> \verbatim
 *>
 *>  If UPLO = 'U', then A = U*D*U**T, where
 *>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
 *>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
 *>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
 *>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
 *>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
 *>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
 *>
 *>             (   I    v    0   )   k-s
 *>     U(k) =  (   0    I    0   )   s
 *>             (   0    0    I   )   n-k
 *>                k-s   s   n-k
 *>
 *>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
 *>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
 *>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
 *>
 *>  If UPLO = 'L', then A = L*D*L**T, where
 *>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
 *>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
 *>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
 *>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
 *>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
 *>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
 *>
 *>             (   I    0     0   )  k-1
 *>     L(k) =  (   0    I     0   )  s
 *>             (   0    v     I   )  n-k-s+1
 *>                k-1   s  n-k-s+1
 *>
 *>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
 *>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
 *>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
 *> \endverbatim
 *
 *> \par Contributors:
 *  ==================
 *>
 *> \verbatim
 *>
 *>   April 2012, Igor Kozachenko,
 *>                  Computer Science Division,
 *>                  University of California, Berkeley
 *>
 *>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
 *>                  School of Mathematics,
 *>                  University of Manchester
 *>
 *> \endverbatim
 */

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
sytrf_rook(const char *uplo, i_type n, V *a, i_type lda, i_type *ipiv, V* work,
      i_type lwork, i_type *info);

BLAS_EXPORT void ssytrf_rook(const char *uplo, i_type n, s_type *a, i_type lda, 
                        i_type *ipiv, s_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void dsytrf_rook(const char *uplo, i_type n, d_type *a, i_type lda, 
                        i_type *ipiv, d_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void csytrf_rook(const char *uplo, i_type n, c_type *a, i_type lda, 
                        i_type *ipiv, c_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void zsytrf_rook(const char *uplo, i_type n, z_type *a, i_type lda, 
                        i_type *ipiv, z_type *work, i_type lwork, i_type *info);

//-----------------------------------------------------------------------
//                          SYTRF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSYTRF computes the factorization of a real symmetric matrix A using
*  the Bunch-Kaufman diagonal pivoting method.  The form of the
*  factorization is
*
*     A = U*D*U**T  or  A = L*D*L**T
*
*  where U (or L) is a product of permutation and unit upper (lower)
*  triangular matrices, and D is symmetric and block diagonal with
*  1-by-1 and 2-by-2 diagonal blocks.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, the block diagonal matrix D and the multipliers used
*          to obtain the factor U or L (see below for further details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D.
*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*          interchanged and D(k,k) is a 1-by-1 diagonal block.
*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of WORK.  LWORK >=1.  For best performance
*          LWORK >= N*NB, where NB is the block size returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
*                has been completed, but the block diagonal matrix D is
*                exactly singular, and division by zero will occur if it
*                is used to solve a system of equations.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', then A = U*D*U**T, where
*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    v    0   )   k-s
*     U(k) =  (   0    I    0   )   s
*             (   0    0    I   )   n-k
*                k-s   s   n-k
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
*
*  If UPLO = 'L', then A = L*D*L**T, where
*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    0     0   )  k-1
*     L(k) =  (   0    I     0   )  s
*             (   0    v     I   )  n-k-s+1
*                k-1   s  n-k-s+1
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
sytrf(const char *uplo, i_type n, V *a, i_type lda, i_type *ipiv, V* work,
      i_type lwork, i_type *info);

BLAS_EXPORT void ssytrf(const char *uplo, i_type n, s_type *a, i_type lda, 
                        i_type *ipiv, s_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void dsytrf(const char *uplo, i_type n, d_type *a, i_type lda, 
                        i_type *ipiv, d_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void csytrf(const char *uplo, i_type n, c_type *a, i_type lda, 
                        i_type *ipiv, c_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void zsytrf(const char *uplo, i_type n, z_type *a, i_type lda, 
                        i_type *ipiv, z_type *work, i_type lwork, i_type *info);

//-----------------------------------------------------------------------
//                          HETRF
//-----------------------------------------------------------------------
/*
*> \par Purpose:
 *  =============
 *>
 *> \verbatim
 *>
 *> ZHETRF_ROOK computes the factorization of a complex Hermitian matrix A
 *> using the bounded Bunch-Kaufman ("rook") diagonal pivoting method.
 *> The form of the factorization is
 *>
 *>    A = U*D*U**T  or  A = L*D*L**T
 *>
 *> where U (or L) is a product of permutation and unit upper (lower)
 *> triangular matrices, and D is Hermitian and block diagonal with
 *> 1-by-1 and 2-by-2 diagonal blocks.
 *>
 *> This is the blocked version of the algorithm, calling Level 3 BLAS.
 *> \endverbatim
 *
 *  Arguments:
 *  ==========
 *
 *> \param[in] UPLO
 *> \verbatim
 *>          UPLO is CHARACTER*1
 *>          = 'U':  Upper triangle of A is stored;
 *>          = 'L':  Lower triangle of A is stored.
 *> \endverbatim
 *>
 *> \param[in] N
 *> \verbatim
 *>          N is INTEGER
 *>          The order of the matrix A.  N >= 0.
 *> \endverbatim
 *>
 *> \param[in,out] A
 *> \verbatim
 *>          A is COMPLEX*16 array, dimension (LDA,N)
 *>          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
 *>          N-by-N upper triangular part of A contains the upper
 *>          triangular part of the matrix A, and the strictly lower
 *>          triangular part of A is not referenced.  If UPLO = 'L', the
 *>          leading N-by-N lower triangular part of A contains the lower
 *>          triangular part of the matrix A, and the strictly upper
 *>          triangular part of A is not referenced.
 *>
 *>          On exit, the block diagonal matrix D and the multipliers used
 *>          to obtain the factor U or L (see below for further details).
 *> \endverbatim
 *>
 *> \param[in] LDA
 *> \verbatim
 *>          LDA is INTEGER
 *>          The leading dimension of the array A.  LDA >= max(1,N).
 *> \endverbatim
 *>
 *> \param[out] IPIV
 *> \verbatim
 *>          IPIV is INTEGER array, dimension (N)
 *>          Details of the interchanges and the block structure of D.
 *>
 *>          If UPLO = 'U':
 *>             Only the last KB elements of IPIV are set.
 *>
 *>             If IPIV(k) > 0, then rows and columns k and IPIV(k) were
 *>             interchanged and D(k,k) is a 1-by-1 diagonal block.
 *>
 *>             If IPIV(k) < 0 and IPIV(k-1) < 0, then rows and
 *>             columns k and -IPIV(k) were interchanged and rows and
 *>             columns k-1 and -IPIV(k-1) were inerchaged,
 *>             D(k-1:k,k-1:k) is a 2-by-2 diagonal block.
 *>
 *>          If UPLO = 'L':
 *>             Only the first KB elements of IPIV are set.
 *>
 *>             If IPIV(k) > 0, then rows and columns k and IPIV(k)
 *>             were interchanged and D(k,k) is a 1-by-1 diagonal block.
 *>
 *>             If IPIV(k) < 0 and IPIV(k+1) < 0, then rows and
 *>             columns k and -IPIV(k) were interchanged and rows and
 *>             columns k+1 and -IPIV(k+1) were inerchaged,
 *>             D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
 *> \endverbatim
 *>
 *> \param[out] WORK
 *> \verbatim
 *>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK)).
 *>          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
 *> \endverbatim
 *>
 *> \param[in] LWORK
 *> \verbatim
 *>          LWORK is INTEGER
 *>          The length of WORK.  LWORK >=1.  For best performance
 *>          LWORK >= N*NB, where NB is the block size returned by ILAENV.
 *>
 *>          If LWORK = -1, then a workspace query is assumed; the routine
 *>          only calculates the optimal size of the WORK array, returns
 *>          this value as the first entry of the WORK array, and no error
 *>          message related to LWORK is issued by XERBLA.
 *> \endverbatim
 *>
 *> \param[out] INFO
 *> \verbatim
 *>          INFO is INTEGER
 *>          = 0:  successful exit
 *>          < 0:  if INFO = -i, the i-th argument had an illegal value
 *>          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
 *>                has been completed, but the block diagonal matrix D is
 *>                exactly singular, and division by zero will occur if it
 *>                is used to solve a system of equations.
 *> \endverbatim
 *
 *  Authors:
 *  ========
 *
 *> \author Univ. of Tennessee
 *> \author Univ. of California Berkeley
 *> \author Univ. of Colorado Denver
 *> \author NAG Ltd.
 *
 *> \date June 2016
 *
 *> \ingroup complex16HEcomputational
 *
 *> \par Further Details:
 *  =====================
 *>
 *> \verbatim
 *>
 *>  If UPLO = 'U', then A = U*D*U**T, where
 *>     U = P(n)*U(n)* ... *P(k)U(k)* ...,
 *>  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
 *>  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
 *>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
 *>  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
 *>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
 *>
 *>             (   I    v    0   )   k-s
 *>     U(k) =  (   0    I    0   )   s
 *>             (   0    0    I   )   n-k
 *>                k-s   s   n-k
 *>
 *>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
 *>  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
 *>  and A(k,k), and v overwrites A(1:k-2,k-1:k).
 *>
 *>  If UPLO = 'L', then A = L*D*L**T, where
 *>     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
 *>  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
 *>  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
 *>  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
 *>  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
 *>  that if the diagonal block D(k) is of order s (s = 1 or 2), then
 *>
 *>             (   I    0     0   )  k-1
 *>     L(k) =  (   0    I     0   )  s
 *>             (   0    v     I   )  n-k-s+1
 *>                k-1   s  n-k-s+1
 *>
 *>  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
 *>  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
 *>  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
 *> \endverbatim
 *
 *> \par Contributors:
 *  ==================
 *>
 *> \verbatim
 *>
 *>  June 2016,  Igor Kozachenko,
 *>                  Computer Science Division,
 *>                  University of California, Berkeley
 *>
 *>  September 2007, Sven Hammarling, Nicholas J. Higham, Craig Lucas,
 *>                  School of Mathematics,
 *>                  University of Manchester
 *>
 *> \endverbatim
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
hetrf_rook(const char *uplo, i_type n, V *a, i_type lda, i_type *ipiv, V* work, 
      i_type lwork, i_type *info);

BLAS_EXPORT void shetrf_rook(const char *uplo, i_type n, s_type *a, i_type lda, 
                        i_type *ipiv, s_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void dhetrf_rook(const char *uplo, i_type n, d_type *a, i_type lda,
                        i_type *ipiv, d_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void chetrf_rook(const char *uplo, i_type n, c_type *a, i_type lda,
                        i_type *ipiv, c_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void zhetrf_rook(const char *uplo, i_type n, z_type *a, i_type lda,
                        i_type *ipiv, z_type *work, i_type lwork, i_type *info);

//-----------------------------------------------------------------------
//                          HETRF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHETRF computes the factorization of a complex Hermitian matrix A
*  using the Bunch-Kaufman diagonal pivoting method.  The form of the
*  factorization is
*
*     A = U*D*U**H  or  A = L*D*L**H
*
*  where U (or L) is a product of permutation and unit upper (lower)
*  triangular matrices, and D is Hermitian and block diagonal with
*  1-by-1 and 2-by-2 diagonal blocks.
*
*  This is the blocked version of the algorithm, calling Level 3 BLAS.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX*16 array, dimension (LDA,N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the leading
*          N-by-N upper triangular part of A contains the upper
*          triangular part of the matrix A, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of A contains the lower
*          triangular part of the matrix A, and the strictly upper
*          triangular part of A is not referenced.
*
*          On exit, the block diagonal matrix D and the multipliers used
*          to obtain the factor U or L (see below for further details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (output) INTEGER array, dimension (N)
*          Details of the interchanges and the block structure of D.
*          If IPIV(k) > 0, then rows and columns k and IPIV(k) were
*          interchanged and D(k,k) is a 1-by-1 diagonal block.
*          If UPLO = 'U' and IPIV(k) = IPIV(k-1) < 0, then rows and
*          columns k-1 and -IPIV(k) were interchanged and D(k-1:k,k-1:k)
*          is a 2-by-2 diagonal block.  If UPLO = 'L' and IPIV(k) =
*          IPIV(k+1) < 0, then rows and columns k+1 and -IPIV(k) were
*          interchanged and D(k:k+1,k:k+1) is a 2-by-2 diagonal block.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of WORK.  LWORK >=1.  For best performance
*          LWORK >= N*NB, where NB is the block size returned by ILAENV.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, D(i,i) is exactly zero.  The factorization
*                has been completed, but the block diagonal matrix D is
*                exactly singular, and division by zero will occur if it
*                is used to solve a system of equations.
*
*  Further Details
*  ===============
*
*  If UPLO = 'U', then A = U*D*U', where
*     U = P(n)*U(n)* ... *P(k)U(k)* ...,
*  i.e., U is a product of terms P(k)*U(k), where k decreases from n to
*  1 in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and U(k) is a unit upper triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    v    0   )   k-s
*     U(k) =  (   0    I    0   )   s
*             (   0    0    I   )   n-k
*                k-s   s   n-k
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(1:k-1,k).
*  If s = 2, the upper triangle of D(k) overwrites A(k-1,k-1), A(k-1,k),
*  and A(k,k), and v overwrites A(1:k-2,k-1:k).
*
*  If UPLO = 'L', then A = L*D*L', where
*     L = P(1)*L(1)* ... *P(k)*L(k)* ...,
*  i.e., L is a product of terms P(k)*L(k), where k increases from 1 to
*  n in steps of 1 or 2, and D is a block diagonal matrix with 1-by-1
*  and 2-by-2 diagonal blocks D(k).  P(k) is a permutation matrix as
*  defined by IPIV(k), and L(k) is a unit lower triangular matrix, such
*  that if the diagonal block D(k) is of order s (s = 1 or 2), then
*
*             (   I    0     0   )  k-1
*     L(k) =  (   0    I     0   )  s
*             (   0    v     I   )  n-k-s+1
*                k-1   s  n-k-s+1
*
*  If s = 1, D(k) overwrites A(k,k), and v overwrites A(k+1:n,k).
*  If s = 2, the lower triangle of D(k) overwrites A(k,k), A(k+1,k),
*  and A(k+1,k+1), and v overwrites A(k+2:n,k:k+1).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
hetrf(const char *uplo, i_type n, V *a, i_type lda, i_type *ipiv, V* work, 
      i_type lwork, i_type *info);

BLAS_EXPORT void shetrf(const char *uplo, i_type n, s_type *a, i_type lda, 
                        i_type *ipiv, s_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void dhetrf(const char *uplo, i_type n, d_type *a, i_type lda,
                        i_type *ipiv, d_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void chetrf(const char *uplo, i_type n, c_type *a, i_type lda,
                        i_type *ipiv, c_type *work, i_type lwork, i_type *info);
BLAS_EXPORT void zhetrf(const char *uplo, i_type n, z_type *a, i_type lda,
                        i_type *ipiv, z_type *work, i_type lwork, i_type *info);

//-----------------------------------------------------------------------
//                          DLALN2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLALN2 solves a system of the form  (ca A - w D ) X = s B
*  or (ca A' - w D) X = s B   with possible scaling ("s") and
*  perturbation of A.  (A' means A-transpose.)
*
*  A is an NA x NA real matrix, ca is a real scalar, D is an NA x NA
*  real diagonal matrix, w is a real or complex value, and X and B are
*  NA x 1 matrices -- real if w is real, complex if w is complex.  NA
*  may be 1 or 2.
*
*  If w is complex, X and B are represented as NA x 2 matrices,
*  the first column of each being the real part and the second
*  being the imaginary part.
*
*  "s" is a scaling factor (.LE. 1), computed by DLALN2, which is
*  so chosen that X can be computed without overflow.  X is further
*  scaled if necessary to assure that norm(ca A - w D)*norm(X) is less
*  than overflow.
*
*  If both singular values of (ca A - w D) are less than SMIN,
*  SMIN*identity will be used instead of (ca A - w D).  If only one
*  singular value is less than SMIN, one element of (ca A - w D) will be
*  perturbed enough to make the smallest singular value roughly SMIN.
*  If both singular values are at least SMIN, (ca A - w D) will not be
*  perturbed.  In any case, the perturbation will be at most some small
*  multiple of max( SMIN, ulp*norm(ca A - w D) ).  The singular values
*  are computed by infinity-norm approximations, and thus will only be
*  correct to a factor of 2 or so.
*
*  Note: all input quantities are assumed to be smaller than overflow
*  by a reasonable factor.  (See BIGNUM.)
*
*  Arguments
*  ==========
*
*  LTRANS  (input) LOGICAL
*          =.TRUE.:  A-transpose will be used.
*          =.FALSE.: A will be used (not transposed.)
*
*  NA      (input) INTEGER
*          The size of the matrix A.  It may (only) be 1 or 2.
*
*  NW      (input) INTEGER
*          1 if "w" is real, 2 if "w" is complex.  It may only be 1
*          or 2.
*
*  SMIN    (input) DOUBLE PRECISION
*          The desired lower bound on the singular values of A.  This
*          should be a safe distance away from underflow or overflow,
*          say, between (underflow/machine precision) and  (machine
*          precision * overflow ).  (See BIGNUM and ULP.)
*
*  CA      (input) DOUBLE PRECISION
*          The coefficient c, which A is multiplied by.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,NA)
*          The NA x NA matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of A.  It must be at least NA.
*
*  D1      (input) DOUBLE PRECISION
*          The 1,1 element in the diagonal matrix D.
*
*  D2      (input) DOUBLE PRECISION
*          The 2,2 element in the diagonal matrix D.  Not used if NW=1.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,NW)
*          The NA x NW matrix B (right-hand side).  If NW=2 ("w" is
*          complex), column 1 contains the real part of B and column 2
*          contains the imaginary part.
*
*  LDB     (input) INTEGER
*          The leading dimension of B.  It must be at least NA.
*
*  WR      (input) DOUBLE PRECISION
*          The real part of the scalar "w".
*
*  WI      (input) DOUBLE PRECISION
*          The imaginary part of the scalar "w".  Not used if NW=1.
*
*  X       (output) DOUBLE PRECISION array, dimension (LDX,NW)
*          The NA x NW matrix X (unknowns), as computed by DLALN2.
*          If NW=2 ("w" is complex), on exit, column 1 will contain
*          the real part of X and column 2 will contain the imaginary
*          part.
*
*  LDX     (input) INTEGER
*          The leading dimension of X.  It must be at least NA.
*
*  SCALE   (output) DOUBLE PRECISION
*          The scale factor that B must be multiplied by to insure
*          that overflow does not occur when computing X.  Thus,
*          (ca A - w D) X  will be SCALE*B, not B (ignoring
*          perturbations of A.)  It will be at most 1.
*
*  XNORM   (output) DOUBLE PRECISION
*          The infinity-norm of X, when X is regarded as an NA x NW
*          real matrix.
*
*  INFO    (output) INTEGER
*          An error flag.  It will be set to zero if no error occurs,
*          a negative number if an argument is in error, or a positive
*          number if  ca A - w D  had to be perturbed.
*          The possible values are:
*          = 0: No error occurred, and (ca A - w D) did not have to be
*                 perturbed.
*          = 1: (ca A - w D) had to be perturbed to make its smallest
*               (or only) singular value greater than SMIN.
*          NOTE: In the interests of speed, this routine does not
*                check the inputs for errors.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
laln2(bool LTRANS, i_type NA, i_type NW, V SMIN, V CA, const V* A, i_type LDA,
      V D1, V D2, const V* B, i_type LDB, V WR, V WI, V* X, i_type LDX, V* SCALE, 
      V* XNORM, i_type* INFO );

BLAS_EXPORT void slaln2(bool LTRANS, i_type NA, i_type NW, s_type SMIN, s_type CA,
                        const s_type* A, i_type LDA, s_type D1, s_type D2, 
                        const s_type* B, i_type LDB, s_type WR, s_type WI, s_type* X,
                        i_type LDX, s_type* SCALE, s_type* XNORM, i_type* INFO );
BLAS_EXPORT void dlaln2(bool LTRANS, i_type NA, i_type NW, d_type SMIN, d_type CA,
                        const d_type* A, i_type LDA, d_type D1, d_type D2, 
                        const d_type* B, i_type LDB, d_type WR, d_type WI, d_type* X,
                        i_type LDX, d_type* SCALE, d_type* XNORM, i_type* INFO );

//-----------------------------------------------------------------------
//                          DLAGV2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLAGV2 computes the Generalized Schur factorization of a real 2-by-2
*  matrix pencil (A,B) where B is upper triangular. This routine
*  computes orthogonal (rotation) matrices given by CSL, SNL and CSR,
*  SNR such that
*
*  1) if the pencil (A,B) has two real eigenvalues (include 0/0 or 1/0
*     types), then
*
*     [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
*     [  0  a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
*
*     [ b11 b12 ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
*     [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ],
*
*  2) if the pencil (A,B) has a pair of complex conjugate eigenvalues,
*     then
*
*     [ a11 a12 ] := [  CSL  SNL ] [ a11 a12 ] [  CSR -SNR ]
*     [ a21 a22 ]    [ -SNL  CSL ] [ a21 a22 ] [  SNR  CSR ]
*
*     [ b11  0  ] := [  CSL  SNL ] [ b11 b12 ] [  CSR -SNR ]
*     [  0  b22 ]    [ -SNL  CSL ] [  0  b22 ] [  SNR  CSR ]
*
*     where b11 >= b22 > 0.
*
*
*  Arguments
*  =========
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, 2)
*          On entry, the 2 x 2 matrix A.
*          On exit, A is overwritten by the ``A-part'' of the
*          generalized Schur form.
*
*  LDA     (input) INTEGER
*          THe leading dimension of the array A.  LDA >= 2.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, 2)
*          On entry, the upper triangular 2 x 2 matrix B.
*          On exit, B is overwritten by the ``B-part'' of the
*          generalized Schur form.
*
*  LDB     (input) INTEGER
*          THe leading dimension of the array B.  LDB >= 2.
*
*  ALPHAR  (output) DOUBLE PRECISION array, dimension (2)
*  ALPHAI  (output) DOUBLE PRECISION array, dimension (2)
*  BETA    (output) DOUBLE PRECISION array, dimension (2)
*          (ALPHAR(k)+i*ALPHAI(k))/BETA(k) are the eigenvalues of the
*          pencil (A,B), k=1,2, i = sqrt(-1).  Note that BETA(k) may
*          be zero.
*
*  CSL     (output) DOUBLE PRECISION
*          The cosine of the left rotation matrix.
*
*  SNL     (output) DOUBLE PRECISION
*          The sine of the left rotation matrix.
*
*  CSR     (output) DOUBLE PRECISION
*          The cosine of the right rotation matrix.
*
*  SNR     (output) DOUBLE PRECISION
*          The sine of the right rotation matrix.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
lagv2(V* A, i_type LDA, V* B, i_type LDB, V* ALPHAR, V* ALPHAI, V* BETA, 
      V& CSL, V& SNL, V& CSR, V& SNR );

//-----------------------------------------------------------------------
//                          DSTEVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSTEVD computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric tridiagonal matrix. If eigenvectors are desired, it
*  uses a divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix
*          A.
*          On exit, if INFO = 0, the eigenvalues in ascending order.
*
*  E       (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix A, stored in elements 1 to N-1 of E.
*          On exit, the contents of E are destroyed.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with D(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*                                         dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If JOBZ  = 'N' or N <= 1 then LWORK must be at least 1.
*          If JOBZ  = 'V' and N > 1 then LWORK must be at least
*                         ( 1 + 4*N + N**2 ).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK and IWORK
*          arrays, returns these values as the first entries of the WORK
*          and IWORK arrays, and no error message related to LWORK or
*          LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If JOBZ  = 'N' or N <= 1 then LIWORK must be at least 1.
*          If JOBZ  = 'V' and N > 1 then LIWORK must be at least 3+5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK and IWORK arrays, and no error message related to
*          LWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of E did not converge to zero.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
stevd(const char * JOBZ, i_type N, V* D, V* E, V* Z, i_type LDZ, V* WORK, 
      i_type LWORK, i_type* IWORK, i_type LIWORK, i_type* INFO );

BLAS_EXPORT void sstevd(const char * JOBZ, i_type N, s_type* D, s_type* E, 
                        s_type* Z, i_type LDZ, s_type* WORK, i_type LWORK, 
                        i_type* IWORK, i_type LIWORK, i_type* INFO );

BLAS_EXPORT void dstevd(const char * JOBZ, i_type N, d_type* D, d_type* E,
                        d_type* Z, i_type LDZ, d_type* WORK, i_type LWORK,
                        i_type* IWORK, i_type LIWORK, i_type* INFO );

//-----------------------------------------------------------------------
//                          DLARFT
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLARFT forms the triangular factor T of a real block reflector H
*  of order n, which is defined as a product of k elementary reflectors.
*
*  If DIRECT = 'F', H = H(1) H(2) . . . H(k) and T is upper triangular;
*
*  If DIRECT = 'B', H = H(k) . . . H(2) H(1) and T is lower triangular.
*
*  If STOREV = 'C', the vector which defines the elementary reflector
*  H(i) is stored in the i-th column of the array V, and
*
*     H  =  I - V * T * V'
*
*  If STOREV = 'R', the vector which defines the elementary reflector
*  H(i) is stored in the i-th row of the array V, and
*
*     H  =  I - V' * T * V
*
*  Arguments
*  =========
*
*  DIRECT  (input) CHARACTER*1
*          Specifies the order in which the elementary reflectors are
*          multiplied to form the block reflector:
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Specifies how the vectors which define the elementary
*          reflectors are stored (see also Further Details):
*          = 'C': columnwise
*          = 'R': rowwise
*
*  N       (input) INTEGER
*          The order of the block reflector H. N >= 0.
*
*  K       (input) INTEGER
*          The order of the triangular factor T (= the number of
*          elementary reflectors). K >= 1.
*
*  V       (input/output) DOUBLE PRECISION array, dimension
*                               (LDV,K) if STOREV = 'C'
*                               (LDV,N) if STOREV = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C', LDV >= max(1,N); if STOREV = 'R', LDV >= K.
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i).
*
*  T       (output) DOUBLE PRECISION array, dimension (LDT,K)
*          The k by k triangular factor T of the block reflector.
*          If DIRECT = 'F', T is upper triangular; if DIRECT = 'B', T is
*          lower triangular. The rest of the array is not used.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  Further Details
*  ===============
*
*  The shape of the matrix V and the storage of the vectors which define
*  the H(i) is best illustrated by the following example with n = 5 and
*  k = 3. The elements equal to 1 are not stored; the corresponding
*  array elements are modified but restored on exit. The rest of the
*  array is not used.
*
*  DIRECT = 'F' and STOREV = 'C':         DIRECT = 'F' and STOREV = 'R':
*
*               V = (  1       )                 V = (  1 v1 v1 v1 v1 )
*                   ( v1  1    )                     (     1 v2 v2 v2 )
*                   ( v1 v2  1 )                     (        1 v3 v3 )
*                   ( v1 v2 v3 )
*                   ( v1 v2 v3 )
*
*  DIRECT = 'B' and STOREV = 'C':         DIRECT = 'B' and STOREV = 'R':
*
*               V = ( v1 v2 v3 )                 V = ( v1 v1  1       )
*                   ( v1 v2 v3 )                     ( v2 v2 v2  1    )
*                   (  1 v2 v3 )                     ( v3 v3 v3 v3  1 )
*                   (     1 v3 )
*                   (        1 )
*
*  =====================================================================
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
larft(const char *direct, const char *storev, i_type n, i_type k, const V *v,
      i_type ldv, const V* tau, V* t, i_type ldt);

BLAS_EXPORT void clarft(const char *direct, const char *storev, i_type n, 
                        i_type k, const c_type *v, i_type ldv, const c_type* tau,
                        c_type* t, i_type ldt);

BLAS_EXPORT void dlarft(const char *direct, const char *storev, i_type n, 
                        i_type k, const d_type *v, i_type ldv, const d_type* tau,
                        d_type* t, i_type ldt);

BLAS_EXPORT void slarft(const char *direct, const char *storev, i_type n, 
                        i_type k, const s_type *v, i_type ldv, const s_type *tau,
                        s_type *t, i_type ldt);

BLAS_EXPORT void zlarft(const char *direct, const char *storev, i_type n, 
                        i_type k, const z_type *v, i_type ldv, const z_type *tau,
                        z_type * t, i_type ldt);

//-----------------------------------------------------------------------
//                          DLARFB
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLARFB applies a real block reflector H or its transpose H' to a
*  real m by n matrix C, from either the left or the right.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply H or H' from the Left
*          = 'R': apply H or H' from the Right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': apply H (No transpose)
*          = 'T': apply H' (Transpose)
*
*  DIRECT  (input) CHARACTER*1
*          Indicates how H is formed from a product of elementary
*          reflectors
*          = 'F': H = H(1) H(2) . . . H(k) (Forward)
*          = 'B': H = H(k) . . . H(2) H(1) (Backward)
*
*  STOREV  (input) CHARACTER*1
*          Indicates how the vectors which define the elementary
*          reflectors are stored:
*          = 'C': Columnwise
*          = 'R': Rowwise
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  K       (input) INTEGER
*          The order of the matrix T (= the number of elementary
*          reflectors whose product defines the block reflector).
*
*  V       (input) DOUBLE PRECISION array, dimension
*                                (LDV,K) if STOREV = 'C'
*                                (LDV,M) if STOREV = 'R' and SIDE = 'L'
*                                (LDV,N) if STOREV = 'R' and SIDE = 'R'
*          The matrix V. See further details.
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.
*          If STOREV = 'C' and SIDE = 'L', LDV >= max(1,M);
*          if STOREV = 'C' and SIDE = 'R', LDV >= max(1,N);
*          if STOREV = 'R', LDV >= K.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,K)
*          The triangular k by k matrix T in the representation of the
*          block reflector.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= K.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by H*C or H'*C or C*H or C*H'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDA >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,K)
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          If SIDE = 'L', LDWORK >= max(1,N);
*          if SIDE = 'R', LDWORK >= max(1,M).
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
larfb(const char *side, const char *trans, const char *direct, const char * storev, 
    i_type m, i_type n, i_type k, const V* v, i_type ldv,  const V* t, i_type ldt, 
    V* c__, i_type ldc, V* work, i_type ldwork);

BLAS_EXPORT void clarfb(const char *side, const char *trans, const char *direct, 
                        const char * storev, i_type m, i_type n, i_type k,
                        const c_type* v, i_type ldv, const c_type* t, i_type ldt, 
                        c_type* c__, i_type ldc, c_type* work, i_type ldwork);

BLAS_EXPORT void dlarfb(const char *side, const char *trans, const char *direct,
                        const char * storev, i_type m, i_type n, i_type k, 
                        const d_type* v, i_type ldv, const d_type* t, i_type ldt, 
                        d_type* c__, i_type ldc, d_type* work, i_type ldwork);

BLAS_EXPORT void slarfb(const char *side, const char *trans, const char *direct,
                        const char * storev, i_type m, i_type n, i_type k, 
                        const s_type* v, i_type ldv, const s_type* t, i_type ldt, 
                        s_type* c__, i_type ldc, s_type* work, i_type ldwork);

BLAS_EXPORT void zlarfb(const char *side, const char *trans, const char *direct, 
                        const char * storev, i_type m, i_type n, i_type k, 
                        const z_type* v, i_type  ldv, const z_type* t, i_type ldt, 
                        z_type* c__, i_type ldc, z_type* work, i_type ldwork);

//-----------------------------------------------------------------------
//                          DGEQR2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGEQR2 computes a QR factorization of a real m by n matrix A:
*  A = Q * R.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the m by n matrix A.
*          On exit, the elements on and above the diagonal of the array
*          contain the min(m,n) by n upper trapezoidal matrix R (R is
*          upper triangular if m >= n); the elements below the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(1) H(2) . . . H(k), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
*  and tau in TAU(i).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
geqr2(i_type m, i_type n, V* a, i_type lda, V* tau, V* work, i_type& info);

BLAS_EXPORT void cgeqr2(i_type m, i_type n, c_type* a, i_type lda, 
                        c_type* tau, c_type* work, i_type& info);

BLAS_EXPORT void dgeqr2(i_type m, i_type n, d_type* a, i_type lda, d_type* tau,
                        d_type* work, i_type& info);

BLAS_EXPORT void sgeqr2(i_type m, i_type n, s_type* a, i_type lda, s_type* tau, 
                        s_type* work, i_type& info);

BLAS_EXPORT void zgeqr2(i_type m, i_type n, d_type* a, i_type lda, d_type* tau,
                        d_type* work, i_type& info);

//-----------------------------------------------------------------------
//                          DLARFG
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLARFG generates a real elementary reflector H of order n, such
*  that
*
*        H * ( alpha ) = ( beta ),   H' * H = I.
*            (   x   )   (   0  )
*
*  where alpha and beta are scalars, and x is an (n-1)-element real
*  vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a real scalar and v is a real (n-1)-element
*  vector.
*
*  If the elements of x are all zero, then tau = 0 and H is taken to be
*  the unit matrix.
*
*  Otherwise  1 <= tau <= 2.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) DOUBLE PRECISION
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) DOUBLE PRECISION array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX <> 0.
*
*  TAU     (output) DOUBLE PRECISION
*          The value tau.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
larfg(i_type n, V* alpha, V* x, i_type incx, V* tau);

BLAS_EXPORT void clarfg(i_type n, c_type* alpha, c_type* x, i_type incx, 
                        c_type* tau);

BLAS_EXPORT void dlarfg(i_type n, d_type* alpha, d_type* x, i_type incx,
                        d_type* tau);

BLAS_EXPORT void slarfg(i_type n, s_type* alpha, s_type* x, i_type incx,
                        s_type* tau);

BLAS_EXPORT void zlarfg(i_type n, d_type* alpha, d_type* x, i_type incx,
                        z_type* tau);

//-----------------------------------------------------------------------
//                          DLARF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLARF applies a real elementary reflector H to a real m by n matrix
*  C, from either the left or the right. H is represented in the form
*
*        H = I - tau * v * v'
*
*  where tau is a real scalar and v is a real vector.
*
*  If tau = 0, then H is taken to be the unit matrix.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': form  H * C
*          = 'R': form  C * H
*
*  M       (input) INTEGER
*          The number of rows of the matrix C.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C.
*
*  V       (input) DOUBLE PRECISION array, dimension
*                     (1 + (M-1)*abs(INCV)) if SIDE = 'L'
*                  or (1 + (N-1)*abs(INCV)) if SIDE = 'R'
*          The vector v in the representation of H. V is not used if
*          TAU = 0.
*
*  INCV    (input) INTEGER
*          The increment between elements of v. INCV <> 0.
*
*  TAU     (input) DOUBLE PRECISION
*          The value tau in the representation of H.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the m by n matrix C.
*          On exit, C is overwritten by the matrix H * C if SIDE = 'L',
*          or C * H if SIDE = 'R'.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension
*                         (N) if SIDE = 'L'
*                      or (M) if SIDE = 'R'
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
larf(const char *side, i_type m, i_type n, V* v, i_type incv, const V* tau,
     V* c__, i_type ldc, V* work);

BLAS_EXPORT void clarf(const char *side, i_type m, i_type n, c_type* v, 
                       i_type incv, const c_type* tau, c_type* c__, i_type ldc,
                       c_type* work);

BLAS_EXPORT void dlarf(const char *side, i_type m, i_type n, d_type* v, 
	                   i_type incv, const d_type* tau, d_type* c__, i_type ldc, 
                       d_type* work);

BLAS_EXPORT void slarf(const char *side, i_type m, i_type n, s_type* v, 
	                   i_type incv, const s_type* tau, s_type* c__, i_type ldc, 
                       s_type* work);

BLAS_EXPORT void zlarf(const char *side, i_type m, i_type n, d_type* v, 
                       i_type incv, const d_type* tau, d_type* c__, i_type ldc, 
                       d_type* work);

//-----------------------------------------------------------------------
//                          DGEBRD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGEBRD reduces a general real M-by-N matrix A to upper or lower
*  bidiagonal form B by an orthogonal transformation: Q**T * A * P = B.
*
*  If m >= n, B is upper bidiagonal; if m < n, B is lower bidiagonal.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows in the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns in the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N general matrix to be reduced.
*          On exit,
*          if m >= n, the diagonal and the first superdiagonal are
*            overwritten with the upper bidiagonal matrix B; the
*            elements below the diagonal, with the array TAUQ, represent
*            the orthogonal matrix Q as a product of elementary
*            reflectors, and the elements above the first superdiagonal,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors;
*          if m < n, the diagonal and the first subdiagonal are
*            overwritten with the lower bidiagonal matrix B; the
*            elements below the first subdiagonal, with the array TAUQ,
*            represent the orthogonal matrix Q as a product of
*            elementary reflectors, and the elements above the diagonal,
*            with the array TAUP, represent the orthogonal matrix P as
*            a product of elementary reflectors.
*          See Further Details.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  D       (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The diagonal elements of the bidiagonal matrix B:
*          D(i) = A(i,i).
*
*  E       (output) DOUBLE PRECISION array, dimension (min(M,N)-1)
*          The off-diagonal elements of the bidiagonal matrix B:
*          if m >= n, E(i) = A(i,i+1) for i = 1,2,...,n-1;
*          if m < n, E(i) = A(i+1,i) for i = 1,2,...,m-1.
*
*  TAUQ    (output) DOUBLE PRECISION array dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix Q. See Further Details.
*
*  TAUP    (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors which
*          represent the orthogonal matrix P. See Further Details.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.  LWORK >= max(1,M,N).
*          For optimum performance LWORK >= (M+N)*NB, where NB
*          is the optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  The matrices Q and P are represented as products of elementary
*  reflectors:
*
*  If m >= n,
*
*     Q = H(1) H(2) . . . H(n)  and  P = G(1) G(2) . . . G(n-1)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i-1) = 0, v(i) = 1, and v(i+1:m) is stored on exit in A(i+1:m,i);
*  u(1:i) = 0, u(i+1) = 1, and u(i+2:n) is stored on exit in A(i,i+2:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  If m < n,
*
*     Q = H(1) H(2) . . . H(m-1)  and  P = G(1) G(2) . . . G(m)
*
*  Each H(i) and G(i) has the form:
*
*     H(i) = I - tauq * v * v'  and G(i) = I - taup * u * u'
*
*  where tauq and taup are real scalars, and v and u are real vectors;
*  v(1:i) = 0, v(i+1) = 1, and v(i+2:m) is stored on exit in A(i+2:m,i);
*  u(1:i-1) = 0, u(i) = 1, and u(i+1:n) is stored on exit in A(i,i+1:n);
*  tauq is stored in TAUQ(i) and taup in TAUP(i).
*
*  The contents of A on exit are illustrated by the following examples:
*
*  m = 6 and n = 5 (m > n):          m = 5 and n = 6 (m < n):
*
*    (  d   e   u1  u1  u1 )           (  d   u1  u1  u1  u1  u1 )
*    (  v1  d   e   u2  u2 )           (  e   d   u2  u2  u2  u2 )
*    (  v1  v2  d   e   u3 )           (  v1  e   d   u3  u3  u3 )
*    (  v1  v2  v3  d   e  )           (  v1  v2  e   d   u4  u4 )
*    (  v1  v2  v3  v4  d  )           (  v1  v2  v3  e   d   u5 )
*    (  v1  v2  v3  v4  v5 )
*
*  where d and e denote diagonal and off-diagonal elements of B, vi
*  denotes an element of the vector defining H(i), and ui an element of
*  the vector defining G(i).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gebrd(i_type M, i_type N, V* A, i_type LDA, typename details::real_type<V>::type* D, 
      typename details::real_type<V>::type* E, V* tauq, V* taup, V* WORK, 
      i_type lwork, i_type& info);

BLAS_EXPORT void 
dgebrd(i_type M, i_type N, d_type* A, i_type LDA, d_type* D, d_type* E, d_type* tauq, 
       d_type* taup, d_type* WORK, i_type lwork, i_type& info);

BLAS_EXPORT void 
sgebrd(i_type M, i_type N, s_type* A, i_type LDA, s_type* D, s_type* E, s_type* tauq, 
       s_type* taup, s_type* WORK, i_type lwork, i_type& info);

BLAS_EXPORT void 
cgebrd(i_type M, i_type N, c_type* A, i_type LDA, s_type* D, s_type* E, c_type* tauq, 
       c_type* taup, c_type* WORK, i_type lwork, i_type& info);

BLAS_EXPORT void 
zgebrd(i_type M, i_type N, z_type* A, i_type LDA, d_type* D, d_type* E, z_type* tauq, 
       z_type* taup, z_type* WORK, i_type lwork, i_type& info);

//-----------------------------------------------------------------------
//                          DGGHRD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGGHRD reduces a pair of real matrices (A,B) to generalized upper
*  Hessenberg form using orthogonal transformations, where A is a
*  general matrix and B is upper triangular.  The form of the
*  generalized eigenvalue problem is
*     A*x = lambda*B*x,
*  and B is typically made upper triangular by computing its QR
*  factorization and moving the orthogonal matrix Q to the left side
*  of the equation.
*
*  This subroutine simultaneously reduces A to a Hessenberg matrix H:
*     Q**T*A*Z = H
*  and transforms B to another upper triangular matrix T:
*     Q**T*B*Z = T
*  in order to reduce the problem to its standard form
*     H*y = lambda*T*y
*  where y = Z**T*x.
*
*  The orthogonal matrices Q and Z are determined as products of Givens
*  rotations.  They may either be formed explicitly, or they may be
*  postmultiplied into input matrices Q1 and Z1, so that
*
*       Q1 * A * Z1**T = (Q1*Q) * H * (Z1*Z)**T
*
*       Q1 * B * Z1**T = (Q1*Q) * T * (Z1*Z)**T
*
*  If Q1 is the orthogonal matrix from the QR factorization of B in the
*  original equation A*x = lambda*B*x, then DGGHRD reduces the original
*  problem to generalized Hessenberg form.
*
*  Arguments
*  =========
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not compute Q;
*          = 'I': Q is initialized to the unit matrix, and the
*                 orthogonal matrix Q is returned;
*          = 'V': Q must contain an orthogonal matrix Q1 on entry,
*                 and the product Q1*Q is returned.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not compute Z;
*          = 'I': Z is initialized to the unit matrix, and the
*                 orthogonal matrix Z is returned;
*          = 'V': Z must contain an orthogonal matrix Z1 on entry,
*                 and the product Z1*Z is returned.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          ILO and IHI mark the rows and columns of A which are to be
*          reduced.  It is assumed that A is already upper triangular
*          in rows and columns 1:ILO-1 and IHI+1:N.  ILO and IHI are
*          normally set by a previous call to SGGBAL; otherwise they
*          should be set to 1 and N respectively.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the N-by-N general matrix to be reduced.
*          On exit, the upper triangle and the first subdiagonal of A
*          are overwritten with the upper Hessenberg matrix H, and the
*          rest is set to zero.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB, N)
*          On entry, the N-by-N upper triangular matrix B.
*          On exit, the upper triangular matrix T = Q**T B Z.  The
*          elements below the diagonal are set to zero.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*          On entry, if COMPQ = 'V', the orthogonal matrix Q1,
*          typically from the QR factorization of B.
*          On exit, if COMPQ='I', the orthogonal matrix Q, and if
*          COMPQ = 'V', the product Q1*Q.
*          Not referenced if COMPQ='N'.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= N if COMPQ='V' or 'I'; LDQ >= 1 otherwise.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
*          On entry, if COMPZ = 'V', the orthogonal matrix Z1.
*          On exit, if COMPZ='I', the orthogonal matrix Z, and if
*          COMPZ = 'V', the product Z1*Z.
*          Not referenced if COMPZ='N'.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.
*          LDZ >= N if COMPZ='V' or 'I'; LDZ >= 1 otherwise.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  This routine reduces A to Hessenberg and B to triangular form by
*  an unblocked reduction, as described in _Matrix_Computations_,
*  by Golub and Van Loan (Johns Hopkins Press.)
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gghrd(const char* COMPQ, const char* COMPZ, i_type N, i_type ILO, i_type IHI, V* A, 
      i_type LDA, V* B, i_type LDB, V* Q, i_type LDQ, V* Z, i_type LDZ, i_type& INFO );

//-----------------------------------------------------------------------
//                          DSBTRD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSBTRD reduces a real symmetric band matrix A to symmetric
*  tridiagonal form T by an orthogonal similarity transformation:
*  Q**T * A * Q = T.
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'N':  do not form Q;
*          = 'V':  form Q;
*          = 'U':  update a matrix X, by forming X*Q.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB,N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*          On exit, the diagonal elements of AB are overwritten by the
*          diagonal elements of the tridiagonal matrix T; if KD > 0, the
*          elements on the first superdiagonal (if UPLO = 'U') or the
*          first subdiagonal (if UPLO = 'L') are overwritten by the
*          off-diagonal elements of T; the rest of AB is overwritten by
*          values generated during the reduction.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD+1.
*
*  D       (output) DOUBLE PRECISION array, dimension (N)
*          The diagonal elements of the tridiagonal matrix T.
*
*  E       (output) DOUBLE PRECISION array, dimension (N-1)
*          The off-diagonal elements of the tridiagonal matrix T:
*          E(i) = T(i,i+1) if UPLO = 'U'; E(i) = T(i+1,i) if UPLO = 'L'.
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if VECT = 'U', then Q must contain an N-by-N
*          matrix X; if VECT = 'N' or 'V', then Q need not be set.
*
*          On exit:
*          if VECT = 'V', Q contains the N-by-N orthogonal matrix Q;
*          if VECT = 'U', Q contains the product X*Q;
*          if VECT = 'N', the array Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1, and LDQ >= N if VECT = 'V' or 'U'.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
sbtrd(const char *VECT, const char *UPLO, i_type N, i_type KD, V *AB, i_type LDAB, 
      typename details::real_type<V>::type* D, typename details::real_type<V>::type* E, 
      V* Q, i_type LDQ, V* work, i_type& info);

//-----------------------------------------------------------------------
//                          DSYEVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSYEVD computes all eigenvalues and, optionally, eigenvectors of a
*  real symmetric matrix A. If eigenvectors are desired, it uses a
*  divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Because of large use of BLAS of level 3, DSYEVD needs N**2 more
*  workspace than DSYEVX.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA, N)
*          On entry, the symmetric matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*                                         dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N <= 1,               LWORK must be at least 1.
*          If JOBZ = 'N' and N > 1, LWORK must be at least 2*N+1.
*          If JOBZ = 'V' and N > 1, LWORK must be at least
*                                                1 + 6*N + 2*N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK and IWORK
*          arrays, returns these values as the first entries of the WORK
*          and IWORK arrays, and no error message related to LWORK or
*          LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If N <= 1,                LIWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
*          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK and IWORK arrays, and no error message related to
*          LWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
*                to converge; i off-diagonal elements of an intermediate
*                tridiagonal form did not converge to zero;
*                if INFO = i and JOBZ = 'V', then the algorithm failed
*                to compute an eigenvalue while working on the submatrix
*                lying in rows and columns INFO/(N+1) through
*                mod(INFO,N+1).
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
syevd(const char *jobv,const char *uplo, i_type n, V *a,i_type lda, 
      typename details::real_type<V>::type* w, V *work,i_type lwork, i_type* iwork,
      i_type liwork, i_type *info);

//-----------------------------------------------------------------------
//                          ZHEEVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  CHEEVD computes all eigenvalues and, optionally, eigenvectors of a
*  complex Hermitian matrix A.  If eigenvectors are desired, it uses a
*  divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) COMPLEX array, dimension (LDA, N)
*          On entry, the Hermitian matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of A contains the
*          upper triangular part of the matrix A.  If UPLO = 'L',
*          the leading N-by-N lower triangular part of A contains
*          the lower triangular part of the matrix A.
*          On exit, if JOBZ = 'V', then if INFO = 0, A contains the
*          orthonormal eigenvectors of the matrix A.
*          If JOBZ = 'N', then on exit the lower triangle (if UPLO='L')
*          or the upper triangle (if UPLO='U') of A, including the
*          diagonal, is destroyed.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  W       (output) REAL array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  WORK    (workspace/output) COMPLEX array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The length of the array WORK.
*          If N <= 1,                LWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LWORK must be at least N + 1.
*          If JOBZ  = 'V' and N > 1, LWORK must be at least 2*N + N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK, RWORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) REAL array,
*                                         dimension (LRWORK)
*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*
*  LRWORK  (input) INTEGER
*          The dimension of the array RWORK.
*          If N <= 1,                LRWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LRWORK must be at least N.
*          If JOBZ  = 'V' and N > 1, LRWORK must be at least
*                         1 + 5*N + 2*N**2.
*
*          If LRWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If N <= 1,                LIWORK must be at least 1.
*          If JOBZ  = 'N' and N > 1, LIWORK must be at least 1.
*          If JOBZ  = 'V' and N > 1, LIWORK must be at least 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i and JOBZ = 'N', then the algorithm failed
*                to converge; i off-diagonal elements of an intermediate
*                tridiagonal form did not converge to zero;
*                if INFO = i and JOBZ = 'V', then the algorithm failed
*                to compute an eigenvalue while working on the submatrix
*                lying in rows and columns INFO/(N+1) through
*                mod(INFO,N+1).
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
heevd(const char *jobv,const char *uplo, i_type n, V *a,i_type lda, 
      typename details::real_type<V>::type* w, V *work,i_type lwork, 
      typename details::real_type<V>::type* rwork, i_type lrwork, 
      i_type* iwork, i_type liwork, i_type *info);

//-----------------------------------------------------------------------
//                          DSBEV
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSBEV computes all the eigenvalues and, optionally, eigenvectors of
*  a real symmetric band matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, AB is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the first
*          superdiagonal and the diagonal of the tridiagonal matrix T
*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
*          the diagonal and first subdiagonal of T are returned in the
*          first two rows of AB.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD + 1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (max(1,3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
sbev(const char *jobv,const char *uplo, i_type n, i_type LD, V *a,i_type lda, 
     typename details::real_type<V>::type* w, V* Z, i_type LDZ, V *work, 
     i_type *info);

//-----------------------------------------------------------------------
//                          ZHBEV
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHBEV computes all the eigenvalues and, optionally, eigenvectors of
*  a complex Hermitian band matrix A.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX*16 array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the Hermitian band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, AB is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the first
*          superdiagonal and the diagonal of the tridiagonal matrix T
*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
*          the diagonal and first subdiagonal of T are returned in the
*          first two rows of AB.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD + 1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (max(1,3*N-2))
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
hbev(const char *jobv,const char *uplo, i_type n, i_type LD, V *a,i_type lda, 
     typename details::real_type<V>::type* w, V* Z, i_type LDZ, V *work, 
     typename details::real_type<V>::type* RWORK, i_type *info);

//-----------------------------------------------------------------------
//                          DSBEVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSBEVD computes all the eigenvalues and, optionally, eigenvectors of
*  a real symmetric band matrix A. If eigenvectors are desired, it uses
*  a divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, AB is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the first
*          superdiagonal and the diagonal of the tridiagonal matrix T
*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
*          the diagonal and first subdiagonal of T are returned in the
*          first two rows of AB.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD + 1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*                                         dimension (LWORK)
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          IF N <= 1,                LWORK must be at least 1.
*          If JOBZ  = 'N' and N > 2, LWORK must be at least 2*N.
*          If JOBZ  = 'V' and N > 2, LWORK must be at least
*                         ( 1 + 5*N + 2*N**2 ).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK and IWORK
*          arrays, returns these values as the first entries of the WORK
*          and IWORK arrays, and no error message related to LWORK or
*          LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array LIWORK.
*          If JOBZ  = 'N' or N <= 1, LIWORK must be at least 1.
*          If JOBZ  = 'V' and N > 2, LIWORK must be at least 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK and IWORK arrays, and no error message related to
*          LWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
sbevd(const char *jobv,const char *uplo, i_type n, i_type LD, V *a,i_type lda, 
      typename details::real_type<V>::type* w, V* Z, i_type LDZ, V *work, 
      i_type lwork, i_type* iwork, i_type liwork, i_type *info);

//-----------------------------------------------------------------------
//                          ZHBEVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHBEVD computes all the eigenvalues and, optionally, eigenvectors of
*  a complex Hermitian band matrix A.  If eigenvectors are desired, it
*  uses a divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangle of A is stored;
*          = 'L':  Lower triangle of A is stored.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  KD      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KD >= 0.
*
*  AB      (input/output) COMPLEX*16 array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the Hermitian band
*          matrix A, stored in the first KD+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(kd+1+i-j,j) = A(i,j) for max(1,j-kd)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+kd).
*
*          On exit, AB is overwritten by values generated during the
*          reduction to tridiagonal form.  If UPLO = 'U', the first
*          superdiagonal and the diagonal of the tridiagonal matrix T
*          are returned in rows KD and KD+1 of AB, and if UPLO = 'L',
*          the diagonal and first subdiagonal of T are returned in the
*          first two rows of AB.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KD + 1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) COMPLEX*16 array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the orthonormal
*          eigenvectors of the matrix A, with the i-th column of Z
*          holding the eigenvector associated with W(i).
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N <= 1,               LWORK must be at least 1.
*          If JOBZ = 'N' and N > 1, LWORK must be at least N.
*          If JOBZ = 'V' and N > 1, LWORK must be at least 2*N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK, RWORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  RWORK   (workspace/output) DOUBLE PRECISION array,
*                                         dimension (LRWORK)
*          On exit, if INFO = 0, RWORK(1) returns the optimal LRWORK.
*
*  LRWORK  (input) INTEGER
*          The dimension of array RWORK.
*          If N <= 1,               LRWORK must be at least 1.
*          If JOBZ = 'N' and N > 1, LRWORK must be at least N.
*          If JOBZ = 'V' and N > 1, LRWORK must be at least
*                        1 + 5*N + 2*N**2.
*
*          If LRWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of array IWORK.
*          If JOBZ = 'N' or N <= 1, LIWORK must be at least 1.
*          If JOBZ = 'V' and N > 1, LIWORK must be at least 3 + 5*N .
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK, RWORK
*          and IWORK arrays, returns these values as the first entries
*          of the WORK, RWORK and IWORK arrays, and no error message
*          related to LWORK or LRWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = i, the algorithm failed to converge; i
*                off-diagonal elements of an intermediate tridiagonal
*                form did not converge to zero.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
hbevd(const char *jobv,const char *uplo, i_type n, i_type LD, V *a,i_type lda, 
      typename details::real_type<V>::type* w, V* Z, i_type LDZ, V *work,
      i_type lwork, typename details::real_type<V>::type *rwork, i_type lrwork,
      i_type* iwork, i_type liwork, i_type *info);

//-----------------------------------------------------------------------
//                          DTREVC
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTREVC computes some or all of the right and/or left eigenvectors of
*  a real upper quasi-triangular matrix T.
*  Matrices of this type are produced by the Schur factorization of
*  a real general matrix:  A = Q*T*Q**T, as computed by DHSEQR.
*  
*  The right eigenvector x and the left eigenvector y of T corresponding
*  to an eigenvalue w are defined by:
*  
*     T*x = w*x,     (y**H)*T = w*(y**H)
*  
*  where y**H denotes the conjugate transpose of y.
*  The eigenvalues are not input to this routine, but are read directly
*  from the diagonal blocks of T.
*  
*  This routine returns the matrices X and/or Y of right and left
*  eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
*  input matrix.  If Q is the orthogonal factor that reduces a matrix
*  A to Schur form T, then Q*X and Q*Y are the matrices of right and
*  left eigenvectors of A.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'R':  compute right eigenvectors only;
*          = 'L':  compute left eigenvectors only;
*          = 'B':  compute both right and left eigenvectors.
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A':  compute all right and/or left eigenvectors;
*          = 'B':  compute all right and/or left eigenvectors,
*                  backtransformed by the matrices in VR and/or VL;
*          = 'S':  compute selected right and/or left eigenvectors,
*                  as indicated by the logical array SELECT.
*
*  SELECT  (input/output) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
*          computed.
*          If w(j) is a real eigenvalue, the corresponding real
*          eigenvector is computed if SELECT(j) is .TRUE..
*          If w(j) and w(j+1) are the real and imaginary parts of a
*          complex eigenvalue, the corresponding complex eigenvector is
*          computed if either SELECT(j) or SELECT(j+1) is .TRUE., and
*          on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is set to
*          .FALSE..
*          Not referenced if HOWMNY = 'A' or 'B'.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
*          The upper quasi-triangular matrix T in Schur canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*          of Schur vectors returned by DHSEQR).
*          On exit, if SIDE = 'L' or 'B', VL contains:
*          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*Y;
*          if HOWMNY = 'S', the left eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VL, in the same order as their
*                           eigenvalues.
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part, and the second the imaginary part.
*          Not referenced if SIDE = 'R'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.  LDVL >= 1, and if
*          SIDE = 'L' or 'B', LDVL >= N.
*
*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*          of Schur vectors returned by DHSEQR).
*          On exit, if SIDE = 'R' or 'B', VR contains:
*          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*          if HOWMNY = 'B', the matrix Q*X;
*          if HOWMNY = 'S', the right eigenvectors of T specified by
*                           SELECT, stored consecutively in the columns
*                           of VR, in the same order as their
*                           eigenvalues.
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part and the second the imaginary part.
*          Not referenced if SIDE = 'L'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.  LDVR >= 1, and if
*          SIDE = 'R' or 'B', LDVR >= N.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR actually
*          used to store the eigenvectors.
*          If HOWMNY = 'A' or 'B', M is set to N.
*          Each selected real eigenvector occupies one column and each
*          selected complex eigenvector occupies two columns.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The algorithm used in this program is basically backward (forward)
*  substitution, with scaling to make the the code robust against
*  possible overflow.
*
*  Each eigenvector is normalized so that the element of largest
*  magnitude has magnitude 1; here the magnitude of a complex number
*  (x,y) is taken to be |x| + |y|.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trevc(const char* side, const char* howmny, i_type* select, i_type n, V * t,
      i_type ldt, V *vl, i_type ldvl, V* vr, i_type ldvr, i_type mm, i_type m,
      V* work, i_type* info);

//-----------------------------------------------------------------------
//                          DTRSNA
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTRSNA estimates reciprocal condition numbers for specified
*  eigenvalues and/or right eigenvectors of a real upper
*  quasi-triangular matrix T (or of any matrix Q*T*Q**T with Q
*  orthogonal).
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elements equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for
*          eigenvalues (S) or eigenvectors (SEP):
*          = 'E': for eigenvalues only (S);
*          = 'V': for eigenvectors only (SEP);
*          = 'B': for both eigenvalues and eigenvectors (S and SEP).
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A': compute condition numbers for all eigenpairs;
*          = 'S': compute condition numbers for selected eigenpairs
*                 specified by the array SELECT.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
*          condition numbers are required. To select condition numbers
*          for the eigenpair corresponding to a real eigenvalue w(j),
*          SELECT(j) must be set to .TRUE.. To select condition numbers
*          corresponding to a complex conjugate pair of eigenvalues w(j)
*          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
*          set to .TRUE..
*          If HOWMNY = 'A', SELECT is not referenced.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
*          The upper quasi-triangular matrix T, in Schur canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  VL      (input) DOUBLE PRECISION array, dimension (LDVL,M)
*          If JOB = 'E' or 'B', VL must contain left eigenvectors of T
*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
*          must be stored in consecutive columns of VL, as returned by
*          DHSEIN or DTREVC.
*          If JOB = 'V', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.
*          LDVL >= 1; and if JOB = 'E' or 'B', LDVL >= N.
*
*  VR      (input) DOUBLE PRECISION array, dimension (LDVR,M)
*          If JOB = 'E' or 'B', VR must contain right eigenvectors of T
*          (or of any Q*T*Q**T with Q orthogonal), corresponding to the
*          eigenpairs specified by HOWMNY and SELECT. The eigenvectors
*          must be stored in consecutive columns of VR, as returned by
*          DHSEIN or DTREVC.
*          If JOB = 'V', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.
*          LDVR >= 1; and if JOB = 'E' or 'B', LDVR >= N.
*
*  S       (output) DOUBLE PRECISION array, dimension (MM)
*          If JOB = 'E' or 'B', the reciprocal condition numbers of the
*          selected eigenvalues, stored in consecutive elements of the
*          array. For a complex conjugate pair of eigenvalues two
*          consecutive elements of S are set to the same value. Thus
*          S(j), SEP(j), and the j-th columns of VL and VR all
*          correspond to the same eigenpair (but not in general the
*          j-th eigenpair, unless all eigenpairs are selected).
*          If JOB = 'V', S is not referenced.
*
*  SEP     (output) DOUBLE PRECISION array, dimension (MM)
*          If JOB = 'V' or 'B', the estimated reciprocal condition
*          numbers of the selected eigenvectors, stored in consecutive
*          elements of the array. For a complex eigenvector two
*          consecutive elements of SEP are set to the same value. If
*          the eigenvalues cannot be reordered to compute SEP(j), SEP(j)
*          is set to 0; this can only occur when the true value would be
*          very small anyway.
*          If JOB = 'E', SEP is not referenced.
*
*  MM      (input) INTEGER
*          The number of elements in the arrays S (if JOB = 'E' or 'B')
*           and/or SEP (if JOB = 'V' or 'B'). MM >= M.
*
*  M       (output) INTEGER
*          The number of elements of the arrays S and/or SEP actually
*          used to store the estimated condition numbers.
*          If HOWMNY = 'A', M is set to N.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LDWORK,N+6)
*          If JOB = 'E', WORK is not referenced.
*
*  LDWORK  (input) INTEGER
*          The leading dimension of the array WORK.
*          LDWORK >= 1; and if JOB = 'V' or 'B', LDWORK >= N.
*
*  IWORK   (workspace) INTEGER array, dimension (2*(N-1))
*          If JOB = 'E' or type V is complex, IWORK is not referenced.
*
*  RWORK   (workspace) INTEGER array, dimension (N)
*          If JOB = 'E' or V is real, RWORK is not referenced.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The reciprocal of the condition number of an eigenvalue lambda is
*  defined as
*
*          S(lambda) = |v'*u| / (norm(u)*norm(v))
*
*  where u and v are the right and left eigenvectors of T corresponding
*  to lambda; v' denotes the conjugate-transpose of v, and norm(u)
*  denotes the Euclidean norm. These reciprocal condition numbers always
*  lie between zero (very badly conditioned) and one (very well
*  conditioned). If n = 1, S(lambda) is defined to be 1.
*
*  An approximate error bound for a computed eigenvalue W(i) is given by
*
*                      EPS * norm(T) / S(i)
*
*  where EPS is the machine precision.
*
*  The reciprocal of the condition number of the right eigenvector u
*  corresponding to lambda is defined as follows. Suppose
*
*              T = ( lambda  c  )
*                  (   0    T22 )
*
*  Then the reciprocal condition number is
*
*          SEP( lambda, T22 ) = sigma-min( T22 - lambda*I )
*
*  where sigma-min denotes the smallest singular value. We approximate
*  the smallest singular value by the reciprocal of an estimate of the
*  one-norm of the inverse of T22 - lambda*I. If n = 1, SEP(1) is
*  defined to be abs(T(1,1)).
*
*  An approximate error bound for a computed right eigenvector VR(i)
*  is given by
*
*                      EPS * norm(T) / SEP(i)
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trsna(const char *job, const char *howmny, const i_type *select, i_type n,
      const V* t, i_type ldt, const V* vl, i_type ldvl, const V* vr, i_type ldvr,
      typename details::real_type<V>::type* s, 
      typename details::real_type<V>::type * sep, i_type mm, i_type m, V* work, 
      i_type ldwork, i_type* iwork, typename details::real_type<V>::type* rwork,
      i_type* info);

//-----------------------------------------------------------------------
//                          DLACN2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLACN2 estimates the 1-norm of a square, real matrix A.
*  Reverse communication is used for evaluating matrix-vector products.
*
*  Arguments
*  =========
*
*  N      (input) INTEGER
*         The order of the matrix.  N >= 1.
*
*  V      (workspace) DOUBLE PRECISION array, dimension (N)
*         On the final return, V = A*W,  where  EST = norm(V)/norm(W)
*         (W is not returned).
*
*  X      (input/output) DOUBLE PRECISION array, dimension (N)
*         On an intermediate return, X should be overwritten by
*               A * X,   if KASE=1,
*               A' * X,  if KASE=2,
*         and DLACN2 must be re-called with all the other parameters
*         unchanged.
*
*  ISGN   (workspace) INTEGER array, dimension (N)
*
*  EST    (input/output) DOUBLE PRECISION
*         On entry with KASE = 1 or 2 and ISAVE(1) = 3, EST should be
*         unchanged from the previous call to DLACN2.
*         On exit, EST is an estimate (a lower bound) for norm(A). 
*
*  KASE   (input/output) INTEGER
*         On the initial call to DLACN2, KASE should be 0.
*         On an intermediate return, KASE will be 1 or 2, indicating
*         whether X should be overwritten by A * X  or A' * X.
*         On the final return from DLACN2, KASE will again be 0.
*
*  ISAVE  (input/output) INTEGER array, dimension (3)
*         ISAVE is used to save variables between calls to DLACN2
*
*  Further Details
*  ======= =======
*
*  Contributed by Nick Higham, University of Manchester.
*  Originally named SONEST, dated March 16, 1988.
*
*  Reference: N.J. Higham, "FORTRAN codes for estimating the one-norm of
*  a real or complex matrix, with applications to condition estimation",
*  ACM Trans. Math. Soft., vol. 14, no. 4, pp. 381-396, December 1988.
*
*  This is a thread safe version of DLACON, which uses the array ISAVE
*  in place of a SAVE statement, as follows:
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
lacn2(i_type n, V* v, V* x, i_type* isgn, typename details::real_type<V>::type* est, 
      i_type* kase, i_type* isave);

//-----------------------------------------------------------------------
//                          HSEIN
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DHSEIN uses inverse iteration to find specified right and/or left
*  eigenvectors of a real upper Hessenberg matrix H.
*
*  The right eigenvector x and the left eigenvector y of the matrix H
*  corresponding to an eigenvalue w are defined by:
*
*               H * x = w * x,     y**h * H = w * y**h
*
*  where y**h denotes the conjugate transpose of the vector y.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'R': compute right eigenvectors only;
*          = 'L': compute left eigenvectors only;
*          = 'B': compute both right and left eigenvectors.
*
*  EIGSRC  (input) CHARACTER*1
*          Specifies the source of eigenvalues supplied in (WR,WI):
*          = 'Q': the eigenvalues were found using DHSEQR; thus, if
*                 H has zero subdiagonal elements, and so is
*                 block-triangular, then the j-th eigenvalue can be
*                 assumed to be an eigenvalue of the block containing
*                 the j-th row/column.  This property allows DHSEIN to
*                 perform inverse iteration on just one diagonal block.
*          = 'N': no assumptions are made on the correspondence
*                 between eigenvalues and diagonal blocks.  In this
*                 case, DHSEIN must always perform inverse iteration
*                 using the whole matrix H.
*
*  INITV   (input) CHARACTER*1
*          = 'N': no initial vectors are supplied;
*          = 'U': user-supplied initial vectors are stored in the arrays
*                 VL and/or VR.
*
*  SELECT  (input/output) LOGICAL array, dimension (N)
*          Specifies the eigenvectors to be computed. To select the
*          real eigenvector corresponding to a real eigenvalue WR(j),
*          SELECT(j) must be set to .TRUE.. To select the complex
*          eigenvector corresponding to a complex eigenvalue
*          (WR(j),WI(j)), with complex conjugate (WR(j+1),WI(j+1)),
*          either SELECT(j) or SELECT(j+1) or both must be set to
*          .TRUE.; then on exit SELECT(j) is .TRUE. and SELECT(j+1) is
*          .FALSE..
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  H       (input) DOUBLE PRECISION array, dimension (LDH,N)
*          The upper Hessenberg matrix H.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max(1,N).
*
*  WR      (input/output) DOUBLE PRECISION array, dimension (N)
*  WI      (input) DOUBLE PRECISION array, dimension (N)
*          On entry, the real and imaginary parts of the eigenvalues of
*          H; a complex conjugate pair of eigenvalues must be stored in
*          consecutive elements of WR and WI.
*          On exit, WR may have been altered since close eigenvalues
*          are perturbed slightly in searching for independent
*          eigenvectors.
*
*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
*          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must
*          contain starting vectors for the inverse iteration for the
*          left eigenvectors; the starting vector for each eigenvector
*          must be in the same column(s) in which the eigenvector will
*          be stored.
*          On exit, if SIDE = 'L' or 'B', the left eigenvectors
*          specified by SELECT will be stored consecutively in the
*          columns of VL, in the same order as their eigenvalues. A
*          complex eigenvector corresponding to a complex eigenvalue is
*          stored in two consecutive columns, the first holding the real
*          part and the second the imaginary part.
*          If SIDE = 'R', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.
*          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
*
*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
*          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must
*          contain starting vectors for the inverse iteration for the
*          right eigenvectors; the starting vector for each eigenvector
*          must be in the same column(s) in which the eigenvector will
*          be stored.
*          On exit, if SIDE = 'R' or 'B', the right eigenvectors
*          specified by SELECT will be stored consecutively in the
*          columns of VR, in the same order as their eigenvalues. A
*          complex eigenvector corresponding to a complex eigenvalue is
*          stored in two consecutive columns, the first holding the real
*          part and the second the imaginary part.
*          If SIDE = 'L', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.
*          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR required to
*          store the eigenvectors; each selected real eigenvector
*          occupies one column and each selected complex eigenvector
*          occupies two columns.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension ((N+2)*N)
*
*  IFAILL  (output) INTEGER array, dimension (MM)
*          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left
*          eigenvector in the i-th column of VL (corresponding to the
*          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the
*          eigenvector converged satisfactorily. If the i-th and (i+1)th
*          columns of VL hold a complex eigenvector, then IFAILL(i) and
*          IFAILL(i+1) are set to the same value.
*          If SIDE = 'R', IFAILL is not referenced.
*
*  IFAILR  (output) INTEGER array, dimension (MM)
*          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right
*          eigenvector in the i-th column of VR (corresponding to the
*          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the
*          eigenvector converged satisfactorily. If the i-th and (i+1)th
*          columns of VR hold a complex eigenvector, then IFAILR(i) and
*          IFAILR(i+1) are set to the same value.
*          If SIDE = 'L', IFAILR is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, i is the number of eigenvectors which
*                failed to converge; see IFAILL and IFAILR for further
*                details.
*
*  Further Details
*  ===============
*
*  Each eigenvector is normalized so that the element of largest
*  magnitude has magnitude 1; here the magnitude of a complex number
*  (x,y) is taken to be |x|+|y|.
*
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
hsein(const char *side, const char *eigsrc, const char *initv, i_type* select,
      i_type n, const V* h, i_type ldh, typename details::real_type<V>::type* wr,
      typename details::real_type<V>::type* wi, V* vl, i_type ldvl, V* vr, 
      i_type ldvr, i_type mm, i_type m, V* work, i_type* ifaill, i_type* ifailr,
      i_type* info);

//-----------------------------------------------------------------------
//                          ZHSEIN
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHSEIN uses inverse iteration to find specified right and/or left
*  eigenvectors of a complex upper Hessenberg matrix H.
*
*  The right eigenvector x and the left eigenvector y of the matrix H
*  corresponding to an eigenvalue w are defined by:
*
*               H * x = w * x,     y**h * H = w * y**h
*
*  where y**h denotes the conjugate transpose of the vector y.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'R': compute right eigenvectors only;
*          = 'L': compute left eigenvectors only;
*          = 'B': compute both right and left eigenvectors.
*
*  EIGSRC  (input) CHARACTER*1
*          Specifies the source of eigenvalues supplied in W:
*          = 'Q': the eigenvalues were found using ZHSEQR; thus, if
*                 H has zero subdiagonal elements, and so is
*                 block-triangular, then the j-th eigenvalue can be
*                 assumed to be an eigenvalue of the block containing
*                 the j-th row/column.  This property allows ZHSEIN to
*                 perform inverse iteration on just one diagonal block.
*          = 'N': no assumptions are made on the correspondence
*                 between eigenvalues and diagonal blocks.  In this
*                 case, ZHSEIN must always perform inverse iteration
*                 using the whole matrix H.
*
*  INITV   (input) CHARACTER*1
*          = 'N': no initial vectors are supplied;
*          = 'U': user-supplied initial vectors are stored in the arrays
*                 VL and/or VR.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          Specifies the eigenvectors to be computed. To select the
*          eigenvector corresponding to the eigenvalue W(j),
*          SELECT(j) must be set to .TRUE..
*
*  N       (input) INTEGER
*          The order of the matrix H.  N >= 0.
*
*  H       (input) COMPLEX*16 array, dimension (LDH,N)
*          The upper Hessenberg matrix H.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max(1,N).
*
*  W       (input/output) COMPLEX*16 array, dimension (N)
*          On entry, the eigenvalues of H.
*          On exit, the real parts of W may have been altered since
*          close eigenvalues are perturbed slightly in searching for
*          independent eigenvectors.
*
*  VL      (input/output) COMPLEX*16 array, dimension (LDVL,MM)
*          On entry, if INITV = 'U' and SIDE = 'L' or 'B', VL must
*          contain starting vectors for the inverse iteration for the
*          left eigenvectors; the starting vector for each eigenvector
*          must be in the same column in which the eigenvector will be
*          stored.
*          On exit, if SIDE = 'L' or 'B', the left eigenvectors
*          specified by SELECT will be stored consecutively in the
*          columns of VL, in the same order as their eigenvalues.
*          If SIDE = 'R', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL.
*          LDVL >= max(1,N) if SIDE = 'L' or 'B'; LDVL >= 1 otherwise.
*
*  VR      (input/output) COMPLEX*16 array, dimension (LDVR,MM)
*          On entry, if INITV = 'U' and SIDE = 'R' or 'B', VR must
*          contain starting vectors for the inverse iteration for the
*          right eigenvectors; the starting vector for each eigenvector
*          must be in the same column in which the eigenvector will be
*          stored.
*          On exit, if SIDE = 'R' or 'B', the right eigenvectors
*          specified by SELECT will be stored consecutively in the
*          columns of VR, in the same order as their eigenvalues.
*          If SIDE = 'L', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.
*          LDVR >= max(1,N) if SIDE = 'R' or 'B'; LDVR >= 1 otherwise.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR required to
*          store the eigenvectors (= the number of .TRUE. elements in
*          SELECT).
*
*  WORK    (workspace) COMPLEX*16 array, dimension (N*N)
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  IFAILL  (output) INTEGER array, dimension (MM)
*          If SIDE = 'L' or 'B', IFAILL(i) = j > 0 if the left
*          eigenvector in the i-th column of VL (corresponding to the
*          eigenvalue w(j)) failed to converge; IFAILL(i) = 0 if the
*          eigenvector converged satisfactorily.
*          If SIDE = 'R', IFAILL is not referenced.
*
*  IFAILR  (output) INTEGER array, dimension (MM)
*          If SIDE = 'R' or 'B', IFAILR(i) = j > 0 if the right
*          eigenvector in the i-th column of VR (corresponding to the
*          eigenvalue w(j)) failed to converge; IFAILR(i) = 0 if the
*          eigenvector converged satisfactorily.
*          If SIDE = 'L', IFAILR is not referenced.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, i is the number of eigenvectors which
*                failed to converge; see IFAILL and IFAILR for further
*                details.
*
*  Further Details
*  ===============
*
*  Each eigenvector is normalized so that the element of largest
*  magnitude has magnitude 1; here the magnitude of a complex number
*  (x,y) is taken to be |x|+|y|.
*
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid_complex<void,V>::type
hsein(const char *side, const char *eigsrc, const char *initv, i_type* select,
      i_type n, const V* h, i_type ldh, V* w, V* vl, i_type ldvl, V* vr, 
      i_type ldvr, i_type mm, i_type m, V* work, 
      typename details::real_type<V>::type* rwork, i_type* ifaill, 
      i_type* ifailr, i_type* info);

//-----------------------------------------------------------------------
//                          DSTEIN
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSTEIN computes the eigenvectors of a real symmetric tridiagonal
*  matrix T corresponding to specified eigenvalues, using inverse
*  iteration.
*
*  The maximum number of iterations allowed for each eigenvector is
*  specified by an internal parameter MAXITS (currently set to 5).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix.  N >= 0.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal elements of the tridiagonal matrix
*          T, in elements 1 to N-1.
*
*  M       (input) INTEGER
*          The number of eigenvectors to be found.  0 <= M <= N.
*
*  W       (input) DOUBLE PRECISION array, dimension (N)
*          The first M elements of W contain the eigenvalues for
*          which eigenvectors are to be computed.  The eigenvalues
*          should be grouped by split-off block and ordered from
*          smallest to largest within the block.  ( The output array
*          W from DSTEBZ with ORDER = 'B' is expected here. )
*
*  IBLOCK  (input) INTEGER array, dimension (N)
*          The submatrix indices associated with the corresponding
*          eigenvalues in W; IBLOCK(i)=1 if eigenvalue W(i) belongs to
*          the first submatrix from the top, =2 if W(i) belongs to
*          the second submatrix, etc.  ( The output array IBLOCK
*          from DSTEBZ is expected here. )
*
*  ISPLIT  (input) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to
*          ISPLIT( 1 ), the second of rows/columns ISPLIT( 1 )+1
*          through ISPLIT( 2 ), etc.
*          ( The output array ISPLIT from DSTEBZ is expected here. )
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, M)
*          The computed eigenvectors.  The eigenvector associated
*          with the eigenvalue W(i) is stored in the i-th column of
*          Z.  Any vector which fails to converge is set to its current
*          iterate after MAXITS iterations.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= max(1,N).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (5*N)
*
*  IWORK   (workspace) INTEGER array, dimension (N)
*
*  IFAIL   (output) INTEGER array, dimension (M)
*          On normal exit, all elements of IFAIL are zero.
*          If one or more eigenvectors fail to converge after
*          MAXITS iterations, then their indices are stored in
*          array IFAIL.
*
*  INFO    (output) INTEGER
*          = 0: successful exit.
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, then i eigenvectors failed to converge
*               in MAXITS iterations.  Their indices are stored in
*               array IFAIL.
*
*  Internal Parameters
*  ===================
*
*  MAXITS  INTEGER, default = 5
*          The maximum number of iterations performed.
*
*  EXTRA   INTEGER, default = 2
*          The number of iterations performed after norm growth
*          criterion is satisfied, should be at least 1.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
stein(i_type n, const typename details::real_type<V>::type* d, 
      const typename details::real_type<V>::type* e, i_type m, 
      const typename details::real_type<V>::type* w, const i_type* iblock, 
      const i_type* isplit, V* z, i_type ldz, 
      typename details::real_type<V>::type * work, i_type* iwork, 
      i_type* ifail, i_type* info);

//-----------------------------------------------------------------------
//                          DTGEVC
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTGEVC computes some or all of the right and/or left eigenvectors of
*  a pair of real matrices (S,P), where S is a quasi-triangular matrix
*  and P is upper triangular.  Matrix pairs of this type are produced by
*  the generalized Schur factorization of a matrix pair (A,B):
*
*     A = Q*S*Z**T,  B = Q*P*Z**T
*
*  as computed by DGGHRD + DHGEQZ.
*
*  The right eigenvector x and the left eigenvector y of (S,P)
*  corresponding to an eigenvalue w are defined by:
*  
*     S*x = w*P*x,  (y**H)*S = w*(y**H)*P,
*  
*  where y**H denotes the conjugate tranpose of y.
*  The eigenvalues are not input to this routine, but are computed
*  directly from the diagonal blocks of S and P.
*  
*  This routine returns the matrices X and/or Y of right and left
*  eigenvectors of (S,P), or the products Z*X and/or Q*Y,
*  where Z and Q are input matrices.
*  If Q and Z are the orthogonal factors from the generalized Schur
*  factorization of a matrix pair (A,B), then Z*X and Q*Y
*  are the matrices of right and left eigenvectors of (A,B).
* 
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'R': compute right eigenvectors only;
*          = 'L': compute left eigenvectors only;
*          = 'B': compute both right and left eigenvectors.
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A': compute all right and/or left eigenvectors;
*          = 'B': compute all right and/or left eigenvectors,
*                 backtransformed by the matrices in VR and/or VL;
*          = 'S': compute selected right and/or left eigenvectors,
*                 specified by the logical array SELECT.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If HOWMNY='S', SELECT specifies the eigenvectors to be
*          computed.  If w(j) is a real eigenvalue, the corresponding
*          real eigenvector is computed if SELECT(j) is .TRUE..
*          If w(j) and w(j+1) are the real and imaginary parts of a
*          complex eigenvalue, the corresponding complex eigenvector
*          is computed if either SELECT(j) or SELECT(j+1) is .TRUE.,
*          and on exit SELECT(j) is set to .TRUE. and SELECT(j+1) is
*          set to .FALSE..
*          Not referenced if HOWMNY = 'A' or 'B'.
*
*  N       (input) INTEGER
*          The order of the matrices S and P.  N >= 0.
*
*  S       (input) DOUBLE PRECISION array, dimension (LDS,N)
*          The upper quasi-triangular matrix S from a generalized Schur
*          factorization, as computed by DHGEQZ.
*
*  LDS     (input) INTEGER
*          The leading dimension of array S.  LDS >= max(1,N).
*
*  P       (input) DOUBLE PRECISION array, dimension (LDP,N)
*          The upper triangular matrix P from a generalized Schur
*          factorization, as computed by DHGEQZ.
*          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks
*          of S must be in positive diagonal form.
*
*  LDP     (input) INTEGER
*          The leading dimension of array P.  LDP >= max(1,N).
*
*  VL      (input/output) DOUBLE PRECISION array, dimension (LDVL,MM)
*          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*          contain an N-by-N matrix Q (usually the orthogonal matrix Q
*          of left Schur vectors returned by DHGEQZ).
*          On exit, if SIDE = 'L' or 'B', VL contains:
*          if HOWMNY = 'A', the matrix Y of left eigenvectors of (S,P);
*          if HOWMNY = 'B', the matrix Q*Y;
*          if HOWMNY = 'S', the left eigenvectors of (S,P) specified by
*                      SELECT, stored consecutively in the columns of
*                      VL, in the same order as their eigenvalues.
*
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part, and the second the imaginary part.
*
*          Not referenced if SIDE = 'R'.
*
*  LDVL    (input) INTEGER
*          The leading dimension of array VL.  LDVL >= 1, and if
*          SIDE = 'L' or 'B', LDVL >= N.
*
*  VR      (input/output) DOUBLE PRECISION array, dimension (LDVR,MM)
*          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*          contain an N-by-N matrix Z (usually the orthogonal matrix Z
*          of right Schur vectors returned by DHGEQZ).
*
*          On exit, if SIDE = 'R' or 'B', VR contains:
*          if HOWMNY = 'A', the matrix X of right eigenvectors of (S,P);
*          if HOWMNY = 'B' or 'b', the matrix Z*X;
*          if HOWMNY = 'S' or 's', the right eigenvectors of (S,P)
*                      specified by SELECT, stored consecutively in the
*                      columns of VR, in the same order as their
*                      eigenvalues.
*
*          A complex eigenvector corresponding to a complex eigenvalue
*          is stored in two consecutive columns, the first holding the
*          real part and the second the imaginary part.
*          
*          Not referenced if SIDE = 'L'.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR.  LDVR >= 1, and if
*          SIDE = 'R' or 'B', LDVR >= N.
*
*  MM      (input) INTEGER
*          The number of columns in the arrays VL and/or VR. MM >= M.
*
*  M       (output) INTEGER
*          The number of columns in the arrays VL and/or VR actually
*          used to store the eigenvectors.  If HOWMNY = 'A' or 'B', M
*          is set to N.  Each selected real eigenvector occupies one
*          column and each selected complex eigenvector occupies two
*          columns.
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (6*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  the 2-by-2 block (INFO:INFO+1) does not have a complex
*                eigenvalue.
*
*  Further Details
*  ===============
*
*  Allocation of workspace:
*  ---------- -- ---------
*
*     WORK( j ) = 1-norm of j-th column of A, above the diagonal
*     WORK( N+j ) = 1-norm of j-th column of B, above the diagonal
*     WORK( 2*N+1:3*N ) = real part of eigenvector
*     WORK( 3*N+1:4*N ) = imaginary part of eigenvector
*     WORK( 4*N+1:5*N ) = real part of back-transformed eigenvector
*     WORK( 5*N+1:6*N ) = imaginary part of back-transformed eigenvector
*
*  Rowwise vs. columnwise solution methods:
*  ------- --  ---------- -------- -------
*
*  Finding a generalized eigenvector consists basically of solving the
*  singular triangular system
*
*   (A - w B) x = 0     (for right) or:   (A - w B)**H y = 0  (for left)
*
*  Consider finding the i-th right eigenvector (assume all eigenvalues
*  are real). The equation to be solved is:
*       n                   i
*  0 = sum  C(j,k) v(k)  = sum  C(j,k) v(k)     for j = i,. . .,1
*      k=j                 k=j
*
*  where  C = (A - w B)  (The components v(i+1:n) are 0.)
*
*  The "rowwise" method is:
*
*  (1)  v(i) := 1
*  for j = i-1,. . .,1:
*                          i
*      (2) compute  s = - sum C(j,k) v(k)   and
*                        k=j+1
*
*      (3) v(j) := s / C(j,j)
*
*  Step 2 is sometimes called the "dot product" step, since it is an
*  inner product between the j-th row and the portion of the eigenvector
*  that has been computed so far.
*
*  The "columnwise" method consists basically in doing the sums
*  for all the rows in parallel.  As each v(j) is computed, the
*  contribution of v(j) times the j-th column of C is added to the
*  partial sums.  Since FORTRAN arrays are stored columnwise, this has
*  the advantage that at each step, the elements of C that are accessed
*  are adjacent to one another, whereas with the rowwise method, the
*  elements accessed at a step are spaced LDS (and LDP) words apart.
*
*  When finding left eigenvectors, the matrix in question is the
*  transpose of the one in storage, so the rowwise method then
*  actually accesses columns of A and B at each step, and so is the
*  preferred method.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
tgevc(const char *side, const char *howmny, i_type *select, i_type n, 
      const V* s, i_type lds, const V* p, i_type ldp, V* vl, i_type ldvl,
      V* vr, i_type ldvr, i_type mm, i_type m, V* work, i_type* info);

//-----------------------------------------------------------------------
//                          DTGSNA
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTGSNA estimates reciprocal condition numbers for specified
*  eigenvalues and/or eigenvectors of a matrix pair (A, B) in
*  generalized real Schur canonical form (or of any matrix pair
*  (Q*A*Z', Q*B*Z') with orthogonal matrices Q and Z, where
*  Z' denotes the transpose of Z.
*
*  (A, B) must be in generalized real Schur form (as returned by DGGES),
*  i.e. A is block upper triangular with 1-by-1 and 2-by-2 diagonal
*  blocks. B is upper triangular.
*
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for
*          eigenvalues (S) or eigenvectors (DIF):
*          = 'E': for eigenvalues only (S);
*          = 'V': for eigenvectors only (DIF);
*          = 'B': for both eigenvalues and eigenvectors (S and DIF).
*
*  HOWMNY  (input) CHARACTER*1
*          = 'A': compute condition numbers for all eigenpairs;
*          = 'S': compute condition numbers for selected eigenpairs
*                 specified by the array SELECT.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          If HOWMNY = 'S', SELECT specifies the eigenpairs for which
*          condition numbers are required. To select condition numbers
*          for the eigenpair corresponding to a real eigenvalue w(j),
*          SELECT(j) must be set to .TRUE.. To select condition numbers
*          corresponding to a complex conjugate pair of eigenvalues w(j)
*          and w(j+1), either SELECT(j) or SELECT(j+1) or both, must be
*          set to .TRUE..
*          If HOWMNY = 'A', SELECT is not referenced.
*
*  N       (input) INTEGER
*          The order of the square matrix pair (A, B). N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The upper quasi-triangular matrix A in the pair (A,B).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB,N)
*          The upper triangular matrix B in the pair (A,B).
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  VL      (input) DOUBLE PRECISION array, dimension (LDVL,M)
*          If JOB = 'E' or 'B', VL must contain left eigenvectors of
*          (A, B), corresponding to the eigenpairs specified by HOWMNY
*          and SELECT. The eigenvectors must be stored in consecutive
*          columns of VL, as returned by DTGEVC.
*          If JOB = 'V', VL is not referenced.
*
*  LDVL    (input) INTEGER
*          The leading dimension of the array VL. LDVL >= 1.
*          If JOB = 'E' or 'B', LDVL >= N.
*
*  VR      (input) DOUBLE PRECISION array, dimension (LDVR,M)
*          If JOB = 'E' or 'B', VR must contain right eigenvectors of
*          (A, B), corresponding to the eigenpairs specified by HOWMNY
*          and SELECT. The eigenvectors must be stored in consecutive
*          columns ov VR, as returned by DTGEVC.
*          If JOB = 'V', VR is not referenced.
*
*  LDVR    (input) INTEGER
*          The leading dimension of the array VR. LDVR >= 1.
*          If JOB = 'E' or 'B', LDVR >= N.
*
*  S       (output) DOUBLE PRECISION array, dimension (MM)
*          If JOB = 'E' or 'B', the reciprocal condition numbers of the
*          selected eigenvalues, stored in consecutive elements of the
*          array. For a complex conjugate pair of eigenvalues two
*          consecutive elements of S are set to the same value. Thus
*          S(j), DIF(j), and the j-th columns of VL and VR all
*          correspond to the same eigenpair (but not in general the
*          j-th eigenpair, unless all eigenpairs are selected).
*          If JOB = 'V', S is not referenced.
*
*  DIF     (output) DOUBLE PRECISION array, dimension (MM)
*          If JOB = 'V' or 'B', the estimated reciprocal condition
*          numbers of the selected eigenvectors, stored in consecutive
*          elements of the array. For a complex eigenvector two
*          consecutive elements of DIF are set to the same value. If
*          the eigenvalues cannot be reordered to compute DIF(j), DIF(j)
*          is set to 0; this can only occur when the true value would be
*          very small anyway.
*          If JOB = 'E', DIF is not referenced.
*
*  MM      (input) INTEGER
*          The number of elements in the arrays S and DIF. MM >= M.
*
*  M       (output) INTEGER
*          The number of elements of the arrays S and DIF used to store
*          the specified condition numbers; for each selected real
*          eigenvalue one element is used, and for each selected complex
*          conjugate pair of eigenvalues, two elements are used.
*          If HOWMNY = 'A', M is set to N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK. LWORK >= max(1,N).
*          If JOB = 'V' or 'B' LWORK >= 2*N*(N+2)+16.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (N + 6) or (N + 2) if V is 
*          complex.
*          If JOB = 'E', IWORK is not referenced.
*
*  INFO    (output) INTEGER
*          =0: Successful exit
*          <0: If INFO = -i, the i-th argument had an illegal value
*
*
*  Further Details
*  ===============
*
*  The reciprocal of the condition number of a generalized eigenvalue
*  w = (a, b) is defined as
*
*       S(w) = (|u'Av|**2 + |u'Bv|**2)**(1/2) / (norm(u)*norm(v))
*
*  where u and v are the left and right eigenvectors of (A, B)
*  corresponding to w; |z| denotes the absolute value of the complex
*  number, and norm(u) denotes the 2-norm of the vector u.
*  The pair (a, b) corresponds to an eigenvalue w = a/b (= u'Av/u'Bv)
*  of the matrix pair (A, B). If both a and b equal zero, then (A B) is
*  singular and S(I) = -1 is returned.
*
*  An approximate error bound on the chordal distance between the i-th
*  computed generalized eigenvalue w and the corresponding exact
*  eigenvalue lambda is
*
*       chord(w, lambda) <= EPS * norm(A, B) / S(I)
*
*  where EPS is the machine precision.
*
*  The reciprocal of the condition number DIF(i) of right eigenvector u
*  and left eigenvector v corresponding to the generalized eigenvalue w
*  is defined as follows:
*
*  a) If the i-th eigenvalue w = (a,b) is real
*
*     Suppose U and V are orthogonal transformations such that
*
*                U'*(A, B)*V  = (S, T) = ( a   * ) ( b  * )  1
*                                        ( 0  S22 ),( 0 T22 )  n-1
*                                          1  n-1     1 n-1
*
*     Then the reciprocal condition number DIF(i) is
*
*                Difl((a, b), (S22, T22)) = sigma-min( Zl ),
*
*     where sigma-min(Zl) denotes the smallest singular value of the
*     2(n-1)-by-2(n-1) matrix
*
*         Zl = [ kron(a, In-1)  -kron(1, S22) ]
*              [ kron(b, In-1)  -kron(1, T22) ] .
*
*     Here In-1 is the identity matrix of size n-1. kron(X, Y) is the
*     Kronecker product between the matrices X and Y.
*
*     Note that if the default method for computing DIF(i) is wanted
*     (see DLATDF), then the parameter DIFDRI (see below) should be
*     changed from 3 to 4 (routine DLATDF(IJOB = 2 will be used)).
*     See DTGSYL for more details.
*
*  b) If the i-th and (i+1)-th eigenvalues are complex conjugate pair,
*
*     Suppose U and V are orthogonal transformations such that
*
*                U'*(A, B)*V = (S, T) = ( S11  *  ) ( T11  * )  2
*                                       ( 0    S22 ),( 0    T22) n-2
*                                         2    n-2     2    n-2
*
*     and (S11, T11) corresponds to the complex conjugate eigenvalue
*     pair (w, conjg(w)). There exist unitary matrices U1 and V1 such
*     that
*
*         U1'*S11*V1 = ( s11 s12 )   and U1'*T11*V1 = ( t11 t12 )
*                      (  0  s22 )                    (  0  t22 )
*
*     where the generalized eigenvalues w = s11/t11 and
*     conjg(w) = s22/t22.
*
*     Then the reciprocal condition number DIF(i) is bounded by
*
*         min( d1, max( 1, |real(s11)/real(s22)| )*d2 )
*
*     where, d1 = Difl((s11, t11), (s22, t22)) = sigma-min(Z1), where
*     Z1 is the complex 2-by-2 matrix
*
*              Z1 =  [ s11  -s22 ]
*                    [ t11  -t22 ],
*
*     This is done by computing (using real arithmetic) the
*     roots of the characteristical polynomial det(Z1' * Z1 - lambda I),
*     where Z1' denotes the conjugate transpose of Z1 and det(X) denotes
*     the determinant of X.
*
*     and d2 is an upper bound on Difl((S11, T11), (S22, T22)), i.e. an
*     upper bound on sigma-min(Z2), where Z2 is (2n-2)-by-(2n-2)
*
*              Z2 = [ kron(S11', In-2)  -kron(I2, S22) ]
*                   [ kron(T11', In-2)  -kron(I2, T22) ]
*
*     Note that if the default method for computing DIF is wanted (see
*     DLATDF), then the parameter DIFDRI (see below) should be changed
*     from 3 to 4 (routine DLATDF(IJOB = 2 will be used)). See DTGSYL
*     for more details.
*
*  For each eigenvalue/vector specified by SELECT, DIF stores a
*  Frobenius norm-based estimate of Difl.
*
*  An approximate error bound for the i-th computed eigenvector VL(i) or
*  VR(i) is given by
*
*             EPS * norm(A, B) / DIF(i).
*
*  See ref. [2-3] for more details and further references.
*
*  Based on contributions by
*     Bo Kagstrom and Peter Poromaa, Department of Computing Science,
*     Umea University, S-901 87 Umea, Sweden.
*
*  References
*  ==========
*
*  [1] B. Kagstrom; A Direct Method for Reordering Eigenvalues in the
*      Generalized Real Schur Form of a Regular Matrix Pair (A, B), in
*      M.S. Moonen et al (eds), Linear Algebra for Large Scale and
*      Real-Time Applications, Kluwer Academic Publ. 1993, pp 195-218.
*
*  [2] B. Kagstrom and P. Poromaa; Computing Eigenspaces with Specified
*      Eigenvalues of a Regular Matrix Pair (A, B) and Condition
*      Estimation: Theory, Algorithms and Software,
*      Report UMINF - 94.04, Department of Computing Science, Umea
*      University, S-901 87 Umea, Sweden, 1994. Also as LAPACK Working
*      Note 87. To appear in Numerical Algorithms, 1996.
*
*  [3] B. Kagstrom and P. Poromaa, LAPACK-Style Algorithms and Software
*      for Solving the Generalized Sylvester Equation and Estimating the
*      Separation between Regular Matrix Pairs, Report UMINF - 93.23,
*      Department of Computing Science, Umea University, S-901 87 Umea,
*      Sweden, December 1993, Revised April 1994, Also as LAPACK Working
*      Note 75.  To appear in ACM Trans. on Math. Software, Vol 22,
*      No 1, 1996.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
tgsna(const char *job, const char *howmny, const i_type *select, i_type n, 
      const V* a, i_type lda, const V* b, i_type ldb, V* vl, i_type ldvl, 
      V* vr, i_type ldvr, typename details::real_type<V>::type* s, 
      typename details::real_type<V>::type* dif, i_type mm, i_type m, 
      V* work, i_type lwork, i_type* iwork, i_type *info);

//-----------------------------------------------------------------------
//                          DLASSQ
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLASSQ  returns the values  scl  and  smsq  such that
*
*     ( scl**2 )*smsq = x( 1 )**2 +...+ x( n )**2 + ( scale**2 )*sumsq,
*
*  where  x( i ) = X( 1 + ( i - 1 )*INCX ). The value of  sumsq  is
*  assumed to be non-negative and  scl  returns the value
*
*     scl = max( scale, abs( x( i ) ) ).
*
*  scale and sumsq must be supplied in SCALE and SUMSQ and
*  scl and smsq are overwritten on SCALE and SUMSQ respectively.
*
*  The routine makes only one pass through the vector x.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of elements to be used from the vector X.
*
*  X       (input) DOUBLE PRECISION array, dimension (N)
*          The vector for which a scaled sum of squares is computed.
*             x( i )  = X( 1 + ( i - 1 )*INCX ), 1 <= i <= n.
*
*  INCX    (input) INTEGER
*          The increment between successive values of the vector X.
*          INCX > 0.
*
*  SCALE   (input/output) DOUBLE PRECISION
*          On entry, the value  scale  in the equation above.
*          On exit, SCALE is overwritten with  scl , the scaling factor
*          for the sum of squares.
*
*  SUMSQ   (input/output) DOUBLE PRECISION
*          On entry, the value  sumsq  in the equation above.
*          On exit, SUMSQ is overwritten with  smsq , the basic sum of
*          squares from which  scl  has been factored out.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
lassq(i_type n, const V* x, i_type incx, typename details::real_type<V>::type& scale, 
      typename details::real_type<V>::type& sumsq);

//-----------------------------------------------------------------------
//                          DLABAD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLABAD takes as input the values computed by DLAMCH for underflow and
*  overflow, and returns the square root of each of these values if the
*  log of LARGE is sufficiently large.  This subroutine is intended to
*  identify machines with a large exponent range, such as the Crays, and
*  redefine the underflow and overflow limits to be the square roots of
*  the values computed by DLAMCH.  This subroutine is needed because
*  DLAMCH does not compensate for poor arithmetic in the upper half of
*  the exponent range, as is found on a Cray.
*
*  Arguments
*  =========
*
*  SMALL   (input/output) DOUBLE PRECISION
*          On entry, the underflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of SMALL, otherwise unchanged.
*
*  LARGE   (input/output) DOUBLE PRECISION
*          On entry, the overflow threshold as computed by DLAMCH.
*          On exit, if LOG10(LARGE) is sufficiently large, the square
*          root of LARGE, otherwise unchanged.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
labad(V& small, V& large);

//-----------------------------------------------------------------------
//                          DLASCL
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLASCL multiplies the M by N real matrix A by the real scalar
*  CTO/CFROM.  This is done without over/underflow as long as the final
*  result CTO*A(I,J)/CFROM does not over/underflow. TYPE specifies that
*  A may be full, upper triangular, lower triangular, upper Hessenberg,
*  or banded.
*
*  Arguments
*  =========
*
*  TYPE    (input) CHARACTER*1
*          TYPE indices the storage type of the input matrix.
*          = 'G':  A is a full matrix.
*          = 'L':  A is a lower triangular matrix.
*          = 'U':  A is an upper triangular matrix.
*          = 'H':  A is an upper Hessenberg matrix.
*          = 'B':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the lower
*                  half stored.
*          = 'Q':  A is a symmetric band matrix with lower bandwidth KL
*                  and upper bandwidth KU and with the only the upper
*                  half stored.
*          = 'Z':  A is a band matrix with lower bandwidth KL and upper
*                  bandwidth KU.
*
*  KL      (input) INTEGER
*          The lower bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  KU      (input) INTEGER
*          The upper bandwidth of A.  Referenced only if TYPE = 'B',
*          'Q' or 'Z'.
*
*  CFROM   (input) DOUBLE PRECISION
*  CTO     (input) DOUBLE PRECISION
*          The matrix A is multiplied by CTO/CFROM. A(I,J) is computed
*          without over/underflow if the final result CTO*A(I,J)/CFROM
*          can be represented without over/underflow.  CFROM must be
*          nonzero.
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          The matrix to be multiplied by CTO/CFROM.  See TYPE for the
*          storage type.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  INFO    (output) INTEGER
*          0  - successful exit
*          <0 - if INFO = -i, the i-th argument had an illegal value.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
lascl(const char *type, i_type kl, i_type ku, 
      typename details::real_type<V>::type cfrom, 
      typename details::real_type<V>::type cto, i_type m, i_type n, V* a, 
      i_type lda, i_type& info);

//-----------------------------------------------------------------------
//                          DGGBAL
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGGBAL balances a pair of general real matrices (A,B).  This
*  involves, first, permuting A and B by similarity transformations to
*  isolate eigenvalues in the first 1 to ILO-1 and last IHI+1 to N
*  elements on the diagonal; and second, applying a diagonal similarity
*  transformation to rows and columns ILO to IHI to make the rows
*  and columns as close in norm as possible. Both steps are optional.
*
*  Balancing may reduce the 1-norm of the matrices, and improve the
*  accuracy of the computed eigenvalues and/or eigenvectors in the
*  generalized eigenvalue problem A*x = lambda*B*x.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies the operations to be performed on A and B:
*          = 'N':  none:  simply set ILO = 1, IHI = N, LSCALE(I) = 1.0
*                  and RSCALE(I) = 1.0 for i = 1,...,N.
*          = 'P':  permute only;
*          = 'S':  scale only;
*          = 'B':  both permute and scale.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the input matrix A.
*          On exit,  A is overwritten by the balanced matrix.
*          If JOB = 'N', A is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the input matrix B.
*          On exit,  B is overwritten by the balanced matrix.
*          If JOB = 'N', B is not referenced.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
*
*  ILO     (output) INTEGER
*  IHI     (output) INTEGER
*          ILO and IHI are set to integers such that on exit
*          A(i,j) = 0 and B(i,j) = 0 if i > j and
*          j = 1,...,ILO-1 or i = IHI+1,...,N.
*          If JOB = 'N' or 'S', ILO = 1 and IHI = N.
*
*  LSCALE  (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied
*          to the left side of A and B.  If P(j) is the index of the
*          row interchanged with row j, and D(j)
*          is the scaling factor applied to row j, then
*            LSCALE(j) = P(j)    for J = 1,...,ILO-1
*                      = D(j)    for J = ILO,...,IHI
*                      = P(j)    for J = IHI+1,...,N.
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  RSCALE  (output) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and scaling factors applied
*          to the right side of A and B.  If P(j) is the index of the
*          column interchanged with column j, and D(j)
*          is the scaling factor applied to column j, then
*            LSCALE(j) = P(j)    for J = 1,...,ILO-1
*                      = D(j)    for J = ILO,...,IHI
*                      = P(j)    for J = IHI+1,...,N.
*          The order in which the interchanges are made is N to IHI+1,
*          then 1 to ILO-1.
*
*  WORK    (workspace) REAL array, dimension (lwork)
*          lwork must be at least max(1,6*N) when JOB = 'S' or 'B', and
*          at least 1 when JOB = 'N' or 'P'.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  See R.C. WARD, Balancing the generalized eigenvalue problem,
*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ggbal(const char *job, i_type n, V* a, i_type lda, V* b, i_type ldb, 
      i_type& ilo, i_type& ihi, typename details::real_type<V>::type* lscale,
      typename details::real_type<V>::type* rscale, 
      typename details::real_type<V>::type* work, i_type& info);

//-----------------------------------------------------------------------
//                          DORMQR
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DORMQR overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by DGEQRF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,K)
*          The i-th column must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGEQRF in the first k columns of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If SIDE = 'L', LDA >= max(1,M);
*          if SIDE = 'R', LDA >= max(1,N).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGEQRF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ormqr(const char* side, const char* trans, i_type m, i_type n, i_type k, V* a, 
      i_type lda, const V* tau, V* c, i_type ldc, V* work, i_type lwork, 
      i_type& info);

//-----------------------------------------------------------------------
//                          DHGEQZ
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DHGEQZ computes the eigenvalues of a real matrix pair (H,T),
*  where H is an upper Hessenberg matrix and T is upper triangular,
*  using the double-shift QZ method.
*  Matrix pairs of this type are produced by the reduction to
*  generalized upper Hessenberg form of a real matrix pair (A,B):
*
*     A = Q1*H*Z1**T,  B = Q1*T*Z1**T,
*
*  as computed by DGGHRD.
*
*  If JOB='S', then the Hessenberg-triangular pair (H,T) is
*  also reduced to generalized Schur form,
*  
*     H = Q*S*Z**T,  T = Q*P*Z**T,
*  
*  where Q and Z are orthogonal matrices, P is an upper triangular
*  matrix, and S is a quasi-triangular matrix with 1-by-1 and 2-by-2
*  diagonal blocks.
*
*  The 1-by-1 blocks correspond to real eigenvalues of the matrix pair
*  (H,T) and the 2-by-2 blocks correspond to complex conjugate pairs of
*  eigenvalues.
*
*  Additionally, the 2-by-2 upper triangular diagonal blocks of P
*  corresponding to 2-by-2 blocks of S are reduced to positive diagonal
*  form, i.e., if S(j+1,j) is non-zero, then P(j+1,j) = P(j,j+1) = 0,
*  P(j,j) > 0, and P(j+1,j+1) > 0.
*
*  Optionally, the orthogonal matrix Q from the generalized Schur
*  factorization may be postmultiplied into an input matrix Q1, and the
*  orthogonal matrix Z may be postmultiplied into an input matrix Z1.
*  If Q1 and Z1 are the orthogonal matrices from DGGHRD that reduced
*  the matrix pair (A,B) to generalized upper Hessenberg form, then the
*  output matrices Q1*Q and Z1*Z are the orthogonal factors from the
*  generalized Schur factorization of (A,B):
*
*     A = (Q1*Q)*S*(Z1*Z)**T,  B = (Q1*Q)*P*(Z1*Z)**T.
*  
*  To avoid overflow, eigenvalues of the matrix pair (H,T) (equivalently,
*  of (A,B)) are computed as a pair of values (alpha,beta), where alpha is
*  complex and beta real.
*  If beta is nonzero, lambda = alpha / beta is an eigenvalue of the
*  generalized nonsymmetric eigenvalue problem (GNEP)
*     A*x = lambda*B*x
*  and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
*  alternate form of the GNEP
*     mu*A*y = B*y.
*  Real eigenvalues can be read directly from the generalized Schur
*  form: 
*    alpha = S(i,i), beta = P(i,i).
*
*  Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*       Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*       pp. 241--256.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'E': Compute eigenvalues only;
*          = 'S': Compute eigenvalues and the Schur form. 
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': Left Schur vectors (Q) are not computed;
*          = 'I': Q is initialized to the unit matrix and the matrix Q
*                 of left Schur vectors of (H,T) is returned;
*          = 'V': Q must contain an orthogonal matrix Q1 on entry and
*                 the product Q1*Q is returned.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': Right Schur vectors (Z) are not computed;
*          = 'I': Z is initialized to the unit matrix and the matrix Z
*                 of right Schur vectors of (H,T) is returned;
*          = 'V': Z must contain an orthogonal matrix Z1 on entry and
*                 the product Z1*Z is returned.
*
*  N       (input) INTEGER
*          The order of the matrices H, T, Q, and Z.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          ILO and IHI mark the rows and columns of H which are in
*          Hessenberg form.  It is assumed that A is already upper
*          triangular in rows and columns 1:ILO-1 and IHI+1:N.
*          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.
*
*  H       (input/output) DOUBLE PRECISION array, dimension (LDH, N)
*          On entry, the N-by-N upper Hessenberg matrix H.
*          On exit, if JOB = 'S', H contains the upper quasi-triangular
*          matrix S from the generalized Schur factorization;
*          2-by-2 diagonal blocks (corresponding to complex conjugate
*          pairs of eigenvalues) are returned in standard form, with
*          H(i,i) = H(i+1,i+1) and H(i+1,i)*H(i,i+1) < 0.
*          If JOB = 'E', the diagonal blocks of H match those of S, but
*          the rest of H is unspecified.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max( 1, N ).
*
*  T       (input/output) DOUBLE PRECISION array, dimension (LDT, N)
*          On entry, the N-by-N upper triangular matrix T.
*          On exit, if JOB = 'S', T contains the upper triangular
*          matrix P from the generalized Schur factorization;
*          2-by-2 diagonal blocks of P corresponding to 2-by-2 blocks of S
*          are reduced to positive diagonal form, i.e., if H(j+1,j) is
*          non-zero, then T(j+1,j) = T(j,j+1) = 0, T(j,j) > 0, and
*          T(j+1,j+1) > 0.
*          If JOB = 'E', the diagonal blocks of T match those of P, but
*          the rest of T is unspecified.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T.  LDT >= max( 1, N ).
*
*  ALPHAR  (output) DOUBLE PRECISION array, dimension (N)
*          The real parts of each scalar alpha defining an eigenvalue
*          of GNEP.
*
*  ALPHAI  (output) DOUBLE PRECISION array, dimension (N)
*          The imaginary parts of each scalar alpha defining an
*          eigenvalue of GNEP.
*          If ALPHAI(j) is zero, then the j-th eigenvalue is real; if
*          positive, then the j-th and (j+1)-st eigenvalues are a
*          complex conjugate pair, with ALPHAI(j+1) = -ALPHAI(j).
*
*  BETA    (output) DOUBLE PRECISION array, dimension (N)
*          The scalars beta that define the eigenvalues of GNEP.
*          Together, the quantities alpha = (ALPHAR(j),ALPHAI(j)) and
*          beta = BETA(j) represent the j-th eigenvalue of the matrix
*          pair (A,B), in one of the forms lambda = alpha/beta or
*          mu = beta/alpha.  Since either lambda or mu may overflow,
*          they should not, in general, be computed.
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ, N)
*          On entry, if COMPZ = 'V', the orthogonal matrix Q1 used in
*          the reduction of (A,B) to generalized Hessenberg form.
*          On exit, if COMPZ = 'I', the orthogonal matrix of left Schur
*          vectors of (H,T), and if COMPZ = 'V', the orthogonal matrix
*          of left Schur vectors of (A,B).
*          Not referenced if COMPZ = 'N'.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= 1.
*          If COMPQ='V' or 'I', then LDQ >= N.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ, N)
*          On entry, if COMPZ = 'V', the orthogonal matrix Z1 used in
*          the reduction of (A,B) to generalized Hessenberg form.
*          On exit, if COMPZ = 'I', the orthogonal matrix of
*          right Schur vectors of (H,T), and if COMPZ = 'V', the
*          orthogonal matrix of right Schur vectors of (A,B).
*          Not referenced if COMPZ = 'N'.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If COMPZ='V' or 'I', then LDZ >= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
*                     in Schur form, but ALPHAR(i), ALPHAI(i), and
*                     BETA(i), i=INFO+1,...,N should be correct.
*          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
*                     in Schur form, but ALPHAR(i), ALPHAI(i), and
*                     BETA(i), i=INFO-N+1,...,N should be correct.
*
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
hgeqz(const char* job, const char* compq, const char* compz, i_type n, 
      i_type ilo, i_type ihi, V* h, i_type ldh, V* t, i_type ldt, V* alphar, 
      V* alphai, V* beta, V* q, i_type ldq, V* z, i_type ldz, V* work, 
      i_type lwork, i_type& info);

//-----------------------------------------------------------------------
//                          ZHGEQZ
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZHGEQZ computes the eigenvalues of a complex matrix pair (H,T),
*  where H is an upper Hessenberg matrix and T is upper triangular,
*  using the single-shift QZ method.
*  Matrix pairs of this type are produced by the reduction to
*  generalized upper Hessenberg form of a complex matrix pair (A,B):
*  
*     A = Q1*H*Z1**H,  B = Q1*T*Z1**H,
*  
*  as computed by ZGGHRD.
*  
*  If JOB='S', then the Hessenberg-triangular pair (H,T) is
*  also reduced to generalized Schur form,
*  
*     H = Q*S*Z**H,  T = Q*P*Z**H,
*  
*  where Q and Z are unitary matrices and S and P are upper triangular.
*  
*  Optionally, the unitary matrix Q from the generalized Schur
*  factorization may be postmultiplied into an input matrix Q1, and the
*  unitary matrix Z may be postmultiplied into an input matrix Z1.
*  If Q1 and Z1 are the unitary matrices from ZGGHRD that reduced
*  the matrix pair (A,B) to generalized Hessenberg form, then the output
*  matrices Q1*Q and Z1*Z are the unitary factors from the generalized
*  Schur factorization of (A,B):
*  
*     A = (Q1*Q)*S*(Z1*Z)**H,  B = (Q1*Q)*P*(Z1*Z)**H.
*  
*  To avoid overflow, eigenvalues of the matrix pair (H,T)
*  (equivalently, of (A,B)) are computed as a pair of complex values
*  (alpha,beta).  If beta is nonzero, lambda = alpha / beta is an
*  eigenvalue of the generalized nonsymmetric eigenvalue problem (GNEP)
*     A*x = lambda*B*x
*  and if alpha is nonzero, mu = beta / alpha is an eigenvalue of the
*  alternate form of the GNEP
*     mu*A*y = B*y.
*  The values of alpha and beta for the i-th eigenvalue can be read
*  directly from the generalized Schur form:  alpha = S(i,i),
*  beta = P(i,i).
*
*  Ref: C.B. Moler & G.W. Stewart, "An Algorithm for Generalized Matrix
*       Eigenvalue Problems", SIAM J. Numer. Anal., 10(1973),
*       pp. 241--256.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          = 'E': Compute eigenvalues only;
*          = 'S': Computer eigenvalues and the Schur form.
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': Left Schur vectors (Q) are not computed;
*          = 'I': Q is initialized to the unit matrix and the matrix Q
*                 of left Schur vectors of (H,T) is returned;
*          = 'V': Q must contain a unitary matrix Q1 on entry and
*                 the product Q1*Q is returned.
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': Right Schur vectors (Z) are not computed;
*          = 'I': Q is initialized to the unit matrix and the matrix Z
*                 of right Schur vectors of (H,T) is returned;
*          = 'V': Z must contain a unitary matrix Z1 on entry and
*                 the product Z1*Z is returned.
*
*  N       (input) INTEGER
*          The order of the matrices H, T, Q, and Z.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          ILO and IHI mark the rows and columns of H which are in
*          Hessenberg form.  It is assumed that A is already upper
*          triangular in rows and columns 1:ILO-1 and IHI+1:N.
*          If N > 0, 1 <= ILO <= IHI <= N; if N = 0, ILO=1 and IHI=0.
*
*  H       (input/output) COMPLEX*16 array, dimension (LDH, N)
*          On entry, the N-by-N upper Hessenberg matrix H.
*          On exit, if JOB = 'S', H contains the upper triangular
*          matrix S from the generalized Schur factorization.
*          If JOB = 'E', the diagonal of H matches that of S, but
*          the rest of H is unspecified.
*
*  LDH     (input) INTEGER
*          The leading dimension of the array H.  LDH >= max( 1, N ).
*
*  T       (input/output) COMPLEX*16 array, dimension (LDT, N)
*          On entry, the N-by-N upper triangular matrix T.
*          On exit, if JOB = 'S', T contains the upper triangular
*          matrix P from the generalized Schur factorization.
*          If JOB = 'E', the diagonal of T matches that of P, but
*          the rest of T is unspecified.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T.  LDT >= max( 1, N ).
*
*  ALPHA   (output) COMPLEX*16 array, dimension (N)
*          The complex scalars alpha that define the eigenvalues of
*          GNEP.  ALPHA(i) = S(i,i) in the generalized Schur
*          factorization.
*
*  BETA    (output) COMPLEX*16 array, dimension (N)
*          The real non-negative scalars beta that define the
*          eigenvalues of GNEP.  BETA(i) = P(i,i) in the generalized
*          Schur factorization.
*
*          Together, the quantities alpha = ALPHA(j) and beta = BETA(j)
*          represent the j-th eigenvalue of the matrix pair (A,B), in
*          one of the forms lambda = alpha/beta or mu = beta/alpha.
*          Since either lambda or mu may overflow, they should not,
*          in general, be computed.
*
*  Q       (input/output) COMPLEX*16 array, dimension (LDQ, N)
*          On entry, if COMPZ = 'V', the unitary matrix Q1 used in the
*          reduction of (A,B) to generalized Hessenberg form.
*          On exit, if COMPZ = 'I', the unitary matrix of left Schur
*          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
*          left Schur vectors of (A,B).
*          Not referenced if COMPZ = 'N'.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.  LDQ >= 1.
*          If COMPQ='V' or 'I', then LDQ >= N.
*
*  Z       (input/output) COMPLEX*16 array, dimension (LDZ, N)
*          On entry, if COMPZ = 'V', the unitary matrix Z1 used in the
*          reduction of (A,B) to generalized Hessenberg form.
*          On exit, if COMPZ = 'I', the unitary matrix of right Schur
*          vectors of (H,T), and if COMPZ = 'V', the unitary matrix of
*          right Schur vectors of (A,B).
*          Not referenced if COMPZ = 'N'.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1.
*          If COMPZ='V' or 'I', then LDZ >= N.
*
*  WORK    (workspace/output) COMPLEX*16 array, dimension (MAX(1,LWORK))
*          On exit, if INFO >= 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  RWORK   (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1,...,N: the QZ iteration did not converge.  (H,T) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO+1,...,N should be correct.
*          = N+1,...,2*N: the shift calculation failed.  (H,T) is not
*                     in Schur form, but ALPHA(i) and BETA(i),
*                     i=INFO-N+1,...,N should be correct.
*
*  Further Details
*  ===============
*
*  We assume that complex ABS works as long as its value is less than
*  overflow.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid_complex<void,V>::type
hgeqz(const char* job, const char* compq, const char* compz, i_type n,
      i_type ilo, i_type ihi, V* h, i_type ldh, V* t, i_type ldt, V* alpha, 
      V* beta, V* q, i_type ldq, V* z, i_type ldz, V* work, i_type lwork,
      typename details::real_type<V>::type* rwork, i_type& info);

//-----------------------------------------------------------------------
//                          DGGBAK
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGGBAK forms the right or left eigenvectors of a real generalized
*  eigenvalue problem A*x = lambda*B*x, by backward transformation on
*  the computed eigenvectors of the balanced pair of matrices output by
*  DGGBAL.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies the type of backward transformation required:
*          = 'N':  do nothing, return immediately;
*          = 'P':  do backward transformation for permutation only;
*          = 'S':  do backward transformation for scaling only;
*          = 'B':  do backward transformations for both permutation and
*                  scaling.
*          JOB must be the same as the argument JOB supplied to DGGBAL.
*
*  SIDE    (input) CHARACTER*1
*          = 'R':  V contains right eigenvectors;
*          = 'L':  V contains left eigenvectors.
*
*  N       (input) INTEGER
*          The number of rows of the matrix V.  N >= 0.
*
*  ILO     (input) INTEGER
*  IHI     (input) INTEGER
*          The integers ILO and IHI determined by DGGBAL.
*          1 <= ILO <= IHI <= N, if N > 0; ILO=1 and IHI=0, if N=0.
*
*  LSCALE  (input) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and/or scaling factors applied
*          to the left side of A and B, as returned by DGGBAL.
*
*  RSCALE  (input) DOUBLE PRECISION array, dimension (N)
*          Details of the permutations and/or scaling factors applied
*          to the right side of A and B, as returned by DGGBAL.
*
*  M       (input) INTEGER
*          The number of columns of the matrix V.  M >= 0.
*
*  V       (input/output) DOUBLE PRECISION array, dimension (LDV,M)
*          On entry, the matrix of right or left eigenvectors to be
*          transformed, as returned by DTGEVC.
*          On exit, V is overwritten by the transformed eigenvectors.
*
*  LDV     (input) INTEGER
*          The leading dimension of the matrix V. LDV >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*
*  See R.C. Ward, Balancing the generalized eigenvalue problem,
*                 SIAM J. Sci. Stat. Comp. 2 (1981), 141-152.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ggbak(const char* job, const char* side, i_type n, i_type ilo, i_type ihi, 
      const typename details::real_type<V>::type* lscale, 
      const typename details::real_type<V>::type* rscale, i_type m, V* v, 
      i_type ldv, i_type& info);

//-----------------------------------------------------------------------
//                          DLANHS
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLANHS  returns the value of the one norm,  or the Frobenius norm, or
*  the  infinity norm,  or the  element of  largest absolute value  of a
*  Hessenberg matrix A.
*
*  Description
*  ===========
*
*  DLANHS returns the value
*
*     DLANHS = ( max(abs(A(i,j))), NORM = 'M' or 'm'
*              (
*              ( norm1(A),         NORM = '1', 'O' or 'o'
*              (
*              ( normI(A),         NORM = 'I' or 'i'
*              (
*              ( normF(A),         NORM = 'F', 'f', 'E' or 'e'
*
*  where  norm1  denotes the  one norm of a matrix (maximum column sum),
*  normI  denotes the  infinity norm  of a matrix  (maximum row sum) and
*  normF  denotes the  Frobenius norm of a matrix (square root of sum of
*  squares).  Note that  max(abs(A(i,j)))  is not a  matrix norm.
*
*  Arguments
*  =========
*
*  NORM    (input) CHARACTER*1
*          Specifies the value to be returned in DLANHS as described
*          above.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.  When N = 0, DLANHS is
*          set to zero.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,N)
*          The n by n upper Hessenberg matrix A; the part of A below the
*          first sub-diagonal is not referenced.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(N,1).
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (LWORK),
*          where LWORK >= N when NORM = 'I'; otherwise, WORK is not
*          referenced.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<typename details::real_type<V>::type,V>::type
lanhs(const char* NORM, i_type N, const V* A, i_type LDA, V* WORK );

//-----------------------------------------------------------------------
//                          DLAG2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLAG2 computes the eigenvalues of a 2 x 2 generalized eigenvalue
*  problem  A - w B, with scaling as necessary to avoid over-/underflow.
*
*  The scaling factor "s" results in a modified eigenvalue equation
*
*      s A - w B
*
*  where  s  is a non-negative scaling factor chosen so that  w,  w B,
*  and  s A  do not overflow and, if possible, do not underflow, either.
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA, 2)
*          On entry, the 2 x 2 matrix A.  It is assumed that its 1-norm
*          is less than 1/SAFMIN.  Entries less than
*          sqrt(SAFMIN)*norm(A) are subject to being treated as zero.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= 2.
*
*  B       (input) DOUBLE PRECISION array, dimension (LDB, 2)
*          On entry, the 2 x 2 upper triangular matrix B.  It is
*          assumed that the one-norm of B is less than 1/SAFMIN.  The
*          diagonals should be at least sqrt(SAFMIN) times the largest
*          element of B (in absolute value); if a diagonal is smaller
*          than that, then  +/- sqrt(SAFMIN) will be used instead of
*          that diagonal.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= 2.
*
*  SAFMIN  (input) DOUBLE PRECISION
*          The smallest positive number s.t. 1/SAFMIN does not
*          overflow.  (This should always be DLAMCH('S') -- it is an
*          argument in order to avoid having to call DLAMCH frequently.)
*
*  SCALE1  (output) DOUBLE PRECISION
*          A scaling factor used to avoid over-/underflow in the
*          eigenvalue equation which defines the first eigenvalue.  If
*          the eigenvalues are complex, then the eigenvalues are
*          ( WR1  +/-  WI i ) / SCALE1  (which may lie outside the
*          exponent range of the machine), SCALE1=SCALE2, and SCALE1
*          will always be positive.  If the eigenvalues are real, then
*          the first (real) eigenvalue is  WR1 / SCALE1 , but this may
*          overflow or underflow, and in fact, SCALE1 may be zero or
*          less than the underflow threshhold if the exact eigenvalue
*          is sufficiently large.
*
*  SCALE2  (output) DOUBLE PRECISION
*          A scaling factor used to avoid over-/underflow in the
*          eigenvalue equation which defines the second eigenvalue.  If
*          the eigenvalues are complex, then SCALE2=SCALE1.  If the
*          eigenvalues are real, then the second (real) eigenvalue is
*          WR2 / SCALE2 , but this may overflow or underflow, and in
*          fact, SCALE2 may be zero or less than the underflow
*          threshhold if the exact eigenvalue is sufficiently large.
*
*  WR1     (output) DOUBLE PRECISION
*          If the eigenvalue is real, then WR1 is SCALE1 times the
*          eigenvalue closest to the (2,2) element of A B**(-1).  If the
*          eigenvalue is complex, then WR1=WR2 is SCALE1 times the real
*          part of the eigenvalues.
*
*  WR2     (output) DOUBLE PRECISION
*          If the eigenvalue is real, then WR2 is SCALE2 times the
*          other eigenvalue.  If the eigenvalue is complex, then
*          WR1=WR2 is SCALE1 times the real part of the eigenvalues.
*
*  WI      (output) DOUBLE PRECISION
*          If the eigenvalue is real, then WI is zero.  If the
*          eigenvalue is complex, then WI is SCALE1 times the imaginary
*          part of the eigenvalues.  WI will always be non-negative.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
lag2(const V* A, i_type LDA, const V* B, i_type LDB, V SAFMIN, V& SCALE1,
     V& SCALE2, V& WR1, V& WR2, V& WI );

//-----------------------------------------------------------------------
//                          DLAPY3
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLAPY3 returns sqrt(x**2+y**2+z**2), taking care not to cause
*  unnecessary overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*  Z       (input) DOUBLE PRECISION
*          X, Y and Z specify the values x, y and z.
*/
template<class V>
typename details::enable_if_valid_real<V,V>::type
lapy3(const V& X, const V& Y, const V& Z);

//-----------------------------------------------------------------------
//                          DLAPY2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*/
template<class V>
typename details::enable_if_valid_real<V,V>::type
lapy2(const V& X, const V& Y);

//-----------------------------------------------------------------------
//                          DSBGVD
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSBGVD computes all the eigenvalues, and optionally, the eigenvectors
*  of a real generalized symmetric-definite banded eigenproblem, of the
*  form A*x=(lambda)*B*x.  Here A and B are assumed to be symmetric and
*  banded, and B is also positive definite.  If eigenvectors are
*  desired, it uses a divide and conquer algorithm.
*
*  The divide and conquer algorithm makes very mild assumptions about
*  floating point arithmetic. It will work on machines with a guard
*  digit in add/subtract, or on those binary machines without guard
*  digits which subtract like the Cray X-MP, Cray Y-MP, Cray C-90, or
*  Cray-2. It could conceivably fail on hexadecimal or decimal machines
*  without guard digits, but we know of none.
*
*  Arguments
*  =========
*
*  JOBZ    (input) CHARACTER*1
*          = 'N':  Compute eigenvalues only;
*          = 'V':  Compute eigenvalues and eigenvectors.
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  Upper triangles of A and B are stored;
*          = 'L':  Lower triangles of A and B are stored.
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
*
*  KA      (input) INTEGER
*          The number of superdiagonals of the matrix A if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KA >= 0.
*
*  KB      (input) INTEGER
*          The number of superdiagonals of the matrix B if UPLO = 'U',
*          or the number of subdiagonals if UPLO = 'L'.  KB >= 0.
*
*  AB      (input/output) DOUBLE PRECISION array, dimension (LDAB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix A, stored in the first ka+1 rows of the array.  The
*          j-th column of A is stored in the j-th column of the array AB
*          as follows:
*          if UPLO = 'U', AB(ka+1+i-j,j) = A(i,j) for max(1,j-ka)<=i<=j;
*          if UPLO = 'L', AB(1+i-j,j)    = A(i,j) for j<=i<=min(n,j+ka).
*
*          On exit, the contents of AB are destroyed.
*
*  LDAB    (input) INTEGER
*          The leading dimension of the array AB.  LDAB >= KA+1.
*
*  BB      (input/output) DOUBLE PRECISION array, dimension (LDBB, N)
*          On entry, the upper or lower triangle of the symmetric band
*          matrix B, stored in the first kb+1 rows of the array.  The
*          j-th column of B is stored in the j-th column of the array BB
*          as follows:
*          if UPLO = 'U', BB(ka+1+i-j,j) = B(i,j) for max(1,j-kb)<=i<=j;
*          if UPLO = 'L', BB(1+i-j,j)    = B(i,j) for j<=i<=min(n,j+kb).
*
*          On exit, the factor S from the split Cholesky factorization
*          B = S**T*S, as returned by DPBSTF.
*
*  LDBB    (input) INTEGER
*          The leading dimension of the array BB.  LDBB >= KB+1.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          If INFO = 0, the eigenvalues in ascending order.
*
*  Z       (output) DOUBLE PRECISION array, dimension (LDZ, N)
*          If JOBZ = 'V', then if INFO = 0, Z contains the matrix Z of
*          eigenvectors, with the i-th column of Z holding the
*          eigenvector associated with W(i).  The eigenvectors are
*          normalized so Z**T*B*Z = I.
*          If JOBZ = 'N', then Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.  LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(1,N).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If N <= 1,               LWORK >= 1.
*          If JOBZ = 'N' and N > 1, LWORK >= 3*N.
*          If JOBZ = 'V' and N > 1, LWORK >= 1 + 5*N + 2*N**2.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal sizes of the WORK and IWORK
*          arrays, returns these values as the first entries of the WORK
*          and IWORK arrays, and no error message related to LWORK or
*          LIWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if LIWORK > 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If JOBZ  = 'N' or N <= 1, LIWORK >= 1.
*          If JOBZ  = 'V' and N > 1, LIWORK >= 3 + 5*N.
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal sizes of the WORK and
*          IWORK arrays, returns these values as the first entries of
*          the WORK and IWORK arrays, and no error message related to
*          LWORK or LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, and i is:
*             <= N:  the algorithm failed to converge:
*                    i off-diagonal elements of an intermediate
*                    tridiagonal form did not converge to zero;
*             > N:   if INFO = N + i, for 1 <= i <= N, then DPBSTF
*                    returned INFO = i: B is not positive definite.
*                    The factorization of B could not be completed and
*                    no eigenvalues or eigenvectors were computed.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
hbgvd(const char* JOBZ, const char* UPLO, i_type N, i_type KA, i_type KB,
      V* AB, i_type LDAB, V* BB, i_type LDBB, 
      typename details::real_type<V>::type* W, V* Z, i_type LDZ, V* WORK, 
      i_type LWORK, typename details::real_type<V>::type* RWORK, i_type LRWORK,
      i_type* IWORK, i_type LIWORK, i_type& INFO );

//-----------------------------------------------------------------------
//                          DNRM2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DNRM2 returns the euclidean norm of a vector via the function
*  name, so that
*
*     DNRM2 := sqrt( x'*x )
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<typename details::real_type<V>::type,V>::type
nrm2(i_type N, const V* X, i_type INCX);

//-----------------------------------------------------------------------
//                          DLAS2
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DLAS2  computes the singular values of the 2-by-2 matrix
*     [  F   G  ]
*     [  0   H  ].
*  On return, SSMIN is the smaller singular value and SSMAX is the
*  larger singular value.
*
*  Arguments
*  =========
*
*  F       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  G       (input) DOUBLE PRECISION
*          The (1,2) element of the 2-by-2 matrix.
*
*  H       (input) DOUBLE PRECISION
*          The (2,2) element of the 2-by-2 matrix.
*
*  SSMIN   (output) DOUBLE PRECISION
*          The smaller singular value.
*
*  SSMAX   (output) DOUBLE PRECISION
*          The larger singular value.
*
*  Further Details
*  ===============
*
*  Barring over/underflow, all output quantities are correct to within
*  a few units in the last place (ulps), even in the absence of a guard
*  digit in addition/subtraction.
*
*  In IEEE arithmetic, the code works correctly if one matrix element is
*  infinite.
*
*  Overflow will not occur unless the largest singular value itself
*  overflows, or is within a few ulps of overflow. (On machines with
*  partial overflow, like the Cray, overflow may occur if the largest
*  singular value is within a factor of 2 of overflow.)
*
*  Underflow is harmless if underflow is gradual. Otherwise, results
*  may correspond to a matrix modified by perturbations of size near
*  the underflow threshold.
*/
template<class V>
typename details::enable_if_valid_real<void,V>::type
las2(V F, V G, V H, V& SSMIN, V& SSMAX );

//-----------------------------------------------------------------------
//                          ZPTTRF
//-----------------------------------------------------------------------
/*
*  ZPTTRF computes the L*D*L' factorization of a complex Hermitian
*  positive definite tridiagonal matrix A.  The factorization may also
*  be regarded as having the form A = U'*D*U.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, the n diagonal elements of the tridiagonal matrix
*          A.  On exit, the n diagonal elements of the diagonal matrix
*          D from the L*D*L' factorization of A.
*
*  E       (input/output) COMPLEX*16 array, dimension (N-1)
*          On entry, the (n-1) subdiagonal elements of the tridiagonal
*          matrix A.  On exit, the (n-1) subdiagonal elements of the
*          unit bidiagonal factor L from the L*D*L' factorization of A.
*          E can also be regarded as the superdiagonal of the unit
*          bidiagonal factor U from the U'*D*U factorization of A.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*          > 0: if INFO = k, the leading minor of order k is not
*               positive definite; if k < N, the factorization could not
*               be completed, while if k = N, the factorization was
*               completed, but D(N) <= 0.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
pttrf(i_type N, typename details::real_type<V>::type* D, V* E, i_type& info);

//-----------------------------------------------------------------------
//                          PTTRS
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ZPTTRS solves a tridiagonal system of the form
*     A * X = B
*  using the factorization A = U'*D*U or A = L*D*L' computed by ZPTTRF.
*  D is a diagonal matrix specified in the vector D, U (or L) is a unit
*  bidiagonal matrix whose superdiagonal (subdiagonal) is specified in
*  the vector E, and X and B are N by NRHS matrices.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          Specifies the form of the factorization and whether the
*          vector E is the superdiagonal of the upper bidiagonal factor
*          U or the subdiagonal of the lower bidiagonal factor L.
*          = 'U':  A = U'*D*U, E is the superdiagonal of U
*          = 'L':  A = L*D*L', E is the subdiagonal of L
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the diagonal matrix D from the
*          L*D*L' factorization of A.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) subdiagonal elements of the unit bidiagonal factor
*          L from the L*D*L' factorization of A.  E can also be regarded
*          as the superdiagonal of the unit bidiagonal factor U from the
*          factorization A = U'*D*U.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the right hand side vectors B for the system of
*          linear equations.
*          On exit, the solution vectors, X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
pttrs(const char* UPLO, i_type N, i_type NRHS, 
      const typename details::real_type<V>::type* D, const V* E, V* B, 
      i_type LDB, i_type& INFO );

//-----------------------------------------------------------------------
//                          DTRTRI
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTRTRI computes the inverse of a real upper or lower triangular
*  matrix A.
*
*  This is the Level 3 BLAS version of the algorithm.
*
*  Arguments
*  =========
*
*  UPLO    (input) CHARACTER*1
*          = 'U':  A is upper triangular;
*          = 'L':  A is lower triangular.
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the triangular matrix A.  If UPLO = 'U', the
*          leading N-by-N upper triangular part of the array A contains
*          the upper triangular matrix, and the strictly lower
*          triangular part of A is not referenced.  If UPLO = 'L', the
*          leading N-by-N lower triangular part of the array A contains
*          the lower triangular matrix, and the strictly upper
*          triangular part of A is not referenced.  If DIAG = 'U', the
*          diagonal elements of A are also not referenced and are
*          assumed to be 1.
*          On exit, the (triangular) inverse of the original matrix, in
*          the same storage format.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, A(i,i) is exactly zero.  The triangular
*               matrix is singular and its inverse can not be computed.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
trtri(const char* UPLO, const char* DIAG, i_type N, V* A, i_type LDA, 
      i_type& INFO);

//-----------------------------------------------------------------------
//                          DGETRI
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGETRI computes the inverse of a matrix using the LU factorization
*  computed by DGETRF.
*
*  This method inverts U and then computes inv(A) by solving the system
*  inv(A)*L = inv(U) for inv(A).
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the factors L and U from the factorization
*          A = P*L*U as computed by DGETRF.
*          On exit, if INFO = 0, the inverse of the original matrix A.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,N).
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices from DGETRF; for 1<=i<=N, row i of the
*          matrix was interchanged with row IPIV(i).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO=0, then WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,N).
*          For optimal performance LWORK >= N*NB, where NB is
*          the optimal blocksize returned by ILAENV.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, U(i,i) is exactly zero; the matrix is
*                singular and its inverse could not be computed.
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
getri(i_type N, V* A, i_type LDA, const i_type* IPIV, V* WORK, i_type LWORK,
      i_type& INFO );

//-----------------------------------------------------------------------
//                          DGTTRF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGTTRF computes an LU factorization of a real tridiagonal matrix A
*  using elimination with partial pivoting and row interchanges.
*
*  The factorization has the form
*     A = L * U
*  where L is a product of permutation and unit lower bidiagonal
*  matrices and U is upper triangular with nonzeros in only the main
*  diagonal and first two superdiagonals.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  DL      (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, DL must contain the (n-1) sub-diagonal elements of
*          A.
*
*          On exit, DL is overwritten by the (n-1) multipliers that
*          define the matrix L from the LU factorization of A.
*
*  D       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry, D must contain the diagonal elements of A.
*
*          On exit, D is overwritten by the n diagonal elements of the
*          upper triangular matrix U from the LU factorization of A.
*
*  DU      (input/output) DOUBLE PRECISION array, dimension (N-1)
*          On entry, DU must contain the (n-1) super-diagonal elements
*          of A.
*
*          On exit, DU is overwritten by the (n-1) elements of the first
*          super-diagonal of U.
*
*  DU2     (output) DOUBLE PRECISION array, dimension (N-2)
*          On exit, DU2 is overwritten by the (n-2) elements of the
*          second super-diagonal of U.
*
*  IPIV    (output) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= n, row i of the matrix was
*          interchanged with row IPIV(i).  IPIV(i) will always be either
*          i or i+1; IPIV(i) = i indicates a row interchange was not
*          required.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -k, the k-th argument had an illegal value
*          > 0:  if INFO = k, U(k,k) is exactly zero. The factorization
*                has been completed, but the factor U is exactly
*                singular, and division by zero will occur if it is used
*                to solve a system of equations.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gttrf(i_type N, V* DL, V* D, V* DU, V* DU2, i_type* IPIV, i_type& INFO);

//-----------------------------------------------------------------------
//                          DGTTRS
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGTTRS solves one of the systems of equations
*     A*X = B  or  A'*X = B,
*  with a tridiagonal matrix A using the LU factorization computed
*  by DGTTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  A * X = B  (No transpose)
*          = 'T':  A'* X = B  (Transpose)
*          = 'C':  A'* X = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of columns
*          of the matrix B.  NRHS >= 0.
*
*  DL      (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) multipliers that define the matrix L from the
*          LU factorization of A.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the upper triangular matrix U from
*          the LU factorization of A.
*
*  DU      (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) elements of the first super-diagonal of U.
*
*  DU2     (input) DOUBLE PRECISION array, dimension (N-2)
*          The (n-2) elements of the second super-diagonal of U.
*
*  IPIV    (input) INTEGER array, dimension (N)
*          The pivot indices; for 1 <= i <= n, row i of the matrix was
*          interchanged with row IPIV(i).  IPIV(i) will always be either
*          i or i+1; IPIV(i) = i indicates a row interchange was not
*          required.
*
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,NRHS)
*          On entry, the matrix of right hand side vectors B.
*          On exit, B is overwritten by the solution vectors X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,N).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gttrs(const char* TRANS, i_type N, i_type NRHS, const V* DL, const V* D, 
      const V* DU, const V* DU2, const i_type* IPIV, V* B, i_type LDB, 
      i_type& INFO);

//-----------------------------------------------------------------------
//                          DSTEBZ
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DSTEBZ computes the eigenvalues of a symmetric tridiagonal
*  matrix T.  The user may ask for all eigenvalues, all eigenvalues
*  in the half-open interval (VL, VU], or the IL-th through IU-th
*  eigenvalues.
*
*  To avoid overflow, the matrix must be scaled so that its
*  largest element is no greater than overflow**(1/2) *
*  underflow**(1/4) in absolute value, and for greatest
*  accuracy, it should not be much smaller than that.
*
*  See W. Kahan "Accurate Eigenvalues of a Symmetric Tridiagonal
*  Matrix", Report CS41, Computer Science Dept., Stanford
*  University, July 21, 1966.
*
*  Arguments
*  =========
*
*  RANGE   (input) CHARACTER*1
*          = 'A': ("All")   all eigenvalues will be found.
*          = 'V': ("Value") all eigenvalues in the half-open interval
*                           (VL, VU] will be found.
*          = 'I': ("Index") the IL-th through IU-th eigenvalues (of the
*                           entire matrix) will be found.
*
*  ORDER   (input) CHARACTER*1
*          = 'B': ("By Block") the eigenvalues will be grouped by
*                              split-off block (see IBLOCK, ISPLIT) and
*                              ordered from smallest to largest within
*                              the block.
*          = 'E': ("Entire matrix")
*                              the eigenvalues for the entire matrix
*                              will be ordered from smallest to
*                              largest.
*
*  N       (input) INTEGER
*          The order of the tridiagonal matrix T.  N >= 0.
*
*  VL      (input) DOUBLE PRECISION
*  VU      (input) DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for eigenvalues.  Eigenvalues less than or equal
*          to VL, or greater than VU, will not be returned.  VL < VU.
*          Not referenced if RANGE = 'A' or 'I'.
*
*  IL      (input) INTEGER
*  IU      (input) INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest eigenvalues to be returned.
*          1 <= IL <= IU <= N, if N > 0; IL = 1 and IU = 0 if N = 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
*  ABSTOL  (input) DOUBLE PRECISION
*          The absolute tolerance for the eigenvalues.  An eigenvalue
*          (or cluster) is considered to be located if it has been
*          determined to lie in an interval whose width is ABSTOL or
*          less.  If ABSTOL is less than or equal to zero, then ULP*|T|
*          will be used, where |T| means the 1-norm of T.
*
*          Eigenvalues will be computed most accurately when ABSTOL is
*          set to twice the underflow threshold 2*DLAMCH('S'), not zero.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the tridiagonal matrix T.
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) off-diagonal elements of the tridiagonal matrix T.
*
*  M       (output) INTEGER
*          The actual number of eigenvalues found. 0 <= M <= N.
*          (See also the description of INFO=2,3.)
*
*  NSPLIT  (output) INTEGER
*          The number of diagonal blocks in the matrix T.
*          1 <= NSPLIT <= N.
*
*  W       (output) DOUBLE PRECISION array, dimension (N)
*          On exit, the first M elements of W will contain the
*          eigenvalues.  (DSTEBZ may use the remaining N-M elements as
*          workspace.)
*
*  IBLOCK  (output) INTEGER array, dimension (N)
*          At each row/column j where E(j) is zero or small, the
*          matrix T is considered to split into a block diagonal
*          matrix.  On exit, if INFO = 0, IBLOCK(i) specifies to which
*          block (from 1 to the number of blocks) the eigenvalue W(i)
*          belongs.  (DSTEBZ may use the remaining N-M elements as
*          workspace.)
*
*  ISPLIT  (output) INTEGER array, dimension (N)
*          The splitting points, at which T breaks up into submatrices.
*          The first submatrix consists of rows/columns 1 to ISPLIT(1),
*          the second of rows/columns ISPLIT(1)+1 through ISPLIT(2),
*          etc., and the NSPLIT-th consists of rows/columns
*          ISPLIT(NSPLIT-1)+1 through ISPLIT(NSPLIT)=N.
*          (Only the first NSPLIT elements will actually be used, but
*          since the user cannot know a priori what value NSPLIT will
*          have, N words must be reserved for ISPLIT.)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (4*N)
*
*  IWORK   (workspace) INTEGER array, dimension (3*N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  some or all of the eigenvalues failed to converge or
*                were not computed:
*                =1 or 3: Bisection failed to converge for some
*                        eigenvalues; these eigenvalues are flagged by a
*                        negative block number.  The effect is that the
*                        eigenvalues may not be as accurate as the
*                        absolute and relative tolerances.  This is
*                        generally caused by unexpectedly inaccurate
*                        arithmetic.
*                =2 or 3: RANGE='I' only: Not all of the eigenvalues
*                        IL:IU were found.
*                        Effect: M < IU+1-IL
*                        Cause:  non-monotonic arithmetic, causing the
*                                Sturm sequence to be non-monotonic.
*                        Cure:   recalculate, using RANGE='A', and pick
*                                out eigenvalues IL:IU.  In some cases,
*                                increasing the PARAMETER "FUDGE" may
*                                make things work.
*                = 4:    RANGE='I', and the Gershgorin interval
*                        initially used was too small.  No eigenvalues
*                        were computed.
*                        Probable cause: your machine has sloppy
*                                        floating-point arithmetic.
*                        Cure: Increase the PARAMETER "FUDGE",
*                              recompile, and try again.
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid_real<void,V>::type
stebz(const char* RANGE, const char* ORDER, i_type N, V VL, V VU, i_type IL,
      i_type IU, V ABSTOL, const V* D, const V* E, i_type& M, i_type& NSPLIT, 
      V* W, i_type* IBLOCK, i_type* ISPLIT, V* WORK, i_type* IWORK, i_type& INFO );

//-----------------------------------------------------------------------
//                          DGGSVD3
//-----------------------------------------------------------------------
/*
* Purpose:
*  =============
*
* DGGSVD3 computes the generalized singular value decomposition (GSVD)
* of an M-by-N real matrix A and P-by-N real matrix B:
*
*       U**T*A*Q = D1*( 0 R ),    V**T*B*Q = D2*( 0 R )
*
* where U, V and Q are orthogonal matrices.
* Let K+L = the effective numerical rank of the matrix (A**T,B**T)**T,
* then R is a K+L-by-K+L nonsingular upper triangular matrix, D1 and
* D2 are M-by-(K+L) and P-by-(K+L) "diagonal" matrices and of the
* following structures, respectively:
*
* If M-K-L >= 0,
*
*                     K  L
*        D1 =     K ( I  0 )
*                 L ( 0  C )
*             M-K-L ( 0  0 )
*
*                   K  L
*        D2 =   L ( 0  S )
*             P-L ( 0  0 )
*
*                 N-K-L  K    L
*   ( 0 R ) = K (  0   R11  R12 )
*             L (  0    0   R22 )
*
* where
*
*   C = diag( ALPHA(K+1), ... , ALPHA(K+L) ),
*   S = diag( BETA(K+1),  ... , BETA(K+L) ),
*   C**2 + S**2 = I.
*
*   R is stored in A(1:K+L,N-K-L+1:N) on exit.
*
* If M-K-L < 0,
*
*                   K M-K K+L-M
*        D1 =   K ( I  0    0   )
*             M-K ( 0  C    0   )
*
*                     K M-K K+L-M
*        D2 =   M-K ( 0  S    0  )
*             K+L-M ( 0  0    I  )
*               P-L ( 0  0    0  )
*
*                    N-K-L  K   M-K  K+L-M
*   ( 0 R ) =     K ( 0    R11  R12  R13  )
*               M-K ( 0     0   R22  R23  )
*             K+L-M ( 0     0    0   R33  )
*
* where
*
*   C = diag( ALPHA(K+1), ... , ALPHA(M) ),
*   S = diag( BETA(K+1),  ... , BETA(M) ),
*   C**2 + S**2 = I.
*
*   (R11 R12 R13 ) is stored in A(1:M, N-K-L+1:N), and R33 is stored
*   ( 0  R22 R23 )
*   in B(M-K+1:L,N+M-K-L+1:N) on exit.
*
* The routine computes C, S, R, and optionally the orthogonal
* transformation matrices U, V and Q.
*
* In particular, if B is an N-by-N nonsingular matrix, then the GSVD of
* A and B implicitly gives the SVD of A*inv(B):
*                      A*inv(B) = U*(D1*inv(D2))*V**T.
* If ( A**T,B**T)**T  has orthonormal columns, then the GSVD of A and B is
* also equal to the CS decomposition of A and B. Furthermore, the GSVD
* can be used to derive the solution of the eigenvalue problem:
*                      A**T*A x = lambda* B**T*B x.
* In some literature, the GSVD of A and B is presented in the form
*                  U**T*A*X = ( 0 D1 ),   V**T*B*X = ( 0 D2 )
* where U and V are orthogonal and X is nonsingular, D1 and D2 are
* ``diagonal''.  The former GSVD form can be converted to the latter
* form by taking the nonsingular matrix X as
*
*                      X = Q*( I   0    )
*                            ( 0 inv(R) ).
*
*  Arguments:
*  ==========
*
* [in] JOBU
*          JOBU is CHARACTER*1
*          = 'U':  Orthogonal matrix U is computed;
*          = 'N':  U is not computed.
*
* [in] JOBV
*          JOBV is CHARACTER*1
*          = 'V':  Orthogonal matrix V is computed;
*          = 'N':  V is not computed.
*
* [in] JOBQ
*          JOBQ is CHARACTER*1
*          = 'Q':  Orthogonal matrix Q is computed;
*          = 'N':  Q is not computed.
*
* [in] M
*          M is INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
* [in] N
*          N is INTEGER
*          The number of columns of the matrices A and B.  N >= 0.
*
* [in] P
*          P is INTEGER
*          The number of rows of the matrix B.  P >= 0.
*
* [out] K
*          K is INTEGER
*
* [out] L
*          L is INTEGER
*
*          On exit, K and L specify the dimension of the subblocks
*          described in Purpose.
*          K + L = effective numerical rank of (A**T,B**T)**T.
*
* [in,out] A
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, A contains the triangular matrix R, or part of R.
*          See Purpose for details.
*
* [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
* [in,out] B
*          B is DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the P-by-N matrix B.
*          On exit, B contains the triangular matrix R if M-K-L < 0.
*          See Purpose for details.
*
* [in] LDB
*          LDB is INTEGER
*          The leading dimension of the array B. LDB >= max(1,P).
*
* [out] ALPHA
*          ALPHA is DOUBLE PRECISION array, dimension (N)
*
* [out] BETA
*          BETA is DOUBLE PRECISION array, dimension (N)
*
*          On exit, ALPHA and BETA contain the generalized singular
*          value pairs of A and B;
*            ALPHA(1:K) = 1,
*            BETA(1:K)  = 0,
*          and if M-K-L >= 0,
*            ALPHA(K+1:K+L) = C,
*            BETA(K+1:K+L)  = S,
*          or if M-K-L < 0,
*            ALPHA(K+1:M)=C, ALPHA(M+1:K+L)=0
*            BETA(K+1:M) =S, BETA(M+1:K+L) =1
*          and
*            ALPHA(K+L+1:N) = 0
*            BETA(K+L+1:N)  = 0
*
* [out] U
*          U is DOUBLE PRECISION array, dimension (LDU,M)
*          If JOBU = 'U', U contains the M-by-M orthogonal matrix U.
*          If JOBU = 'N', U is not referenced.
*
* [in] LDU
*          LDU is INTEGER
*          The leading dimension of the array U. LDU >= max(1,M) if
*          JOBU = 'U'; LDU >= 1 otherwise.
*
* [out] V
*          V is DOUBLE PRECISION array, dimension (LDV,P)
*          If JOBV = 'V', V contains the P-by-P orthogonal matrix V.
*          If JOBV = 'N', V is not referenced.
*
* [in] LDV
*          LDV is INTEGER
*          The leading dimension of the array V. LDV >= max(1,P) if
*          JOBV = 'V'; LDV >= 1 otherwise.
*
* [out] Q
*          Q is DOUBLE PRECISION array, dimension (LDQ,N)
*          If JOBQ = 'Q', Q contains the N-by-N orthogonal matrix Q.
*          If JOBQ = 'N', Q is not referenced.
*
* [in] LDQ
*          LDQ is INTEGER
*          The leading dimension of the array Q. LDQ >= max(1,N) if
*          JOBQ = 'Q'; LDQ >= 1 otherwise.
*
* [out] WORK
*          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
* [in] LWORK
*          LWORK is INTEGER
*          The dimension of the array WORK.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
* [out]	RWORK	
*          RWORK is DOUBLE PRECISION array, dimension (2*N), referrend only
*          for complex value types
*
* [out] IWORK
*          IWORK is INTEGER array, dimension (N)
*          On exit, IWORK stores the sorting information. More
*          precisely, the following loop will sort ALPHA
*             for I = K+1, min(M,K+L)
*                 swap ALPHA(I) and ALPHA(IWORK(I))
*             endfor
*          such that ALPHA(1) >= ALPHA(2) >= ... >= ALPHA(N).
*
* [out] INFO
*          INFO is INTEGER
*          = 0:  successful exit.
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          > 0:  if INFO = 1, the Jacobi-type procedure failed to
*                converge.  For further details, see subroutine DTGSJA.
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ggsvd3(const char* JOBU, const char* JOBV, const char* JOBQ, i_type M, i_type N,
       i_type P, i_type& K, i_type& L, V* ptr_A, i_type LDA, V* ptr_B, i_type LDB,
       typename details::real_type<V>::type* ptr_alpha, 
       typename details::real_type<V>::type* ptr_beta, V* ptr_U, i_type LDU, V* ptr_V, 
       i_type LDV, V* ptr_Q, i_type LDQ, V* work, i_type lwork, 
       typename details::real_type<V>::type* rwork, i_type* iwork, i_type& info);

//-----------------------------------------------------------------------
//                          DGESVDX
//-----------------------------------------------------------------------
/*
* Purpose:
*  =============
*
*  DGESVDX computes the singular value decomposition (SVD) of a real
*  M-by-N matrix A, optionally computing the left and/or right singular
*  vectors. The SVD is written
* 
*      A = U * SIGMA * transpose(V)
* 
*  where SIGMA is an M-by-N matrix which is zero except for its
*  min(m,n) diagonal elements, U is an M-by-M orthogonal matrix, and
*  V is an N-by-N orthogonal matrix.  The diagonal elements of SIGMA
*  are the singular values of A; they are real and non-negative, and
*  are returned in descending order.  The first min(m,n) columns of
*  U and V are the left and right singular vectors of A.
* 
*  DGESVDX uses an eigenvalue problem for obtaining the SVD, which 
*  allows for the computation of a subset of singular values and 
*  vectors. See DBDSVDX for details.
* 
*  Note that the routine returns V**T, not V.
*   
*  Arguments:
*  ==========
*
* [in] JOBU
*          JOBU is CHARACTER*1
*          Specifies options for computing all or part of the matrix U:
*          = 'V':  the first min(m,n) columns of U (the left singular
*                  vectors) or as specified by RANGE are returned in 
*                  the array U;
*          = 'N':  no columns of U (no left singular vectors) are
*                  computed.
*
* [in] JOBVT
*          JOBVT is CHARACTER*1
*           Specifies options for computing all or part of the matrix
*           V**T:
*           = 'V':  the first min(m,n) rows of V**T (the right singular
*                   vectors) or as specified by RANGE are returned in 
*                   the array VT;
*           = 'N':  no rows of V**T (no right singular vectors) are
*                   computed.
*
* [in] RANGE
*          RANGE is CHARACTER*1
*          = 'A': all singular values will be found.
*          = 'V': all singular values in the half-open interval (VL,VU]
*                 will be found.
*          = 'I': the IL-th through IU-th singular values will be found. 
*
* [in] M
*          M is INTEGER
*          The number of rows of the input matrix A.  M >= 0.
*
* [in] N
*          N is INTEGER
*          The number of columns of the input matrix A.  N >= 0.
*
* [in,out] A
*          A is DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the contents of A are destroyed.
*
* [in] LDA
*          LDA is INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
* [in] VL
*          VL is DOUBLE PRECISION
*          VL >=0.
*
* [in] VU
*          VU is DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for singular values. VU > VL.
*          Not referenced if RANGE = 'A' or 'I'.
*
* [in] IL
*          IL is INTEGER
*
* [in] IU
*          IU is INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest singular values to be returned.
*          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
* [out] NS
*          NS is INTEGER
*          The total number of singular values found,  
*          0 <= NS <= min(M,N).
*          If RANGE = 'A', NS = min(M,N); if RANGE = 'I', NS = IU-IL+1.
*
* [out] S
*          S is DOUBLE PRECISION array, dimension (min(M,N))
*          The singular values of A, sorted so that S(i) >= S(i+1).
*
* [out] U
*          U is DOUBLE PRECISION array, dimension (LDU,UCOL)
*          If JOBU = 'V', U contains columns of U (the left singular 
*          vectors, stored columnwise) as specified by RANGE; if 
*          JOBU = 'N', U is not referenced.
*          Note: The user must ensure that UCOL >= NS; if RANGE = 'V', 
*          the exact value of NS is not known ILQFin advance and an upper 
*          bound must be used.
*
* [in] LDU
*          LDU is INTEGER
*          The leading dimension of the array U.  LDU >= 1; if
*          JOBU = 'V', LDU >= M.
*
* [out] VT
*          VT is DOUBLE PRECISION array, dimension (LDVT,N)
*          If JOBVT = 'V', VT contains the rows of V**T (the right singular 
*          vectors, stored rowwise) as specified by RANGE; if JOBVT = 'N', 
*          VT is not referenced.
*          Note: The user must ensure that LDVT >= NS; if RANGE = 'V', 
*          the exact value of NS is not known in advance and an upper 
*          bound must be used.
*
* [in] LDVT
*          LDVT is INTEGER
*          The leading dimension of the array VT.  LDVT >= 1; if
*          JOBVT = 'V', LDVT >= NS (see above).
*
* [out] WORK
*          WORK is DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK;
*
* [in] LWORK
*          LWORK is INTEGER
*          The dimension of the array WORK.
*          LWORK >= MAX(1,MIN(M,N)*(MIN(M,N)+4)) for the paths (see 
*          comments inside the code):
*             - PATH 1  (M much larger than N) 
*             - PATH 1t (N much larger than M)
*          LWORK >= MAX(1,MIN(M,N)*2+MAX(M,N)) for the other paths.
*          For good performance, LWORK should generally be larger.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
* [out] RWORK
*          RWORK is DOUBLE PRECISION array, dimension (MAX(1,LRWORK))
*          LRWORK >= MIN(M,N)*(MIN(M,N)*2+15*MIN(M,N)).
*          referenced only for complex values
*
* [out] IWORK
*          IWORK is INTEGER array, dimension (12*MIN(M,N))
*          If INFO = 0, the first NS elements of IWORK are zero. If INFO > 0, 
*          then IWORK contains the indices of the eigenvectors that failed 
*          to converge in DBDSVDX/DSTEVX.
*
* [out] INFO
*     INFO is INTEGER
*           = 0:  successful exit
*           < 0:  if INFO = -i, the i-th argument had an illegal value
*           > 0:  if INFO = i, then i eigenvectors failed to converge
*                 in DBDSVDX/DSTEVX.
*                 if INFO = N*2 + 1, an internal error occurred in
*                 DBDSVDX
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gesvdx(const char* JOBU, const char* JOBVT, const char* RANGE, i_type M, 
       i_type N, V* A, i_type LDA, 
       typename details::real_type<V>::type VL, typename details::real_type<V>::type VU, 
       i_type IL, i_type IU, i_type& NS, typename details::real_type<V>::type* S, V* U,
       i_type LDU, V* VT, i_type LDVT, V* WORK, i_type LWORK, 
       typename details::real_type<V>::type* RWORK, i_type* IWORK, i_type& INFO);

//-----------------------------------------------------------------------
//                          DBDSVDX
//-----------------------------------------------------------------------
/*
* Purpose:
* =============
*
*  DBDSVDX computes the singular value decomposition (SVD) of a real
*  N-by-N (upper or lower) bidiagonal matrix B, B = U * S * VT, 
*  where S is a diagonal matrix with non-negative diagonal elements 
*  (the singular values of B), and U and VT are orthogonal matrices 
*  of left and right singular vectors, respectively.
*
*  Given an upper bidiagonal B with diagonal D = [ d_1 d_2 ... d_N ] 
*  and superdiagonal E = [ e_1 e_2 ... e_N-1 ], DBDSVDX computes the 
*  singular value decompositon of B through the eigenvalues and 
*  eigenvectors of the N*2-by-N*2 tridiagonal matrix
*             
*        |  0  d_1                |      
*        | d_1  0  e_1            |         
*  TGK = |     e_1  0  d_2        | 
*        |         d_2  .   .     |      
*        |              .   .   . |
*
*  If (s,u,v) is a singular triplet of B with ||u|| = ||v|| = 1, then 
*  (+/-s,q), ||q|| = 1, are eigenpairs of TGK, with q = P * ( u' +/-v' ) / 
*  sqrt(2) = ( v_1 u_1 v_2 u_2 ... v_n u_n ) / sqrt(2), and 
*  P = [ e_{n+1} e_{1} e_{n+2} e_{2} ... ]. 
*
*  Given a TGK matrix, one can either a) compute -s,-v and change signs 
*  so that the singular values (and corresponding vectors) are already in 
*  descending order (as in DGESVD/DGESDD) or b) compute s,v and reorder 
*  the values (and corresponding vectors). DBDSVDX implements a) by 
*  calling DSTEVX (bisection plus inverse iteration, to be replaced 
*  with a version of the Multiple Relative Robust Representation 
*  algorithm. (See P. Willems and B. Lang, A framework for the MR^3 
*  algorithm: theory and implementation, SIAM J. Sci. Comput., 
*  35:740-766, 2013.)
*
* Arguments:
* ==========
*
* [in] UPLO
*          UPLO is CHARACTER*1
*          = 'U':  B is upper bidiagonal;
*          = 'L':  B is lower bidiagonal.
*
* [in] JOBZ
*          JOBZ is CHARACTER*1
*          = 'N':  Compute singular values only;
*          = 'V':  Compute singular values and singular vectors.
*
* [in] RANGE
*          RANGE is CHARACTER*1
*          = 'A': all singular values will be found.
*          = 'V': all singular values in the half-open interval [VL,VU)
*                 will be found.
*          = 'I': the IL-th through IU-th singular values will be found.
*
* [in] N
*          N is INTEGER
*          The order of the bidiagonal matrix.  N >= 0.
* 
* [in] D
*          D is DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the bidiagonal matrix B.
* 
* [in] E
*          E is DOUBLE PRECISION array, dimension (max(1,N-1))
*          The (n-1) superdiagonal elements of the bidiagonal matrix
*          B in elements 1 to N-1.
*
* [in] VL
*          VL is DOUBLE PRECISION
*          VL >=0.
*
* [in] VU
*         VU is DOUBLE PRECISION
*          If RANGE='V', the lower and upper bounds of the interval to
*          be searched for singular values. VU > VL.
*          Not referenced if RANGE = 'A' or 'I'.
*
* [in] IL
*          IL is INTEGER
*
* [in] IU
*          IU is INTEGER
*          If RANGE='I', the indices (in ascending order) of the
*          smallest and largest singular values to be returned.
*          1 <= IL <= IU <= min(M,N), if min(M,N) > 0.
*          Not referenced if RANGE = 'A' or 'V'.
*
* [out] NS
*          NS is INTEGER
*          The total number of singular values found.  0 <= NS <= N.
*          If RANGE = 'A', NS = N, and if RANGE = 'I', NS = IU-IL+1.
*
* [out] S
*          S is DOUBLE PRECISION array, dimension (N)
*          The first NS elements contain the selected singular values in
*          ascending order.
*
* [out] Z
*          Z is DOUBLE PRECISION array, dimension (2*N,K) )
*          If JOBZ = 'V', then if INFO = 0 the first NS columns of Z
*          contain the singular vectors of the matrix B corresponding to 
*          the selected singular values, with U in rows 1 to N and V
*          in rows N+1 to N*2, i.e.
*          Z = [ U ] 
*              [ V ]
*          If JOBZ = 'N', then Z is not referenced.    
*          Note: The user must ensure that at least K = NS+1 columns are 
*          supplied in the array Z; if RANGE = 'V', the exact value of 
*          NS is not known in advance and an upper bound must be used.
*
* [in] LDZ
*          LDZ is INTEGER
*          The leading dimension of the array Z. LDZ >= 1, and if
*          JOBZ = 'V', LDZ >= max(2,N*2).
*
* [out] WORK
*          WORK is DOUBLE PRECISION array, dimension (14*N)
*
* [out] IWORK
*          IWORK is INTEGER array, dimension (12*N)
*          If JOBZ = 'V', then if INFO = 0, the first NS elements of
*          IWORK are zero. If INFO > 0, then IWORK contains the indices 
*          of the eigenvectors that failed to converge in DSTEVX.
*
* [out] IWORK
*          INFO is INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, then i eigenvectors failed to converge
*                   in DSTEVX. The indices of the eigenvectors
*                   (as returned by DSTEVX) are stored in the
*                   array IWORK.
*                if INFO = N*2 + 1, an internal error occurred.
*
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
bdsvdx(const char* UPLO, const char* JOBZ, const char* RANGE, i_type N,
       const V* D, const V* E, V VL, V VU, i_type IL, i_type IU, i_type& NS,
       V* S, V* Z, i_type LDZ, V* WORK, i_type* IWORK, i_type& INFO );

//-----------------------------------------------------------------------
//                          DORMBR
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  If VECT = 'Q', DORMBR overwrites the general real M-by-N matrix C
*  with
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  If VECT = 'P', DORMBR overwrites the general real M-by-N matrix C
*  with
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      P * C          C * P
*  TRANS = 'T':      P**T * C       C * P**T
*
*  Here Q and P**T are the orthogonal matrices determined by DGEBRD when
*  reducing a real matrix A to bidiagonal form: A = Q * B * P**T. Q and
*  P**T are defined as products of elementary reflectors H(i) and G(i)
*  respectively.
*
*  Let nq = m if SIDE = 'L' and nq = n if SIDE = 'R'. Thus nq is the
*  order of the orthogonal matrix Q or P**T that is applied.
*
*  If VECT = 'Q', A is assumed to have been an NQ-by-K matrix:
*  if nq >= k, Q = H(1) H(2) . . . H(k);
*  if nq < k, Q = H(1) H(2) . . . H(nq-1).
*
*  If VECT = 'P', A is assumed to have been a K-by-NQ matrix:
*  if k < nq, P = G(1) G(2) . . . G(k);
*  if k >= nq, P = G(1) G(2) . . . G(nq-1).
*
*  Arguments
*  =========
*
*  VECT    (input) CHARACTER*1
*          = 'Q': apply Q or Q**T;
*          = 'P': apply P or P**T.
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q, Q**T, P or P**T from the Left;
*          = 'R': apply Q, Q**T, P or P**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q  or P;
*          = 'T':  Transpose, apply Q**T or P**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          If VECT = 'Q', the number of columns in the original
*          matrix reduced by DGEBRD.
*          If VECT = 'P', the number of rows in the original
*          matrix reduced by DGEBRD.
*          K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension
*                                (LDA,min(nq,K)) if VECT = 'Q'
*                                (LDA,nq)        if VECT = 'P'
*          The vectors which define the elementary reflectors H(i) and
*          G(i), whose products determine the matrices Q and P, as
*          returned by DGEBRD.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*          If VECT = 'Q', LDA >= max(1,nq);
*          if VECT = 'P', LDA >= max(1,min(nq,K)).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (min(nq,K))
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i) or G(i) which determines Q or P, as returned
*          by DGEBRD in the array argument TAUQ or TAUP.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q
*          or P*C or P**T*C or C*P or C*P**T.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ormbr(const char* VECT, const char* SIDE, const char* TRANS, i_type M, 
      i_type N, i_type K, const V* A, i_type LDA, const V* TAU, V* C, 
      i_type LDC, V* WORK, i_type LWORK, i_type& INFO );

//-----------------------------------------------------------------------
//                          DGELQF
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGELQF computes an LQ factorization of a real M-by-N matrix A:
*  A = L * Q.
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          The number of rows of the matrix A.  M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix A.  N >= 0.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the M-by-N matrix A.
*          On exit, the elements on and below the diagonal of the array
*          contain the m-by-min(m,n) lower trapezoidal matrix L (L is
*          lower triangular if m <= n); the elements above the diagonal,
*          with the array TAU, represent the orthogonal matrix Q as a
*          product of elementary reflectors (see Further Details).
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.  LDA >= max(1,M).
*
*  TAU     (output) DOUBLE PRECISION array, dimension (min(M,N))
*          The scalar factors of the elementary reflectors (see Further
*          Details).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.  LWORK >= max(1,M).
*          For optimum performance LWORK >= M*NB, where NB is the
*          optimal blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  The matrix Q is represented as a product of elementary reflectors
*
*     Q = H(k) . . . H(2) H(1), where k = min(m,n).
*
*  Each H(i) has the form
*
*     H(i) = I - tau * v * v'
*
*  where tau is a real scalar, and v is a real vector with
*  v(1:i-1) = 0 and v(i) = 1; v(i+1:n) is stored on exit in A(i,i+1:n),
*  and tau in TAU(i).
*
*/

template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
gelqf(i_type M, i_type N, V* A, i_type LDA, V* TAU, V* WORK, i_type LWORK,
      i_type& INFO );

//-----------------------------------------------------------------------
//                          DORMLQ
//-----------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DORMLQ overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'     SIDE = 'R'
*  TRANS = 'N':      Q * C          C * Q
*  TRANS = 'T':      Q**T * C       C * Q**T
*
*  where Q is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(k) . . . H(2) H(1)
*
*  as returned by DGELQF. Q is of order M if SIDE = 'L' and of order N
*  if SIDE = 'R'.
*
*  Arguments
*  =========
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply Q or Q**T from the Left;
*          = 'R': apply Q or Q**T from the Right.
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply Q;
*          = 'T':  Transpose, apply Q**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension
*                               (LDA,M) if SIDE = 'L',
*                               (LDA,N) if SIDE = 'R'
*          The i-th row must contain the vector which defines the
*          elementary reflector H(i), for i = 1,2,...,k, as returned by
*          DGELQF in the first k rows of its array argument A.
*          A is modified by the routine but restored on exit.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,K).
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by DGELQF.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by Q*C or Q**T*C or C*Q**T or C*Q.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If SIDE = 'L', LWORK >= max(1,N);
*          if SIDE = 'R', LWORK >= max(1,M).
*          For optimum performance LWORK >= N*NB if SIDE = 'L', and
*          LWORK >= M*NB if SIDE = 'R', where NB is the optimal
*          blocksize.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*/
template<class V> BLAS_EXPORT
typename details::enable_if_valid<void,V>::type
ormlq(const char* SIDE, const char* TRANS, i_type M, i_type N, i_type K, V* A, 
      i_type LDA, const V* TAU, V* C, i_type LDC, V* WORK, i_type LWORK, i_type& INFO);

};};

#include "matcl-blas-lapack/lapack/details/lapack_impl.inl"