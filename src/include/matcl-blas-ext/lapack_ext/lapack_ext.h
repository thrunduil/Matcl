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

#include "matcl-blas-ext/config_blas_ext.h"
#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"

#include <cstring>

namespace matcl { namespace lapack
{

// additional functions with lapack-type interface

//-----------------------------------------------------------------------
//                      UTILS
//-----------------------------------------------------------------------
//  SYNTAX:
//
//      void perm2int(int_type M,int_type *IPIV, int_type *WORK, bool ZB)
//
//  PURPOSE:
//  -----------------------------------------------------------------------
//
//  transforms permutation vector IPIV to series of interchanges
//
//  ARGUMENTS:
//  -----------------------------------------------------------------------
//  M       (input) INTEGER
//          The number of rows of the vector IPIV.  M >= 0.
//
//  IPIV    (input/output) INTEGER array, dimension M
//          On entry IPIV is a vector of permutations of numbers 0:1:M-1 (1:1:M
//          if ZB is false)
//          On exit IPIV is replaced by vector of interchanges of a vector
//          0:1:M-1 (1:1:M) that produces given vector of permutations
//
//  WORK    (input) INTEGER array, dimension at least M
//  ZB      if true then IPIV contains zero-based indices, otherwise IPIV contains
//          1-based indices
BLAS_EXT_EXPORT void perm2int(i_type M,i_type *IPIV, i_type* WORK, bool zero_based);

// transforms series of interchanges to permutation vector IPIV 
BLAS_EXT_EXPORT void int2perm(i_type M,i_type *IPIV, i_type* WORK, bool zero_based);

//--------------------------------------------------------------------------------
//                              laswpc
//--------------------------------------------------------------------------------
/* Purpose
*  =======
*
*  DLASWPC performs a series of column interchanges on the matrix A.
*  One column interchange is initiated for each of columns K1 through K2 of A.
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The number of rows of the matrix A.
*
*  A       (input/output) DOUBLE PRECISION array, dimension (N,K)
*          On entry, the matrix of row dimension N to which the column
*          interchanges will be applied.
*          On exit, the permuted matrix.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A.
*
*  K1      (input) INTEGER
*          The first element of IPIV for which a column interchange will
*          be done.
*
*  K2      (input) INTEGER
*          The last element of IPIV for which a column interchange will
*          be done.
*
*  IPIV    (input) INTEGER array, dimension (K2*abs(INCX))
*          The vector of pivot indices.  Only the elements in positions
*          K1 through K2 of IPIV are accessed.
*          IPIV(K) = L implies columns K and L are to be interchanged.
*
*  INCX    (input) INTEGER
*          The increment between successive values of IPIV.  If IPIV
*          is negative, the pivots are applied in reverse order.
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
laswpc(i_type n, V *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx);

BLAS_EXT_EXPORT void    claswpc(i_type n, c_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx);
BLAS_EXT_EXPORT void    slaswpc(i_type n, s_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx);
BLAS_EXT_EXPORT void    dlaswpc(i_type n, d_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx);
BLAS_EXT_EXPORT void    zlaswpc(i_type n, z_type *a,i_type lda,i_type k1,i_type k2,const i_type *ipiv,i_type incx);

//--------------------------------------------------------------------------------
//                              GETF2R
//--------------------------------------------------------------------------------
// Purpose
// =======
//
// GETF2R computes an LU factorization of a general m-by-n matrix A
// using rook pivoting with row and column interchanges.
//
// The factorization has the form
//    A = P * L * U * Q
// where P, Q are permutation matrices, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
//
// This is the Level 2 BLAS version of the algorithm.
//
// Arguments
// =========
//
// M       (input) INTEGER
//         The number of rows of the matrix A.  M >= 0.
//
// N       (input/output) INTEGER
//         The number of columns of the matrix A.  N >= 0.
//         On exit, number of nonzero columns
//
// JB      (input) INTEGER
//         The number of columns of the matrix A in the working block.  
//         JB >= 0.
//
// J0      (input) INTEGER
//         Position of the current block in the large matrix.
//
// A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
//         On entry, the m by n matrix to be factored.
//         On exit, the factors L and U from the factorization
//         A = P*L*U; the unit diagonal elements of L are not stored.
//
// LDA     (input) INTEGER
//         The leading dimension of the array A.  LDA >= max(1,M).
//
// IPIV    (output) INTEGER array, dimension (min(M,N))
//         The pivot indices; for 1 <= i <= min(M,N), row i of the
//         matrix was interchanged with row IPIV(i).
//
// IQIV    (input/output) INTEGER array, dimension (N)
//         The column permutation vector; for 1 <= i <= N, column i 
//         is given by IQIV(i) column of the matrix A, on input IQIV
//         gives starting permutation vector, for example 1:N
//
// TOLC    (input) DOUBLE PRECISION
//         tolerace, given column is selected as a pivot if element 
//         is greater than 1/TOLC * max element if given row, 
//         0 <= TOLC <= 1
//
// TOLR    (input) DOUBLE PRECISION
//         tolerace, given row is selected as a pivot if element 
//         is greater than 1/TOLR * max element if given column, 
//         0 <= TOLR <= 1
//
// TOLV    (input) DOUBLE PRECISION
//         tolerace, if absolute value of pivot is less than TOLV,
//         then all elements in given column are changed to zero, column 
//         is moved to the right end of the matrix and marked as singular
//         If TOLV is negative, then TOLV = 10 * U * |A|_F / sqrt(min(M,N)).
//
// WORK    (input/output) working array of size in bytes WORK_SIZE
//         WORK_SIZE is returned by this routine is is called with 
//         INFO[0] = -1
//
// LWORK   (input) INTEGER
//         If LWORK = -1, then a workspace query is assumed; the routine
//         only calculates the optimal size of the WORK array, returns
//         this value as the first entry of the WORK array.
//
// INFO    (output) INTEGER
//         = 0: successful exit
//         < 0: if INFO = -k, the k-th argument had an illegal value
//
template<class T> 
void getf2r(i_type M, i_type* N, i_type JB, i_type J, T *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
            typename lapack::details::real_type<T>::type TOLC,
            typename lapack::details::real_type<T>::type TOLR, 
            typename lapack::details::real_type<T>::type TOLV,
            T* WORK, i_type LWORK, i_type *INFO );

//--------------------------------------------------------------------------------
//                              GETRFR
//--------------------------------------------------------------------------------
// Purpose
// =======
//
// DGETRFR computes an LU factorization of a general M-by-N matrix A
// using rook pivoting
//
// The factorization has the form
//    A = P * L * U * Q
// where P, Q are permutation matrices, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
//
// This is the right-looking Level 3 BLAS version of the algorithm.
//
// Arguments
// =========
//
// M       (input) INTEGER
//         The number of rows of the matrix A.  M >= 0.
//
// N       (input) INTEGER
//         The number of columns of the matrix A.  N >= 0.
//
// A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
//         On entry, the M-by-N matrix to be factored.
//         On exit, the factors L and U from the factorization
//         A = P*L*U; the unit diagonal elements of L are not stored.
//
// LDA     (input) INTEGER
//         The leading dimension of the array A.  LDA >= max(1,M).
//
// IPIV    (output) INTEGER array, dimension (min(M,N))
//         The row pivot indices; for 1 <= i <= min(M,N), row i of the
//         matrix was interchanged with row IPIV(i).
//
// IQIV    (input/output) INTEGER array, dimension (N)
//         The column permutation vector; for 1 <= i <= N, column i 
//         is given by IQIV(i) column of the matrix A. On entry IQIV gives
//         initial column permutation, for example 1:N.
//
// TOLC    (input) DOUBLE PRECISION
//         tolerace, given column is selected as a pivot if element 
//         is greater than 1/TOLC * max element if given row, 
//         0 <= TOLC <= 1
//
// TOLR    (input) DOUBLE PRECISION
//         tolerace, given row is selected as a pivot if element 
//         is greater than 1/TOLR * max element if given column, 
//         0 <= TOLR <= 1
//
// TOLV    (input) DOUBLE PRECISION
//         tolerace, if absolute value of pivot is less than TOLV,
//         then all elements in given column are changed to zero, column 
//         is moved to the right end of the matrix and marked as singular
//         If TOLV is negative, then TOLV = 10 * U * |A|_F / sqrt(min(M,N)).
//
// WORK    (workspace/output) DOUBLE PRECISION array,
//         dimension (MAX(1,LWORK)) 
//         On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
//
// LWORK   (input) INTEGER
//         If LWORK = -1, then a workspace query is assumed; the routine
//         only calculates the optimal size of the WORK array, returns
//         this value as the first entry of the WORK array.
//
// INFO    (output) INTEGER
//         >= 0: successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value
//         >=0:  INFO is estimated rank.
template<class T> BLAS_EXT_EXPORT
void getrfr(i_type M, i_type N, T *A, i_type LDA, i_type *IPIV, i_type *IQIV, 
            typename lapack::details::real_type<T>::type TOLC,
            typename lapack::details::real_type<T>::type TOLR, 
            typename lapack::details::real_type<T>::type TOLV,
            T* WORK, i_type LWORK, i_type *INFO);

//--------------------------------------------------------------------------------
//                              GETRFC
//--------------------------------------------------------------------------------
// Purpose
// =======
//
// DGETRFC computes an LU factorization of a general M-by-N matrix A
// using complete pivoting
//
// The factorization has the form
//    A = P * L * U * Q
// where P, Q are permutation matrices, L is lower triangular with unit
// diagonal elements (lower trapezoidal if m > n), and U is upper
// triangular (upper trapezoidal if m < n).
//
// Arguments
// =========
//
// M       (input) INTEGER
//         The number of rows of the matrix A.  M >= 0.
//
// N       (input) INTEGER
//         The number of columns of the matrix A.  N >= 0.
//
// A       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
//         On entry, the M-by-N matrix to be factored.
//         On exit, the factors L and U from the factorization
//         A = P*L*U; the unit diagonal elements of L are not stored.
//
// LDA     (input) INTEGER
//         The leading dimension of the array A.  LDA >= max(1,M).
//
// IPIV    (output) INTEGER array, dimension (min(M,N))
//         The row pivot indices; for 1 <= i <= min(M,N), row i of the
//         matrix was interchanged with row IPIV(i).
//
// IQIV    (input/output) INTEGER array, dimension (N)
//         The column permutation vector; for 1 <= i <= N, column i 
//         is given by IQIV(i) column of the matrix A. On entry IQIV gives
//         initial column permutation, for example 1:N.
//
// TOLV    (input) DOUBLE PRECISION
//         tolerace, if absolute value of pivot is less than TOLV,
//         then all elements in given column are changed to zero, column 
//         is moved to the right end of the matrix and marked as singular
//         If TOLV is negative, then TOLV = 10 * U * |A|_F / sqrt(min(M,N)).
//
// INFO    (output) INTEGER
//         >= 0: successful exit
//         < 0:  if INFO = -i, the i-th argument had an illegal value
//         >=0:  INFO is estimated rank.
template<class T> BLAS_EXT_EXPORT
void getrfc(i_type M, i_type N, T *A, i_type LDA, i_type *IPIV, i_type *IQIV, 
            typename lapack::details::real_type<T>::type TOLV, i_type *INFO);

// GETF2 computes an LU factorization of a general m-by-n matrix A
// using partial pivoting with row interchanges. Level 2 algorithm
template<class T> BLAS_EXT_EXPORT
void getf2(i_type M, i_type N, T *A, i_type LDA, i_type* IPIV, i_type *INFO );


// getrf_rl computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges. Level 3 right-looking
// algorithm
template<class T> BLAS_EXT_EXPORT
void getrf_rl(i_type M, i_type N, T *A, i_type LDA, i_type *IPIV, i_type *INFO);

// getrf_rl computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges. Level 3 left-looking
// algorithm
template<class T> BLAS_EXT_EXPORT
void getrf_ll(i_type M,i_type N, T* A,i_type LDA,i_type* IPIV,i_type* INFO);

// getrf_rl computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges. Level 3 Crout algorithm
template<class T> BLAS_EXT_EXPORT
void getrf_cr(i_type M, i_type N, T* A,i_type LDA,i_type* IPIV,i_type* INFO);

// getrf_rl computes an LU factorization of a general M-by-N matrix A
// using partial pivoting with row interchanges. Level 3 recursive algorithm
template<class T> BLAS_EXT_EXPORT
void getrf_rec(i_type M, i_type N, T* A,i_type LDA,i_type* IPIV,i_type* INFO );

/*
* Purpose
* =======
*
* potfp2 computes the Cholesky factorization with complete
* pivoting of a real symmetric positive semi-definite matrix A.
*
* The factorization has the form
* P' * A * P = U' * U , if UPLO = 'U',
* P' * A * P = L * L', if UPLO = 'L',
* where U is an upper triangular matrix and L is lower triangular, and
* P is stored as vector PIV.
*
* This algorithm does not attempt to check that A is positive
* semi-definite. This version of the algorithm calls level 2 BLAS.
*
* Arguments
* =========
*
* UPLO (input) CHARACTER*1
* Specifies whether the upper or lower triangular part of the
* symmetric matrix A is stored.
* = 'U': Upper triangular
* = 'L': Lower triangular
*
* N (input) INTEGER
* The order of the matrix A. N >= 0.
*
* A (input/output) DOUBLE PRECISION array, dimension (LDA,N)
* On entry, the symmetric matrix A. If UPLO = 'U', the leading
* n by n upper triangular part of A contains the upper
* triangular part of the matrix A, and the strictly lower
* triangular part of A is not referenced. If UPLO = 'L', the
* leading n by n lower triangular part of A contains the lower
* triangular part of the matrix A, and the strictly upper
* triangular part of A is not referenced.
*
* On exit, if INFO = 0, the factor U or L from the Cholesky
* factorization as above.
*
* PIV (output) INTEGER array, dimension (N)
* PIV is such that the nonzero entries are P( PIV( K ), K ) = 1.
*
* RANK (output) INTEGER
* The rank of A given by the number of steps the algorithm
* completed.
*
* TOL (input) DOUBLE PRECISION
* User defined tolerance. If TOL < 0, then 10 * U * |A|_F / sqrt(N)
* will be used, otherwise TOL. The algorithm 
* terminates at the (k-1)st step if the pivot <= TOL.
*
* NORM (input) DOUBLE PRECISION
* Upper bound on the second norm of the matrix A. If NORM < 1, then
* NORM = |A|_F / sqrt(N)
*
* LDA (input) INTEGER
* The leading dimension of the array A. LDA >= max(1,N).
*
* WORK DOUBLE PRECISION array, dimension (2*N)
* Work space.
*
* INFO (output) INTEGER
* < 0:  if INFO = -k, the k-th argument had an illegal value
* = 0   algorithm completed successfully.
*/
template<class T> BLAS_EXT_EXPORT
void potfp2(const char * uplo, i_type N, T *A, i_type LDA, i_type* PIV, 
                    i_type& RANK, typename details::real_type<T>::type TOL, 
                    typename details::real_type<T>::type* WORK, i_type& INFO);

/*
*
* Purpose
* =======
*
* potfp3 computes the Cholesky factorization with complete
* pivoting of a real symmetric positive semi-definite matrix A.
*
* The factorization has the form
* P' * A * P = U' * U , if UPLO = 'U',
* P' * A * P = L * L', if UPLO = 'L',
* where U is an upper triangular matrix and L is lower triangular, and
* P is stored as vector PIV.
*
* This algorithm does not attempt to check that A is positive
* semi-definite. This version of the algorithm calls level 3 BLAS.
*
* Arguments
* =========
*
* UPLO (input) CHARACTER*1
* Specifies whether the upper or lower triangular part of the
* symmetric matrix A is stored.
* = 'U': Upper triangular
* = 'L': Lower triangular
*
* N (input) INTEGER
* The order of the matrix A. N >= 0.
*
* A (input/output) DOUBLE PRECISION array, dimension (LDA,N)
* On entry, the symmetric matrix A. If UPLO = 'U', the leading
* n by n upper triangular part of A contains the upper
* triangular part of the matrix A, and the strictly lower
* triangular part of A is not referenced. If UPLO = 'L', the
* leading n by n lower triangular part of A contains the lower
* triangular part of the matrix A, and the strictly upper
* triangular part of A is not referenced.
*
* On exit, if INFO = 0, the factor U or L from the Cholesky
* factorization as above.
*
* PIV (output) INTEGER array, dimension (N)
* PIV is such that the nonzero entries are P( PIV( K ), K ) = 1.
*
* RANK (output) INTEGER
* The rank of A given by the number of steps the algorithm
* completed.
*
* TOL (input) DOUBLE PRECISION
* User defined tolerance. If TOL < 0, then 10 * U * NORM_F
* will be used, otherwise TOL. The algorithm 
* terminates at the (k-1)st step if the pivot <= TOL.
*
* LDA (input) INTEGER
* The leading dimension of the array A. LDA >= max(1,N).
*
* WORK DOUBLE PRECISION array, dimension (2*N)
* Work space.
*
* INFO (output) INTEGER
* < 0:  if INFO = -k, the k-th argument had an illegal value
* = 0   algorithm completed successfully.
*
*/
template<class T> BLAS_EXT_EXPORT
void potfp3(const char * uplo, i_type N, T *A, i_type LDA, i_type* PIV, 
                    i_type& RANK, typename details::real_type<T>::type TOL, 
                    typename details::real_type<T>::type* WORK, i_type& INFO);

//--------------------------------------------------------------------------------
//                              GETRS_REV
//--------------------------------------------------------------------------------

/*
*  Purpose
*  =======
*
*  DGETRS solves a system of linear equations
*     X * A = B  or  X * A' = B
*  with a general N-by-N matrix A using the LU factorization computed
*  by DGETRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations:
*          = 'N':  X * A  = B  (No transpose)
*          = 'T':  X * A' = B  (Transpose)
*          = 'C':  X * A' = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of rows
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
*  B       (input/output) DOUBLE PRECISION array, dimension (NRHS, LDB)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,nrhs).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
getrs_rev(const char *trans,i_type n,i_type nrhs,const V *a,i_type lda,const i_type *ipiv,V *b,i_type ldb,
                        i_type *info);

BLAS_EXT_EXPORT void    cgetrs_rev(const char *trans,i_type n,i_type nrhs,const c_type *a,i_type lda,
                            const i_type *ipiv,c_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    sgetrs_rev(const char *trans,i_type n,i_type nrhs,const s_type *a,i_type lda,
                           const i_type *ipiv,s_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    dgetrs_rev(const char *trans,i_type n,i_type nrhs,const d_type *a,i_type lda,
                           const i_type *ipiv,d_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    zgetrs_rev(const char *trans,i_type n,i_type nrhs,const z_type *a,i_type lda,
                           const i_type *ipiv,z_type *b,i_type ldb,i_type *info);

//--------------------------------------------------------------------------------
//                              GBTRS_REV
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGBTRS solves a system of linear equations
*     X * A = B  or  X * A' = B
*  with a general band matrix A using the LU factorization computed
*  by DGBTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  X * A  = B  (No transpose)
*          = 'T':  X * A' = B  (Transpose)
*          = 'C':  X * A' = B  (Conjugate transpose = Transpose)
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
*          The number of right hand sides, i.e., the number of rows
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
*  B       (input/output) DOUBLE PRECISION array, dimension (NRHS, LDB)
*          On entry, the right hand side matrix B.
*          On exit, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,nrhs).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const V *ab,
                           i_type ldab,const i_type *ipiv,V *b,i_type ldb,i_type *info);

BLAS_EXT_EXPORT void    cgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const c_type *ab,
                           i_type ldab,const i_type *ipiv,c_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    sgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const s_type *ab,
                           i_type ldab,const i_type *ipiv,s_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    dgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const d_type *ab,
                           i_type ldab,const i_type *ipiv,d_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    zgbtrs_rev(const char *trans,i_type n,i_type kl,i_type ku,i_type nrhs,const z_type *ab,
                           i_type ldab,const i_type *ipiv,z_type *b,i_type ldb,i_type *info);

//--------------------------------------------------------------------------------
//                              TRTRS_REV
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTRTRS solves a triangular system of the form
*
*     X * A = B  or  X * A**T = B,
*
*  where A is a triangular matrix of order N, and B is an NRHS-by-N
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
*          = 'N':  X * A    = B  (No transpose)
*          = 'T':  X * A**T = B  (Transpose)
*          = 'C':  X* A**H  = B  (Conjugate transpose = Transpose)
*
*  DIAG    (input) CHARACTER*1
*          = 'N':  A is non-unit triangular;
*          = 'U':  A is unit triangular.
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of rows
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
*  B       (input/output) DOUBLE PRECISION array, dimension (NRHS, LDB)
*          On entry, the right hand side matrix B.
*          On exit, if INFO = 0, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,NRHS).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0: if INFO = i, the i-th diagonal element of A is zero,
*               indicating that the matrix is singular and the solutions
*               X have not been computed.
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
trtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const V *a,i_type lda,V *b,i_type ldb,i_type *info);

BLAS_EXT_EXPORT void    dtrtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const d_type *a,i_type lda,d_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    ctrtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const c_type *a,i_type lda,c_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    strtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const s_type *a,i_type lda,s_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    ztrtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type nrhs,
                           const z_type *a,i_type lda,z_type *b,i_type ldb,i_type *info);

//--------------------------------------------------------------------------------
//                              TBTRS_REV
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTBTRS solves a triangular system of the form
*
*     X * A = B  or  X * A**T = B,
*
*  where A is a triangular band matrix of order N, and B is an
*  NRHS-by-N matrix.  A check is made to verify that A is nonsingular.
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
*          The number of right hand sides, i.e., the number of rows
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
*  B       (input/output) DOUBLE PRECISION array, dimension (NRHS,LDB)
*          On entry, the right hand side matrix B.
*          On exit, if INFO = 0, the solution matrix X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,NRHS).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*          > 0:  if INFO = i, the i-th diagonal element of A is zero,
*                indicating that the matrix is singular and the
*                solutions X have not been computed.
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
tbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const V *ab,i_type ldab,V *b,i_type ldb,i_type *info);

BLAS_EXT_EXPORT void    ctbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs, 
                           const c_type *ab,i_type ldab,c_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    stbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const s_type *ab,i_type ldab,s_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    dtbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs,
                           const d_type *ab,i_type ldab,d_type *b,i_type ldb,i_type *info);
BLAS_EXT_EXPORT void    ztbtrs_rev(const char *uplo,const char *trans,const char *diag,i_type n,i_type kd,i_type nrhs, 
                           const z_type *ab,i_type ldab,z_type *b,i_type ldb,i_type *info);

//--------------------------------------------------------------------------------
//                              DTGSEN3
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTGSEN3 reorders the generalized real Schur decomposition of a real
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
*  DTGSEN3 also computes the generalized eigenvalues
*
*              w(j) = (ALPHAR(j) + i*ALPHAI(j))/BETA(j)
*
*  of the reordered matrix pair (A, B).
*
*  Arguments
*  =========
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
*  K       (input) INTEGER
*          Number of rows of the matrix Q and Z. K >= 0
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
*          and if WANTQ = .TRUE., LDQ >= K.
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
*          If WANTZ = .TRUE., LDZ >= K.
*
*  M       (output) INTEGER
*          The dimension of the specified pair of left and right eigen-
*          spaces (deflating subspaces). 0 <= M <= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array,
*          dimension (MAX(1,LWORK)) 
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
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
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
tgsen3(i_type wantq,i_type wantz,const i_type *select,i_type n, i_type k,
        V *a,i_type lda, V *b,i_type ldb, typename details::complex_type<V>::type *alpha, V *beta, 
        V *q,i_type ldq, V *z, i_type ldz,i_type *m, V *work,i_type lwork, i_type *info);

//--------------------------------------------------------------------------------
//                              DTRSEN3
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DTRSEN3 reorders the real Schur factorization of a real matrix
*  A * Q = Q * T, so that a selected cluster of eigenvalues appears in
*  the leading diagonal blocks of the upper quasi-triangular matrix T,
*  and the leading columns of Q form an orthonormal basis of the
*  corresponding right invariant subspace.
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elemnts equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  COMPQ   (input) CHARACTER*1
*          = 'V': update the matrix Q of Schur vectors;
*          = 'N': do not update Q.
*
*  SELECT  (input) LOGICAL array, dimension (N)
*          SELECT specifies the eigenvalues in the selected cluster. To
*          select a real eigenvalue w(j), SELECT(j) must be set to
*          .TRUE.. To select a complex conjugate pair of eigenvalues
*          w(j) and w(j+1), corresponding to a 2-by-2 diagonal block,
*          either SELECT(j) or SELECT(j+1) or both must be set to
*          .TRUE.; a complex conjugate pair of eigenvalues must be
*          either both included in the cluster or both excluded.
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  K       (input) INTEGER
*          Number of rows of the matrix Q. K >= 0
*
*  T       (input/output) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          canonical form.
*          On exit, T is overwritten by the reordered matrix T, again in
*          Schur canonical form, with the selected eigenvalues in the
*          leading diagonal blocks.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if COMPQ = 'V', the matrix Q of Schur vectors.
*          On exit, if COMPQ = 'V', Q has been postmultiplied by the
*          orthogonal transformation matrix which reorders T; the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If COMPQ = 'N', Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1; and if COMPQ = 'V', LDQ >= max(1,K).
*
*  W       (output) DOUBLE PRECISION COMPLEX array, dimension (N)
*          The real WR and imaginary WI parts, respectively, of the reordered
*          eigenvalues of T. The eigenvalues are stored in the same
*          order as on the diagonal of T, with WR(i) = T(i,i) and, if
*          T(i:i+1,i:i+1) is a 2-by-2 diagonal block, WI(i) > 0 and
*          WI(i+1) = -WI(i). Note that if a complex eigenvalue is
*          sufficiently ill-conditioned, then its value may differ
*          significantly from its value before reordering.
*
*  M       (output) INTEGER
*          The dimension of the specified invariant subspace.
*          0 < = M <= N.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: reordering of T failed because some eigenvalues are too
*               close to separate (the problem is very ill-conditioned);
*               T may have been partially reordered, and WR and WI
*               contain the eigenvalues in the same order as in T; S and
*               SEP (if requested) are set to zero.
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
trsen3(const char* compq, const i_type *select, i_type n, i_type k, V* t, i_type ldt, 
       V* q, i_type ldq, typename details::complex_type<V>::type * w, i_type* m, V* work, 
       i_type lwork, i_type* info);

//--------------------------------------------------------------------------------
//                              qtriu2schur
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  qtriu2schur transform upper quasi-triangular real matrix to upper
*  quasi-triangular matrix in the canonical schur form, and do nothing
*  for complex matrices
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input/output) array, dimension (LDT,N). On entry, the upper 
*          quasi-triangular matrix T. On exit, T is overwritten by upper
*          quasi-traingular matrix in the schur canonical form
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*  wantq   (input) boolean
*          = true: update the matrix Q of Schur vectors;
*          = false: do not update Q.
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if wantq = true, the matrix Q of Schur vectors.
*          On exit, if COMPQ = true, Q has been postmultiplied by the
*          orthogonal transformation matrix which reorders T; the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If wantq = false, Q is not referenced.
*
*  QR      (input) INTEGER
*          Number of rows of Q
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1; and if wantq = wantq, LDQ >= QR.
*
*/

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
qtriu2schur(i_type n, V* T, i_type ldt, bool wantq, V* Q, i_type QR, i_type ldq);

//--------------------------------------------------------------------------------
//                              qtriu2gschur
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  qtriu2gschur transform a matrix pair (T,S) where T is upper quasi-uppertriangular 
*  real matrix and S is upper triangular matrix to a matrix pair in canonical 
*  generalized  schur form, and do nothing for complex matrices
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  M       (input) INTEGER
*          Number of rows of Q and Z factors
*
*  T       (input/output) array, dimension (LDT,N). On entry, the upper 
*          quasi-triangular matrix T. On exit, T is overwritten by upper
*          quasi-traingular matrix in the generalized schur canonical form
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  S       (input/output) array, dimension (LDS,N). On entry, the upper 
*          triangular matrix T. On exit, S is overwritten by upper
*          traingular matrix in the generalized schur canonical form
*
*  LDS     (input) INTEGER
*          The leading dimension of the array S. LDS >= max(1,N).
*
*  wantq   (input) boolean
*          = true: update the matrix Q of Schur vectors;
*          = false: do not update Q.
*
*  wantz   (input) boolean
*          = true: update the matrix Z of Schur vectors;
*          = false: do not update Z.
*
*  Q       (input/output) DOUBLE PRECISION array, dimension (LDQ,N)
*          On entry, if wantq = true, the matrix Q of Schur vectors.
*          On exit, if wantq = true, Q has been postmultiplied by the
*          orthogonal transformation matrix which reorders (T,S); the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If wantq = false, Q is not referenced.
*
*  LDQ     (input) INTEGER
*          The leading dimension of the array Q.
*          LDQ >= 1; and if wantq = wantq, LDQ >= M.
*
*  Z       (input/output) DOUBLE PRECISION array, dimension (LDZ,N)
*          On entry, if wantz = true, the matrix Z of Schur vectors.
*          On exit, if wantz = true, Z has been postmultiplied by the
*          orthogonal transformation matrix which reorders (T,S); the
*          leading M columns of Q form an orthonormal basis for the
*          specified invariant subspace.
*          If wantz = false, Z is not referenced.
*
*  LDZ     (input) INTEGER
*          The leading dimension of the array Z.
*          LDZ >= 1; and if wantq = wantz, LDZ >= M.
*
*/

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
qtriu2gschur(i_type n, i_type m, V* T, i_type ldt, V* S, i_type lds, bool wantq, bool wantz, 
             V* Q, i_type ldq, V* Z, i_type ldz);

//--------------------------------------------------------------------------------
//                              QTRTRS
//--------------------------------------------------------------------------------
/*
*
*  Purpose
*  =======
*
*  QTRTRS solves a triangular system of the form:
*
*     op(A)*X  = C
*
*  where op(A) = A or A**T, and  A is upper quasi-triangular, M-by-M; 
*  the right hand side C and the solution X are M-by-N.
*
*  Arguments
*  =========
*
*  TRANA   (input) CHARACTER*1
*          Specifies the option op(A):
*          = 'N': op(A) = A    (No transpose)
*          = 'T': op(A) = A**T (Transpose)
*          = 'C': op(A) = A**H (Conjugate transpose = Transpose)
*
*  M       (input) INTEGER
*          The order of the matrix A, and the number of rows in the
*          matrices X and C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns in the matrices X and C. N >= 0.
*
*  A       (input) DOUBLE PRECISION array, dimension (LDA,M)
*          The upper quasi-triangular matrix A, in Schur canonical form.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,M).
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N right hand side matrix C.
*          On exit, C is overwritten by the solution matrix X.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M)
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          > 0, if INFO = k, then  k-th element on diagonal of A is close
*               to zero
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
qtrtrs(const char* TRANA, i_type M, i_type N, const V* A, i_type LDA, 
            V* C, i_type LDC, i_type* INFO);

//--------------------------------------------------------------------------------
//                              CHUDL
//--------------------------------------------------------------------------------
/* Purpose
*  =======
*     Given a cholesky decomposition A=R R^T of some matrix A, compute
*     the decomposition of B = A + sigma * X * X^T as a function of R and the
*     vector X.
*
*  Arguments:
*  =========
*      R     double precision(ldr,p), denoting the lower triangular
*            cholesky decomposition of positive definite matrix A. After 
*			 completing this function call the matrix will be overwritten with 
*            the updated cholesky.
*
*      LDR   integer, leading dimension of R.
*
*      P     integer, order of the matrix R.
*
*      SIGMA real scalar
*
*      X     double precision(p), vector used for the rank one update.
*
*      INCX  specifies the increment for the elements of X. INCX must not be zero.
*
*      INFO  (output) INTEGER
*            = 0: successful exit
*            < 0: if INFO = -i, the i-th argument had an illegal value
*            > 0: if INFO = +i, R(i,i) becomes zero, updated matrix is not 
*                 positive definite
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
chudl(V* R, i_type LDR, i_type p, const typename details::real_type<V>::type& SIGMA, V* X, i_type INCX, 
      i_type* INFO);

BLAS_EXT_EXPORT void schudl(s_type* R, i_type LDR, i_type p, const s_type& SIGMA, s_type* X, i_type INCX, 
                        i_type* INFO);
BLAS_EXT_EXPORT void dchudl(d_type* R, i_type LDR, i_type p, const d_type& SIGMA, d_type* X, i_type INCX, 
                        i_type* INFO);
BLAS_EXT_EXPORT void cchudl(c_type* R, i_type LDR, i_type p, const s_type& SIGMA, c_type* X, i_type INCX, 
                        i_type* INFO);
BLAS_EXT_EXPORT void zchudl(z_type* R, i_type LDR, i_type p, const d_type& SIGMA, z_type* X, i_type INCX, 
                        i_type* INFO);

//--------------------------------------------------------------------------------
//                              CHUDU
//--------------------------------------------------------------------------------
/* Purpose
*  =======
*     Given a cholesky decomposition A=R^T R of some matrix A, compute
*     the decomposition of B = A + sigma * X * X^T as a function of R and the
*     vector x.
*
*  Arguments:
*  =========
*      R     double precision(ldr,p), denoting the upper triangular
*            cholesky decomposition of the matrix A. After completing
*            this function call the matrix will be overwritten with
*            the updated cholesky.
*
*      LDR   integer, leading dimension of R.
*
*      P     integer, order of the matrix R.
*
*      SIGMA real scalar
*
*      X     double precision(p), vector used for the rank one update.
*
*      INCX  specifies the increment for the elements of X. INCX must not be zero.
*
*      WORK  double precision array of length 2 * p.
*
*      INFO  (output) INTEGER
*            = 0: successful exit
*            < 0: if INFO = -i, the i-th argument had an illegal value
*            > 0: if INFO = +i, R(i,i) becomes zero, updated matrix is not 
*                 positive definite
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
chudu(V* R, i_type LDR, i_type p, const typename details::real_type<V>::type& SIGMA, V* X, i_type INCX, 
    V* WORK, i_type* INFO);

BLAS_EXT_EXPORT void schudu(s_type* R, i_type LDR, i_type p, const s_type& SIGMA, s_type* X, i_type INCX, 
                        s_type* WORK, i_type* INFO);
BLAS_EXT_EXPORT void dchudu(d_type* R, i_type LDR, i_type p, const d_type& SIGMA, d_type* X, i_type INCX, 
                        d_type* WORK, i_type* INFO);
BLAS_EXT_EXPORT void cchudu(c_type* R, i_type LDR, i_type p, const s_type& SIGMA, c_type* X, i_type INCX, 
                        c_type* WORK, i_type* INFO);
BLAS_EXT_EXPORT void zchudu(z_type* R, i_type LDR, i_type p, const d_type& SIGMA, z_type* X, i_type INCX, 
                        z_type* WORK, i_type* INFO);

//--------------------------------------------------------------------------------
//                              trsen_cond
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  trsen_cond computes the reciprocal condition numbers of
*  the cluster of eigenvalues and/or the invariant subspace.
*
*  T must be in Schur canonical form (as returned by DHSEQR), that is,
*  block upper triangular with 1-by-1 and 2-by-2 diagonal blocks; each
*  2-by-2 diagonal block has its diagonal elemnts equal and its
*  off-diagonal elements of opposite sign.
*
*  Arguments
*  =========
*
*  JOB     (input) CHARACTER*1
*          Specifies whether condition numbers are required for the
*          cluster of eigenvalues (S) or the invariant subspace (SEP):
*          = 'E': for eigenvalues only (S);
*          = 'V': for invariant subspace only (SEP);
*          = 'B': for both eigenvalues and invariant subspace (S and
*                 SEP).
*
*  M       (input) INTEGER
*          Number of eigenvalues in invariant subspace. If M-th eigenvalue
*          is complex, then the invariant subspace contains M+1 eigenvalues
*          0 <= M <= N
*
*  N       (input) INTEGER
*          The order of the matrix T. N >= 0.
*
*  T       (input) DOUBLE PRECISION array, dimension (LDT,N)
*          On entry, the upper quasi-triangular matrix T, in Schur
*          canonical form.
*
*  LDT     (input) INTEGER
*          The leading dimension of the array T. LDT >= max(1,N).
*
*  S       (output) DOUBLE PRECISION
*          If JOB = 'E' or 'B', S is a lower bound on the reciprocal
*          condition number for the selected cluster of eigenvalues.
*          S cannot underestimate the true reciprocal condition number
*          by more than a factor of sqrt(N). If M = 0 or N, S = 1.
*          If JOB = 'V', S is not referenced.
*
*  SEP     (output) DOUBLE PRECISION
*          If JOB = 'V' or 'B', SEP is the estimated reciprocal
*          condition number of the specified invariant subspace. If
*          N = 0, SEP = norm(T).
*          If JOB = 'E', SEP is not referenced.
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          if JOB = 'E', LWORK >= max(1,M*(N-M));
*          if JOB = 'V' or 'B', LWORK >= max(1,2*M*(N-M)).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          if JOB = 'V' or 'B', LIWORK >= max(1,M*(N-M)).
*
*          If LIWORK = -1, then a workspace query is assumed; the
*          routine only calculates the optimal size of the IWORK array,
*          returns this value as the first entry of the IWORK array, and
*          no error message related to LIWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
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
*
*  =====================================================================
*
* Pawel Kowal:
*   This is a direct rewrite of fragments Lapack's TRSEN function responsible
*   for estimating condition numbers
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
trsen_cond(const char* JOB, i_type N, i_type M, const V* T, i_type LDT, 
           typename details::real_type<V>::type& S, typename details::real_type<V>::type& SEP, 
           V* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO );

//--------------------------------------------------------------------------------
//                              tgsen_cond
//--------------------------------------------------------------------------------
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
*           =1: Reciprocal of norms of "projections" onto left and right
*               eigenspaces w.r.t. the selected cluster (PL and PR).
*           =2: Upper bounds on Difu and Difl. F-norm-based estimate
*               (DIF(1:2)).
*           =3: Estimate of Difu and Difl. 1-norm-based estimate
*               (DIF(1:2)).
*               About 5 times as expensive as IJOB = 2.
*           =4: Compute PL, PR and DIF (i.e. 1 and 2 above): Economic
*               version to get it all.
*           =5: Compute PL, PR and DIF (i.e. 1 and 3 above)
*
*  N       (input) INTEGER
*          The order of the matrices A and B. N >= 0.
*
*  M       (input) INTEGER
*          Number of eigenvalues in invariant subspace. If M-th eigenvalue
*          is complex, then the invariant subspace contains M+1 eigenvalues
*          0 <= M <= N
*
*  A       (input) DOUBLE PRECISION array, dimension(LDA,N)
*          On entry, the upper quasi-triangular matrix A, with (A, B) in
*          generalized real Schur canonical form.
*
*  LDA     (input) INTEGER
*          The leading dimension of the array A. LDA >= max(1,N).
*
*  B       (input) DOUBLE PRECISION array, dimension(LDB,N)
*          On entry, the upper triangular matrix B, with (A, B) in
*          generalized real Schur canonical form.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B. LDB >= max(1,N).
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
*          The dimension of the array WORK. 
*          If IJOB = 1, 2 or 4, LWORK >= MAX(4*N+16, 2*M*(N-M)).
*          If IJOB = 3 or 5, LWORK >= MAX(4*N+16, 4*M*(N-M)).
*
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
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
*
*  Further Details
*  ===============
*
*  DTGSEN first collects the selected eigenvalues by computing
*  orthogonal U and W that move them to the top left corner of (A, B).
*  In other words, the selected eigenvalues are the eigenvalues of
*  (A11, B11) in:
*
*                U'*(A, B)*W = (A11 A12) (B11 B12) n1
*                              ( 0  A22),( 0  B22) n2
*                                n1  n2    n1  n2
*
*  where N = n1+n2 and U' means the transpose of U. The first n1 columns
*  of U and W span the specified pair of left and right eigenspaces
*  (deflating subspaces) of (A, B).
*
*  If (A, B) has been obtained from the generalized real Schur
*  decomposition of a matrix pair (C, D) = Q*(A, B)*Z', then the
*  reordered generalized real Schur form of (C, D) is given by
*
*           (C, D) = (Q*U)*(U'*(A, B)*W)*(Z*W)',
*
*  and the first n1 columns of Q*U and Z*W span the corresponding
*  deflating subspaces of (C, D) (Q and Z store Q*U and Z*W, resp.).
*
*  Note that if the selected eigenvalue is sufficiently ill-conditioned,
*  then its value may differ significantly from its value before
*  reordering.
*
*  The reciprocal condition numbers of the left and right eigenspaces
*  spanned by the first n1 columns of U and W (or Q*U and Z*W) may
*  be returned in DIF(1:2), corresponding to Difu and Difl, resp.
*
*  The Difu and Difl are defined as:
*
*       Difu[(A11, B11), (A22, B22)] = sigma-min( Zu )
*  and
*       Difl[(A11, B11), (A22, B22)] = Difu[(A22, B22), (A11, B11)],
*
*  where sigma-min(Zu) is the smallest singular value of the
*  (2*n1*n2)-by-(2*n1*n2) matrix
*
*       Zu = [ kron(In2, A11)  -kron(A22', In1) ]
*            [ kron(In2, B11)  -kron(B22', In1) ].
*
*  Here, Inx is the identity matrix of size nx and A22' is the
*  transpose of A22. kron(X, Y) is the Kronecker product between
*  the matrices X and Y.
*
*  When DIF(2) is small, small changes in (A, B) can cause large changes
*  in the deflating subspace. An approximate (asymptotic) bound on the
*  maximum angular error in the computed deflating subspaces is
*
*       EPS * norm((A, B)) / DIF(2),
*
*  where EPS is the machine precision.
*
*  The reciprocal norm of the projectors on the left and right
*  eigenspaces associated with (A11, B11) may be returned in PL and PR.
*  They are computed as follows. First we compute L and R so that
*  P*(A, B)*Q is block diagonal, where
*
*       P = ( I -L ) n1           Q = ( I R ) n1
*           ( 0  I ) n2    and        ( 0 I ) n2
*             n1 n2                    n1 n2
*
*  and (L, R) is the solution to the generalized Sylvester equation
*
*       A11*R - L*A22 = -A12
*       B11*R - L*B22 = -B12
*
*  Then PL = (F-norm(L)**2+1)**(-1/2) and PR = (F-norm(R)**2+1)**(-1/2).
*  An approximate (asymptotic) bound on the average absolute error of
*  the selected eigenvalues is
*
*       EPS * norm((A, B)) / PL.
*
*  There are also global error bounds which valid for perturbations up
*  to a certain restriction:  A lower bound (x) on the smallest
*  F-norm(E,F) for which an eigenvalue of (A11, B11) may move and
*  coalesce with an eigenvalue of (A22, B22) under perturbation (E,F),
*  (i.e. (A + E, B + F), is
*
*   x = min(Difu,Difl)/((1/(PL*PL)+1/(PR*PR))**(1/2)+2*max(1/PL,1/PR)).
*
*  An approximate bound on x can be computed from DIF(1:2), PL and PR.
*
*  If y = ( F-norm(E,F) / x) <= 1, the angles between the perturbed
*  (L', R') and unperturbed (L, R) left and right deflating subspaces
*  associated with the selected cluster in the (1,1)-blocks can be
*  bounded as
*
*   max-angle(L, L') <= arctan( y * PL / (1 - y * (1 - PL * PL)**(1/2))
*   max-angle(R, R') <= arctan( y * PR / (1 - y * (1 - PR * PR)**(1/2))
*
*  See LAPACK User's Guide section 4.11 or the following references
*  for more information.
*
*  Note that if the default method for computing the Frobenius-norm-
*  based estimate DIF is not wanted (see DLATDF), then the parameter
*  IDIFJB (see below) should be changed from 3 to 4 (routine DLATDF
*  (IJOB = 2 will be used)). See DTGSYL for more details.
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
*      Note 75. To appear in ACM Trans. on Math. Software, Vol 22, No 1,
*      1996.
*
*  =====================================================================
*
* Pawel Kowal:
*   This is a direct rewrite of fragments Lapack's TGSEN function responsible
*   for estimating condition numbers
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
tgsen_cond(i_type ijob, i_type n, i_type m, const V *a,i_type lda, const V *b, i_type ldb, 
        typename details::real_type<V>::type& pl, typename details::real_type<V>::type& pr,
        typename details::real_type<V>::type* dif, V *work, i_type lwork, i_type *iwork, i_type liwork,
        i_type& info);

//--------------------------------------------------------------------------------
//                              GGES2
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  GGES2 computes for a pair of N-by-N real nonsymmetric matrices (A,B),
*  the generalized eigenvalues, the generalized real Schur form (S,T),
*  optionally, the left and/or right matrices of Schur vectors (VSL and
*  VSR). This gives the generalized Schur factorization
*
*           (A,B) = ( (VSL) * S * (VSR)**T, (VSL) * T * (VSR)**T )
*
*  The leading columns of VSL and VSR then form an orthonormal basis for the
*  corresponding left and right eigenspaces (deflating subspaces).
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
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the IWORK array, returns
*          this value as the first entry of the IWORK array, and no error
*          message related to LWORK is issued by XERBLA.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*          = 1,...,N:
*                The QZ iteration failed.  (A,B) are not in Schur
*                form, but ALPHAR(j), ALPHAI(j), and BETA(j) should
*                be correct for j=INFO+1,...,N.
*          > N:  =N+1: other than QZ iteration failed in DHGEQZ.
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gges2(const char *jobvsl, const char *jobvsr, i_type n, V *a, i_type lda, 
    V *b, i_type ldb, typename details::complex_type<V>::type *alpha, V *beta, V *vsl, 
    i_type ldvsl, V *vsr, i_type ldvsr, V *work, i_type lwork, i_type* iwork, i_type liwork, 
    i_type& info);

//--------------------------------------------------------------------------------
//                              DGGHRD
//--------------------------------------------------------------------------------
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
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If LIWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the IWORK array, returns
*          this value as the first entry of the IWORK array, and no error
*          message related to LIWORK is issued.
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
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gghrd2(const char* COMPQ, const char* COMPZ, i_type N, i_type ILO, i_type IHI, V* A, 
      i_type LDA, V* B, i_type LDB, V* Q, i_type LDQ, V* Z, i_type LDZ, 
      V* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, i_type& INFO );

//--------------------------------------------------------------------------------
//                              DGGHRDF
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  DGGHRD reduces a pair of real matrices (A,B) to generalized upper
*  Hessenberg form using orthogonal transformations, where A and B are 
*  general matrices. The form of the generalized eigenvalue problem is
*     A*x = lambda*B*x,
*  and moving the orthogonal matrix Q to the left side of the equation.
*
*  This subroutine simultaneously reduces A to a Hessenberg matrix H:
*     Q**T*A*Z = H
*  and transforms B to another upper triangular matrix T:
*     Q**T*B*Z = T
*  in order to reduce the problem to its standard form
*     H*y = lambda*T*y
*  where y = Z**T*x.
*
*  Arguments
*  =========
*
*  COMPQ   (input) CHARACTER*1
*          = 'N': do not compute Q;
*          = 'V': compute Q
*
*  COMPZ   (input) CHARACTER*1
*          = 'N': do not compute Z;
*          = 'V': compute Z
*
*  N       (input) INTEGER
*          The order of the matrices A and B.  N >= 0.
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
*          On entry, the N-by-N matrix B.
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
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          On exit, if INFO = 0, WORK(1) returns the optimal LWORK.
*
*  LWORK   (input) INTEGER
*          The dimension of the array WORK.
*          If LWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the WORK array, returns
*          this value as the first entry of the WORK array, and no error
*          message related to LWORK is issued.
*
*  IWORK   (workspace/output) INTEGER array, dimension (MAX(1,LIWORK))
*          On exit, if INFO = 0, IWORK(1) returns the optimal LIWORK.
*
*  LIWORK  (input) INTEGER
*          The dimension of the array IWORK.
*          If LIWORK = -1, then a workspace query is assumed; the routine
*          only calculates the optimal size of the IWORK array, returns
*          this value as the first entry of the IWORK array, and no error
*          message related to LIWORK is issued.
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
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gghrdf(const char* COMPQ, const char* COMPZ, i_type N, V* A, i_type LDA, V* B, i_type LDB, 
       V* Q, i_type LDQ, V* Z, i_type LDZ, V* WORK, i_type LWORK, i_type* IWORK, i_type LIWORK, 
       i_type& INFO );

//--------------------------------------------------------------------------------
//                              DGGBAL
//--------------------------------------------------------------------------------
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
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*  this version ggbal uses scaling technique proposed by Lemonnier and Van Dooren
*
*  Lemonnier, Damien, and Paul Van Dooren. "Balancing regular matrix pencils." 
*  SIAM journal on matrix analysis and applications 28.1 (2006): 253-263.
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
ggbal2(const char *job, i_type n, V* a, i_type lda, V* b, i_type ldb, i_type& ilo, i_type& ihi, 
        typename details::real_type<V>::type* lscale, typename details::real_type<V>::type* rscale, 
        i_type& info);

//--------------------------------------------------------------------------------
//                              ROTSEQ
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  apply the plane rotation L-sequence or R-sequence defined by arrays C, S, IND1, IND2
*  of length K.
*
*  The i-th elements of the left (L-) or right (R-) sequence (C, S, IND1, IND2) indicates,
*  that in the i-th step of construction of the sequence rows or  columns given by IND1 and
*  IND2 of a matrix are multiplied by the plane rotation G_i constructed from C(i) and
*  S(i), i.e. the following multiplication was performed:
*
*    [X(:,IND1) X(:,IND2)] := [X(:,IND1) X(:,IND2)] * GR_i   for R-sequences
*    [X(IND1,:);X(IND2,:)] := GL_i * [X(IND1,:);X(IND2,:)]   for L-sequences 
*
*  where
*       GL_i = [C(i)   S(i)],   GR_i = [C(i)  -S'(i)]
*              [-S(i)' C(i)]           [S(i)  C(i)]
*
*  in the following order:
*       GL_K * ... * GL_2 * GL_1 * X    for L-sequences
*       X  * GR_1 * GR_2 * ... * GR_K   for R-sequences
*
*  Left or right sequences can be applied from the left or from the right to a new matrix Y
*  forming:
*
*       op[GL_K * ... * GL_1] * Y   when L-sequence is applied from the left
*       Y * op[GR_1 * ... * GR_K]   when R-sequence is applied from the right
*       Y * op[GL_K * ... * GL_1]   when L-sequence is applied from the right
*       op[GR_1 * ... * GR_K] * X   when R-sequence is applied from the left
*
*  where op is one of the sequence operators
*
*       op[G_1 * ... * G_k] = G_1 * ... * G_k               if op = "N" (no trans)
*       op[G_1 * ... * G_k] = conj(G_1) * ... * conj(G_k)   if op = "Z" (conjugation)
*       op[G_1 * ... * G_k] = trans(G_k) * ... * trans(G_1) if op = "T" (transposition)
*       op[G_1 * ... * G_k] = G_k' * ... * G_1'             if op = "C" (conjugate transposition)
*
*  Arguments
*  =========
*
*  Type    (input) CHARACTER*1
*          = 'L': the sequence was constructed as the left sequence
*          = 'R': the sequence was constructed as the right sequence
*
*  SIDE    (input) CHARACTER*1
*          = 'L': apply rotation seqence from the left
*          = 'R': apply rotation seqence from the right
*
*  TRANS   (input) CHARACTER*1
*          = 'N': no transpose
*          = 'Z': conjugation
*          = 'T': transpose
*          = 'C': conjugate transpose
*
*  K       (input) INTEGER
*          Order of the plane rotation sequence (K >= 0)
*
*  C       (input) DOUBLE PRECISION array, dimension (max(K))
*          Cosines of the plane rotations.
*
*  S       (input) DOUBLE PRECISION array, dimension (max(K))
*          Sines of the plane rotations.
*
*  IND1    (input) INTEGER array, dimension (max(K))
*          First indices of rotation sequence.
*
*  IND2    (input) INTEGER array, dimension (max(0,N-1))
*          Second indices of rotation sequence.
*
*  M       (input) INTEGER
*          Number of rows of the matrix X, (M >= 0)
*
*  N       (input) INTEGER
*          Number of columns of the matrix X, (N >= 0)
*
*  LDIAGS  (input/output) INTEGER
*          On entry upper bound on number of subdiagonals of X, (LDIAGS >= 0)
*          On exit returns new bound on number of subdiagonals of the modified
*          matrix X
*
*  UDIAGS  (input/output) INTEGER
*          On entry upper bound on number of superdiagonals of X, (UDIAGS >= 0)
*          On exit returns new bound on number of superdiagonals of the modified
*          matrix X
*
*  X       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the input matrix of size M x N, 
*          On exit, X is overwritten by the product of X and the sequence of plane 
*          rotations
*
*  LDX     (input) INTEGER
*          The leading dimension of the array X. LDR >= max(1,M).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*/
template<class V1, class V2> BLAS_EXT_EXPORT
typename details::enable_if_valid2<void,V1, V2>::type
rotseq(const char* TYPE, const char* SIDE, const char* TRANS, i_type K, 
       const typename details::real_type<V1>::type* C, const V1* S, const i_type* IND1, const i_type* IND2,
       i_type M, i_type N, i_type& LDIAGS, i_type& UDIAGS, V2* X, i_type LDX, i_type& INFO);

//--------------------------------------------------------------------------------
//                              QRUPHR
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  QRUPHR generates plane rotations that reduce rank-1 update of an upper Hessenberg 
*  square matrix R with K subdiagonals in the form
*                           RR = R + sigma * U * V'
*  where U, V are vectors and sigma is a scalar to upper hessenberg matrix using plane
*  rotations applied from the right 
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          Number of rows and columns of the matrix R.
*
*  K       (input) INTEGER
*          Number of subdiagonals of the matrix R (K >= 0).
*
*  R       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the input upper Hessenberg matrix with K subdiagonals, 
*          On exit,  R is overwritten by the upper Hessenberg matrix with K+1
*          subdiagonals;
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,M).
*
*  SIGMA   (input) DOUBLE PRECISION scalar
*
*  U       (input) DOUBLE PRECISION array, dimension (M)
*          Vector used for the rank one update.
*
*  INCU    (input) INTEGER
*          Distance between two consecutive elements of U.
*
*  V       (input) DOUBLE PRECISION array, dimension (M)
*          Vector used for the rank one update.
*
*  INCV    (input) INTEGER
*          Distance between two consecutive elements of V.
*
*  C       (output) DOUBLE PRECISION array, dimension (max(0,N-2))
*          On exit returns cosines of the plane rotations.
*
*  S       (output) DOUBLE PRECISION array, dimension (max(0,N-2))
*          On exit returns sines of the plane rotations.
*
*  IND1    (output) INTEGER array, dimension (max(0,N-2))
*          On exit returns first indices of rotation sequence.
*
*  IND2    (output) INTEGER array, dimension (max(0,N-2))
*          On exit returns second indices of rotation sequence.
*
*  SLEN    (output) INTEGER
*          Length of arrays C, S, IND1, IND2.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*  vectors C and S have size N - 2; i-th elements of C, S, IND1, IND2 indicates,
*  that in the i-th step columns IND1 and IND2 of RR were multiplied by the plane
*  rotation G_i constructed from C(i) and S(i), i.e. the following multiplication
*  was done:
*       [RR(:,IND1) RR(:,IND2)] := [RR(:,IND1) RR(:,IND2)] * [C(i) -S(i)']
*                                                            [S(i) C(i)]
*  forming the multiplication sequence with matrices G_1, G_2, ..., G_{N - 2}
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
qruphr(i_type N, i_type K, V* R, i_type LDR, const V& SIGMA, const V* U, i_type inc_u, const V* W, 
       i_type inc_w, typename details::real_type<V>::type* C, V* S, i_type* IND1,
       i_type* IND2, i_type& SLEN, i_type& INFO);

//--------------------------------------------------------------------------------
//                              QRUPHL
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  QRUPHL generates plane rotations that reduce rank-1 update of a Hessenberg matrix R
*  with K subdiagonals in the form
*                           RR = R + sigma * U * V'
*  where U, V are vectors and sigma is a scalar to upper hessenberg matrix with K + 1
*  subdiagonals using plane rotations applied from the left
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          Number of rows of the upper triangular matrix R.
*
*  N       (input) INTEGER
*          Number of columns of the upper triangular matrix R.
*
*  K       (input) INTEGER
*          Number of subdiagonals of the matrix R (K >= 0).
*
*  R       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the input upper Hessenberg matrix with K subdiagonals, 
*          On exit,  R is overwritten by the upper Hessenberg matrix with K+1
*          subdiagonals;
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,M).
*
*  SIGMA   (input) DOUBLE PRECISION scalar
*
*  U       (input) DOUBLE PRECISION array, dimension (M)
*          Vector used for the rank one update.
*
*  INCU    (input) INTEGER
*          Distance between two consecutive elements of U.
*
*  V       (input) DOUBLE PRECISION array, dimension (M)
*          Vector used for the rank one update.
*
*  INCV    (input) INTEGER
*          Distance between two consecutive elements of V.
*
*  C       (output) DOUBLE PRECISION array, dimension (max(0,M-1))
*          On exit returns cosines of the plane rotations.
*
*  S       (output) DOUBLE PRECISION array, dimension (max(0,M-1))
*          On exit returns sines of the plane rotations.
*
*  IND1    (output) INTEGER array, dimension (max(0,N-2))
*          On exit returns first indices of rotation sequence.
*
*  IND2    (output) INTEGER array, dimension (max(0,N-2))
*          On exit returns second indices of rotation sequence.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  SLEN    (output) INTEGER
*          Length of arrays C, S, IND1, IND2.
*
*  Further Details
*  ===============
*  vectors C and S have size N - 2; i-th elements of C, S, IND1, IND2 indicates,
*  that in the i-th step columns IND1 and IND2 of RR were multiplied by the plane
*  rotation G_i constructed from C(i) and S(i), i.e. the following multiplication
*  was done:
*       [RR(IND1,:); RR(IND2,:)] := [C(i)   S(i)] * [RR(IND1,:); RR(IND2,:)]
*                                   [-S(i)' C(i)]
*  forming the multiplication sequence with matrices G_1, G_2, ..., G_{N - 2}
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
qruphl(i_type M, i_type N, i_type K, V* R, i_type LDR, const V& SIGMA, const V* U, i_type inc_u, 
       const V* W, i_type inc_w, typename details::real_type<V>::type* C, V* S, i_type* IND1,
       i_type* IND2, i_type& SLEN, i_type& INFO);

//--------------------------------------------------------------------------------
//                              HUUNDR
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  HUUNDR generates plane rotations that anihilate subdiagonals of a square part RS
*  of the matrix R = [R0; RS], with D subdiagonal using plane rotations applied from
*  the right 
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          Number of rows and columns of the matrix R.
*
*  D       (input) INTEGER
*          Number of subdiagonals of the matrix R (D > 0 and D <= N - 1)
*
*  K       (input) INTEGER
*          Number of rows of the matrix R0, (K >= 0)
*
*  R       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix R with D subdiagonals in the squre part RS, 
*          On exit,  R is overwritten by a matrix with zero subdiagonals
*          in the square part.
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,N+K).
*
*  C       (output) DOUBLE PRECISION array, dimension (D * max(0,N-1))
*          On exit returns cosines of the plane rotations.
*
*  S       (output) DOUBLE PRECISION array, dimension (D * max(0,N-1))
*          On exit returns sines of the plane rotations.
*
*  IND1    (output) INTEGER array, dimension (D * max(0,N-1))
*          On exit returns first indices of rotation sequence.
*
*  IND2    (output) INTEGER array, dimension (D * max(0,N-1))
*          On exit returns second indices of rotation sequence.
*
*  SLEN    (output) INTEGER
*          Length of arrays C, S, IND1, IND2.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*  elements of C, S, IND1, IND2 indicates, that in the i-th step columns IND1 and IND2
*  of R were multiplied by the plane rotation G_i constructed from C(i) and S(i), i.e.
*  the following multiplication  was done:
*       [R(:,IND1) R(:,IND2)] := [R(:,IND1) R(:,IND2)] * [C(i) -S(i)']
*                                                        [S(i)  C(i)]
*  forming the multiplication sequence with matrices G_1, G_2, ..., G_{N - 2}
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
huundr(i_type N, i_type D, i_type K, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO);

//--------------------------------------------------------------------------------
//                              HUUNDL
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  HUUNDR generates plane rotations that anihilate subdiagonals of a matrix R
*  with D subdiagonal using plane rotations applied from the left 
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          Number of rows of the matrix R.
*
*  N       (input) INTEGER
*          Number of columns of the matrix R.
*
*  D       (input) INTEGER
*          Number of subdiagonals of the matrix R (D > 0 and D <= M - 1)
*
*  R       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix R with D subdiagonals, 
*          On exit,  R is overwritten by an upper triangular matrix.
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,M).
*
*  C       (output) DOUBLE PRECISION array, dimension (D * max(0,min(N,M-1)) )
*          On exit returns cosines of the plane rotations.
*
*  S       (output) DOUBLE PRECISION array, dimension (D * max(0,min(N,M-1)) )
*          On exit returns sines of the plane rotations.
*
*  IND1    (output) INTEGER array, dimension (D * max(0,min(N,M-1)) )
*          On exit returns first indices of rotation sequence.
*
*  IND2    (output) INTEGER array, dimension (D * max(0,min(N,M-1)) )
*          On exit returns second indices of rotation sequence.
*
*  SLEN    (output) INTEGER
*          Length of arrays C, S, IND1, IND2.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*  elements of C, S, IND1, IND2 indicates, that in the i-th step rows IND1 and IND2
*  of R were multiplied by the plane rotation G_i constructed from C(i) and S(i), i.e.
*  the following multiplication  was done:
*       [R(IND1,:); R(IND2,:)] := [C(i)   S(i)] * [R(IND1,:); R(IND2,:)]
*                                 [-S(i)' C(i)]
*  forming the multiplication sequence with matrices G_1, G_2, ..., G_{N - 2}
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
huundl(i_type M, i_type N, i_type D, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO);

//--------------------------------------------------------------------------------
//                              HLUNDL
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  HLUNDL generates plane rotations that anihilate superdiagonals of a square matrix R
*  with D superdiagonals using plane rotations applied from the left 
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          Number of rows and columns of the matrix R.
*
*  D       (input) INTEGER
*          Number of superdiagonals of the matrix R (D > 0 and D <= N - 1)
*
*  R       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix R with D superdiagonals, 
*          On exit,  R is overwritten by a lower triangular matrix.
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,M).
*
*  C       (output) DOUBLE PRECISION array, dimension (D * max(0,N-1))
*          On exit returns cosines of the plane rotations.
*
*  S       (output) DOUBLE PRECISION array, dimension (D * max(0,N-1))
*          On exit returns sines of the plane rotations.
*
*  IND1    (output) INTEGER array, dimension (D * max(0,N-1))
*          On exit returns first indices of rotation sequence.
*
*  IND2    (output) INTEGER array, dimension (D * max(0,N-1))
*          On exit returns second indices of rotation sequence.
*
*  SLEN    (output) INTEGER
*          Length of arrays C, S, IND1, IND2.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*  elements of C, S, IND1, IND2 indicates, that in the i-th step rows IND1 and IND2
*  of R were multiplied by the plane rotation G_i constructed from C(i) and S(i), i.e.
*  the following multiplication  was done:
*       [R(IND1,:); R(IND2,:)] := [C(i)   S(i)] * [R(IND1,:); R(IND2,:)]
*                                 [-S(i)' C(i)]
*  forming the multiplication sequence with matrices G_1, G_2, ..., G_{N - 2}
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
hlundl(i_type N, i_type D, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO);

//--------------------------------------------------------------------------------
//                              HLUNDR
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  HLUNDR generates plane rotations that anihilate superdiagonals of a matrix R
*  with D superdiagonals using plane rotations applied from the right
*
*  Arguments
*  =========
*
*  M       (input) INTEGER
*          Number of rows of the matrix R.
*
*  N       (input) INTEGER
*          Number of columns of the matrix R.
*
*  D       (input) INTEGER
*          Number of superdiagonals of the matrix R (D > 0 and D <= N - 1)
*
*  R       (input/output) DOUBLE PRECISION array, dimension (LDA,N)
*          On entry, the matrix R with D subdiagonals, 
*          On exit,  R is overwritten by a lower triangular matrix.
*
*  LDR     (input) INTEGER
*          The leading dimension of the array R. LDR >= max(1,M).
*
*  C       (output) DOUBLE PRECISION array, dimension (D * max(0,min(N-1,M)) )
*          On exit returns cosines of the plane rotations.
*
*  S       (output) DOUBLE PRECISION array, dimension (D * max(0,min(N-1,M)) )
*          On exit returns sines of the plane rotations.
*
*  IND1    (output) INTEGER array, dimension (D * max(0,min(N-1,M)) )
*          On exit returns first indices of rotation sequence.
*
*  IND2    (output) INTEGER array, dimension (D * max(0,min(N-1,M)) )
*          On exit returns second indices of rotation sequence.
*
*  SLEN    (output) INTEGER
*          Length of arrays C, S, IND1, IND2.
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value.
*
*  Further Details
*  ===============
*  elements of C, S, IND1, IND2 indicates, that in the i-th step columns IND1 and IND2
*  of R were multiplied by the plane rotation G_i constructed from C(i) and S(i), i.e.
*  the following multiplication  was done:
*       [R(:,IND1) R(:,IND2)] := [C(i)   S(i)] * [R(:,IND1) R(:,IND2)]
*                                [-S(i)' C(i)]
*  forming the multiplication sequence with matrices G_1, G_2, ..., G_{N - 2}
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
hlundr(i_type M, i_type N, i_type D, V* R, i_type LDR, typename details::real_type<V>::type* C, 
        V* S, i_type* IND1, i_type* IND2, i_type& SLEN, i_type& INFO);

//--------------------------------------------------------------------------------
//                              DHGEQZ
//--------------------------------------------------------------------------------
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

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
hgeqz2(const char* job, const char* compq, const char* compz, i_type n, i_type ilo, i_type ihi,
    V* h, i_type ldh, V* t, i_type ldt, V* alphar, V* alphai, V* beta, V* q, i_type ldq, V* z, 
    i_type ldz, V* work, i_type lwork, i_type& info);

//--------------------------------------------------------------------------------
//                              gs_ort
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  gs_ort orthogonalize the vector x againts the matrix Q with orthonormal columns, i.e.
*  Q' * Q = I, using Gram-Schmidt with iterative refinement as in Daniel et al (1976)
*  such that
*       x = Q * s + r, with Q'*r = 0    (1)
*  the matrix Q is given by
*       Q = V           if trans =  'N'
*       Q = V'          if trans != 'N'
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          = 'N': The matrix Q is given by V
*          !='N': The matrix Q is given by V'
*
*  N       (input) INTEGER
*          The length of vector x and number of rows of Q,  N >= 0.
*
*  K       (input) INTEGER
*          The number of columns of Q,  K >= 0, K <= N
*
*  V       (input) DOUBLE PRECISION array, dimension (LDV, K)
*          The matrix with K orthonormal columns of length N. 
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.  LDQ >= 1.
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry the vector to orthogonalize. 
*          On exit, X contains the orthogonalized vector r. if 
*          If INFO > 0, then X is zero to 0 and XNRM is set to 0.
*
*  INCX    (input) INTEGER
*          The distance between two consequtive elements of X. INCX != 0.
*
*  S       (output) DOUBLE PRECISION array, dimension (K)
*          On exit the vector s as in (1)
*
*  INCS    (input) INTEGER
*          The distance between two consequtive elements of S. INCS != 0.
*
*  XRNM    (output) DOUBLE PRECISION
*          On exit the norm of orthogonalized vector r stored in X
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,K))
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: orthogonalization failed, X is possibly a linear combination
*               of columns of Q
*
* Ref: 
*  =========
*
*   Daniel, J. W., W. B. Gragg, L. Kaufman, and G. W. Stewart (1976). 
*   Reorthogonalization and Stable Algorithms for Updating the Gram–Schmidt QR Factorization.
*   Math. Comp., 30(136):772–795.
*  
*   Reichel, L. and W. B. Gragg (1990). FORTRAN Subroutines for Updating the QR Decomposition.
*   ACM Trans. Math. Softw., 16:369–377
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gs_ort(const char* TRANS, i_type N, i_type K, const V* Q, i_type LDQ, V* X, i_type incx, 
            V* s, i_type incs, typename details::real_type<V>::type& x_norm, V* work, i_type& info);

//--------------------------------------------------------------------------------
//                              gsg_ort
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  gsg_ort orthogonalize the vector x againts the matrix Q with orthonormal columns
*  with respect to inner product given by the matrix B, i.e. Q' * B * Q = I, 
*  using Gram-Schmidt with iterative refinement as in Daniel et al (1976) such that
*       x = Q * s + r, with Q' * B * r = 0    (1)
*  the matrix Q is given by
*       Q = V           if trans =  'N'
*       Q = V'          if trans != 'N'
*  B must be positive semi-definite hermitian matrix and is given by a function that
*  evaluates product B * x
*
*  Arguments
*  =========
*
*  CALLB   (input) pointer to function that evaluates B * X
*               void CALLB(CTS, N, X, WB);
*           CTX : context passed to the function gsg_ort
*           N   : integer N passed to the function gsg_ort
*           X   : pointer to a vector of length B, always the same
*                 as pointer X passed to the function gsg_ort
*           WB  : array to store the result of B * X of length N
*                 always the same as WORK_B passed to the function gsg_ort
*
*  CTX     (input) additional context information that is passed to 
*          the function CALLB
*
*  TRANS   (input) CHARACTER*1
*          = 'N': The matrix Q is given by V
*          !='N': The matrix Q is given by V'
*
*  N       (input) INTEGER
*          The length of vector x and number of rows of Q,  N >= 0.
*
*  K       (input) INTEGER
*          The number of columns of Q,  K >= 0, K <= N
*
*  V       (input) DOUBLE PRECISION array, dimension (LDV, K)
*          The matrix with K orthonormal columns of length N. 
*
*  LDV     (input) INTEGER
*          The leading dimension of the array V.  LDQ >= 1.
*
*  X       (input/output) DOUBLE PRECISION array, dimension (N)
*          On entry the vector to orthogonalize. 
*          On exit, X contains the orthogonalized vector r. if 
*          If INFO > 0, then X is zero to 0 and XNRM is set to 0.
*
*  INCX    (input) INTEGER
*          The distance between two consequtive elements of X. INCX != 0.
*
*  S       (output) DOUBLE PRECISION array, dimension (K)
*          On exit the vector s as in (1)
*
*  INCS    (input) INTEGER
*          The distance between two consequtive elements of S. INCS != 0.
*
*  XRNM    (output) DOUBLE PRECISION
*          On exit the norm of orthogonalized vector r stored in X
*
*  WORK_B  (workspace/output) DOUBLE PRECISION array, dimension (N)
*          array that stores result of B * X
*
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,K))
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -i, the i-th argument had an illegal value
*          = 1: orthogonalization failed, X is possibly a linear combination
*               of columns of Q
*
* Ref: 
*  =========
*
*   Daniel, J. W., W. B. Gragg, L. Kaufman, and G. W. Stewart (1976). 
*   Reorthogonalization and Stable Algorithms for Updating the Gram–Schmidt QR Factorization.
*   Math. Comp., 30(136):772–795.
*  
*   Reichel, L. and W. B. Gragg (1990). FORTRAN Subroutines for Updating the QR Decomposition.
*   ACM Trans. Math. Softw., 16:369–377
*/

template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gsg_ort(void (*CALLBACK)(void* CTX, i_type N, V* X, V* WORK), void* CTX, 
        const char* TRANS, i_type N, i_type K, const V* Q, i_type LDQ, V* X, i_type INCX, 
        V* S, i_type INCS, typename details::real_type<V>::type& X_NORM, 
        V* WORK_B, V* WORK, i_type& INFO);

//--------------------------------------------------------------------------------
//                              ptinorm
//--------------------------------------------------------------------------------
/*
*  Purpose
*  =======
*
*  ptinorm computes the the 1-norm of an inverse of symmetric positive definite 
*  tridiagonal matrix from the factorization A = L*D*L**T. The norm is computed
*  by a direct method
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the matrix A.  N >= 0.
*
*  D       (input) DOUBLE PRECISION array, dimension (N)
*          The n diagonal elements of the diagonal matrix D from the
*          factorization of A
*
*  E       (input) DOUBLE PRECISION array, dimension (N-1)
*          The (n-1) off-diagonal elements of the unit bidiagonal factor
*          U or L from the factorization of A
*
*  NORM    (output) DOUBLE PRECISION
*          The 1-norm of inv(A)
*
*  WORK    (workspace) DOUBLE PRECISION array, dimension (N)
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*
*  Further Details
*  ===============
*
*  This is a modification of DPTCON function from Lapack.
*
*  The method used is described in Nicholas J. Higham, "Efficient
*  Algorithms for Computing the Condition Number of a Tridiagonal
*  Matrix", SIAM J. Sci. Stat. Comput., Vol. 7, No. 1, January 1986.
*  
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
ptinorm(i_type N, const typename details::real_type<V>::type* D, const V* E, 
        typename details::real_type<V>::type& NORM, typename details::real_type<V>::type* WORK, 
        i_type& INFO );

//=======================   PTTRS_REV    =================================
/*
*  Purpose
*  =======
*
*  PTTRS_REV solves a tridiagonal system of the form
*     X * A = B
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
*          The number of right hand sides, i.e., the number of rows
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
*  B       (input/output) DOUBLE PRECISION array, dimension (LDB,N)
*          On entry, the right hand side vectors B for the system of
*          linear equations.
*          On exit, the solution vectors, X.
*
*  LDB     (input) INTEGER
*          The leading dimension of the array B.  LDB >= max(1,NRHS).
*
*  INFO    (output) INTEGER
*          = 0: successful exit
*          < 0: if INFO = -k, the k-th argument had an illegal value
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
pttrs_rev(const char* UPLO, i_type N, i_type NRHS, const typename details::real_type<V>::type* D, 
      const V* E, V* B, i_type LDB, i_type& INFO );

//=======================   gttrs_rev    =================================
/*
*  Purpose
*  =======
*
*  gttrs_rev solves one of the systems of equations
*     X * A = B  or  X * A' = B,
*  with a tridiagonal matrix A using the LU factorization computed
*  by DGTTRF.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          Specifies the form of the system of equations.
*          = 'N':  X * A = B  (No transpose)
*          = 'T':  X * A' = B  (Transpose)
*          = 'C':  X * A' = B  (Conjugate transpose = Transpose)
*
*  N       (input) INTEGER
*          The order of the matrix A.
*
*  NRHS    (input) INTEGER
*          The number of right hand sides, i.e., the number of rows
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
*          The leading dimension of the array B.  LDB >= max(1,NRHS).
*
*  INFO    (output) INTEGER
*          = 0:  successful exit
*          < 0:  if INFO = -i, the i-th argument had an illegal value
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
gttrs_rev(const char* TRANS, i_type N, i_type NRHS, const V* DL, const V* D, const V* DU, const V* DU2, 
      const i_type* IPIV, V* B, i_type LDB, i_type& INFO);

/*
*  Purpose
*  =======
*
*  ZLANV2 computes the Schur factorization of a complex 2-by-2
*  nonhermitian matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [  0  DD ] [-SN  CS ]
*
*  Arguments
*  =========
*
*  A       (input/output) COMPLEX*16
*  B       (input/output) COMPLEX*16
*  C       (input/output) COMPLEX*16
*  D       (input/output) COMPLEX*16
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1     (output) COMPLEX*16
*  RT2     (output) COMPLEX*16
*          The two eigenvalues.
*
*  CS      (output) DOUBLE PRECISION
*  SN      (output) COMPLEX*16
*          Parameters of the rotation matrix.
*
*  Further Details
*  ===============
*
*  Implemented by Mark R. Fahey, May 28, 1999
*  this is C++ version of SCALAPACK function ZLANV2
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_complex<void,V>::type
lanv2z(V& A, V& B, V& C, V& D, V& RT1, V& RT2, typename details::real_type<V>::type& CS, V& SN);

//=======================   DLAEV2    =====================================================      
/*
*  Purpose
*  =======
*
*  DLAEV2 computes the eigendecomposition of a 2-by-2 symmetric matrix
*     [  A   B  ]
*     [  B   C  ].
*  On return, RT1 is the eigenvalue of larger absolute value, RT2 is the
*  eigenvalue of smaller absolute value, and (CS1,SN1) is the unit right
*  eigenvector for RT1, giving the decomposition
*
*     [ CS1  SN1 ] [  A   B  ] [ CS1 -SN1 ]  =  [ RT1  0  ]
*     [-SN1  CS1 ] [  B   C  ] [ SN1  CS1 ]     [  0  RT2 ].
*
*  Arguments
*  =========
*
*  A       (input) DOUBLE PRECISION
*          The (1,1) element of the 2-by-2 matrix.
*
*  B       (input) DOUBLE PRECISION
*          The (1,2) element and the conjugate of the (2,1) element of
*          the 2-by-2 matrix.
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
*  CS1     (output) DOUBLE PRECISION
*  SN1     (output) DOUBLE PRECISION
*          The vector (CS1, SN1) is a unit right eigenvector for RT1.
*
* this is C++ version of Lapack function dlaev2
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
laev2(const V& A, const V& B, const V& C, V& RT1, V& RT2, V& CS1, V& SN1 );

//=======================   lanv2   =====================================================
/*
*  Purpose
*  =======
*
*  DLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
*  matrix in standard form:
*
*       [ A  B ] = [ CS -SN ] [ AA  BB ] [ CS  SN ]
*       [ C  D ]   [ SN  CS ] [ CC  DD ] [-SN  CS ]
*
*  where either
*  1) CC = 0 so that AA and DD are real eigenvalues of the matrix, or
*  2) AA = DD and BB*CC < 0, so that AA + or - sqrt(BB*CC) are complex
*  conjugate eigenvalues.
*
*  Arguments
*  =========
*
*  A       (input/output) DOUBLE PRECISION
*  B       (input/output) DOUBLE PRECISION
*  C       (input/output) DOUBLE PRECISION
*  D       (input/output) DOUBLE PRECISION
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1R    (output) DOUBLE PRECISION
*  RT1I    (output) DOUBLE PRECISION
*  RT2R    (output) DOUBLE PRECISION
*  RT2I    (output) DOUBLE PRECISION
*          The real and imaginary parts of the eigenvalues. If the
*          eigenvalues are a complex conjugate pair, RT1I > 0.
*
*  CS      (output) DOUBLE PRECISION
*  SN      (output) DOUBLE PRECISION
*          Parameters of the rotation matrix.
*
*  Further Details
*  ===============
*
*  Modified by V. Sima, Research Institute for Informatics, Bucharest,
*  Romania, to reduce the risk of cancellation errors,
*  when computing real eigenvalues, and to ensure, if possible, that
*  abs(RT1R) >= abs(RT2R).
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
lanv2(V& A, V& B, V& C, V& D, V& RT1R, V& RT1I, V& RT2R, V& RT2I, V& CS, V& SN);

//=======================   DLASV2    =====================================================
/*      
*  Purpose
*  =======
*
*  DLASV2 computes the singular value decomposition of a 2-by-2
*  triangular matrix
*     [  F   G  ]
*     [  0   H  ].
*  On return, abs(SSMAX) is the larger singular value, abs(SSMIN) is the
*  smaller singular value, and (CSL,SNL) and (CSR,SNR) are the left and
*  right singular vectors for abs(SSMAX), giving the decomposition
*
*     [ CSL  SNL ] [  F   G  ] [ CSR -SNR ]  =  [ SSMAX   0   ]
*     [-SNL  CSL ] [  0   H  ] [ SNR  CSR ]     [  0    SSMIN ].
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
*          abs(SSMIN) is the smaller singular value.
*
*  SSMAX   (output) DOUBLE PRECISION
*          abs(SSMAX) is the larger singular value.
*
*  SNL     (output) DOUBLE PRECISION
*  CSL     (output) DOUBLE PRECISION
*          The vector (CSL, SNL) is a unit left singular vector for the
*          singular value abs(SSMAX).
*
*  SNR     (output) DOUBLE PRECISION
*  CSR     (output) DOUBLE PRECISION
*          The vector (CSR, SNR) is a unit right singular vector for the
*          singular value abs(SSMAX).
*
*  Further Details
*  ===============
*
*  Any input parameter may be aliased with any output parameter.
*
*  Barring over/underflow and assuming a guard digit in subtraction, all
*  output quantities are correct to within a few units in the last
*  place (ulps).
*
*  In IEEE arithmetic, the code works correctly if one matrix element is
*  infinite.
*
*  Overflow will not occur unless the largest singular value itself
*  overflows or is within a few ulps of overflow. (On machines with
*  partial overflow, like the Cray, overflow may occur if the largest
*  singular value is within a factor of 2 of overflow.)
*
*  Underflow is harmless if underflow is gradual. Otherwise, results
*  may correspond to a matrix modified by perturbations of size near
*  the underflow threshold.
*
*/
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid_real<void,V>::type
lasv2(V F, V G, V H, V& SSMIN, V& SSMAX, V& SNR, V& CSR, V& SNL, V& CSL);

};};
