      SUBROUTINE AR_DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     February 29, 1992
*
*     .. Scalar Arguments ..
      CHARACTER          UPLO
      INTEGER            LDA, LDB, M, N
*     ..
*     .. Array Arguments ..
      DOUBLE PRECISION   A( LDA, * ), B( LDB, * )
*     ..
*
*  Purpose
*  =======
*
*  AR_DLACPY copies all or part of a two-dimensional matrix A to another
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
*
*  =====================================================================
*
*     .. External Functions ..
      EXTERNAL           DLACPY
*     ..
*     ..
*     .. Executable Statements ..
      CALL DLACPY( UPLO, M, N, A, LDA, B, LDB )
*
*     End of AR_DLACPY
*
      END