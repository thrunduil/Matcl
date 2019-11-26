      REAL             FUNCTION AR_SLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      REAL               X, Y
*     ..
*
*  Purpose
*  =======
*
*  AR_SLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) REAL
*  Y       (input) REAL
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. External Functions ..
      REAL               SLAPY2
      EXTERNAL           SLAPY2
*     ..
*     ..
*     .. Executable Statements ..
      AR_SLAPY2 = SLAPY2(X, Y)
*
*     End of AR_SLAPY2
*
      END