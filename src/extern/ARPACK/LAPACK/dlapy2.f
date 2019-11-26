      DOUBLE PRECISION FUNCTION AR_DLAPY2( X, Y )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992
*
*     .. Scalar Arguments ..
      DOUBLE PRECISION   X, Y
*     ..
*
*  Purpose
*  =======
*
*  AR_DLAPY2 returns sqrt(x**2+y**2), taking care not to cause unnecessary
*  overflow.
*
*  Arguments
*  =========
*
*  X       (input) DOUBLE PRECISION
*  Y       (input) DOUBLE PRECISION
*          X and Y specify the values x and y.
*
*  =====================================================================
*
*     .. External Functions ..
      DOUBLE PRECISION   DLAPY2
      EXTERNAL           DLAPY2
*     ..
*     ..
*     .. Executable Statements ..
      AR_DLAPY2 = DLAPY2(X, Y)
*
*     End of AR_DLAPY2
*
      END