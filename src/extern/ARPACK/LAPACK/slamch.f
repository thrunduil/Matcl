      REAL             FUNCTION AR_SLAMCH( CMACH )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     October 31, 1992 
*
*     .. Scalar Arguments ..
      CHARACTER          CMACH
*     ..
*
*  Purpose
*  =======
*
*  AR_SLAMCH determines single precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by AR_SLAMCH:
*          = 'E' or 'e',   AR_SLAMCH := eps
*          = 'S' or 's ,   AR_SLAMCH := sfmin
*          = 'B' or 'b',   AR_SLAMCH := base
*          = 'P' or 'p',   AR_SLAMCH := eps*base
*          = 'N' or 'n',   AR_SLAMCH := t
*          = 'R' or 'r',   AR_SLAMCH := rnd
*          = 'M' or 'm',   AR_SLAMCH := emin
*          = 'U' or 'u',   AR_SLAMCH := rmin
*          = 'L' or 'l',   AR_SLAMCH := emax
*          = 'O' or 'o',   AR_SLAMCH := rmax
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
*
* =====================================================================
*
*     .. External Functions ..
      REAL               SLAMCH
      EXTERNAL           SLAMCH
*     ..
*     ..
*     .. Executable Statements ..
      AR_SLAMCH = SLAMCH(CMACH)
*
*     End of AR_SLAMCH
*
      END