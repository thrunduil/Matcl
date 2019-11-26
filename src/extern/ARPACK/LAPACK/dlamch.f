      DOUBLE PRECISION FUNCTION AR_DLAMCH( CMACH )
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
*  AR_DLAMCH determines double precision machine parameters.
*
*  Arguments
*  =========
*
*  CMACH   (input) CHARACTER*1
*          Specifies the value to be returned by AR_DLAMCH:
*          = 'E' or 'e',   AR_DLAMCH := eps
*          = 'S' or 's ,   AR_DLAMCH := sfmin
*          = 'B' or 'b',   AR_DLAMCH := base
*          = 'P' or 'p',   AR_DLAMCH := eps*base
*          = 'N' or 'n',   AR_DLAMCH := t
*          = 'R' or 'r',   AR_DLAMCH := rnd
*          = 'M' or 'm',   AR_DLAMCH := emin
*          = 'U' or 'u',   AR_DLAMCH := rmin
*          = 'L' or 'l',   AR_DLAMCH := emax
*          = 'O' or 'o',   AR_DLAMCH := rmax
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
      DOUBLE PRECISION   DLAMCH
      EXTERNAL           DLAMCH
*     ..
*     ..
*     .. Executable Statements ..
      AR_DLAMCH = DLAMCH(CMACH)
*
*     End of AR_DLAMCH
*
      END