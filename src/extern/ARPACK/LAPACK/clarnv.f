      SUBROUTINE AR_CLARNV( IDIST, ISEED, N, X )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            IDIST, N
*     ..
*     .. Array Arguments ..
      INTEGER            ISEED( 4 )
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  AR_CLARNV returns a vector of n random complex numbers from a uniform or
*  normal distribution.
*
*  Arguments
*  =========
*
*  IDIST   (input) INTEGER
*          Specifies the distribution of the random numbers:
*          = 1:  real and imaginary parts each uniform (0,1)
*          = 2:  real and imaginary parts each uniform (-1,1)
*          = 3:  real and imaginary parts each normal (0,1)
*          = 4:  uniformly distributed on the disc abs(z) < 1
*          = 5:  uniformly distributed on the circle abs(z) = 1
*
*  ISEED   (input/output) INTEGER array, dimension (4)
*          On entry, the seed of the random number generator; the array
*          elements must be between 0 and 4095, and ISEED(4) must be
*          odd.
*          On exit, the seed is updated.
*
*  N       (input) INTEGER
*          The number of random numbers to be generated.
*
*  X       (output) COMPLEX array, dimension (N)
*          The generated random numbers.
*
*  Further Details
*  ===============
*
*  This routine calls the auxiliary routine AR_SLARUV to generate random
*  real numbers from a uniform (0,1) distribution, in batches of up to
*  128 using vectorisable code. The Box-Muller method is used to
*  transform numbers from a uniform to a normal distribution.
*
*  =====================================================================
*
*     .. External Functions ..
      EXTERNAL           CLARNV
*     ..
*     ..
*     .. Executable Statements ..
      CALL CLARNV( IDIST, ISEED, N, X )
*
*     End of AR_CLARNV
*
      END
