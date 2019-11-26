      SUBROUTINE AR_CLARFG( N, ALPHA, X, INCX, TAU )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      INTEGER            INCX, N
      COMPLEX            ALPHA, TAU
*     ..
*     .. Array Arguments ..
      COMPLEX            X( * )
*     ..
*
*  Purpose
*  =======
*
*  AR_CLARFG generates a complex elementary reflector H of order n, such
*  that
*
*        H' * ( alpha ) = ( beta ),   H' * H = I.
*             (   x   )   (   0  )
*
*  where alpha and beta are scalars, with beta real, and x is an
*  (n-1)-element complex vector. H is represented in the form
*
*        H = I - tau * ( 1 ) * ( 1 v' ) ,
*                      ( v )
*
*  where tau is a complex scalar and v is a complex (n-1)-element
*  vector. Note that H is not hermitian.
*
*  If the elements of x are all zero and alpha is real, then tau = 0
*  and H is taken to be the unit matrix.
*
*  Otherwise  1 <= real(tau) <= 2  and  abs(tau-1) <= 1 .
*
*  Arguments
*  =========
*
*  N       (input) INTEGER
*          The order of the elementary reflector.
*
*  ALPHA   (input/output) COMPLEX
*          On entry, the value alpha.
*          On exit, it is overwritten with the value beta.
*
*  X       (input/output) COMPLEX array, dimension
*                         (1+(N-2)*abs(INCX))
*          On entry, the vector x.
*          On exit, it is overwritten with the vector v.
*
*  INCX    (input) INTEGER
*          The increment between elements of X. INCX > 0.
*
*  TAU     (output) COMPLEX
*          The value tau.
*
*  =====================================================================
*
*     .. External Functions ..
      EXTERNAL           CLARFG
*     ..
*     ..
*     .. Executable Statements ..
      CALL CLARFG( N, ALPHA, X, INCX, TAU )
*
*     End of AR_CLARFG
*
      END
