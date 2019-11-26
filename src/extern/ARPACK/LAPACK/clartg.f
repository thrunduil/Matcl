      SUBROUTINE AR_CLARTG( F, G, CS, SN, R )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      REAL               CS
      COMPLEX            F, G, R, SN
*     ..
*
*  Purpose
*  =======
*
*  AR_CLARTG generates a plane rotation so that
*
*     [  CS  SN  ]     [ F ]     [ R ]
*     [  __      ]  .  [   ]  =  [   ]   where CS**2 + |SN|**2 = 1.
*     [ -SN  CS  ]     [ G ]     [ 0 ]
*
*  This is a faster version of the BLAS1 routine CROTG, except for
*  the following differences:
*     F and G are unchanged on return.
*     If G=0, then CS=1 and SN=0.
*     If F=0 and (G .ne. 0), then CS=0 and SN=1 without doing any
*        floating point operations.
*
*  Arguments
*  =========
*
*  F       (input) COMPLEX
*          The first component of vector to be rotated.
*
*  G       (input) COMPLEX
*          The second component of vector to be rotated.
*
*  CS      (output) REAL
*          The cosine of the rotation.
*
*  SN      (output) COMPLEX
*          The sine of the rotation.
*
*  R       (output) COMPLEX
*          The nonzero component of the rotated vector.
*
*  =====================================================================
*
*     .. External Functions ..
      EXTERNAL           CLARTG
*     ..
*     ..
*     .. Executable Statements ..
      CALL CLARTG( F, G, CS, SN, R )
*
*     End of AR_CLARTG
*
      END