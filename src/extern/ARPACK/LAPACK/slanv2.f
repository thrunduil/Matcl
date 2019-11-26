      SUBROUTINE AR_SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
*
*  -- LAPACK auxiliary routine (version 2.0) --
*     Univ. of Tennessee, Univ. of California Berkeley, NAG Ltd.,
*     Courant Institute, Argonne National Lab, and Rice University
*     September 30, 1994
*
*     .. Scalar Arguments ..
      REAL               A, B, C, CS, D, RT1I, RT1R, RT2I, RT2R, SN
*     ..
*
*  Purpose
*  =======
*
*  AR_SLANV2 computes the Schur factorization of a real 2-by-2 nonsymmetric
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
*  A       (input/output) REAL
*  B       (input/output) REAL
*  C       (input/output) REAL
*  D       (input/output) REAL
*          On entry, the elements of the input matrix.
*          On exit, they are overwritten by the elements of the
*          standardised Schur form.
*
*  RT1R    (output) REAL
*  RT1I    (output) REAL
*  RT2R    (output) REAL
*  RT2I    (output) REAL
*          The real and imaginary parts of the eigenvalues. If the
*          eigenvalues are both real, abs(RT1R) >= abs(RT2R); if the
*          eigenvalues are a complex conjugate pair, RT1I > 0.
*
*  CS      (output) REAL
*  SN      (output) REAL
*          Parameters of the rotation matrix.
*
*  =====================================================================
*
*     .. External Functions ..
      EXTERNAL           SLANV2
*     ..
*     ..
*     .. Executable Statements ..
      CALL SLANV2( A, B, C, D, RT1R, RT1I, RT2R, RT2I, CS, SN )
*
*     End of AR_SLANV2
*
      END