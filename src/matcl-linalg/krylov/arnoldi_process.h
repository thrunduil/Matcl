/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#pragma once

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl { namespace raw
{

namespace md = matcl::details;

template<class Val>
struct callback_aitr
{    
    //if type = 1 eval  Y = OP * X
    //if type = 2 eval  Y = B * X
    //X is stored in array X
    //Y must be stored in array Y
    //if type = 1 and bmat == 'G', then B * X is stored in BX
    // X, Y and BX are arrays of length N
    virtual void eval(Integer type, Integer N, Val* X, Val* Y, Val* BX) = 0;
};

//---------------------------------------------------------------------------
//                          aitr
//---------------------------------------------------------------------------
/*
*  Apply NP additional steps to a K step Arnoldi factorization.
*
*  Input:  OP*U_{k} - U_{k}*H = r_{k}*e_{k}^T
*
*          with U_{k}^T * B  * U_{k} = I, U_{k}^T * B * r_{k} = 0.
*
*  Output: OP*U_{k+p} - U_{k+p}*H = r_{k+p}*e_{k+p}^T
*
*          with U_{k+p}^T * B * U_{k+p} = I, U_{k+p}^T * B * r_{k+p} = 0.
*
*          U_{k+p} spans the Krylov subspace {x, OP*x, ..., OP^(k+p-1)*x},
*          where x is the initial vector.
*          H is upper hessenberg matrix and H is hermitian if sym == true
*
*  If sym == true, then a linear operator OP must be symmetric with respect 
*  to a positive semi-definite symmetric matrix B, i.e.
*                   
*       B*OP = (OP')*B, or equivalently < x,OPy > = < OPx,y >  
*
*  where < z,w > = z'Bw. Then H is real tridiagonal symmetric matrix.
*
*  Arguments
*  =========
*
*  SYM     Logical (Input)
*          = true   : assume that OP is symmetric with respect to B-product
*          = false  : use nonsymmetric version
*
*  GEN     Logical  (INPUT)
*          BMAT specifies the type of matrix B that defines the
*          semi-inner product for the operator OP.
*          B = false: B = I, i.e. the standard product
*          B = true : B is non unit.
*
*  N       Integer.  (INPUT)
*          Dimension of the eigenproblem. (N >= 0)
*
*  K       Integer.  (INPUT)
*          Current order of H and the number of columns of V. (K >= 0 && K <= N)
*
*  NP      Integer.  (INPUT)
*          Number of additional Arnoldi steps to take. (NP >= 0 && K + NP <= N)
*
*  TOL     Double (INPUT)
*          Stop iterations if |RESID|_2 < TOL * |H|_2
*
*  RNORM   Double precision scalar.  (INPUT/OUTPUT)
*          On INPUT the B-norm of r_{k}.
*          On OUTPUT the B-norm of the updated residual r_{k+p}.
*
*  U       Double precision N by K+NP array.  (INPUT/OUTPUT)
*          On INPUT:  U contains the Arnoldi vectors in the first K 
*          columns.
*          On OUTPUT: U contains the new NP Arnoldi vectors in the next
*          NP columns.  The first K columns are unchanged.
*
*  LDU     Integer.  (INPUT)
*          Leading dimension of U. (LDU >= N)
*
*  H_NS    Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
*          H is used to store the generated upper Hessenberg matrix.
*          Not referred if SYM = true.
*
*  H_S     Double precision (K+NP) by 2 array.  (INPUT/OUTPUT)
*          H_S is used to store the generated symmetric tridiagonal matrix
*          with the subdiagonal in the first column starting at H(2,1)
*          and the main diagonal in the second column. Not referred if 
*          SYM = false.
*
*  LDH     Integer.  (INPUT)
*          Leading dimension of H_NS or H_S. (LDH >= K + NP)
*
*  SCALE   Double (INPUT/OUTPUT)
*  SUMSQ   Double (INPUT/OUTPUT)
*          SCALE and SUMSQ defines the Frobenius norm of the matrix H_NS or H_S given by
*          |H|_F = SCALE * sqrt(SUMSQ)
*
*  WORK    Double precision work array of length 4*N.
*          On INPUT:  WORK(1:N) contains the residual vector r_{k}.
*          residual vector must be orthogonal againts U. One may set r_k different
*          from output from previous call to this function if invariant subspace
*          of OP was found.
*          On OUTPUT: WORK(1:N) contains the residual vector r_{k+p}.
*
*  INFO    Integer.  (OUTPUT)
*          < 0: -INFO argument has invalid value
*          = 0: Normal exit.
*          > 0: Size of an invariant subspace of OP is found that is
*               less than K + NP.
* 
*  References:
*  =========
*
*  1. D.C. Sorensen, "Implicit Application of Polynomial Filters in
*     a k-Step Arnoldi Method", SIAM J. Matr. Anal. Apps., 13 (1992),
*     pp 357-385.
*  2. R.B. Lehoucq, "Analysis and Implementation of an Implicitly 
*     Restarted Arnoldi Iteration", Rice University Technical Report
*     TR95-13, Department of Computational and Applied Mathematics.
*
*/
template<class V>
void aitr(bool SYM, bool GEN, Integer N, Integer K, Integer NP, 
    typename md::real_type<V>::type TOL, typename md::real_type<V>::type& RNORM, 
    V* U, Integer LDU, V* H_NS, typename md::real_type<V>::type* H_S, Integer LDH, 
    typename md::real_type<V>::type& SCALE, typename md::real_type<V>::type& SUMSQ,
    V* WORK, callback_aitr<V>* callback, Integer& INFO);

//---------------------------------------------------------------------------
//                          baitr
//---------------------------------------------------------------------------
template<class Val>
struct callback_baitr
{    
    //eval Y = OP * X, where X and Y are MxN matrices
    virtual void eval(Integer M, Integer N, Val* X, Integer LDX, Val* Y, Integer LDY) = 0;
};

/*
*  generate additional p Arnoldi vectors to a K vector Block Arnoldi factorization.
*
*  Input:  OP*U_{k} - U_{k}*H = r_{k}*e_{k}^T
*
*          with U_{k}^T * U_{k} = I, U_{k}^T * r_{k} = 0.
*
*  Output: OP*U_{k+p} - U_{k+p}*H = r_{k+p}*e_{k+p}^T
*
*          with U_{k+p}^T * U_{k+p} = I, U_{k+p}^T * r_{k+p} = 0.
*
*          U_{k+p} spans the Krylov subspace {x, OP*x, ..., OP^((k+p)/kb - 1)*x}
*          where x is the initial matrix of size N x KB
*
*          H is block upper hessenberg matrix with at most 2 * KB - 1 nonzero subdiagonals
*
*  Arguments
*  =========
*
*  SYM     Logical (Input)
*          = true   : assume that OP is symmetric
*          = false  : use nonsymmetric version
*
*  N       Integer.  (INPUT)
*          Dimension of the eigenproblem. (N >= 0)
*
*  K       Integer.  (INPUT)
*          Current order of H and the number of columns of V. (K >= 0 && K <= N)
*
*  NP      Integer.  (INPUT)
*          Number of additional Arnoldi vectors to generate. (NP >= 0 && K + NP <= N)
*
*  KB      Integer. (INPUT/OUTPUT)
*          On INPUT:  Number of residual vectors. (KB >= 1)
*          On OUTPUT: Number of residual vectors.
*           
*  TOL     Double (INPUT)
*          Stop iterations if |RESID|_2 < TOL * |H|_2
*
*  X       Double precision N by KB array.  (INPUT/OUTPUT)
*          On INPUT:  X(1:N,1:KB) contains the initial vector x_{k}.
*          One may set X different from output from previous call to this function
*          if invariant subspace of OP was found.
*          On OUTPUT: X(1:N,1:KB) contains the initial vector for next interations.
*
*  LDX     Integer.  (INPUT)
*          Leading dimension of X. (LDX >= N)
*
*  INIT    Logical. (Input)
*          Set INIT = true on first call and any subsequent call if X is different 
*          that the result from the previos call. Otherwise set false.
*
*  RESID   Double precision N by KB array.  (OUTPUT)
*          On OUTPUT: contains the residual vector.
*
*  LDR     Integer.  (INPUT)
*          Leading dimension of RESID. (LDR >= N)
*
*  NORMR   Double (OUTPUT)
*          On OUTPUT: the Frobenius norm of RESID
*
*  U       Double precision N by K+NP array.  (INPUT/OUTPUT)
*          On INPUT:  U contains the Arnoldi vectors in the first K 
*          columns.
*          On OUTPUT: U contains the new NP Arnoldi vectors in the next
*          NP columns.  The first K columns are unchanged.
*          Arnoldi vectors are represented as elementary reflectors.
*
*  LDU     Integer.  (INPUT)
*          Leading dimension of U. (U >= N)
*
*  TAU     Double precision K + NP array.  (INPUT/OUTPUT) 
*          On INPUT: elements 1:K stores the scalar factors of the elementary
*          reflectors.
*          On OUTPUT: elements 1:K+NP stores the scalar factors of the elementary
*          reflectors.
*
*  H_NS    Double precision (K+NP) by (K+NP) array.  (INPUT/OUTPUT)
*          H is used to store the generated block upper Hessenberg matrix with 
*          at most 2 * KB - 1 nonzero subdiagonals.
*          Not referred if SYM = true.
*
*  H_S     Double precision min(K+NK,2*KB) * (K+NP) array.  (INPUT/OUTPUT)
*          H_S is used to store the generated hermitian matrix with at most min(K+NK,2*KB)
*          subdiagonals. H_S(1,:) stores the main diagonal, H_S(i,:) stores (i-1)-th
*          subdiagonal. First element in given diagonal is stored in the first column.
*          Not referred if SYM = false.
*
*  LDH     Integer.  (INPUT)
*          Leading dimension of H_NS or H_S, (LDH >= K + NP if SYM == false and 
*          LDH >= min(K+NK,2*KB) if SYM == true).
*
*  SCALE   Double (INPUT/OUTPUT)
*  SUMSQ   Double (INPUT/OUTPUT)
*          SCALE and SUMSQ defines the Frobenius norm of the matrix H_NS or H_S given by
*          |H|_F = SCALE * sqrt(SUMSQ)
*
*  T       Double precision (LDT, NTB) array.  (INPUT/OUTPUT) 
*          NTB = (K + NP)/NB. NB is returned by baitr_block_size()
*          On INPUT: accumulated NT triangular factors of size NB*NB of the block reflectors.
*          On OUTPUT: accumulated NT triangular factors of size NB*NB of the block reflectors.
*
*  LDT     Integer.  (INPUT)
*          Leading dimension of T, LDT >= NB*NB. NB is returned by baitr_block_size().
*
*  NT      Integer. (INPUT/OUTPUT) 
*          On INPUT: Number of accumulated triangular factors, should be set to 0 on the first
*          call.
*          On OUTPUT: Number of accumulated triangular factors.
*       
*  WORK    Double precision work array of length max(1,LWORK)
*          On exit, if INFO=0, WORK(1) returns the optimal LWORK.
*
*  LWORK   Integer. (INPUT)
*          The dimension of the array WORK. If LWORK = -1, then a workspace 
*          query is assumed; the routine only calculates the optimal size of
*          the WORK array, returns this value as the first entry of the WORK 
*          array.
*
*  INFO    Integer.  (OUTPUT)
*          < 0: -INFO argument has invalid value
*          = 0: Normal exit.
*          > 0: Size of an invariant subspace of OP is found that is
*               less than K + NP.
* 
*  References:
*  =========
*
*/

template<class V>
void baitr(bool SYM, Integer N, Integer K, Integer NP, Integer& KB, typename md::real_type<V>::type TOL, 
    V* X, Integer LDX, bool INIT, V* RESID, Integer LDR, typename md::real_type<V>::type& NORMR, V* U, 
    Integer LDU, V* TAU, V* H_NS, V* H_S, Integer LDH, typename md::real_type<V>::type& SCALE, 
    typename md::real_type<V>::type& SUMSQ, V* T, Integer LDT, Integer& NT, V* WORK, Integer LWORK, 
    callback_baitr<V>* callback, Integer& INFO);

// return block size used by baitr<V>
template<class V>
Integer baitr_block_size();

/*
*  Purpose
*  =======
*
*  apply_U overwrites the general real M-by-N matrix C with
*
*                  SIDE = 'L'
*  TRANS = 'N':      U * C   
*  TRANS = 'T':      U**T * C
*
*  where U is a real orthogonal matrix defined as the product of k
*  elementary reflectors
*
*        Q = H(1) H(2) . . . H(k)
*
*  as returned by baitr. U is of order M.
*
*  Arguments
*  =========
*
*  TRANS   (input) CHARACTER*1
*          = 'N':  No transpose, apply U;
*          = 'T':  Transpose, apply U**T.
*
*  M       (input) INTEGER
*          The number of rows of the matrix C. M >= 0.
*
*  N       (input) INTEGER
*          The number of columns of the matrix C. N >= 0.
*
*  K       (input) INTEGER
*          The number of elementary reflectors whose product defines
*          the matrix Q.
*          If SIDE = 'L', M >= K >= 0;
*          if SIDE = 'R', N >= K >= 0.
*
*  U       (input) DOUBLE PRECISION M by K array.
*          U contains the Arnoldi vectors in the first K columns returned
*          by baitr.
*
*  LDU     (input) INTEGER
*          Leading dimension of U. (U >= M)
*
*  TAU     (input) DOUBLE PRECISION array, dimension (K)
*          TAU(i) must contain the scalar factor of the elementary
*          reflector H(i), as returned by baitr.
*
*  C       (input/output) DOUBLE PRECISION array, dimension (LDC,N)
*          On entry, the M-by-N matrix C.
*          On exit, C is overwritten by U*C or U**T*C.
*
*  LDC     (input) INTEGER
*          The leading dimension of the array C. LDC >= max(1,M).
*
*  T       (input/output) Double precision (LDT, NTB) array.
*          NTB = (K + NP)/NB. NB is returned by baitr_block_size()
*          On INPUT: accumulated NT triangular factors of size NB*NB of the block reflectors.
*          On OUTPUT: accumulated NT triangular factors of size NB*NB of the block reflectors.
*
*  LDT     (input) INTEGER
*          Leading dimension of T, LDT >= NB*NB. NB is returned by baitr_block_size().
*
*  NT      (input/output) INTEGER
*          On INPUT: Number of accumulated triangular factors.
*          On OUTPUT: Number of accumulated triangular factors.
*       
*  WORK    (workspace/output) DOUBLE PRECISION array, dimension (MAX(1,LWORK))
*          Minimum LWORK = N; optimal LWORK = N * NB, where NB = baitr_block_size()
*
*  LDWORK  (input) INTEGER
*          The dimension of the array WORK.
*/
template<class V>
void apply_U(const char* TRANS, Integer M, Integer N, Integer K, V* U, Integer LDU, const V* TAU,
        V* X, Integer LDX, V* T, Integer LDT, Integer& NT, Integer NTB, V* WORK, Integer LDWORK);

}}