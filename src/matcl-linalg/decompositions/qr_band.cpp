/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

#include "matcl-linalg/decompositions/qr_band.h"
#include <algorithm>

namespace matcl { namespace lapack
{

template<class V>   
struct get_trans_type{};

template<>
struct get_trans_type<d_type> { static const char* eval() { return "trans"; }; };

template<>
struct get_trans_type<s_type> { static const char* eval() { return "trans"; }; };

template<>
struct get_trans_type<c_type> { static const char* eval() { return "conj trans"; }; };

template<>
struct get_trans_type<z_type> { static const char* eval() { return "conj trans"; }; };

//=======================   gbbqrf    ================================

/*********************************************************************
c
cc DGBBQRF QR factors an M by N band matrix in GB format, using blocking.
c
c  Discussion:
c
c    DGBBQRF computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a packed matrix: dense A is M by N, whereas
c    band A is (2*ML+MU+NB) by N.
c    Input matrix AB has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML+MU+NB-1 superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, MU = 1, and NB = 2.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   r   0
c      x  x  x  0  0  0       v1  r   r   r   r   r
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    03 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer NB, the block size to use.
c    1 <= NB.
c
c    Input, integer M, the number of rows in the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns in the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision A(LDA,N).
c    A is stored as a packed matrix: dense A is M by N, whereas
c    band A is (2*ML+MU+NB) by N.
c    On entry, the M by N matrix A in band storage in rows ML+NB to
c    2*ML+MU+NB; rows 1 to ML+NB-1 of the array need not be set.
c    This matrix has ML subdiagonals and MU superdiagonals with data.
c    On exit, the elements on and above the diagonal of the array
c    contain the upper band matrix R.
c    The output matrix R has ML+MU+NB-1 superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of A.
c    max ( 1, M ) <= LDA.
c
c    Output, double precision TAU(min(M,N)), the scalar factors of the 
c    elementary reflectors.
c
c    Workspace, double precision WORK(DIMWORK).
c    NB*NB + MIN(N,ML+MU)*NB + MIN(M,ML+NB)*NB <= DIMWORK.
c
c    Output, integer INFO, error flag.
c    0, no errors.
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
c
*/

template<class Val>
void lapack::gbbqrf(i_type nb, i_type m, i_type n, i_type ml, i_type mu, Val* a, 
                i_type lda, Val* tau, Val* work, i_type& info )
{
    //  Test the input arguments.
    info    = 0;

    if ( nb < 1 )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ml < 0 )
        info = -4;
    else if ( mu < 0 )
        info = -5;
    else if ( lda < std::max(1, 2 * ml + mu + nb ) )
        info = -7;

    if ( info != 0 )
        return;

    //  Quick return if possible.
    if ( m == 0 || n == 0 || ml == 0 )
        return;

    //  Adjust matrix sizes: AM, AN, AML, AMU, and ANB.
    i_type am   = std::min(m, n + ml );
    i_type an   = std::min(n, m + mu );
    i_type aml  = std::min(m - 1, ml );
    i_type amu  = std::min(n - 1, mu );
    i_type anb  = nb;
    i_type ii, jj;

    if ( aml < anb )
        anb = aml;

    i_type mn   = std::min(am, an);

    i_type it   = 1;
    i_type iv   = it + anb * anb;
    i_type irwk = iv + std::min (aml + anb, am) * anb;

    i_type ncol = std::min(an, am - aml);

    const char* trans_t     = get_trans_type<Val>::eval();

    //  Factorization of full band matrix A(:,1:ncol).
    for (i_type j = 1; j <= ncol; j += anb)
    {
        i_type jb   = std::min( anb, ncol - j + 1 );

        //  Factorize block A(j:j+aml+jb-1,j:j+jb-1).
        lapack::gebqr2(aml + jb, jb, aml, amu, &a[ml+mu+nb-1 + (j-1)*lda], lda - 1, 
                        &tau[j-1], &work[irwk-1], info );

        if ( info != 0 )
        {
            info = info - 1000;
            return;
        };     

        if ( j + jb <= an )
        {
            //  WORK(IV) := Y (the lower part of Y is padded with zeros).
            for (jj = 1; jj <= jb; ++jj)
            {
                for (ii = 1; ii <= aml + jj; ++ii)
                {
                    work[iv-1+(aml+jb)*(jj-1)+ii-1] = a[ml+mu+nb+ii-jj-1 + (j+jj-1-1)*lda];
                };

                for (ii = aml + jj + 1; ii <= aml + jb; ++ii)
                {
                    work[iv-1+(aml+jb)*(jj-1)+ii-1] = Val(0.0);
                };
            };

            //  Form the triangular factor T of the block reflector.
            lapack::larft("forward", "columnwise", aml + jb, jb, &work[iv-1], aml + jb, 
                  &tau[j-1], &work[it-1], jb );

            //  Apply block reflector to A(j:j+aml+jb-1,j+jb:an) from the left.
            i_type jlc = std::min( an, j + jb - 1 + aml + amu );

            lapack::larfb("left", trans_t, "forward", "columnwise", aml + jb, jlc - j - jb + 1, jb,
                &work[iv-1], aml + jb, &work[it-1], jb, &a[ml+mu+nb-jb-1 +(j+jb-1)*lda], lda - 1,
                &work[irwk-1], jlc - j - jb + 1 );
        }
    };

    //  Factorization of rectangular matrix A(:,ncol+1:mn).
    for (i_type j = ncol + 1; j <= mn; j += anb)
    {
        i_type jb = std::min( anb, mn - j + 1 );

        //  Factorize block A(j:am,j:j+jb-1).
        geqr2( am - j + 1, jb, &a[ml+mu+nb-1 +(j-1)*lda], lda - 1, &tau[j-1], &work[irwk-1], info);

        if ( info != 0 )
        {
            info = info - 1000;
            return;
        };

        if ( j + jb <= an )
        {
            //  Form the triangular factor T of the block reflector.
            lapack::larft("forward", "columnwise", am - j + 1, jb,
                &a[ml+mu+nb-1+(j-1)*lda], lda - 1, &tau[j-1], &work[it-1], jb );

            //  Apply block reflector to A(j:am,j+jb:an) from the left.
            lapack::larfb("left", trans_t, "forward", "columnwise", am - j + 1, an - j - jb + 1, jb,
                &a[ml+mu+nb-1 + (j-1)*lda], lda - 1, &work[it-1], jb, &a[ml+mu+nb-jb-1 + (j+jb-1)*lda], 
                lda - 1, &work[irwk-1], an - j - jb + 1 );
        };
    };

    return;
};

//dgbbqr2
// from band_qr DGBBQR2 QR factors an M by N band matrix in banded format, with no blocking. 
template<class Val>
void gbbqr2(i_type M, i_type N, i_type ML, i_type MU, Val* A, 
                i_type LDA, Val* tau, Val* WORK, i_type& INFO);

template void lapack::gbbqrf(i_type nb, i_type m, i_type n, i_type ml, i_type mu, d_type* a, 
                i_type lda, d_type* tau, d_type* work, i_type& info );
template void lapack::gbbqrf(i_type nb, i_type m, i_type n, i_type ml, i_type mu, s_type* a, 
                i_type lda, s_type* tau, s_type* work, i_type& info );
template void lapack::gbbqrf(i_type nb, i_type m, i_type n, i_type ml, i_type mu, c_type* a, 
                i_type lda, c_type* tau, c_type* work, i_type& info );
template void lapack::gbbqrf(i_type nb, i_type m, i_type n, i_type ml, i_type mu, z_type* a, 
                i_type lda, z_type* tau, z_type* work, i_type& info );

//=======================   DGBBQR2    ================================
/*
c*********************************************************************72
c
cc DGBBQR2 QR factors an M by N band matrix in GB format, without blocking.
c
c  Discussion:
c
c    DGBBQR2 computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a packed band matrix.
c    Input matrix A has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML + MU superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, and MU = 1.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   0   0
c      x  x  x  0  0  0       v1  r   r   r   r   0
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision A(LDA,N).
c    A is stored as a packed matrix: dense A is M by N, whereas
c    band A is (2*ML+MU+1) by N.
c    On entry, the M by N matrix AB in band storage in rows ML+1 to
c    2*ML+MU+1; rows 1 to ML of the array need not be set.
c    This matrix has ML subdiagonals and MU superdiagonals with data.
c    On exit, the elements on and above the diagonal of the array
c    contain the upper band matrix R.
c    The output matrix R has ML+MU superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of A.
c    2*ML+MU+1 <= LDA.
c
c    Output, double precision TAU(min(M,N)), the scalar factors of the 
c    elementary reflectors.
c
c    Workspace, double precision WORK(min(N,MU+ML)).
c
c    Output, integer INFO, error flag.
c    0: successful exit
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
c
*/
template<class Val>
void lapack::gbbqr2(i_type m, i_type n, i_type ml, i_type mu, Val* a, 
                i_type lda, Val* tau, Val* work, i_type& info)
{
    //  Test the input arguments.
    info = 0;

    if ( m < 0 )
        info = - 1;
    else if ( n < 0 )
        info = - 2;
    else if ( ml < 0 )
        info = - 3;
    else if ( mu < 0 )
        info = - 4;
    else if ( lda < 2 * ml + mu + 1 )
        info = - 6;

    if ( info != 0 )
        return;

    //  Quick return if possible.
    i_type mn = std::min( m, n );

    if ( mn == 0 )
        return;

    //  Set rows 1:ML to zero.
    lapack::laset( "all", ml, n, Val(0.0), Val(0.0), a, lda );

    for (i_type j = 1; j <= mn; ++j)
    {
        i_type mh   = std::min( ml + 1, m - j + 1 );

        //  Generate reflector H(j) to annihilate A(j+1:j+mh-1,j).
        lapack::larfg( mh, &a[ml+mu +(j-1)*lda], &a[ml+mu+1+std::min(1,ml)-1 + (j-1)*lda], 1, 
                      &tau[j-1]);

        i_type nh   = std::min( n - j, mu + ml );

        //  Apply reflector H(j) to rest of matrix: A(j:j+mh-1,j+1:j+nh).
        if ( 0 < nh )
        {
            Val diag = a[ml+mu + (j-1)*lda];

            a[ml+mu + (j-1)*lda]  = Val(1.0);

            Val tau_conj            = conj(tau[j-1]);

            lapack::larf("left", mh, nh, &a[ml+mu +(j-1)*lda], 1, &tau_conj,
                &a[ml+mu-1 + j*lda], lda - 1, work);

            a[ml+mu + (j-1)*lda]  = diag;
        };
    };

    return;
};

template void gbbqr2(i_type m, i_type n, i_type ml, i_type mu, c_type* a, 
                i_type lda, c_type* tau, c_type* work, i_type& info);
template void gbbqr2(i_type m, i_type n, i_type ml, i_type mu, d_type* a, 
                i_type lda, d_type* tau, d_type* work, i_type& info);
template void gbbqr2(i_type m, i_type n, i_type ml, i_type mu, s_type* a, 
                i_type lda, s_type* tau, s_type* work, i_type& info);
template void gbbqr2(i_type m, i_type n, i_type ml, i_type mu, z_type* a, 
                i_type lda, z_type* tau, z_type* work, i_type& info);

//=======================   gebqr2    ================================
/*
c
cc DGEBQR2 QR factors an M by N band matrix in GE format, with no blocking.
c
c  Discussion:
c
c    DGEBQR2 computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a general matrix.
c    Input matrix A has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML+MU superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, and MU = 1.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   0   0
c      x  x  x  0  0  0       v1  r   r   r   r   0
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c    The matrix Q is represented as a product of elementary reflectors
c
c      Q = H(1) H(2) . . . H(k), where k = min(m,n).
c
c    Each H(i) has the form
c
c      H(i) = I - tau * v * v'
c
c    where tau is a real scalar, and v is a real vector with
c    v(1:i-1) = 0 and v(i) = 1; v(i+1:m) is stored on exit in A(i+1:m,i),
c    and tau in TAU(i).
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer M, the number of rows of the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision(LDA,N).
c    On entry, the M by N matrix A. It has ML subdiagonals and MU
c    superdiagonals.
c    On exit, the elements on and above the diagonal of the array
c    contain the min(M,N) by N upper band matrix R.
c    The output matrix R has ML+MU superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of the array A.  
c    max ( 1, M ) <= LDA.
c
c    Output, double precision TAU(min(M,N)),
c    the scalar factors of the elementary reflectors.
c
c    Workspace, double precision WORK(min(N,MU+ML)).
c
c    Output, integer INFO, error flag.
c    0, no errors.
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
*/
template<class Val>
void gebqr2(i_type m, i_type n, i_type ml, i_type mu, Val* a, i_type lda, Val* tau, 
            Val* work, i_type& info )
{
    //  Test the input arguments.
    info = 0;

    if ( m < 0 )
        info = -1;
    else if ( n < 0 )
        info = -2;
    else if ( ml < 0 )
        info = -3;
    else if ( mu < 0 )
        info = -4;
    else if ( lda < std::max ( 1, m ) )
        info = -6;

    if ( info != 0 )
        return;

    //  Quick return if possible.
    i_type mn = std::min( m, n );
    if ( mn == 0 )
        return;

    for (i_type j = 1; j <= mn; ++j)
    {
        i_type mh   = std::min( ml + 1, m - j + 1 );

        //  Generate reflector H(j) to annihilate A(j+1:j+mh-1,j).
        lapack::larfg( mh, &a[j-1 + (j-1)*lda], &a[std::min(j+1,m) - 1 + (j-1)*lda], 
                      1, &tau[j-1]);

        i_type nh   = std::min( n - j, mu + ml );

        if ( 0 < nh )
        {
            //  Apply reflector H(j) to rest of matrix: A(j:j+mh-1,j+1:j+nh).
            Val diag = a[j-1 + (j-1)*lda];

            a[j-1 +(j-1)*lda]   = Val(1.0);
            Val tau_conj        = conj(tau[j-1]);

            lapack::larf("left", mh, nh, &a[j-1 +(j-1)*lda], 1, &tau_conj, 
                         &a[j-1 + j * lda], lda, work);

            a[j-1 + (j-1)*lda]  = diag;
        };        
    };

    return;
};

template void gebqr2(i_type, i_type, i_type, i_type, d_type* a, i_type, d_type*, 
                        d_type*, i_type&);
template void gebqr2(i_type, i_type, i_type, i_type, s_type* a, i_type, s_type*, 
                        s_type*, i_type&);
template void gebqr2(i_type, i_type, i_type, i_type, c_type* a, i_type, c_type*, 
                        c_type*, i_type&);
template void gebqr2(i_type, i_type, i_type, i_type, z_type* a, i_type, z_type*, 
                        z_type*, i_type&);

//=======================   gebqrf    ================================
/*
cc DGEBQRF QR factors an M by N band matrix stored in GE format, with blocking.
c
c  Discussion:
c
c    DGEBQRF computes a QR factorization of a real M by N band matrix A
c    with lower band ML and upper band MU: A = Q * R.
c
c    A is stored as a general matrix.
c    Input matrix A has ML subdiagonals and MU superdiagonals.
c    Output matrix R has ML+MU+NB-1 superdiagonals.
c
c    Example of A with M = 8, N = 6, ML = 2, MU = 1, and NB = 2.
c    Left, input matrix; Right, output matrix.
c
c      x  x  0  0  0  0       r   r   r   r   r   0
c      x  x  x  0  0  0       v1  r   r   r   r   r
c      x  x  x  x  0  0       v1  v2  r   r   r   r
c      0  x  x  x  x  0       0   v2  v3  r   r   r
c      0  0  x  x  x  x       0   0   v3  v4  r   r
c      0  0  0  x  x  x       0   0   0   v4  v5  r
c      0  0  0  0  x  x       0   0   0   0   v5  v6
c      0  0  0  0  0  x       0   0   0   0   0   v6
c
c  Licensing:
c
c    Copyright (c) 2007, Universidad Jaume I de Castellon
c    All rights reserved.
c
c    Redistribution and use in source and binary forms, with or without
c    modification, are permitted provided that the following conditions 
c    are met:
c      * Redistributions of source code must retain the above copyright
c        notice, this list of conditions and the following disclaimer.
c      * Redistributions in binary form must reproduce the above copyright
c        notice, this list of conditions and the following disclaimer in the
c        documentation and/or other materials provided with the distribution.
c      * Neither the name of the <organization> nor the
c        names of its contributors may be used to endorse or promote 
c        products derived from this software without specific prior written 
c        permission.
c
c    This software is provided by <copyright holder> ''as is'' and any
c    express or implied warranties, including, but not limited to, the 
c    implied warranties of merchantability and fitness for a particular 
c    purpose are disclaimed.  In no event shall <copyright holder> be liable 
c    for any direct, indirect, incidental, special, exemplary, or 
c    consequential damages (including, but not limited to, procurement of 
c    substitute goods or services; loss of use, data, or profits; or 
c    business interruption) however caused and on any theory of liability, 
c    whether in contract, strict liability, or tort (including negligence 
c    or otherwise) arising in any way out of the use of this software, even 
c    if advised of the possibility of such damage.
c
c  Modified:
c
c    04 April 2010
c
c  Author:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    Dept. de Ingenieria y Ciencia de Computadores, Univ. Jaume I
c    12.080 Castellon, Spain
c    {gquintan,remon,quintana}@icc.uji.es
c
c  Reference:
c
c    Alfredo Remon, Enrique Quintana-Orti, Gregorio Quintana-Orti,
c    LAPACK-Style Codes for the QR Factorization of Banded Matrices,
c    To Appear.
c
c  Parameters:
c
c    Input, integer NB, the block size.
c    1 <= NB.
c
c    Input, integer M, the number of rows of the matrix.
c    0 <= M.
c
c    Input, integer N, the number of columns of the matrix.
c    0 <= N.
c
c    Input, integer ML, the number of nonzero subdiagonals.
c    0 <= ML.
c
c    Input, integer MU, the number of nonzero superdiagonals.
c    0 <= MU.
c
c    Input/output, double precision A(LDA,N).
c    On entry, the M by N matrix A. It has ML subdiagonals and MU
c    superdiagonals.
c    On exit, the elements on and above the diagonal of the array
c    contain the min(M,N) by N upper band matrix R.
c    The output matrix R has ML+MU+NB-1 superdiagonals.
c    The elements within the band below the diagonal, with the array
c    TAU, represent the orthogonal matrix Q as a product of
c    elementary reflectors.
c
c    Input, integer LDA, the leading dimension of the array A.  
c    max ( 1, M ) <= LDA.
c
c    Output, double precision TAU(min(M,N)), the scalar factors of the 
c    elementary reflectors.
c
c    Workspace, double precision WORK(DIMWORK).
c    NB*NB + min(N,ML+MU)*NB + min(M,ML+NB)*NB <= DIMWORK.
c
c    Output, integer INFO, error flag.
c    0, no errors.
c    nonzero, if INFO = -I, the I-th argument had an illegal value.
*/
template<class Val>
void lapack::gebqrf(i_type nb, i_type m, i_type n, i_type ml, i_type mu, Val *a, i_type lda, 
            Val* tau, Val* work, i_type& info )
{
    //  Test the input arguments.
    info = 0;

    if ( nb < 1 )
        info = -1;
    else if ( m < 0 )
        info = -2;
    else if ( n < 0 )
        info = -3;
    else if ( ml < 0 )
        info = -4;
    else if ( mu < 0 )
        info = -5;
    else if ( lda < std::max( 1, m ) )
        info = -7;

    if ( info != 0 )
        return;

    //  Quick return if possible.
    if ( m == 0 || n == 0 || ml == 0 )
        return;

    //  Adjust matrix sizes: AM, AN, AML, AMU, and ANB.
    i_type am   = std::min( m, n + ml );
    i_type an   = std::min( n, m + mu );
    i_type aml  = std::min( m - 1, ml );
    i_type amu  = std::min( n - 1, mu );
    i_type anb  = nb;

    if ( aml < anb )
        anb     = aml;

    i_type mn   = std::min ( am, an );

    i_type it   = 1;
    i_type iv   = it + anb * anb;
    i_type irwk = iv + std::min (aml + anb, am )*anb;

    i_type ncol = std::min( an, am - aml );

    const char* trans_t     = get_trans_type<Val>::eval();

    //  Factorization of full band matrix A(:,1:ncol).
    for (i_type j = 1; j <= ncol; j += anb)
    {
        i_type jb   = std::min( anb, ncol - j + 1 );

        //  Factorize block A(j:j+aml+jb-1,j:j+jb-1).
        lapack::gebqr2( aml + jb, jb, aml, amu, &a[j-1 + (j-1)*lda], lda, &tau[j-1],
                &work[irwk-1], info );

        if ( info != 0 )
        {
            info = info - 1000;
            return;
        };

        if ( j + jb <= an )
        {
            //  WORK(IV) := Y (the lower part of Y is padded with zeros).
            for (i_type jj = 1; jj <= jb; ++jj)
            {
                for (i_type ii = 1; ii <= aml + jj; ++ii)
                {
                    work[iv-1+(aml+jb)*(jj-1)+ii-1] = a[j+ii-1 -1 + (j+jj-1) * lda];
                };

                for (i_type ii = aml + jj + 1; ii <= aml + jb; ++ii)
                {
                    work[iv-1+(aml+jb)*(jj-1)+ii-1] = Val(0.0);
                };
            };

            //  Form the triangular factor T of the block reflector.
            lapack::larft("forward", "columnwise", aml + jb, jb, &work[iv-1], aml + jb, 
                          &tau[j-1], &work[it-1], jb );

            //  Apply block reflector to A(j:j+aml+jb-1,j+jb:an) from the left.
            i_type jlc  = std::min( an, j + jb - 1 + aml + amu );

            lapack::larfb("left", trans_t, "forward", "columnwise", aml + jb, 
                          jlc - j - jb + 1, jb, &work[iv-1], aml + jb, &work[it-1], jb,
                          &a[j-1 +(j+jb-1)*lda], lda, &work[irwk-1], jlc - j - jb + 1 );
        };
    };

    //  Factorization of rectangular matrix A(:,ncol+1:mn).
    for (i_type j = ncol + 1; j <= mn; j += anb)
    {
        i_type jb = std::min( anb, mn - j + 1 );
    
        //  Factorize block A(j:am,j:j+jb-1).
        lapack::geqr2( am - j + 1, jb, &a[j-1 + (j-1)*lda], lda, &tau[j-1], 
                      &work[irwk-1], info );

        if ( info != 0 )
        {
            info = info - 1000;
            return;
        };

        if ( j + jb <= an )
        {
            //  Form the triangular factor T of the block reflector.
            lapack::larft("forward", "columnwise", am - j + 1, jb, &a[j-1 +(j-1)*lda], 
                          lda, &tau[j-1], &work[it-1], jb );

            //  Apply block reflector to A(j:am,j+jb:an) from the left.
            lapack::larfb("left", "transpose", "forward", "columnwise", 
                          am - j + 1, an - j - jb + 1, jb, &a[j-1 +(j-1)*lda], lda, &work[it-1],
                          jb, &a[j-1 + (j+jb-1)*lda], lda, &work[irwk-1], an - j - jb + 1 );

        };
    };

    return;
};
    
template void gebqrf(i_type, i_type, i_type, i_type, i_type, d_type*, i_type, 
                        d_type*, d_type*, i_type&);
template void gebqrf(i_type, i_type, i_type, i_type, i_type, s_type*, i_type, 
                        s_type*, s_type*, i_type&);
template void gebqrf(i_type, i_type, i_type, i_type, i_type, c_type*, i_type, 
                        c_type*, c_type*, i_type&);
template void gebqrf(i_type, i_type, i_type, i_type, i_type, z_type*, i_type, 
                        z_type*, z_type*, i_type&);

};};