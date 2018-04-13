/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"

namespace matcl { namespace lapack
{

template<class T>
void lapack::getf2r(i_type M, i_type* Nptr, i_type JB, i_type J0, T *A, i_type LDA, i_type* IPIV, i_type* IQIV, 
            typename lapack::details::real_type<T>::type TOLC,
            typename lapack::details::real_type<T>::type TOLR, 
            typename lapack::details::real_type<T>::type TOL,
            T* WORK, i_type LWORK, i_type *INFO )
{

    T  ONE     = 1.;
    T  ZERO    = 0.;    
    i_type     J, JP, JQ, JC, JP0, JP2, N = *Nptr;

    using TR    = typename lapack::details::real_type<T>::type;

    --A;
    --IPIV;
    --IQIV;

    bool test = (LWORK == -1);
    
    //  Test the input parameters.
    *INFO = 0;
    if ( M < 0 )
        *INFO = -1;
    else if ( N < 0 )
        * INFO = -2;
    else if ( JB < 0 )
        * INFO = -3;
    else if ( J0 < 0 )
        * INFO = -4;
    else if( LDA < lapack::maximum((i_type)1,M) )
        *INFO = -6;
    else if (TOLC > 1 || TOLC < 0)
        *INFO = -9;
    else if (TOLR > 1 || TOLR < 0)
        *INFO = -10;

    if( *INFO != 0 )
        return;

    i_type LWRK;

    {
        i_type work_size_1	= M + JB * N;
		i_type work_size_2	= M * sizeof(i_type);

		if (work_size_2 % sizeof(T) != 0)
			work_size_2		= work_size_2/sizeof(T) + 1;

        LWRK = work_size_1 + work_size_2;
    };

    if (test)
    {
        WORK[0] = TR(LWRK);
        return;
    };

    if (LWORK < LWRK)
    {
        *INFO = -13;
        return;
    };

    // Quick return if possible
    if( M == 0 || N == 0 || JB == 0)
    {
        *Nptr = 0;
        return;   
    };

    const TR U  = lamch<TR>( "E" );

    TR TOLV     = TOL;

    if (TOLV < 0)
    {
        TR NORM_F       = lapack::lange<T>("F", M, N, A + 1, LDA, nullptr);
        NORM_F          = NORM_F / sqrt(TR(lapack::minimum(M,N)));
        TOLV            = TR(10.0) * U * NORM_F;
    }        

    i_type Nold = N;
    i_type LDR  = N;
    i_type Jmax = lapack::minimum( M, N );
    Jmax        = lapack::minimum( Jmax, JB );

	T* work_C   = (T*)(WORK);
	T* work_R   = work_C + M;    
	i_type* work_I   = (i_type*)(work_R + JB*N); 

    --work_I;
    --work_C;
    --work_R;

    for (J = 1; J <= M; ++J)
        work_I[J] = J;

    for(J = 1; J <= Jmax; ++J)
    {
        //compute J column
        lapack::copy(M,A+1+(J-1)*LDA,1,work_C+1,1);
        if (J > 1)
        {
            //apply deferred interchanges
            lapack::laswp( 1, work_C+1, M, 1, J-1, IPIV+1, 1 );

            lapack::gemv( "No transpose", M-J+1, J-1, -ONE, A+J, LDA, work_R + J, LDR, 
                ONE, work_C+J, 1 );
        };

        // Find pivot
        JC = J;
        JP = J + lapack::amax( M-J+1, work_C+J, 1 ) - 1;

        for(;;)
        {
            //compute JP row
            JP0 = work_I[JP];

            lapack::copy(Nold-J+1,A+JP0+(J-1)*LDA,LDA,work_R+J + (J-1)*LDR,1);
            lapack::gemv( "No Trans", N-J+1, J-1, -ONE, work_R + J, LDR, A+JP, LDA, 
                    ONE, work_R+J + (J-1)*LDR, 1);

            //column pivot candidate
            JQ = J + lapack::amax( N-J+1, work_R+J + (J-1)*LDR, 1 ) - 1;

            if (lapack::abs(work_R[JC+ (J-1)*LDR]) < TOLC * lapack::abs(work_R[JQ+ (J-1)*LDR]))
            {
                // column test failed, check column JQ

                //compute JQ column
                lapack::copy(M,A+1+(JQ-1)*LDA,1,work_C+1,1);

                if (J > 1)
                {
                    //apply deferred interchanges
                    lapack::laswp( 1, work_C+1, M, 1, J-1, IPIV+1, 1 );

                    lapack::gemv( "No transpose", M-J+1, J-1, -ONE, A+J, LDA, work_R + JQ, LDR, 
                                ONE, work_C+J, 1 );
                };

                //row pivot candidate
                JP2 = J + lapack::amax( M-J+1, work_C+J, 1 ) - 1;
                if (lapack::abs(work_C[JP]) < TOLR * lapack::abs(work_C[JP2]))
                {
                    // row test failed, check row JP2
                    JP = JP2;
                    JC = JQ;
                    continue;
                }
                else
                {
                    //pivot (JP, JQ) is acceptable
                    break;
                };
            }
            else
            {
                //pivot (JP, JC) is acceptable
                JQ = JC;
                break;
            };
        };

        //value test
        if (lapack::abs(work_C[JP]) <= TOLV)
        {            
            if (JQ != N)
            {
                //permute upper block
                lapack::swap( J0-1, A+1-(J0-1)+(JQ-1)*LDA, 1, A+1-(J0-1)+(N-1)*LDA, 1 );                

                //permute current block
                lapack::swap( J-1, work_R + JQ, LDR, work_R + N, LDR );
                lapack::swap( M, A+1+(JQ-1)*LDA, 1, A+1+(N-1)*LDA, 1 );                 

                { i_type tmp = IQIV[JQ]; IQIV[JQ] = IQIV[N]; IQIV[N] = tmp; };
            };            
            
            if (J != JQ)
            {
                //permute upper block
                lapack::swap( J0-1, A+1-(J0-1)+(J-1)*LDA, 1, A+1-(J0-1)+(JQ-1)*LDA, 1 );                

                //permute current block
                lapack::swap( J-1, work_R + J, LDR, work_R + JQ, LDR );
                lapack::swap( M, A+1+(J-1)*LDA, 1, A+1+(JQ-1)*LDA, 1 );

                { i_type tmp = IQIV[J]; IQIV[J] = IQIV[JQ]; IQIV[JQ] = tmp; };
            };

            memset(A+1+(N-1)*LDA,0,M*sizeof(T));

            //remove current column
            --J;            
            --N;
            Jmax = lapack::minimum( Jmax, N );
            continue;
        };

        IPIV[J] = JP;
        {
            i_type tmp  = IQIV[J]; 
            IQIV[J]     = IQIV[JQ]; 
            IQIV[JQ]    = tmp;
        };

        i_type tmp  = work_I[J];
        work_I[J]   = work_I[JP];
        work_I[JP]  = tmp;

        // Apply the interchange to rows 1:M.
        if( JQ != J )
        {       
            //update upper block
            lapack::swap( J0-1, A+1-(J0-1)+(J-1)*LDA, 1, A+1-(J0-1)+(JQ-1)*LDA, 1 );

            lapack::swap( J, work_R + J, LDR, work_R + JQ, LDR );
            lapack::copy(M,A+1+(J-1)*LDA,1,A+1+(JQ-1)*LDA,1);
        }

        if (J > 1 || JQ != J)
            lapack::copy(M-J+1,work_C+J,1,A+J+(J-1)*LDA,1);

        if( A[JP+(J-1)*LDA] != ZERO )
        {
            // Apply the interchange to columns 1:J.
            if( JP != J )
                lapack::swap( J, A+J, LDA, A+JP, LDA );

            // Compute elements J+1:M of J-th column.
            if( J < M )
                lapack::scal( M-J, ONE / A[J+(J-1)*LDA], A+J+1+(J-1)*LDA, 1 );
        };
        
        // update A(1:J-1,J)
        lapack::copy(J-1,work_R+J,LDR,A+1 + (J-1)*LDA,1);
    };

    for (; J <= lapack::minimum( M, Nold ); ++J)
        IPIV[J] = J;

    // Apply interchanges to columns Jmax+1:N.
    lapack::laswp( Nold-Jmax, A+1+(Jmax)*LDA, LDA, 1, Jmax, IPIV+1, 1 );

    //update A(1:Jmax,Jmax+1:N)
    for (J = Jmax + 1; J <= Nold; ++J)
        lapack::copy(Jmax,work_R+J,LDR,A+1+(J-1)*LDA,1);

    *Nptr = N;
    return;
};

template void BLAS_EXT_EXPORT
getf2r<s_type>(i_type M, i_type* N, i_type JB, i_type J0, s_type *A, i_type LDA, i_type* IPIV, 
                                     i_type* IQIV, s_type TOLC, s_type TOLR, s_type TOLV,  
                                     s_type* WORK, i_type LWORK, i_type *INFO );

template void BLAS_EXT_EXPORT
getf2r<d_type>(i_type M, i_type* N, i_type JB, i_type J0, d_type *A, i_type LDA, i_type* IPIV, 
                                     i_type* IQIV, d_type TOLC, d_type TOLR, d_type TOLV,  
                                     d_type* WORK, i_type LWORK, i_type *INFO );

template void BLAS_EXT_EXPORT
getf2r<z_type>(i_type M, i_type* N, i_type JB, i_type J0, z_type *A, i_type LDA, i_type* IPIV, 
                                     i_type* IQIV, d_type TOLC, d_type TOLR, d_type TOLV,  
                                     z_type* WORK, i_type LWORK, i_type *INFO );

template void BLAS_EXT_EXPORT
getf2r<c_type>(i_type M, i_type* N, i_type JB, i_type J0, c_type *A, i_type LDA, i_type* IPIV, 
                                     i_type* IQIV, s_type TOLC, s_type TOLR, s_type TOLV,  
                                     c_type* WORK, i_type LWORK, i_type *INFO );
};};
