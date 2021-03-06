/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2021
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
//  SYNTAX:
//
//      void perm2int(int_type M,int_type *IPIV, int_type *WORK, bool ZB)
//
//  PURPOSE:
//  -----------------------------------------------------------------------
//
//  transforms permutation vector IPIV to series of interchanges
//
//  ARGUMENTS:
//  -----------------------------------------------------------------------
//  M       (input) INTEGER
//          The number of rows of the vector IPIV.  M >= 0.
//
//  IPIV    (input/output) INTEGER array, dimension M
//          On entry IPIV is a vector of permutations of numbers 0:1:M-1 (1:1:M
//          if ZB is false)
//          On exit IPIV is replaced by vector of interchanges of a vector
//          0:1:M-1 (1:1:M) that produces given vector of permutations
//
//  WORK    (input) INTEGER array, dimension at least M
//  ZB      if true then IPIV contains zero-based indices, otherwise IPIV contains
//          1-based indices
void perm2int(i_type M,i_type *IPIV, i_type* WORK, bool ZB)
{
    if (M<=0)
        return;
    
    if (ZB == true)
    {
        for(i_type i=0; i<M;i++)
            WORK[IPIV[i]]	= i;

        for (i_type i=0;i<M;i++)
        {
            IPIV[WORK[i]]	= IPIV[i];
            WORK[IPIV[i]]	= WORK[i];
        };
    }
    else
    {
        for(i_type i = 0; i < M; i++)
            WORK[IPIV[i] - 1]	= i + 1;

        for (i_type i = 0; i < M; i++)
        {
            IPIV[WORK[i] - 1]	= IPIV[i];
            WORK[IPIV[i] - 1]	= WORK[i];
        };
    }
};

};};