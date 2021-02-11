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

#include "special_cases.h"

namespace matcl { namespace test
{
void special_cases::add_special_cases(std::vector<Matrix>& matrices)
{    
    // special case for exactly singular 2x2 cage in >100x100 zero matrix
    // designed for LDL' testing
    {
        Matrix special = zeros(100,100);    
        for(Integer index = 1; index < special.rows() - 2; index += 11)
        {
            special(index,index)     = -1;      special(index,index+1)   = -2;
            special(index+1,index)   = -2;      special(index+1,index+1) = -4; 

            matrices.push_back(special);
        }
    }
    // special case - A is tril sparse with zero non-zero element
    // designed for linsolve (triang) testing
    {
        Matrix temp = spzeros(3,3);
        temp(1,3) = 320.;                                       // will be quitely zeroed out...
        sparse_matrix<Real,false> rep(temp);
        rep.ptr_x()[0] = 0;                                     // ...now
        Matrix newlhs = 200. * speye(3, 3) + temp;              // make it nonsingular
        matrices.push_back(newlhs);
    }
    // special case - A - is triu sparse with zero non-zero element
    // designed for linsolve (triang) testing
    {
        Matrix temp = spzeros(3,3);
        temp(3,1) = 320.;                                       // will be quitely zeroed out...
        sparse_matrix<Real,false> rep(temp);
        rep.ptr_x()[0] = 0;                                     // ...now
        temp(1,3) = 340.;                                       // make it non-lowertriangular
        Matrix newlhs = 200. * speye(3, 3) + temp;              // make it nonsingular  
        matrices.push_back(newlhs);
    }
};

}}