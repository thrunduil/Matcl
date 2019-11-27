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
#pragma once

#include "matcl-linalg/decompositions/givens.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/lib_functions/func_binary.h"

namespace matcl { namespace details
{

template<class Val, bool Is_real>
struct construct_givens_impl
{
    static void eval(const Val& a, const Val& b, Val& c, Val& s, Val& r)
    {
        using VR    = typename real_type<Val>::type;

        Val roe     = b;
        VR abs_a    = abs(a);
        VR abs_b    = abs(b);

        if (abs_a > abs_b)
	        roe     = a;
        
        VR scale    = abs_a + abs_b;
        
        if (scale == VR(0)) 
        {
            c       = Val(1);
            s       = Val(0);
            r       = Val(0);
        }
        else
        {
            Val as  = a / scale;
            Val bs  = b / scale;
            r       = scale * sqrt(as * as + bs * bs);
            r       = copysign(r, roe);
            c       = a / r;
            s       = b / r;
        };

        return;
    };
};

template<class Val>
struct construct_givens_impl<Val,false>
{
    using VR    = typename real_type<Val>::type;

    static void eval(const Val& a, const Val& b, VR& c, Val& s, Val& r)
    {        
        Val roe     = b;
        VR abs_a    = matcl::abs(a);
        VR abs_b    = matcl::abs(b);

        if (abs_a == VR())
        {
            c       = VR(0);
            s       = Val(1.0);
	        r       = b;
        }
        else
        {
            VR scale    = abs_a + abs_b;
            Val as      = a / scale;
            Val bs      = b / scale;
            VR norm     = scale * sqrt(abs2(as) + abs2(bs));
            Val alpha   = a/abs_a;
            c           = abs_a / norm;
            s           = alpha * conj(b) / norm;
            r           = alpha * norm;
        }
        
        return;
    };
};

}};

namespace matcl
{

template<class Val, class Enable>
void matcl::construct_givens(const Val& a, const Val& b, typename md::real_type<Val>::type& c, Val& s)
{
    static const bool is_real = details::is_float_real_scalar<Val>::value;

    Val r;
    details::construct_givens_impl<Val,is_real>::eval(a, b, c, s, r);
};

template<class Val, class Enable>
void matcl::construct_givens(const Val& a, const Val& b, typename md::real_type<Val>::type& c, Val& s, Val& r)
{
    static const bool is_real = details::is_float_real_scalar<Val>::value;
    details::construct_givens_impl<Val,is_real>::eval(a, b, c, s, r);
};

};
