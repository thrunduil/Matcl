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
#pragma once

#include "matcl-scalar/lib_functions/func_matrix.h"
#include "matcl-scalar/lib_functions/func_binary.h"

namespace matcl
{

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::kron(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::mmul_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::mmul(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;

    return mrd::mmul_helper<SP1,SP2>::eval(A, B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
operator*(const S1& x, const S2& y)             
{ 
    return mmul<S1, S2, Enable>(x,y); 
};

};