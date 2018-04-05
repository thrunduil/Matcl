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

#pragma once

#include "matcl-specfunc/lib_functions/hankel.h"

namespace matcl { namespace details
{

template<class T>
struct MATCL_SF_EXPORT hankel_helper
{
    using return_type   = typename md::unify_types<T,Float_complex>::type;

    static return_type  eval_cyl_hankel_1(const T& v, const T& x);
    static return_type  eval_cyl_hankel_2(const T& v, const T& x);
    static return_type  eval_sph_hankel_1(const T& v, const T& x);
    static return_type  eval_sph_hankel_2(const T& v, const T& x);
};

}};

namespace matcl
{

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
matcl::cyl_hankel_1(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::hankel_helper<S>::eval_cyl_hankel_1(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
matcl::cyl_hankel_2(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::hankel_helper<S>::eval_cyl_hankel_2(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
matcl::sph_hankel_1(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::hankel_helper<S>::eval_sph_hankel_1(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float_complex>::type
matcl::sph_hankel_2(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::hankel_helper<S>::eval_sph_hankel_2(S(v),S(x));
};

};
