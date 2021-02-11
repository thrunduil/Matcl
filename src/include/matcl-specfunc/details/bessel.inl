/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-specfunc/lib_functions/bessel.h"

namespace matcl { namespace details
{

template<class T>
struct MATCL_SF_EXPORT bessel_helper
{
    using return_type   = typename md::unify_types<T,Float>::type;

    static return_type  eval_cyl_bessel_j(const T& v, const T& x);
    static return_type  eval_cyl_neumann(const T& v, const T& x);    
    static return_type  eval_cyl_bessel_j_zero(const T& v, Integer m); 
    static return_type  eval_cyl_neumann_zero(const T& v, Integer m);    
    static return_type  eval_cyl_bessel_j_dif(const T& v, const T& x);
    static return_type  eval_cyl_neumann_dif(const T& v, const T& x);
    
    static return_type  eval_cyl_bessel_i(const T& v, const T& x);
    static return_type  eval_cyl_bessel_k(const T& v, const T& x);
    static return_type  eval_cyl_bessel_i_dif(const T& v, const T& x);
    static return_type  eval_cyl_bessel_k_dif(const T& v, const T& x);

    static return_type  eval_sph_bessel(Integer n, const T& x);
    static return_type  eval_sph_neumann(Integer n, const T& x);
    static return_type  eval_sph_bessel_dif(Integer n, const T& x);
    static return_type  eval_sph_neumann_dif(Integer n, const T& x);
};

}};

namespace matcl
{

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_bessel_j(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_bessel_j(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_neumann(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_neumann(S(v),S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::cyl_bessel_j_zero(const S1& v, Integer m)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using S     = typename md::unify_types<SP1,Float>::type;
    return details::bessel_helper<S>::eval_cyl_bessel_j_zero(S(v),m);
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::cyl_neumann_zero(const S1& v, Integer m)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using S     = typename md::unify_types<SP1,Float>::type;
    return details::bessel_helper<S>::eval_cyl_neumann_zero(S(v),m);
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_bessel_j_dif(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_bessel_j_dif(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_neumann_dif(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_neumann_dif(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_bessel_i(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_bessel_i(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_bessel_k(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_bessel_k(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_bessel_i_dif(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_bessel_i_dif(S(v),S(x));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::cyl_bessel_k_dif(const S1& v, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::bessel_helper<S>::eval_cyl_bessel_k_dif(S(v),S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::sph_bessel(Integer n, const S1& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using S     = typename md::unify_types<SP1,Float>::type;
    return details::bessel_helper<S>::eval_sph_bessel(n , S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::sph_neumann(Integer n, const S1& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using S     = typename md::unify_types<SP1,Float>::type;
    return details::bessel_helper<S>::eval_sph_neumann(n , S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::sph_bessel_dif(Integer n, const S1& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using S     = typename md::unify_types<SP1,Float>::type;
    return details::bessel_helper<S>::eval_sph_bessel_dif(n , S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::sph_neumann_dif(Integer n, const S1& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using S     = typename md::unify_types<SP1,Float>::type;
    return details::bessel_helper<S>::eval_sph_neumann_dif(n , S(x));
};

};
