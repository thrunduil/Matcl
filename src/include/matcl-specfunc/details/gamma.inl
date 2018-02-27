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

#include "matcl-specfunc/lib_functions/gamma.h"

namespace matcl { namespace details
{

template<class T>
struct MATCL_SF_EXPORT gamma_helper
{
    using return_type   = typename md::unify_types<T,Float>::type;

    static return_type  eval_gamma(const T& arg);    
    static return_type  eval_gamma1pm1(const T& arg);

    static return_type  eval_gammaln(const T& arg);
    static return_type  eval_gammaln(const T& arg, int& sign);

    static return_type  eval_digamma(const T& arg);
    static return_type  eval_trigamma(const T& arg);
    static return_type  eval_polygamma(const T& arg, Integer n);

    static return_type  eval_gamma_ratio(const T& x, const T& y);
    static return_type  eval_gamma_delta_ratio(const T& x, const T& delta);

    static return_type  eval_igamma_lower(const T& a, const T& z);
    static return_type  eval_igamma_upper(const T& a, const T& z);
    static return_type  eval_igamma_lower_norm(const T& a, const T& z);
    static return_type  eval_igamma_upper_norm(const T& a, const T& z);

    static return_type  eval_igamma_lower_inv(const T& a, const T& q);
    static return_type  eval_igamma_upper_inv(const T& a, const T& q);
    static return_type  eval_igamma_lower_inva(const T& x, const T& q);
    static return_type  eval_igamma_upper_inva(const T& x, const T& q);

    static return_type  eval_igamma_lower_dif(const T& a, const T& x);
};

}};

namespace matcl
{

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::gamma(const S1& n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::gamma_helper<S>::eval_gamma(S(n));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::gammaln(const S1& n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::gamma_helper<S>::eval_gammaln(S(n));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::gamma1pm1(const S1& n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::gamma_helper<S>::eval_gamma1pm1(S(n));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::gammaln(const S1& x, int& sign)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::gamma_helper<S>::eval_gammaln(S(x),sign);
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::digamma(const S1& x)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::gamma_helper<S>::eval_digamma(S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::trigamma(const S1& x)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::gamma_helper<S>::eval_trigamma(S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::polygamma(const S1& x, Integer n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::gamma_helper<S>::eval_polygamma(S(x),n);
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::gamma_ratio(const S1& x, const S2& y)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_gamma_ratio(S(x),S(y));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::gamma_delta_ratio(const S1& x, const S2& y)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_gamma_delta_ratio(S(x),S(y));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_lower(const S1& a, const S2& z)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_lower(S(a),S(z));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_upper(const S1& a, const S2& z)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_upper(S(a),S(z));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_lower_norm(const S1& a, const S2& z)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_lower_norm(S(a),S(z));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_upper_norm(const S1& a, const S2& z)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_upper_norm(S(a),S(z));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_lower_inv(const S1& a, const S2& q)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_lower_inv(S(a),S(q));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_upper_inv(const S1& a, const S2& q)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_upper_inv(S(a),S(q));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_lower_inva(const S1& x, const S2& q)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_lower_inva(S(x),S(q));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_upper_inva(const S1& x, const S2& q)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_upper_inva(S(x),S(q));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::igamma_lower_dif(const S1& a, const S2& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types2<SP1,SP2,Float>::type;
    return details::gamma_helper<S>::eval_igamma_lower_dif(S(a),S(x));
};

};
