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

#include "matcl-specfunc/lib_functions/airy.h"

namespace matcl { namespace details
{

template<class T>
struct MATCL_SF_EXPORT airy_helper
{
    using return_type   = typename md::unify_types<T,Float>::type;

    static return_type  eval_airy_ai(const T& arg);
    static return_type  eval_airy_bi(const T& arg);
    static return_type  eval_airy_ai_dif(const T& arg);
    static return_type  eval_airy_bi_dif(const T& arg);
    static return_type  eval_airy_ai_zero(Integer arg);
    static return_type  eval_airy_bi_zero(Integer arg);
};

}};

namespace matcl
{

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::airy_ai(const S1& x)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::airy_helper<S>::eval_airy_ai(S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::airy_bi(const S1& x)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::airy_helper<S>::eval_airy_bi(S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::airy_ai_dif(const S1& x)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::airy_helper<S>::eval_airy_ai_dif(S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::airy_bi_dif(const S1& x)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::airy_helper<S>::eval_airy_bi_dif(S(x));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::airy_ai_zero(Integer m)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::airy_helper<S>::eval_airy_ai_zero(m);
};

inline Real airy_ai_zero(Integer m)
{
    return details::airy_helper<Real>::eval_airy_ai_zero(m);
};

inline Float fairy_ai_zero(Integer m)
{
    return details::airy_helper<Float>::eval_airy_ai_zero(m);
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::airy_bi_zero(Integer m)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::airy_helper<S>::eval_airy_bi_zero(m);
};

inline Real airy_bi_zero(Integer m)
{
    return details::airy_helper<Real>::eval_airy_bi_zero(m);
};

inline Float fairy_bi_zero(Integer m)
{
    return details::airy_helper<Float>::eval_airy_bi_zero(m);
};

};
