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

#include "matcl-specfunc/lib_functions/error.h"

namespace matcl { namespace details
{

template<class T>
struct MATCL_SF_EXPORT error_helper
{
    using return_type   = typename md::unify_types<T,Float>::type;

    static return_type  eval_erf(const T& arg);
    static return_type  eval_erfc(const T& arg);
    static return_type  eval_erf_inv(const T& arg);
    static return_type  eval_erfc_inv(const T& arg);
};

}};

namespace matcl
{

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type erf(const S1& n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::error_helper<S>::eval_erf(S(n));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type erfc(const S1& n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::error_helper<S>::eval_erfc(S(n));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type erf_inv(const S1& n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::error_helper<S>::eval_erf_inv(S(n));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type erfc_inv(const S1& n)
{
    using S = typename md::unify_types_promote<S1,Float>::type;
    return details::error_helper<S>::eval_erfc_inv(S(n));
};

};
