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

#include "matcl-specfunc/lib_functions/elliptic.h"

namespace matcl { namespace details
{

template<class T>
struct MATCL_SF_EXPORT elliptic_helper
{
    using return_type   = typename md::unify_types<T,Float>::type;

    static return_type  eval_ellint_rf(const T& x, const T& y, const T& z);
    static return_type  eval_ellint_rd(const T& x, const T& y, const T& z);
    static return_type  eval_ellint_rj(const T& x, const T& y, const T& z, const T& p);
    static return_type  eval_ellint_rc(const T& x, const T& y);
    static return_type  eval_ellint_rg(const T& x, const T& y, const T& z);
    static return_type  eval_ellint_1(const T& k, const T& phi);
    static return_type  eval_ellint_1(const T& k);
    static return_type  eval_ellint_2(const T& k, const T& phi);
    static return_type  eval_ellint_2(const T& k);
    static return_type  eval_ellint_3(const T& k, const T& n, const T& phi);
    static return_type  eval_ellint_3(const T& k, const T& n);
    static return_type  eval_ellint_d(const T& k, const T& phi);
    static return_type  eval_ellint_d(const T& k);
    static return_type  eval_jacobi_sn(const T& k, const T& u);
    static return_type  eval_jacobi_cn(const T& k, const T& u);
    static return_type  eval_jacobi_dn(const T& k, const T& u);
    static void         eval_jacobi_elliptic(const T& k, const T& u, T& sn, T& cn, T& dn);

};

}};

namespace matcl
{

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ellint_rf(const S1& x, const S2& y, const S3& z)
{
    using S     = typename md::unify_types3_promote<S1,S2,S3,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_rf(S(x), S(y), S(z));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ellint_rd(const S1& x, const S2& y, const S3& z)
{
    using S     = typename md::unify_types3_promote<S1,S2,S3,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_rd(S(x), S(y), S(z));
};

template<class S1, class S2, class S3, class S4, class Enable>
typename md::unify_types4_promote<S1,S2,S3,S4,Float>::type
matcl::ellint_rj(const S1& x, const S2& y, const S3& z, const S4& p)
{
    using S     = typename md::unify_types4_promote<S1,S2,S3,S4,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_rj(S(x), S(y), S(z), S(p));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::ellint_rc(const S1& x, const S2& y)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_rc(S(x), S(y));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ellint_rg(const S1& x, const S2& y, const S3& z)
{
    using S     = typename md::unify_types3_promote<S1,S2,S3,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_rg(S(x), S(y), S(z));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::ellint_1(const S1& k, const S2& phi)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_1(S(k), S(phi));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::ellint_1(const S1& k)
{
    using S     = typename md::unify_types_promote<S1,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_1(S(k));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::ellint_2(const S1& k, const S2& phi)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_2(S(k), S(phi));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::ellint_2(const S1& k)
{
    using S     = typename md::unify_types_promote<S1,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_2(S(k));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ellint_3(const S1& k, const S2& n, const S3& phi)
{
    using S     = typename md::unify_types3_promote<S1,S2,S3,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_3(S(k), S(n), S(phi));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::ellint_3(const S1& k, const S2& n)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_3(S(k), S(n));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::ellint_d(const S1& k, const S2& phi)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_d(S(k), S(phi));
};

template<class S1, class Enable>
typename md::unify_types_promote<S1,Float>::type
matcl::ellint_d(const S1& k)
{
    using S     = typename md::unify_types_promote<S1,Float>::type;
    return md::elliptic_helper<S>::eval_ellint_d(S(k));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::jacobi_sn(const S1& k, const S2& u)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_jacobi_sn(S(k), S(u));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::jacobi_cn(const S1& k, const S2& u)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_jacobi_cn(S(k), S(u));
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::jacobi_dn(const S1& k, const S2& u)
{
    using S     = typename md::unify_types2_promote<S1,S2,Float>::type;
    return md::elliptic_helper<S>::eval_jacobi_dn(S(k), S(u));
};

template<class S1, class S2, class Ret, class Enable>
void matcl::jacobi_elliptic(const S1& k, const S2& u, Ret& sn, Ret& cn, Ret& dn)
{
    using S     = typename md::unify_types3_promote<S1,S2,Ret,Float>::type;
    
    S ret_sn, ret_cn, ret_dn;
    md::elliptic_helper<S>::eval_jacobi_elliptic(S(k), S(u), ret_sn, ret_cn, ret_dn);

    sn = ret_sn;
    cn = ret_cn;
    dn = ret_dn;
};

};
