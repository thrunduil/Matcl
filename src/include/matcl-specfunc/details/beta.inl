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

#include "matcl-specfunc/lib_functions/beta.h"

namespace matcl { namespace details
{

template<class T>
struct MATCL_SF_EXPORT beta_helper
{
    using return_type   = typename md::unify_types<T, Float>::type;

    static return_type  eval_beta(const T& A, const T& B);
    static return_type  eval_ibeta(const T& a, const T& b, const T& x);
    static return_type  eval_ibetac(const T& a, const T& b, const T& x);
    static return_type  eval_ibeta_norm(const T& a, const T& b, const T& x);
    static return_type  eval_ibetac_norm(const T& a, const T& b, const T& x);
    static return_type  eval_ibeta_inv(const T& a, const T& b, const T& p);
    static return_type  eval_ibetac_inv(const T& a, const T& b, const T& q);
    static return_type  eval_ibeta_inva(const T& b, const T& x, const T& p);
    static return_type  eval_ibetac_inva(const T& b, const T& x, const T& q);
    static return_type  eval_ibeta_invb(const T& a, const T& x, const T& p);
    static return_type  eval_ibetac_invb(const T& a, const T& x, const T& q);
    static return_type  eval_ibeta_dif(const T& a, const T& b, const T& x);
};

}};

namespace matcl
{

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
beta(const S1& A,const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using S     = typename md::unify_types<SP1,SP2>::type;
    return details::beta_helper<S>::eval_beta(S(A),S(B));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibeta(const S1& a, const S2& b, const S3& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibeta(S(a),S(b),S(x));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibetac(const S1& a, const S2& b, const S3& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibetac(S(a),S(b),S(x));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibeta_norm(const S1& a, const S2& b, const S3& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibeta_norm(S(a),S(b),S(x));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibetac_norm(const S1& a, const S2& b, const S3& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibetac_norm(S(a),S(b),S(x));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibeta_inv(const S1& a, const S2& b, const S3& p)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibeta_inv(S(a),S(b),S(p));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibetac_inv(const S1& a, const S2& b, const S3& q)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibetac_inv(S(a),S(b),S(q));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibeta_inva(const S1& b, const S2& x, const S3& p)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibeta_inva(S(b),S(x),S(p));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibetac_inva(const S1& b, const S2& x, const S3& q)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibetac_inva(S(b),S(x),S(q));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibeta_invb(const S1& a, const S2& x, const S3& p)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibeta_invb(S(a),S(x),S(p));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibetac_invb(const S1& a, const S2& x, const S3& q)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibetac_invb(S(a),S(x),S(q));
};

template<class S1, class S2, class S3, class Enable>
typename md::unify_types3_promote<S1,S2,S3,Float>::type
matcl::ibeta_dif(const S1& a, const S2& b, const S3& x)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using SP3   = typename md::promote_scalar<S3>::type;
    using S     = typename md::unify_types2<SP1,SP2,SP3>::type;
    return details::beta_helper<S>::eval_ibeta_dif(S(a),S(b),S(x));
};

};
