/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017
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

#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/details/matfunc_helpers.h"
#include "matcl-scalar/details/func_unary.inl"

namespace matcl { namespace details
{

template<class S>
struct MATCL_SCALAR_EXPORT basic_binfunc_helper
{
    using SF = typename unify_types<S,Float>::type;

    static SF   eval_powm1(const S& x, const S& y);    
};

template<class S>
struct copysign_helper
{
    static S eval(const S& x, const S& y)
    {
        return mrd::scal_func::copysign(x,y);
    }
};
template<>
struct copysign_helper<Object>
{
    static Object eval(const Object& x, const Object& y)
    {
        return dynamic::copysign(x,y);
    }
};

template<class S>
struct nextafter_helper
{
    static S eval(const S& x, const S& y)
    {
        return std::nextafter(x,y);
    }
};
template<>
struct nextafter_helper<Object>
{
    static Object eval(const Object& x, const Object& y)
    {
        return dynamic::nextafter(x,y);
    }
};

template<class S>
struct float_distance_helper
{
    static S eval(const S& x, const S& y)
    {
        return mrd::scal_func::float_distance(x,y);
    }
};
template<>
struct float_distance_helper<Object>
{
    static Object eval(const Object& x, const Object& y)
    {
        return dynamic::float_distance(x,y);
    }
};

template<class S>
struct fdim_helper
{
    static S eval(const S& x, const S& y)
    {
        return std::fdim(x,y);
    }
};
template<>
struct fdim_helper<Object>
{
    static Object eval(const Object& x, const Object& y)
    {
        return dynamic::fdim(x,y);
    }
};

};};

namespace matcl
{

namespace mrd = matcl::raw::details;

template<class S1, class S2, class Enable>
bool matcl::op_and(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::op_and_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
bool matcl::op_or(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::op_or_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
bool matcl::op_xor(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::op_xor_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::elem_and(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::elem_and_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::elem_or(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::elem_or_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::elem_xor(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::elem_xor_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::eeq(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::eeq_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::eeq_nan(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::eeq_nan_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::neq(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::neq_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::neq_nan(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::neq_nan_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::lt(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::lt_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::leq(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::leq_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::gt(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::gt_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::bool_or_object2<S1,S2>::type
matcl::geq(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::geq_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::mul(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::elem_mul_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::plus(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return matcl::raw::details::plus_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::minus(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::minus_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::max(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::max_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::min(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::min_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,S2>::type
matcl::mod(const S1& A, const S2& B)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::mod_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,S2>::type
matcl::rem(const S1& A, const S2& B)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::rem_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::real_unify_types2_promote<S1,S2,Float>::type
matcl::atan2(const S1& A,const S2& B)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::atan2_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::real_unify_types2_promote<S1,S2,Float>::type
matcl::hypot(const S1& A,const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::hypot_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::div(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::div_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::div_0(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::div_0_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::div_1(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::div_1_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::unify_types_promote<S1,S2>::type
matcl::idiv(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::idiv_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename mrd::pow_return<S1,S2>::type
matcl::pow(const S1& A, const S2& B)
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    
    return mrd::pow_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename mrd::pow_c_return<S1,S2>::type
matcl::pow_c(const S1& A, const S2& B) 
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    return mrd::pow_c_helper<SP1,SP2>::eval(A,B);
};

template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
matcl::copysign(const S1& x, const S2& y)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using S = typename md::real_unify_types_promote<S1,Float>::type;
    return details::copysign_helper<S>::eval(S(real(x)), S(real(y)));
}

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2, Float>::type
matcl::fdim(const S1& x, const S2& y)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using S = typename md::unify_types2_promote<S1,S2,Float>::type;
    return details::fdim_helper<S>::eval(S(x), S(y));
}

template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
matcl::nextafter(const S1& x, const S2& y)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using S = typename md::real_unify_types_promote<S1,Float>::type;
    return details::nextafter_helper<S>::eval(S(real(x)), S(real(y)));
}

template<class S1, class S2, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
matcl::float_distance(const S1& x, const S2& y)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using S = typename md::real_unify_types_promote<S1,Float>::type;
    return details::float_distance_helper<S>::eval(S(real(x)), S(real(y)));
}

template<class S1, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
matcl::nextabove(const S1& x)
{
    static_assert(!md::is_complex<S1>::value, "not available for complex");

    using S = typename md::real_unify_types_promote<S1,Float>::type;
    return mrd::nextabove_helper<S>::eval(S(real(x)));
}

template<class S1, class Enable>
typename md::real_unify_types_promote<S1,Float>::type
matcl::nextbelow(const S1& x)
{
    static_assert(!md::is_complex<S1>::value, "not available for complex");

    using S = typename md::real_unify_types_promote<S1,Float>::type;
    return mrd::nextbelow_helper<S>::eval(S(real(x)));
}

template<class S1, class S2, class Enable>
typename md::unify_types2_promote<S1,S2,Float>::type
matcl::powm1(const S1& x, const S2& y)
{
    static_assert(!md::is_complex<S1>::value && !md::is_complex<S2>::value, 
                  "not available for complex");

    using S = typename md::unify_types2_promote<S1,S2,Float>::type;
    return details::basic_binfunc_helper<S>::eval_powm1(S(x), S(y));
};

inline Real matcl::fma_f(Real a, Real b, Real c)
{
    return mrd::fma_f_helper<Real>::eval(a,b,c);
};

inline Float matcl::fma_f(Float a, Float b, Float c)
{
    return mrd::fma_f_helper<Float>::eval(a,b,c);
};

inline Object matcl::fma_f(const Object& a, const Object& b, const Object& c)
{
    return mrd::fma_f_helper<Object>::eval(a,b,c);
};

inline Real matcl::fms_f(Real a, Real b, Real c)
{
    return mrd::fms_f_helper<Real>::eval(a,b,c);
};

inline Float matcl::fms_f(Float a, Float b, Float c)
{
    return mrd::fms_f_helper<Float>::eval(a,b,c);
};

inline Object matcl::fms_f(const Object& a, const Object& b, const Object& c)
{
    return mrd::fms_f_helper<Object>::eval(a,b,c);
};

//
inline Real matcl::fma_a(Real a, Real b, Real c)
{
    return mrd::fma_a_helper<Real>::eval(a,b,c);
};

inline Float matcl::fma_a(Float a, Float b, Float c)
{
    return mrd::fma_a_helper<Float>::eval(a,b,c);
};

inline Object matcl::fma_a(const Object& a, const Object& b, const Object& c)
{
    return mrd::fma_a_helper<Object>::eval(a,b,c);
};

inline Real matcl::fms_a(Real a, Real b, Real c)
{
    return mrd::fms_a_helper<Real>::eval(a,b,c);
};

inline Float matcl::fms_a(Float a, Float b, Float c)
{
    return mrd::fms_a_helper<Float>::eval(a,b,c);
};

inline Object matcl::fms_a(const Object& a, const Object& b, const Object& c)
{
    return mrd::fms_a_helper<Object>::eval(a,b,c);
};

inline Real matcl::dot2_a(Real a, Real b, Real c, Real d)
{
    return mrd::dot2_a_helper<Real>::eval(a,b,c,d);
}

inline Float matcl::dot2_a(Float a, Float b, Float c, Float d)
{
    return mrd::dot2_a_helper<Float>::eval(a,b,c,d);
}

inline Object matcl::dot2_a(const Object& a, const Object& b, const Object& c, 
                                    const Object& d)
{
    return mrd::dot2_a_helper<Object>::eval(a,b,c,d);
}

}
