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

#include "matcl-scalar/details/scalfunc_helpers.h"

#pragma warning (push)
#pragma warning (disable:4127) //conditional expression is constant

namespace matcl { namespace raw { namespace details
{

inline Integer gl_minus(Integer a)                          { return -a; };
inline Real    gl_minus(const Real& a)                      { return -a; };
inline Float   gl_minus(const Float& a)                     { return -a; };
inline Object  gl_minus(const Object& a)                    { return dynamic::operator-(a); };

inline Integer gl_mult(Integer a, Integer b)                { return a*b; };
inline Real    gl_mult(const Real& a, const Real& b)        { return a*b; };
inline Real    gl_mult(const Real& a, const Float& b)       { return a*b; };
inline Real    gl_mult(const Float& a, const Real& b)       { return a*b; };
inline Float   gl_mult(const Float& a, const Float& b)      { return a*b; };
inline Object  gl_mult(const Object& a, const Object& b)    { return dynamic::operator*(a,b); };
inline Object  gl_emul(const Object& a, const Object& b)    { return dynamic::elem_mul(a,b); };

inline Real    gl_div(const Real& a, const Real& b)         { return a/b; };
inline Float   gl_div(const Float& a, const Float& b)       { return a/b; };
inline Real    gl_div(const Real& a, const Float& b)        { return a/b; };
inline Real    gl_div(const Float& a, const Real& b)        { return a/b; };
inline Object  gl_div(const Object& a, const Object& b)     { return dynamic::operator/(a,b); };

inline Real    gl_div_0(const Real& a, const Real& b)       { return (a == 0 && b == 0) ? 0 : a/b; };
inline Float   gl_div_0(const Float& a, const Float& b)     { return (a == 0 && b == 0) ? 0 : a/b; };
inline Real    gl_div_0(const Real& a, const Float& b)      { return (a == 0 && b == 0) ? 0 : a/b; };
inline Real    gl_div_0(const Float& a, const Real& b)      { return (a == 0 && b == 0) ? 0 : a/b; };
inline Object  gl_div_0(const Object& a, const Object& b)   { return dynamic::div_0(a,b); };

inline Real    gl_div_1(const Real& a, const Real& b)       { return (a == 0 && b == 0) ? 1 : a/b; };
inline Float   gl_div_1(const Float& a, const Float& b)     { return (a == 0 && b == 0) ? 1 : a/b; };
inline Real    gl_div_1(const Real& a, const Float& b)      { return (a == 0 && b == 0) ? 1 : a/b; };
inline Real    gl_div_1(const Float& a, const Real& b)      { return (a == 0 && b == 0) ? 1 : a/b; };
inline Object  gl_div_1(const Object& a, const Object& b)   { return dynamic::div_1(a,b); };

inline Integer gl_minus(Integer a, Integer b)               { return a-b; };
inline Real    gl_minus(const Real& a, const Real& b)       { return a-b; };
inline Real    gl_minus(const Real& a, const Float& b)      { return a-b; };
inline Real    gl_minus(const Float& a, const Real& b)      { return a-b; };
inline Float   gl_minus(const Float& a, const Float& b)     { return a-b; };
inline Object  gl_minus(const Object& a, const Object& b)   { return dynamic::operator-(a,b); };

inline Integer gl_plus(Integer a, Integer b)                { return a+b; };
inline Real    gl_plus(const Real& a, const Real& b)        { return a+b; };
inline Real    gl_plus(const Real& a, const Float& b)       { return a+b; };
inline Real    gl_plus(const Float& a, const Real& b)       { return a+b; };
inline Float   gl_plus(const Float& a, const Float& b)      { return a+b; };
inline Object  gl_plus(const Object& a, const Object& b)    { return dynamic::operator+(a,b); };

inline Integer gl_idiv(Integer a, Integer b)                { return a/b; };
inline Real    gl_idiv(const Real& a, const Real& b)        { return a/b; };
inline Real    gl_idiv(const Float& a, const Real& b)       { return a/b; };
inline Real    gl_idiv(const Real& a, const Float& b)       { return a/b; };
inline Float   gl_idiv(const Float& a, const Float& b)      { return a/b; };
inline Object  gl_idiv(const Object& a, const Object& b)    { return dynamic::idiv(a,b); };

template<class T> struct promote_int						{ using type = T;};
template<>		  struct promote_int<Integer>				{ using type = Real;};

template<class T1, class T2>
struct test_range_pow
{
    static bool eval(T1 A,T2 B)
    {
        if (A >= 0 || mrd::isint_helper<T2>::eval(B))
            return true;

        return false;
    };
};
template<class T1>
struct test_range_pow<T1,Integer>
{
    static bool eval(T1 ,Integer )
    {
        return true;
    };
};

template<class ret, bool isc, bool iso, class S1, class S2>
struct plus_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_plus(arg1,arg2);
    };
};
template<class ret,class S1,class S2>
struct plus_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::plus_c(arg1,arg2);
    };
};
template<class ret,bool isc, class S1,class S2>
struct plus_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_plus(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso,class S1,class S2>
struct minus_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_minus(arg1,arg2);
    };
};
template<class ret,class S1,class S2>
struct minus_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::minus_c(arg1,arg2);
    };
};
template<class ret,bool isc, class S1,class S2>
struct minus_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_minus(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso, class S1,class S2>
struct div_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_div(arg1,arg2);
    };
};
template<class ret,bool isc,bool iso, class S1,class S2>
struct div_0_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_div_0(arg1,arg2);
    };
};
template<class ret,bool isc,bool iso, class S1,class S2>
struct div_1_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_div_1(arg1,arg2);
    };
};

template<class ret,class S1,class S2>
struct div_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::div_c(arg1,arg2);
    };
};
template<class ret,class S1,class S2>
struct div_0_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::div_0_c(arg1,arg2);
    };
};
template<class ret,class S1,class S2>
struct div_1_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::div_1_c(arg1,arg2);
    };
};

template<class ret,bool isc, class S1,class S2>
struct div_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_div(Object(arg1),Object(arg2));
    };
};
template<class ret,bool isc, class S1,class S2>
struct div_0_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_div_0(Object(arg1),Object(arg2));
    };
};
template<class ret,bool isc, class S1,class S2>
struct div_1_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_div_1(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso,class S1,class S2>
struct elem_mul_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_mult(arg1,arg2);
    };
};
template<class ret,class S1,class S2>
struct elem_mul_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::mul_c(arg1,arg2);
    };
};
template<class ret,bool isc,class S1,class S2>
struct elem_mul_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_emul(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso,class S1,class S2>
struct mmul_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_mult(arg1,arg2);
    };
};
template<class ret,class S1,class S2>
struct mmul_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::mul_c(arg1,arg2);
    };
};
template<class ret,bool isc,class S1,class S2>
struct mmul_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_mult(Object(arg1),Object(arg2));
    };
};

template<class ret, bool isc, bool iso, class S1,class S2>
struct min_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gt_helper<S1,S2>::eval(arg1,arg2) ? ret(arg2) : ret(arg1); 
    };
};
template<class ret,class S1,class S2>
struct min_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::gt_c(arg1,arg2)? ret(arg2) : ret(arg1); 
    };
};
template<class ret,class S1,class S2>
struct min_helper_impl<ret,false,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::min(Object(arg1), Object(arg2));
    };
};

template<class ret,bool isc,bool iso,class S1,class S2>
struct max_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 > arg2)? ret(arg1) : ret(arg2); 
    };
};
template<class ret,class S1,class S2>
struct max_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::gt_c(arg1,arg2)? ret(arg1) : ret(arg2); 
    };
};
template<class ret,bool isc,class S1,class S2>
struct max_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::max(Object(arg1), Object(arg2));
    };
};

template<class ret,bool isc,bool iso,class S1,class S2>
struct eeq_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 == arg2); 
    };
};
template<class ret,class S1,class S2>
struct eeq_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::eeq_c(arg1,arg2);
    };
};
template<class ret,bool isc,class S1,class S2>
struct eeq_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::operator==(Object(arg1),Object(arg2));
    };
};

template<class ret, bool isc, bool isf, bool iso,class S1,class S2>
struct eeq_nan_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 == arg2); 
    };
};
template<class ret, class S1,class S2>
struct eeq_nan_helper_impl<ret, false,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        if (mrd::scal_func::isnan(arg1) == true)
        {
            if (mrd::scal_func::isnan(arg2) == true)
                return true;
            else
                return false;
        }
        else
        {
            return (arg1 == arg2); 
        }
    };
};
template<class ret,bool isf, class S1,class S2>
struct eeq_nan_helper_impl<ret,true,isf,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::eeq_nan_c(arg1,arg2);
    };
};
template<class ret,bool isc,bool isf, class S1,class S2>
struct eeq_nan_helper_impl<ret,isc,isf,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::eeq_nan(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso, class S1,class S2>
struct neq_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 != arg2); 
    };
};
template<class ret,class S1,class S2>
struct neq_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::neq_c(arg1,arg2);
    };
};
template<class ret,bool isc, class S1,class S2>
struct neq_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::operator!=(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool isf, bool iso, class S1,class S2>
struct neq_nan_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 != arg2); 
    };
};
template<class ret, class S1, class S2>
struct neq_nan_helper_impl<ret, false, true, false, S1, S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        if (mrd::scal_func::isnan(arg1) == true)
        {
            if (mrd::scal_func::isnan(arg2) == true)
                return false;
            else
                return true;
        }
        else
        {
            return (arg1 != arg2); 
        }
    };
};
template<class ret, bool isf, class S1,class S2>
struct neq_nan_helper_impl<ret,true,isf,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::neq_nan_c(arg1,arg2);
    };
};
template<class ret,bool isc, bool isf,class S1,class S2>
struct neq_nan_helper_impl<ret,isc,isf,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::neq_nan(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso,class S1,class S2>
struct leq_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 <= arg2); 
    };
};
template<class ret,class S1,class S2>
struct leq_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::leq_c(arg1,arg2);
    };
};
template<class ret,bool isc,class S1,class S2>
struct leq_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::operator<=(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso,class S1,class S2>
struct geq_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 >= arg2); 
    };
};
template<class ret,class S1,class S2>
struct geq_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::geq_c(arg1,arg2);
    };
};
template<class ret,bool isc,class S1,class S2>
struct geq_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::operator>=(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso, class S1,class S2>
struct lt_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 < arg2); 
    };
};
template<class ret,class S1,class S2>
struct lt_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::lt_c(arg1,arg2);
    };
};
template<class ret,bool isc,class S1,class S2>
struct lt_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::operator<(Object(arg1),Object(arg2));
    };
};

template<class ret,bool isc,bool iso, class S1,class S2>
struct gt_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return (arg1 > arg2); 
    };
};
template<class ret,class S1,class S2>
struct gt_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::gt_c(arg1,arg2);
    };
};
template<class ret,bool isc, class S1,class S2>
struct gt_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::operator>(Object(arg1),Object(arg2));
    };
};

template<class T1, class T2>
struct elem_mul_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static const bool isc = md::is_complex<S1>::value 
                            ||md::is_complex<S2>::value;
    static const bool iso = std::is_same<S1,Object>::value 
                            ||std::is_same<S2,Object>::value;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        return elem_mul_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};
template<>
struct elem_mul_helper<Integer,Integer>
{
    using T1        = Integer;
    using T2        = Integer;
    using ret       = Integer ;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return gl_mult(arg1,arg2);
    };
};

template<class T1, class T2>
struct mmul_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static const bool isc = md::is_complex<S1>::value 
                            ||md::is_complex<S2>::value;
    static const bool iso = std::is_same<S1,Object>::value 
                            ||std::is_same<S2,Object>::value;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        return mmul_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};
template<>
struct mmul_helper<Integer,Integer>
{
    using T1        = Integer;
    using T2        = Integer;
    using ret       = Integer ;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return gl_mult(arg1,arg2);
    };
};

template<class T1, class T2>
struct div_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;
        return div_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<>
struct div_helper<Integer,Integer>
{
    using T1    = Real;
    using T2    = Real;
    using ret   = Real;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return gl_div(arg1,arg2); 
    };
};

template<class T1, class T2>
struct div_0_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;
        return div_0_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<>
struct div_0_helper<Integer,Integer>
{
    using T1    = Real;
    using T2    = Real;
    using ret   = Real;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return gl_div_0(arg1,arg2); 
    };
};

template<class T1, class T2>
struct div_1_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;
        return div_1_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<>
struct div_1_helper<Integer,Integer>
{
    using T1    = Real;
    using T2    = Real;
    using ret   = Real;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return gl_div_1(arg1,arg2); 
    };
};

template<class T1, class T2>
struct minus_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
    static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        return minus_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<>
struct minus_helper<Integer,Integer>
{
    using T1    = Integer;
    using T2    = Integer;
    using ret   = Integer;

    static ret eval(const T1& arg1, const T2& arg2)
    {	
        return gl_minus(arg1,arg2);
    };
};

template<class T1, class T2>
struct plus_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
    static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        return plus_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<>
struct plus_helper<Integer,Integer>
{
    using T1    = Integer;
    using T2    = Integer;
    using ret   = Integer;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return gl_plus(arg1,arg2);
    };
};

template<class T1, class T2>
struct max_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;
        return max_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<>
struct max_helper<Integer,Integer>
{
    using T1    = Integer;
    using T2    = Integer;
    using ret   = Integer;

    static ret eval( T1 arg1, T2 arg2)						
    {	
        return (arg1 > arg2) ? arg1 : arg2; 
    };
};

template<class T1, class T2>
struct min_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value || md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value || std::is_same<S2,Object>::value;
        return min_helper_impl<ret, isc, iso, S1, S2>::eval(arg1,arg2);	
    };
};

template<>
struct min_helper<Integer,Integer>
{
    using T1    = Integer;
    using T2    = Integer;
    using ret   = Integer;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return (arg1 < arg2)? arg1 : arg2; 
    };
};

template<class S1, class S2>
struct lt_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;
        return lt_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class S1, class S2>
struct gt_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;
        return gt_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class S1, class S2>
struct eeq_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value || md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value || std::is_same<S2,Object>::value;
        return eeq_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class S1, class S2>
struct eeq_nan_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value || md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value || std::is_same<S2,Object>::value;
        static const bool isf = md::is_float_real_scalar<S1>::value 
                                || md::is_float_real_scalar<S1>::value;
        return eeq_nan_helper_impl<ret,isc,isf,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class S1, class S2>
struct neq_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;

        return neq_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class S1, class S2>
struct neq_nan_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value ||std::is_same<S2,Object>::value;
        static const bool isf = md::is_float_real_scalar<S1>::value 
                                && md::is_float_real_scalar<S1>::value;

        return neq_nan_helper_impl<ret,isc,isf,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class S1, class S2>
struct leq_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value || md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value || std::is_same<S2,Object>::value;
        return leq_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class S1, class S2>
struct geq_helper
{
    using ret = typename md::bool_or_object2<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value || md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value || std::is_same<S2,Object>::value;
        return geq_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};

template<class T1, class T2, bool is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct op_or_helper
{
    using ret = typename md::bool_or_object2<T1,T2>::type;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return cast_bool_helper<T1>::eval(arg1) || cast_bool_helper<T2>::eval(arg2); 
    };
};
template<class T1, class T2>
struct op_or_helper<T1, T2, true>
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return dynamic::op_or(Object(arg1), Object(arg2)); 
    };
};

template<class T1, class T2, bool is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct elem_or_helper
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return cast_bool_helper<T1>::eval(arg1) || cast_bool_helper<T2>::eval(arg2); 
    };
};
template<class T1, class T2>
struct elem_or_helper<T1,T2,true>
{
    using ret = Object;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return dynamic::elem_or(Object(arg1), Object(arg2));
    };
};

template<class T1, class T2, bool is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct op_and_helper
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return cast_bool_helper<T1>::eval(arg1) && cast_bool_helper<T2>::eval(arg2); 
    };
};
template<class T1, class T2>
struct op_and_helper<T1, T2, true>
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return dynamic::op_and(Object(arg1), Object(arg2)); 
    };
};
template<class T1, class T2, bool is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct elem_and_helper
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return cast_bool_helper<T1>::eval(arg1) && cast_bool_helper<T2>::eval(arg2); 
    };
};
template<class T1, class T2>
struct elem_and_helper<T1,T2,true>
{
    using ret = Object;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return dynamic::elem_and(Object(arg1), Object(arg2));
    };
};

template<class T1, class T2, bool is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct op_xor_helper
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return cast_bool_helper<T1>::eval(arg1) ^ cast_bool_helper<T2>::eval(arg2); 
    };
};
template<class T1, class T2>
struct op_xor_helper<T1, T2, true>
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return dynamic::op_xor(Object(arg1), Object(arg2)); 
    };
};
template<class T1, class T2, bool is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct elem_xor_helper
{
    using ret = bool;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return cast_bool_helper<T1>::eval(arg1) ^ cast_bool_helper<T2>::eval(arg2); 
    };
};
template<class T1, class T2>
struct elem_xor_helper<T1,T2,true>
{
    using ret = Object;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        return dynamic::elem_xor(Object(arg1), Object(arg2));
    };
};

template<class T1, class T2,bool iso>
struct atan2_helper_impl
{
    using ret0  = typename md::unify_types2<T1, T2, Float>::type;
    using ret   = typename md::real_type<ret0>::type;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return scal_func::atan2((ret)arg1,(ret)arg2);
    };
};

template<class T1, class T2>
struct atan2_helper_impl<complex<T1>,T2,false>
{
    using ret0  = typename md::unify_types<T1, T2>::type;
    using ret   = typename md::real_type<ret0>::type;

    static ret eval(const complex<T1>& arg1, const T2& arg2)	
    {	
        if (imag(arg1) == 0)
            return atan2_helper<T1,T2>::eval(real(arg1),arg2); 
        else
            throw error::function_not_defined_for_complex("atan2");
    };
};

template<class T1, class T2>
struct atan2_helper_impl<T1,complex<T2>,false>
{
    using ret0  = typename md::unify_types<T1, T2>::type;
    using ret   = typename md::real_type<ret0>::type;

    static ret eval(const T1& arg1,const complex<T2>& arg2)	
    {	
        if (imag(arg2) == 0)
            return atan2_helper<T1,T2>::eval(arg1,real(arg2)); 
        else
            throw error::function_not_defined_for_complex("atan2");
    };
};

template<class T1, class T2>
struct atan2_helper_impl<complex<T1>,complex<T2>,false>
{
    using ret0 = typename md::unify_types<T1, T2>::type;
    using ret   = typename md::real_type<ret0>::type;

    static ret eval(const complex<T1>& arg1, const complex<T2>& arg2)	
    {	
        if (imag(arg1) == 0 && imag(arg2) == 0)
            return atan2_helper_impl<T1,T2,false>::eval(real(arg1),real(arg2)); 
        else
            throw error::function_not_defined_for_complex("atan2");
    };
};

template<class T>
struct atan2_helper_impl<T,Object,true>
{
    using ret = Object;
    static ret eval( const T& arg1,const Object& arg2)	
    {	
        return dynamic::atan2(Object(arg1),arg2);
    };
};
template<class T>
struct atan2_helper_impl<Object,T,true>
{
    using ret = Object;
    static ret eval(const Object& arg1,const T& arg2)	
    {	
        return dynamic::atan2(arg1,Object(arg2));
    };
};
template<>
struct atan2_helper_impl<Object,Object,true>
{
    using ret = Object;
    static ret eval(const Object& arg1,const Object& arg2)	
    {	
        return dynamic::atan2(arg1,arg2);
    };
};
template<class T1, class T2>
struct atan2_helper
{
    using S1 = typename promote_int<T1> :: type;
    using S2 = typename promote_int<T2> :: type;
    static const bool iso = std::is_same<S1,Object>::value 
                            ||std::is_same<S2,Object>::value;
    using ret = typename atan2_helper_impl<S1,S2,iso>::ret;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        return atan2_helper_impl<S1,S2,iso>::eval(arg1,arg2);
    };
};

//hypot
template<class T1, class T2,bool iso>
struct hypot_helper_impl
{
    using ret0  = typename md::unify_types2<T1, T2, Float>::type;
    using ret   = typename md::real_type<ret0>::type;

    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        return scal_func::hypot((ret)arg1,(ret)arg2);
    };
};

template<class T1, class T2>
struct hypot_helper_impl<complex<T1>,T2,false>
{
    using ret0  = typename md::unify_types<T1, T2>::type;
    using ret   = typename md::real_type<ret0>::type;
    using T1P   = typename md::unify_types<T1, ret>::type;
    using T2P   = typename md::unify_types<T2, ret>::type;

    static ret eval(const complex<T1>& arg1, const T2& arg2)	
    {	
        T1P h1  = hypot_helper<T1P,T1P>::eval(real(arg1),imag(arg1));
        return hypot_helper<T1P,T2>::eval(h1,arg2);
    };
};

template<class T1, class T2>
struct hypot_helper_impl<T1,complex<T2>,false>
{
    using ret0  = typename md::unify_types<T1, T2>::type;
    using ret   = typename md::real_type<ret0>::type;
    using T1P   = typename md::unify_types<T1, ret>::type;
    using T2P   = typename md::unify_types<T2, ret>::type;

    static ret eval(const T1& arg1,const complex<T2>& arg2)	
    {	
        T2P h2  = hypot_helper<T2P,T2P>::eval(real(arg2),imag(arg2));
        return hypot_helper<T1,T2P>::eval(arg1,h2);
    };
};

template<class T1, class T2>
struct hypot_helper_impl<complex<T1>,complex<T2>,false>
{
    using ret0  = typename md::unify_types<T1, T2>::type;
    using ret   = typename md::real_type<ret0>::type;
    using T1P   = typename md::unify_types<T1, ret>::type;
    using T2P   = typename md::unify_types<T2, ret>::type;

    static ret eval(const complex<T1>& arg1,const complex<T2>& arg2)	
    {	
        T1P h1  = hypot_helper<T1P,T1P>::eval(real(arg1),imag(arg1));
        T2P h2  = hypot_helper<T2P,T2P>::eval(real(arg2),imag(arg2));
        return hypot_helper<T1P,T2P>::eval(h1,h2);
    };
};

template<class T>
struct hypot_helper_impl<T,Object,true>
{
    using ret = Object;
    static ret eval( const T& arg1,const Object& arg2)	
    {	
        return dynamic::hypot(Object(arg1),arg2);
    };
};
template<class T>
struct hypot_helper_impl<Object,T,true>
{
    using ret = Object;
    static ret eval(const Object& arg1,const T& arg2)	
    {	
        return dynamic::hypot(arg1,Object(arg2));
    };
};
template<>
struct hypot_helper_impl<Object,Object,true>
{
    using ret = Object;
    static ret eval(const Object& arg1,const Object& arg2)	
    {	
        return dynamic::hypot(arg1,arg2);
    };
};

template<class T1, class T2>
struct hypot_helper
{
    using S1 = typename promote_int<T1> :: type;
    using S2 = typename promote_int<T2> :: type;
    static const bool iso = std::is_same<S1,Object>::value 
                            ||std::is_same<S2,Object>::value;
    using ret = typename hypot_helper_impl<S1,S2,iso>::ret;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        return hypot_helper_impl<S1,S2,iso>::eval(arg1,arg2);
    };
};

//mod
template<class ret, bool isc,bool iso,class S1,class S2>
struct mod_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        if (arg2 == S2())
            return arg2;
        if (arg1 == S1())
            return std::copysign(ret(), arg2);

        if (scal_func::isinf(arg2) == true)
        {
            if (scal_func::isinf(arg1) == true)
                return constants::nan<ret>();

            //general formula will produce nan
            if (scal_func::signbit(arg1) == scal_func::signbit(arg2))
                return arg1;
            else
                return arg2;
        }

        ret n   = floor_helper<ret>::eval(gl_div(arg1,arg2));
        ret res = gl_minus(arg1, gl_mult(n, arg2));
        return std::copysign(res, arg2);
    };
};
template<class ret,class S1,class S2>
struct mod_helper_impl<ret,true,false,S1,S2>
{
    using S1R   = typename md::real_type<S1>::type;
    using S2R   = typename md::real_type<S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {
        if (imag(arg1) == 0 && imag(arg2) == 0)
            return mod_helper_impl<ret, false, false, S1R, S2R>::eval(real(arg1), real(arg2));
        else
            throw error::function_not_defined_for_complex("mod");
    };
};
template<class ret,bool isc,class S1,class S2>
struct mod_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::mod(Object(arg1),Object(arg2));
    };
};

template<class T1, class T2>
struct mod_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret0  = typename md::unify_types<S1,S2>::type;
    using ret   = typename md::real_type<ret0>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value  ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value 
                                ||std::is_same<S2,Object>::value;
        return mod_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);			
    };
};
template<>
struct mod_helper<Integer,Integer>
{
    using T1    = Integer;
    using T2    = Integer;
    using ret   = Integer;
    static ret eval(const T1& arg1, const T2& arg2)						
    {	
        if (arg2 == 0 || arg1 == 0)
            return 0;

        Integer ret = arg1 % arg2; 

        //ret must have the same sign as arg2
        if (ret < 0)
            return arg2 > 0 ? ret + arg2 : ret;
        if (ret > 0)
            return arg2 < 0 ? ret + arg2 : ret;

        return ret;
    };
};

//rem
template<class ret,bool isc,bool iso,class S1,class S2>
struct rem_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        if (arg2 == S2() || arg1 == S1())
            return std::copysign(ret(), arg1);

        ret n = trunc_helper<ret>::eval(gl_div(arg1,arg2));

        if (n == 0)
            return arg1;

        ret out = gl_minus(arg1, gl_mult(n,arg2)); 
        return std::copysign(out, arg1);
    };
};
template<class ret,class S1,class S2>
struct rem_helper_impl<ret,true,false,S1,S2>
{
    using S1R   = typename md::real_type<S1>::type;
    using S2R   = typename md::real_type<S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {
        if (imag(arg1) == 0 && imag(arg2) == 0)
            return rem_helper_impl<ret, false, false, S1R, S2R>::eval(real(arg1), real(arg2));
        else
            throw error::function_not_defined_for_complex("rem");
    };
};
template<class ret,bool isc,class S1,class S2>
struct rem_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return dynamic::rem(Object(arg1),Object(arg2));
    };
};

template<class T1, class T2>
struct rem_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret0  = typename md::unify_types<S1,S2>::type;
    using ret   = typename md::real_type<ret0>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value  ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value 
                                ||std::is_same<S2,Object>::value;
        return rem_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);			
    };
};
template<>
struct rem_helper<Integer,Integer>
{
    using T1    = Integer;
    using T2    = Integer;
    using ret   = Integer;
    static ret eval(const T1& arg1, const T2& arg2)	
    {	
        if (arg2 == 0 || arg1 == 0)
            return 0;

        return arg1 % arg2;
    };
};

//idiv
template<class ret,bool isc,bool iso, class S1,class S2>
struct idiv_helper_impl
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_div(arg1,arg2);
    };
};
template<class ret,class S1,class S2>
struct idiv_helper_impl<ret,true,false,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return md::div_c(arg1,arg2);
    };
};
template<class ret,bool isc, class S1,class S2>
struct idiv_helper_impl<ret,isc,true,S1,S2>
{
    static ret eval(const S1& arg1, const S2& arg2)
    {
        return gl_idiv(arg1,arg2);
    };
};

template<class T1, class T2>
struct idiv_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = typename promote_int<T2> :: type;
    using ret   = typename md::unify_types<S1,S2>::type;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        static const bool isc = md::is_complex<S1>::value 
                                ||md::is_complex<S2>::value;
        static const bool iso = std::is_same<S1,Object>::value 
                                ||std::is_same<S2,Object>::value;
        return idiv_helper_impl<ret,isc,iso,S1,S2>::eval(arg1,arg2);	
    };
};
template<>
struct idiv_helper<Integer,Integer>
{
    using T1    = Integer;
    using T2    = Integer;
    using ret   = Integer;
    static ret eval( T1 arg1, T2 arg2)						
    {	
        //avoids unhandled exceptions -> integer overflow
        return	(arg2 == 0)? 0 : gl_idiv(arg1,arg2); 
    };
};

//pow
template<class Ret, class T1, class T2>
struct MATCL_SCALAR_EXPORT pow_complex
{
    static Ret eval(const T1& arg1, const T2& arg2);
};

template<typename T, typename S1, typename S2, bool isc>
struct pow_helper_eval_impl
{
    static T eval(const S1& arg1, const S2& arg2)
    {
        return std::pow(arg1,arg2); 
    };
};

template<typename T, typename S1, typename S2>
struct pow_helper_eval_impl<T, S1, S2, true>
{
    static T eval(const S1& arg1, const S2& arg2)
    {
        return pow_complex<T, S1,S2>::eval(arg1,arg2);
    };
};

// implement power without conversion to complex
template<class T1, class T2,
        bool Is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct pow_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = T2;
    using ret   = typename md::select_if
                <   std::is_same<T2,Integer>::value, 
                    S1, 
                    typename md::unify_types<S1,S2>::type
                >::type;

    static const bool isc = md::is_complex<S1>::value ||md::is_complex<S2>::value;

    static ret eval(const S1& arg1, const S2& arg2)
    {	
        return pow_helper_eval_impl<ret, S1, S2, isc>::eval( arg1, arg2 );
    };
};

template<class T1, class T2>
struct pow_helper<T1, T2, true>
{
    using ret = Object;
    static Object eval(const T1& A, const T2& B)
    {
        Object obj1(A);
        Object obj2(B);
        return dynamic::pow(std::move(obj1),std::move(obj2));
    };
};

// pow return type
template<class S1, class S2>
struct pow_return
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using type  = typename pow_helper<SP1,SP2>::ret;
};

// pow_c return type
template<class S1, class S2>
struct pow_c_return
{
    using SP1   = typename md::promote_scalar<S1>::type;
    using SP2   = typename md::promote_scalar<S2>::type;
    using type0 = typename pow_helper<SP1,SP2>::ret;
    using type  = typename md::unify_types<type0,Float_complex>::type;
};

// implements pow_c
template<class T1, class T2, 
        bool Is_obj = md::is_object<T1>::value || md::is_object<T2>::value>
struct pow_c_helper
{
    using S1    = typename promote_int<T1> :: type;
    using S2    = T2;
    using ret   = typename pow_c_return<T1, T2>::type;

    static ret eval(const S1& A, const S2& B)
    {
        static const bool isc = md::is_complex<ret>::value;

        if (isc || mrd::test_range_pow<S1,S2>::eval(A,B))
            return ret(mrd::pow_helper<S1,S2>::eval(A,B));
        else
            return pow_helper_eval_impl<ret, ret, S2, isc>::eval( (ret)A, B );
    };
};

template<class T1, class T2>
struct pow_c_helper<T1, T2, true>
{
    using ret = Object;

    static Object eval(const T1& A, const T2& B)
    {
        Object obj1(A);
        Object obj2(B);
        return dynamic::pow_c(std::move(obj1),std::move(obj2));
    };
};

}}};

#pragma warning(pop)