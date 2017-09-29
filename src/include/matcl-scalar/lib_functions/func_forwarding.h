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

#include "matcl-scalar/object.h"
#include "matcl-scalar/lib_functions/func_unary.h"
#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-core/general/result_of.h"
#include "matcl-dynamic/object_type_traits.h"

// different functions defined for scalar types that should be forwarded
// to one of existing function; these functions are defined only for scalar
// types defined in other libraries

// following functions can be redefined; template versions are allowed,
// but redefinined functions cannot have variadic template arguments unless
// at least one additional argument is required

namespace matcl
{

//--------------------------------------------------------------
//          UNARY AND BINARY FUNCTIONS
//--------------------------------------------------------------

//missing functions for string
inline bool cast_bool(const std::string& s)  { return s.empty()? false : true; };
inline bool operator!(const std::string& s)  { return s.empty()? true : false; };

// equivalent to (bool)x
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_bool<S>::type
cast_bool(const S& x,Args...)   { return (bool)(x); };

template<class ... Args>
struct get_elem_0{};

template<class T, class ... Args>
struct get_elem_0<T, Args...>{ using type = T;};

// equivalent to arg function
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_arg<S>::type
angle(const S& x,Args...)       { return arg(x); };

// equivalent to exp2 function
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_exp2<S>::type
pow2(const S& x,Args...)        { return exp2(x); };

// equivalent to exp10 function
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_exp10<S>::type
pow10(const S& x,Args...)       { return exp10(x); };

// equivalent to trunc function
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_trunc<S>::type
fix(const S& x,Args...)         { return trunc(x); };

// equivalent to itrunc function
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_itrunc<S>::type
ifix(const S& x,Args...)        { return itrunc(x); };

// equivalent to casting to bool
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type,
    class Enable2= typename result_of::result_of_cast_bool<S>::type>
bool 
is_true(const S& x,Args...)     { return cast_bool(x); };

// equivalent to operator!
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type,
    class Enable2= typename result_of::result_of_cast_bool<S>::type>
bool 
is_false(const S& x,Args...)    { return !cast_bool(x); };

// equivalent to operator!
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type,
    class Enable2= typename result_of::result_of_cast_bool<S>::type>
bool 
neg(const S& x,Args...)         { return !cast_bool(x); };

template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type,
    class Enable2= typename result_of::result_of_cast_bool<S>::type>
bool 
is_zero(const S& x,Args...)    { return matcl::dynamic::object_type_traits<S>::is_zero(x); };

template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type,
    class Enable2= typename result_of::result_of_cast_bool<S>::type>
bool 
is_one(const S& x,Args...)    { return matcl::dynamic::object_type_traits<S>::is_one(x); };

// equivalent to inv
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_inv<S>::type
invs(const S& x,Args...)        { return inv(x); };

// equivalent to unary operator-
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_uminus<S>::type
uminus(const S& x,Args...)      { return -(x); };

// operator+
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
S operator+(const S& x,Args...) { return x; };

// call sqrt for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_sqrt<typename make_complex_type<S>::type>::type
sqrt_c(const S& x,Args...)      { return sqrt(typename make_complex_type<S>::type(x)); };

// call sqrt1pm1 for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_sqrt1pm1<typename make_complex_type<S>::type>::type
sqrt1pm1_c(const S& x,Args...)  { return sqrt1pm1(typename make_complex_type<S>::type(x)); };

// call log for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_log<typename make_complex_type<S>::type>::type
log_c(const S& x,Args...)       { return log(typename make_complex_type<S>::type(x)); };

// call log1p for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_log1p<typename make_complex_type<S>::type>::type
log1p_c(const S& x,Args...)     { return log1p(typename make_complex_type<S>::type(x)); };

// call log2 for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_log2<typename make_complex_type<S>::type>::type
log2_c(const S& x,Args...)      { return log2(typename make_complex_type<S>::type(x)); };

// call log10 for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_log10<typename make_complex_type<S>::type>::type
log10_c(const S& x,Args...)     { return log10(typename make_complex_type<S>::type(x)); };

// call asin for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_asin<typename make_complex_type<S>::type>::type
asin_c(const S& x,Args...)      { return asin(typename make_complex_type<S>::type(x)); };

// call acos for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_acos<typename make_complex_type<S>::type>::type
acos_c(const S& x,Args...)      { return acos(typename make_complex_type<S>::type(x)); };

// call asec for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_asec<typename make_complex_type<S>::type>::type
asec_c(const S& x,Args...)      { return asec(typename make_complex_type<S>::type(x)); };

// call acsc for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_acsc<typename make_complex_type<S>::type>::type
acsc_c(const S& x,Args...)      { return acsc(typename make_complex_type<S>::type(x)); };

// call acosh for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_acosh<typename make_complex_type<S>::type>::type
acosh_c(const S& x,Args...)     { return acosh(typename make_complex_type<S>::type(x)); };

// call atanh for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_atanh<typename make_complex_type<S>::type>::type
atanh_c(const S& x,Args...)     { return atanh(typename make_complex_type<S>::type(x)); };

// call acoth for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_acoth<typename make_complex_type<S>::type>::type
acoth_c(const S& x,Args...)     { return acoth(typename make_complex_type<S>::type(x)); };

// call asech for complex
template<class S, class ... Args, 
    class Enable = typename md::enable_if_external_scalar<S,sizeof...(Args) == 0,void>::type>
typename result_of::result_of_asech<typename make_complex_type<S>::type>::type
asech_c(const S& x,Args...)     { return asech(typename make_complex_type<S>::type(x)); };

// equivalent to x + y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_plus<S1, S2>::type
plus(const S1& x, const S2& y, Args...)     { return x + y; };

// equivalent to x - y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_minus<S1, S2>::type
minus(const S1& x, const S2& y, Args...)    { return x - y; };

// equivalent to x * y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_op_mul<S1, S2>::type
mul(const S1& x, const S2& y, Args...)      { return x * y; };

// equivalent to x * y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_op_mul<S1, S2>::type
mmul(const S1& x, const S2& y, Args...)     { return x * y; };

// equivalent to x / y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_div<S1, S2>::type
div(const S1& x, const S2& y, Args...)      { return x / y; };

// return x / y; 0/0 gives 0
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_div<S1, S2>::type
div_0(const S1& x, const S2& y, Args...)
{
    using ret  = typename result_of::result_of_div<S1, S2>::type;

    using OT1   = object_type<S1>;
    using OT2   = object_type<S2>;

    if (OT1::is_zero(x) && OT2::is_zero(y))
        return ret();
    else
        return x / y; 
};

// return x / y; 0/0 gives 1
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_div<S1, S2>::type
div_1(const S1& x, const S2& y, Args...)
{
    using ret  = typename result_of::result_of_div<S1, S2>::type;

    using OT1   = object_type<S1>;
    using OT2   = object_type<S2>;
    using OTR   = object_type<ret>;

    if (OT1::is_zero(x) && OT2::is_zero(y))
        return OTR::make_one();
    else
        return x / y; 
};

// equivalent to x * y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_op_mul<S1, S2>::type
kron(const S1& x, const S2& y, Args...)     { return x * y; };

// equivalent to x == y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_eeq<S1, S2>::type
eeq(const S1& x, const S2& y, Args...)      { return x == y; };

// equivalent to x == y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_eeq<S1, S2>::type
eeq_nan(const S1& x, const S2& y, Args...)  { return x == y; };

// equivalent to x != y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_neq<S1, S2>::type
neq(const S1& x, const S2& y, Args...)      { return x != y; };

// equivalent to x != y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_neq<S1, S2>::type
neq_nan(const S1& x, const S2& y, Args...)  { return x != y; };

// equivalent to x >= y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_geq<S1, S2>::type
geq(const S1& x, const S2& y, Args...)      { return x >= y; };

// equivalent to x <= y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_leq<S1, S2>::type
leq(const S1& x, const S2& y, Args...)      { return x <= y; };

// equivalent to x > y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_gt<S1, S2>::type
gt(const S1& x, const S2& y, Args...)       { return x > y; };

// equivalent to x < y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_lt<S1, S2>::type
lt(const S1& x, const S2& y, Args...)       { return x < y; };

// equivalent to cast_bool(x) ^ cast_bool(y)
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
op_xor(const S1& x, const S2& y, Args...)   { return (cast_bool(x) ^ cast_bool(y)) ? true : false; };

// equivalent to cast_bool(x) ^ cast_bool(y)
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
elem_xor(const S1& x, const S2& y, Args...) { return (cast_bool(x) ^ cast_bool(y)) ? true : false; };

// equivalent to cast_bool(x) || cast_bool(y)
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
op_or(const S1& x, const S2& y, Args...)    { return cast_bool(x) || cast_bool(y); };

// equivalent to cast_bool(x) || cast_bool(y)
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
elem_or(const S1& x, const S2& y, Args...)  { return cast_bool(x) || cast_bool(y); };

// equivalent to cast_bool(x) && cast_bool(y)
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
op_and(const S1& x, const S2& y, Args...)   { return cast_bool(x) && cast_bool(y); };

// equivalent to cast_bool(x) && cast_bool(y)
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
elem_and(const S1& x, const S2& y, Args...) { return cast_bool(x) && cast_bool(y); };

// equivalent to cast_bool(x) && cast_bool(y)
template<class S1, class ... Args, 
    class Enable1= typename std::enable_if<sizeof...(Args) == 1>::type,
    class S2     = typename std::tuple_element<0, std::tuple<Args...>>::type,
    class Enable2= typename md::enable_if_external_scalar2<S1,S2,true>::type,
    class Enable3= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
operator&&(const S1& x, const Args&... y)   { return op_and(x,y...); };

// equivalent to cast_bool(x) && cast_bool(y)
template<class S1, class ... Args, 
    class Enable1= typename std::enable_if<sizeof...(Args) == 1>::type,
    class S2     = typename std::tuple_element<0, std::tuple<Args...>>::type,
    class Enable2= typename md::enable_if_external_scalar2<S1,S2,true>::type,
    class Enable3= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
operator&(const S1& x, const Args&... y)    { return elem_and(x,y...); };

// equivalent to cast_bool(x) || cast_bool(y)
template<class S1, class ... Args, 
    class Enable1= typename std::enable_if<sizeof...(Args) == 1>::type,
    class S2     = typename std::tuple_element<0, std::tuple<Args...>>::type,
    class Enable2= typename md::enable_if_external_scalar2<S1,S2,true>::type,
    class Enable3= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
operator||(const S1& x, const Args&... y)   { return op_or(x,y...); };

// equivalent to cast_bool(x) || cast_bool(y)
template<class S1, class ... Args, 
    class Enable1= typename std::enable_if<sizeof...(Args) == 1>::type,
    class S2     = typename std::tuple_element<0, std::tuple<Args...>>::type,
    class Enable2= typename md::enable_if_external_scalar2<S1,S2,true>::type,
    class Enable3= typename std::enable_if<result_of::has_cast_bool<S1>::value 
                        && result_of::has_cast_bool<S2>::value>::type>
bool 
operator|(const S1& x, const Args&... y)    { return elem_or(x,y...); };

// equivalent to (x > y) ? x : y
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename result_of::result_of_gt<S1, S2>::type>
typename md::unify_types<S1,S2>::type
max(const S1& x, const S2& y, Args...) 
{
    using ret = typename md::unify_types<S1,S2>::type;
    return (x > y)? ret(x) : ret(y); 
};

// equivalent to (x > y) ? y : x
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type,
    class Enable2= typename result_of::result_of_gt<S1, S2>::type>
typename md::unify_types<S1,S2>::type
min(const S1& x, const S2& y, Args...) 
{
    using ret = typename md::unify_types<S1,S2>::type;
    return (x > y)? ret(y) : ret(x); 
};

// call pow for complex
template<class S1, class S2, class ... Args, 
    class Enable = typename md::enable_if_external_scalar2<S1,S2,sizeof...(Args) == 0>::type>
typename result_of::result_of_pow<typename make_complex_type<S1>::type, S2>::type
pow_c(const S1& x, const S2& y, Args...)      
{ 
    return pow(typename make_complex_type<S1>::type(x), y); 
};

//TODO
/*
//--------------------------------------------------------------
//                  IO FUNCTIONS
//--------------------------------------------------------------

// disp function
template<class S1,
    class Enable = typename md::enable_if_external_scalar<S1,true>::type>
void disp(const S1& x, const disp_stream_ptr& os = default_disp_stream(),
                                     const options& opts = options())  
{ 
    matcl::disp(object_type<S1>(x), os, opts);
};

template<class S1,
    class Enable = typename md::enable_if_external_scalar<S1,true>::type>
void disp_header(const S1& x, const disp_stream_ptr& os = default_disp_stream(),
                                     const options& opts = options())  
{ 
    matcl::disp_header(object_type<S1>(x), os, opts);
};

template<class S1, class ... Args,
    class Enable = typename md::enable_if_external_scalar<S1,sizeof...(Args) == 0>::type>
std::string to_string(const S1& x, Args...)  
{ 
    return matcl::to_string(object_type<S1>(x));
};
*/

};
