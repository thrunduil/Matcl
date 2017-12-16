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

#include "matcl-simd/details/arch/simd_impl.h"
#include "matcl-simd/details/func/simd_func_def.h"
#include "matcl-simd/simd.h"

namespace matcl { namespace simd { namespace details
{

template<class T, int Bits>
struct simd_reverse<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return x;
    }
};

template<class T, int Bits>
struct simd_mult<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data * y.data);
    };
};

template<class T, int Bits>
struct simd_div<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data / y.data);
    };
};

template<class T, int Bits>
struct simd_plus<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data + y.data);
    };
};

template<class T, int Bits>
struct simd_minus<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(x.data - y.data);
    };
};

template<class T, int Bits>
struct simd_uminus<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(-x.data);
    };
};

template<class T, int Bits>
struct simd_sum_all<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static T eval(const simd_type& x)
    {
        return x.first();
    };
};

template<class T, int Bits>
struct simd_fma_f<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fma_f(x.data, y.data, z.data));
    };
};

template<class T, int Bits>
struct simd_fms_f<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fms_f(x.data, y.data, z.data));
    };
};

template<class T, int Bits>
struct simd_fnma_f<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fnma_f(x.data, y.data, z.data));
    };
};

template<class T, int Bits>
struct simd_fnms_f<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fnms_f(x.data, y.data, z.data));
    };
};

template<class T, int Bits>
struct simd_fma_a<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fma_a(x.data, y.data, z.data));
    };
};

template<class T, int Bits>
struct simd_fms_a<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fms_a(x.data, y.data, z.data));
    };
};

//
template<class T, int Bits>
struct simd_fnma_a<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fnma_a(x.data, y.data, z.data));
    };
};

template<class T, int Bits>
struct simd_fnms_a<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        return simd_type(scalar_func::fnms_a(x.data, y.data, z.data));
    };
};

template<class T, int Bits>
struct simd_abs<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(std::abs(x.data));
    };
};

template<class T, int Bits>
struct simd_bitwise_or<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(scalar_func::bitwise_or(x.data, y.data));
    };
};

template<class T, int Bits>
struct simd_bitwise_xor<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(scalar_func::bitwise_xor(x.data, y.data));
    };
};

template<class T, int Bits>
struct simd_bitwise_and<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(scalar_func::bitwise_and(x.data, y.data));
    };
};

template<class T, int Bits>
struct simd_bitwise_andnot<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type(scalar_func::bitwise_andnot(x.data, y.data));
    };
};

template<class T, int Bits>
struct simd_bitwise_not<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(scalar_func::bitwise_not(x.data));
    };
};

template<class T, int Bits>
struct simd_shift_left<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return simd_type(scalar_func::shift_left(x.data, y));
    };
};

template<class T, int Bits>
struct simd_shift_right<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return simd_type(scalar_func::shift_right(x.data, y));
    };
};

template<class T, int Bits>
struct simd_shift_right_arithmetic<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        return simd_type(scalar_func::shift_right_arithmetic(x.data, y));
    };
};

template<class T, int Bits>
struct simd_max<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type((x.data > y.data) ? x.data : y.data);
    };
};

template<class T, int Bits>
struct simd_min<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        return simd_type((x.data < y.data) ? x.data : y.data);
    };
};

template<class T, int Bits>
struct simd_sqrt<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(std::sqrt(x.data));
    };
};

template<class T, int Bits>
struct simd_round<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(matcl::simd::scalar_func::round(x.data));
    };
};

template<class T, int Bits>
struct simd_floor<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(matcl::simd::scalar_func::floor(x.data));
    };
};

template<class T, int Bits>
struct simd_ceil<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(matcl::simd::scalar_func::ceil(x.data));
    };
};

template<class T, int Bits>
struct simd_trunc<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        return simd_type(matcl::simd::scalar_func::trunc(x.data));
    };
};

template<class T, int Bits>
struct simd_eeq<T, Bits, scalar_nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        T val_true  = get_val_true();

        return simd_type((x.data == y.data) ? val_true : val_false);
    };
};

template<class T, int Bits>
struct simd_neq<T, Bits, scalar_nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        T val_true  = get_val_true();
        return simd_type((x.data != y.data) ? val_true : val_false);
    };
};

template<class T, int Bits>
struct simd_lt<T, Bits, scalar_nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        T val_true  = get_val_true();
        return simd_type(scalar_func::lt(x.data, y.data) ? val_true : val_false);
    };
};

template<class T, int Bits>
struct simd_gt<T, Bits, scalar_nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        T val_true  = get_val_true();
        return simd_type(scalar_func::gt(x.data, y.data) ? val_true : val_false);
    };
};

template<class T, int Bits>
struct simd_leq<T, Bits, scalar_nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        T val_true  = get_val_true();
        return simd_type(scalar_func::leq(x.data, y.data) ? val_true : val_false);
    };
};

template<class T, int Bits>
struct simd_geq<T, Bits, scalar_nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        T val_true  = get_val_true();
        return simd_type(scalar_func::geq(x.data, y.data) ? val_true : val_false);
    };
};

template<class T, int Bits>
struct simd_any_nan<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static bool eval(const simd_type& x)
    {
        return matcl::simd::scalar_func::is_nan(x.data);
    };
};

template<class T, int Bits>
struct simd_is_nan<T, Bits, scalar_nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static simd_type eval(const simd_type& x)
    {
        T val_true  = get_val_true();
        return simd_type(matcl::simd::scalar_func::is_nan(x.data) ? val_true : val_false);
    };
};

template<class T, int Bits>
struct simd_any<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static bool eval(const simd_type& x)
    {        
        return !(x.data == T());
    };
};

template<class T, int Bits>
struct simd_all<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;
    
    force_inline
    static bool eval(const simd_type& x)
    {
        return !(x.data == T());
    };
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<class T, int Bits>
struct simd_if_then_else<T, Bits, scalar_nosimd_tag>
{
    using simd_type = simd<T, Bits, scalar_nosimd_tag>;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        return simd_type(scalar_func::if_then_else(test.data, val_true.data, val_false.data));
    };
};

}}}
