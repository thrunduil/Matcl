/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017
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
#include <fenv.h> 

namespace matcl { namespace simd { namespace details
{

template<class T>   struct get_uint_type{};
template<>          struct get_uint_type<float>     {   using type = uint32_t; };
template<>          struct get_uint_type<int32_t>   {   using type = uint32_t; };
template<>          struct get_uint_type<double>    {   using type = uint64_t; };
template<>          struct get_uint_type<int64_t>   {   using type = uint64_t; };

template<class T>   struct get_int_type{};
template<>          struct get_int_type<float>      {   using type = int32_t; };
template<>          struct get_int_type<int32_t>    {   using type = int32_t; };
template<>          struct get_int_type<double>     {   using type = int64_t; };
template<>          struct get_int_type<int64_t>    {   using type = int64_t; };

template<class T>
constexpr T eps_inv() { return T(1) / std::numeric_limits<T>::epsilon(); };

template<class T>
struct simd_bool_base;

template<>
struct simd_bool_base<float>
{
    force_inline
    static float get_val_true()
    {
        return true_value<float>::get();
    };

    static constexpr float val_false    = 0.0f;
};

template<>
struct simd_bool_base<double>
{
    force_inline
    static double get_val_true()
    {
        return true_value<double>::get();
    };

    static constexpr double val_false   = 0.0;
};

template<>
struct simd_bool_base<int32_t>
{
    force_inline
    static int32_t get_val_true()
    {
        return true_value<int32_t>::get();
    };

    static constexpr int32_t val_false   = 0;
};

template<>
struct simd_bool_base<int64_t>
{
    force_inline
    static int64_t get_val_true()
    {
        return true_value<int64_t>::get();
    };

    static constexpr int64_t val_false   = 0;
};

template<class T, int Bits>
struct simd_reverse<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        static const int last = vector_size - 1;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = x.data[last - i];

        return res;
    }
};

template<class T, int Bits>
struct simd_mult<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = x.data[i] * y.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_div<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = x.data[i] / y.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_plus<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = x.data[i] + y.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_minus<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = x.data[i] - y.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_uminus<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = -x.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_sum_all<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static T eval(const simd_type& x)
    {
        T res = x.data[0];

        for (int i = 1; i < vector_size; ++i)
            res += x.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_sub_add<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    static_assert(vector_size % 2 == 0 && vector_size >= 2, "invalid simd type");

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; i += 2)
        {
            res.data[i]     = x.data[i] - y.data[i];
            res.data[i+1]   = x.data[i+1] + y.data[i+1];
        };

        return res;
    };
};

template<class T, int Bits>
struct simd_fma_f<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fma_f(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_fms_f<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fms_f(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

//
template<class T, int Bits>
struct simd_fnma_f<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fnma_f(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_fnms_f<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fnms_f(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_fma_a<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fma_a(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_fms_a<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fms_a(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

//
template<class T, int Bits>
struct simd_fnma_a<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fnma_a(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_fnms_a<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::fnms_a(x.data[i], y.data[i], z.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_abs<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = std::abs(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_bitwise_or<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::bitwise_or(x.data[i], y.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_bitwise_xor<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::bitwise_xor(x.data[i], y.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_bitwise_and<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::bitwise_and(x.data[i], y.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_bitwise_andnot<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::bitwise_andnot(x.data[i], y.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_bitwise_not<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::bitwise_not(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_shift_left<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::shift_left(x.data[i], y);

        return res;
    };
};

template<class T, int Bits>
struct simd_shift_right<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::shift_right(x.data[i], y);

        return res;
    };
};

template<class T, int Bits>
struct simd_shift_right_arithmetic<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, unsigned int y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::shift_right_arithmetic(x.data[i], y);

        return res;
    };
};

template<class T, int Bits>
struct simd_max<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = (x.data[i] > y.data[i]) ? x.data[i] : y.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_min<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = (x.data[i] < y.data[i]) ? x.data[i] : y.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_sqrt<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = std::sqrt(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_round<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;
        
        for (int i = 0; i < vector_size; ++i)
            res.data[i] = matcl::simd::scalar_func::round(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_floor<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::floor(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_ceil<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::ceil(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_trunc<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::trunc(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_eeq<T, Bits, nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        T val_true  = get_val_true();

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = (x.data[i] == y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_neq<T, Bits, nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        T val_true  = get_val_true();

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = (x.data[i] != y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_lt<T, Bits, nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        const T val_true  = get_val_true();

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::lt(x.data[i], y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_gt<T, Bits, nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        const T val_true  = get_val_true();

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::gt(x.data[i], y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_leq<T, Bits, nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        const T val_true  = get_val_true();

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::leq(x.data[i], y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_geq<T, Bits, nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        const T val_true  = get_val_true();

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::geq(x.data[i], y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_any_nan<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static bool eval(const simd_type& x)
    {
        bool res = false;
        for (int i = 0; i < vector_size; ++i)
            res |= matcl::simd::scalar_func::is_nan(x.data[i]);

        return res;
    };
};

template<class T, int Bits>
struct simd_is_nan<T, Bits, nosimd_tag> : simd_bool_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x)
    {
        simd_type res;

        const T val_true  = get_val_true();

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = matcl::simd::scalar_func::is_nan(x.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_any<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static bool eval(const simd_type& x)
    {        
        bool res = true;
        for (int i = 0; i < vector_size; ++i)
            res &= (x.data[i] == T());

        return !res;
    };
};

template<class T, int Bits>
struct simd_all<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static bool eval(const simd_type& x)
    {
        bool res = false;
        for (int i = 0; i < vector_size; ++i)
            res |= (x.data[i] == T());

        return !res;
    };
};

//-----------------------------------------------------------------------
//                   MATHEMATICAL FUNCTIONS
//-----------------------------------------------------------------------
template<class T, int Bits>
struct simd_pow2k<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    using int_type  = typename details::get_int_type<T>::type;
    using simd_int  = simd<int_type, Bits, nosimd_tag>;

    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& k)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = matcl::simd::scalar_func::pow2k(k.data[i]);

        return res;
    };

    force_inline
    static simd_type eval_i(const simd_int& k)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = matcl::simd::scalar_func::pow2ki(k.data[i]);

        return res;
    };
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<class T, int Bits>
struct simd_if_then_else<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;

    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& test, const simd_type& val_true,
                          const simd_type& val_false)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = scalar_func::if_then_else(test.data[i], val_true.data[i], val_false.data[i]);

        return res;
    };
};

}}}
