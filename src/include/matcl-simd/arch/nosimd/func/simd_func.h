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

#include "matcl-simd/arch/simd_impl.h"
#include "matcl-simd/func/simd_func_def.h"
#include <fenv.h> 

namespace matcl { namespace simd
{

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
struct simd_fma<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = x.data[i] * y.data[i] + z.data[i];

        return res;
    };
};

template<class T, int Bits>
struct simd_fms<T, Bits, nosimd_tag>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y, const simd_type& z)
    {
        simd_type res;

        for (int i = 0; i < vector_size; ++i)
            res.data[i] = x.data[i] * y.data[i] - z.data[i];

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

        ::fesetround(FE_TONEAREST);
        
        for (int i = 0; i < vector_size; ++i)
            res.data[i] = std::nearbyint(x.data[i]);

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
            res.data[i] = std::floor(x.data[i]);

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
            res.data[i] = std::ceil(x.data[i]);

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
            res.data[i] = std::trunc(x.data[i]);

        return res;
    };
};

template<class T>
struct simd_cmp_base;

template<>
struct simd_cmp_base<float>
{
    force_inline
    static float get_val_true()
    {
        return true_value<float>::get();
    };

    static constexpr float val_false    = 0.0f;
};

template<>
struct simd_cmp_base<double>
{
    force_inline
    static double get_val_true()
    {
        return true_value<double>::get();
    };

    static constexpr double val_false   = 0.0;
};

template<class T, int Bits>
struct simd_eeq<T, Bits, nosimd_tag> : simd_cmp_base<T>
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
struct simd_neq<T, Bits, nosimd_tag> : simd_cmp_base<T>
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
struct simd_lt<T, Bits, nosimd_tag> : simd_cmp_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        T val_true  = get_val_true();

        // do not use comparison operators; 
        // for example VS likes to change x <= y to !(x > y), which is clearly wrong
        for (int i = 0; i < vector_size; ++i)
            res.data[i] = std::isless(x.data[i], y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_gt<T, Bits, nosimd_tag> : simd_cmp_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        T val_true  = get_val_true();

        // do not use comparison operators; 
        // for example VS likes to change x <= y to !(x > y), which is clearly wrong
        for (int i = 0; i < vector_size; ++i)
            res.data[i] = std::isgreater(x.data[i], y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_leq<T, Bits, nosimd_tag> : simd_cmp_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        T val_true  = get_val_true();

        // do not use comparison operators; 
        // for example VS likes to change x <= y to !(x > y), which is clearly wrong
        for (int i = 0; i < vector_size; ++i)
            res.data[i] = std::islessequal(x.data[i], y.data[i]) ? val_true : val_false;

        return res;
    };
};

template<class T, int Bits>
struct simd_geq<T, Bits, nosimd_tag> : simd_cmp_base<T>
{
    using simd_type = simd<T, Bits, nosimd_tag>;
    
    static const int 
    vector_size     = simd_type::vector_size;

    force_inline
    static simd_type eval(const simd_type& x, const simd_type& y)
    {
        simd_type res;

        T val_true  = get_val_true();

        // do not use comparison operators; 
        // for example VS likes to change x <= y to !(x > y), which is clearly wrong
        for (int i = 0; i < vector_size; ++i)
            res.data[i] = std::isgreaterequal(x.data[i], y.data[i]) ? val_true : val_false;

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
            res |= matcl::is_nan(x.data[i]);

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

}}
