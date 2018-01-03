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

#include "matcl-core/float/fma_dekker.h"
#include "matcl-simd/details/utils.h"
#include <cmath>

namespace matcl { namespace simd { namespace scalar_func
{

template<class T> struct integer_type{};
template<>        struct integer_type<float>    { using type = int32_t;};
template<>        struct integer_type<double>   { using type = int64_t;};

template<class T>
struct is_signed_impl
{
    using int_type  = typename integer_type<T>::type;

    static bool eval(const T& x)
    {
        int_type xi     = *reinterpret_cast<const int_type*>(&x);
        return xi < 0;
    }
};

//-------------------------------------------------------------------
//                         round
//-------------------------------------------------------------------
template<class T>
struct round_impl
{
    force_inline
    static T eval(T x)
    {
        // std::round has wrong tie rule (we require round to even in case
        // of ties); std::nearbyint is very slow;

        const T einv = details::eps_inv<T>();
        
        T xa        = std::abs(x);
        T tmp       = xa + einv;
        tmp         = tmp - einv;
        T resa      = (xa < einv) ? tmp : xa;
        bool sign   = is_signed_impl<T>::eval(x);
        return sign ? -resa : resa;
    };
};

//-------------------------------------------------------------------
//                         trunc
//-------------------------------------------------------------------
template<class T>
struct maximum_int{};

template<>
struct maximum_int<double>
{
    // 2^53
    static constexpr double value = 9007199254740992.0;
};

template<>
struct maximum_int<float>
{
    // 2^24
    static constexpr float value = 16777216.0;
};

template<class T>
struct trunc_impl
{
    force_inline
    static T eval(T x)
    {
        // std::trunc is insanely slow

        using int_type  = typename integer_type<T>::type;
        const T max_int = maximum_int<T>::value;
        return std::abs(x) < max_int ? static_cast<T>(static_cast<int_type>(x)) : x;
    };
};

template<class T>
struct is_finite_impl
{};

template<>
struct is_finite_impl<double>
{
    force_inline
    static double eval(const double& x)
    {
        int64_t xi      = reinterpret_cast<const int64_t&>(x);

        // mask selecting all bits in the exponent
        int64_t mask    = 0x7FF0000000000000ll;

        xi              = xi & mask;

        // return true if at least one bit in the exponent is not set
        bool res        = (xi != mask);
        return (res == true) ? true_value<double>::get() : 0.0;
    }
};

template<>
struct is_finite_impl<float>
{
    force_inline
    static float eval(const float& x)
    {
        int32_t xi      = reinterpret_cast<const int32_t&>(x);

        // mask selecting all bits in the exponent
        int32_t mask    = 0x7F800000;

        xi              = xi & mask;

        // return true if at least one bit in the exponent is not set
        bool res        = (xi != mask);
        return (res == true) ? true_value<float>::get() : 0.0f;
    }
};

//-------------------------------------------------------------------
//                         missing scalar functions
//-------------------------------------------------------------------

force_inline float floor(float f)       { return std::floor(f); };
force_inline float ceil(float f)        { return std::ceil(f); };
force_inline float trunc(float f)       { return trunc_impl<float>::eval(f); };
force_inline float round(float f)       { return round_impl<float>::eval(f); };

force_inline double floor(double f)     { return std::floor(f); };
force_inline double ceil(double f)      { return std::ceil(f); };
force_inline double trunc(double f)     { return trunc_impl<double>::eval(f); };
force_inline double round(double f)     { return round_impl<double>::eval(f); };

force_inline int32_t floor(int32_t f)   { return f; };
force_inline int32_t ceil(int32_t f)    { return f; };
force_inline int32_t trunc(int32_t f)   { return f; };
force_inline int32_t round(int32_t f)   { return f; };

force_inline int64_t floor(int64_t f)   { return f; };
force_inline int64_t ceil(int64_t f)    { return f; };
force_inline int64_t trunc(int64_t f)   { return f; };
force_inline int64_t round(int64_t f)   { return f; };

force_inline bool is_nan(double f)      { return f != f; }
force_inline bool is_nan(float f)       { return f != f; }

force_inline bool is_finite(double f)   { return is_finite_impl<double>::eval(f); };
force_inline bool is_finite(float f)    { return is_finite_impl<float>::eval(f); };

//-------------------------------------------------------------------
//                         convertions
//-------------------------------------------------------------------
force_inline float      convert_double_float(double x)      { return (float)x; };
force_inline int32_t    convert_double_int32(double x)      { return (int32_t)round(x); };
force_inline int64_t    convert_double_int64(double x)      { return (int64_t)round(x); };

force_inline double     convert_float_double(float x)       { return (double)x; };
force_inline int32_t    convert_float_int32(float x)        { return (int32_t)round(x); };
force_inline int64_t    convert_float_int64(float x)        { return (int64_t)round(x); };

force_inline int64_t    convert_int32_int64(int32_t x)      { return (int64_t)x; }
force_inline double     convert_int32_double(int32_t x)     { return (double)x; }
force_inline float      convert_int32_float(int32_t x)      { return (float)x; }

force_inline int32_t    convert_int64_int32(int64_t x)      { return (int32_t)x; }
force_inline double     convert_int64_double(int64_t x)     { return (double)x; }
force_inline float      convert_int64_float(int64_t x)      { return (float)x; }

//-------------------------------------------------------------------
//                         fma
//-------------------------------------------------------------------

force_inline float fma_f(float x, float y, float z)
{
    return x * y + z;
};

force_inline float fms_f(float x, float y, float z)
{
    return x * y - z;
};

force_inline float fnma_f(float x, float y, float z)
{
    return z - x * y;
};

force_inline float fnms_f(float x, float y, float z)
{
    return -(x * y + z);
};

force_inline double fma_f(double x, double y, double z)
{
    return x * y + z;
};

force_inline double fms_f(double x, double y, double z)
{
    return x * y - z;
};

force_inline double fnma_f(double x, double y, double z)
{
    return z - x * y;
};

force_inline double fnms_f(double x, double y, double z)
{
    return -(x * y + z);
};

force_inline float fma_a(float x, float y, float z)
{
    return fma_dekker(x, y, z);
}

force_inline float fms_a(float x, float y, float z)
{
    return fma_dekker(x, y, -z);
}

force_inline float fnma_a(float x, float y, float z)
{
    return fma_dekker(-x, y, z);
}

force_inline float fnms_a(float x, float y, float z)
{
    return -fma_dekker(x, y, z);
}

force_inline double fma_a(double x, double y, double z)
{
    return fma_dekker(x, y, z);
}

force_inline double fms_a(double x, double y, double z)
{
    return fma_dekker(x, y, -z);
}

force_inline double fnma_a(double x, double y, double z)
{
    return fma_dekker(-x, y, z);
}

force_inline double fnms_a(double x, double y, double z)
{
    return -fma_dekker(x, y, z);
}

//-------------------------------------------------------------------
//                         bit manipulation
//-------------------------------------------------------------------
template<class T>
force_inline
T bitwise_or(const T& x, const T& y)
{
    using uint_type     = details::unsigned_integer_type<T>::type;

    const uint_type* xi   = reinterpret_cast<const uint_type*>(&x);
    const uint_type* yi   = reinterpret_cast<const uint_type*>(&y);

    uint_type res   = (*xi) | (*yi);

    return *reinterpret_cast<const T*>(&res);
};

template<class T>
force_inline
T bitwise_xor(const T& x, const T& y)
{
    using uint_type     = details::unsigned_integer_type<T>::type;

    const uint_type* xi   = reinterpret_cast<const uint_type*>(&x);
    const uint_type* yi   = reinterpret_cast<const uint_type*>(&y);

    uint_type res   = (*xi) ^ (*yi);

    return *reinterpret_cast<const T*>(&res);
};

template<class T>
force_inline
T bitwise_and(const T& x, const T& y)
{
    using uint_type     = details::unsigned_integer_type<T>::type;

    const uint_type* xi   = reinterpret_cast<const uint_type*>(&x);
    const uint_type* yi   = reinterpret_cast<const uint_type*>(&y);

    uint_type res   = (*xi) & (*yi);

    return *reinterpret_cast<const T*>(&res);
};

template<class T>
force_inline
T bitwise_andnot(const T& x, const T& y)
{
    using uint_type     = details::unsigned_integer_type<T>::type;

    const uint_type* xi   = reinterpret_cast<const uint_type*>(&x);
    const uint_type* yi   = reinterpret_cast<const uint_type*>(&y);

    uint_type res   = (~(*xi)) & (*yi);

    return *reinterpret_cast<const T*>(&res);
};

template<class T>
force_inline
T bitwise_not(const T& x)
{
    using uint_type     = details::unsigned_integer_type<T>::type;

    const uint_type* xi   = reinterpret_cast<const uint_type*>(&x);

    uint_type res   = ~(*xi);

    return *reinterpret_cast<const T*>(&res);
};

template<class T>
force_inline
T shift_left(const T& x, unsigned int y)
{
    using uint_type     = details::unsigned_integer_type<T>::type;
        
    static const 
    unsigned max_shift  = sizeof(uint_type) * 8;

    const uint_type* xi = reinterpret_cast<const uint_type*>(&x);
    uint_type ires      = (*xi) << y;
    T res               = *reinterpret_cast<const T*>(&ires);

    return (y >= max_shift) ? T() : res;
};

template<class T>
force_inline
T shift_right(const T& x, unsigned int y)
{
    using uint_type     = details::unsigned_integer_type<T>::type;

    static const 
    unsigned max_shift  = sizeof(uint_type) * 8;

    const uint_type* xi = reinterpret_cast<const uint_type*>(&x);
    uint_type ires      = (*xi) >> y;
    T res               = *reinterpret_cast<const T*>(&ires);

    return (y >= max_shift) ? T() : res;
};

template<class T>
force_inline
T shift_right_arithmetic(const T& x, unsigned int y)
{
    using int_type      = details::integer_type<T>::type;

    static const 
    unsigned max_shift  = sizeof(int_type) * 8 - 1;

    const int_type* xi  = reinterpret_cast<const int_type*>(&x);
    int_type res;

    if(y > max_shift)
        res             = (*xi) >= 0? int_type(0) : int_type(-1);
    else
        res             = (*xi) >> y;

    return *reinterpret_cast<const T*>(&res);
};

//-------------------------------------------------------------------
//                         conditional
//-------------------------------------------------------------------
template<class T>
force_inline
T if_then_else(T test, T val_true, T val_false)
{
    return (test == T()) ? val_false : val_true;
};

//-------------------------------------------------------------------
//                         comparison
//-------------------------------------------------------------------
// we cannot inline this code by hand; VS in this case will change
// x <= y to !(x > y), which is clearly wrong; on the other hand
// function std::islessequal is very slow;
// this dummy function seems to force VS to generate correct code,
// but still not optimal

template<class T>
force_inline
bool lt(const T& x, const T& y)
{
    return x < y;
}

template<class T>
force_inline
bool gt(const T& x, const T& y)
{
    return x > y;
}

template<class T>
force_inline
bool leq(const T& x, const T& y)
{
    return x <= y;
}

template<class T>
force_inline
bool geq(const T& x, const T& y)
{
    return x >= y;
}

}}}
