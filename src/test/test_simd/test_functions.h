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

#include "test_simd_config.h"
#include "matcl-core/config.h"
#include "matcl-simd/math_functions.h"
#include "matcl-mp/matcl_mp.h"

namespace test_functions
{

namespace ms = matcl::simd;

template<class T>
struct convert{};

template<int Bits, class Tag>
struct convert<ms::simd<float, Bits, Tag>>
{
    using simd_ret = ms::simd<float, Bits, Tag>;

    template<class Arg>
    static simd_ret eval(const Arg& x)
    {
        return x.convert_to_float();
    }
};

template<int Bits, class Tag>
struct convert<ms::simd<double, Bits, Tag>>
{
    using simd_ret = ms::simd<double, Bits, Tag>;

    template<class Arg>
    static simd_ret eval(const Arg& x)
    {
        return x.convert_to_double();
    }
};

template<class T>
struct cast{};

template<int Bits, class Tag>
struct cast<ms::simd<float, Bits, Tag>>
{
    using simd_ret = ms::simd<float, Bits, Tag>;

    template<class Arg>
    static simd_ret eval(const Arg& x)
    {
        return x.reinterpret_as_float();
    }
};

template<int Bits, class Tag>
struct cast<ms::simd<double, Bits, Tag>>
{
    using simd_ret = ms::simd<double, Bits, Tag>;

    template<class Arg>
    static simd_ret eval(const Arg& x)
    {
        return x.reinterpret_as_double();
    }
};

template<class T>
struct make_zero{};

template<class T, int Bits, class Tag>
struct make_zero<ms::simd<T,Bits, Tag>>
{
    using simd_type = ms::simd<T,Bits, Tag>;

    static simd_type eval()
    {
        return simd_type::zero();
    }
};

// missing scalar functions
template<class T>
inline T reverse(const T& x)    { return x; };

template<class T>
struct test_cast
{};

template<class T>
struct test_scal_cast
{};

template<int Bits, class Tag>
struct test_cast<ms::simd<float, Bits, Tag>>
{
    using simd_type = ms::simd<float, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xd_l = x.convert_low_to_double();
        auto xd_h = x.convert_high_to_double();

        return simd_type(xd_l.convert_to_float(), xd_h.convert_to_float());
    }

    static simd_type eval_int32(const simd_type& x)
    {
        auto xi = x.convert_to_int32();
        return xi.convert_to_float();
    }
    static simd_type eval_int64(const simd_type& x)
    {
        auto xd_l = x.convert_low_to_int64();
        auto xd_h = x.convert_high_to_int64();

        return simd_type(xd_l.convert_to_float(), xd_h.convert_to_float());
    }
};

template<int Bits, class Tag>
struct test_scal_cast<ms::simd<float, Bits, Tag>>
{
    using simd_type = ms::simd<float, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xd = x.convert_to_double();
        return xd.convert_to_float();
    }

    static simd_type eval_int32(const simd_type& x)
    {
        auto xi = x.convert_to_int32();
        return xi.convert_to_float();
    }
    static simd_type eval_int64(const simd_type& x)
    {
        auto xi = x.convert_to_int64();
        return xi.convert_to_float();
    }
};

template<int Bits, class Tag>
struct test_cast<ms::simd<double, Bits, Tag>>
{
    using simd_type = ms::simd<double, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_float = typename simd_type::simd_float;
        simd_float xf = x.cast_to_float();

        return simd_type(xf.cast_low_to_double(), xf.cast_high_to_double());
    };

    static simd_type eval_int32(const simd_type& x)
    {
        auto xi = x.convert_to_int32();
        return xdi.convert_to_double();
    }
    static simd_type eval_int64(const simd_type& x)
    {
        auto xi = x.convert_to_int64();
        return xdi.convert_to_double();
    }
};

template<int Bits, class Tag>
struct test_scal_cast<ms::simd<double, Bits, Tag>>
{
    using simd_type = ms::simd<double, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_float = typename simd_type::simd_float;
        simd_float xf = x.convert_to_float();
        return xf.convert_to_double();
    };

    static simd_type eval_int32(const simd_type& x)
    {
        auto xi = x.convert_to_int32();
        return xi.convert_to_double();
    }
    static simd_type eval_int64(const simd_type& x)
    {
        auto xi = x.convert_to_int64();
        return xi.convert_to_double();
    }
};

template<int Bits>
struct test_scal_cast<ms::simd<double, Bits, ms::sse_tag>>
{
    using simd_type = ms::simd<double, Bits, ms::sse_tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_float = typename simd_type::simd_float;
        simd_float xf = x.convert_to_float();
        return xf.convert_low_to_double();
    };

    static simd_type eval_int32(const simd_type& x)
    {
        auto xi = x.convert_to_int32();
        return xi.convert_low_to_double();
    }
    static simd_type eval_int64(const simd_type& x)
    {
        auto xi = x.convert_to_int64();
        return xi.convert_to_double();
    }
};

template<class Tag>
struct test_cast<ms::simd<double, 128, Tag>>
{
    using simd_type = ms::simd<double, 128, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xf = x.convert_to_float();

        return simd_type(xf.convert_low_to_double());
    };

    static simd_type eval_int32(const simd_type& x)
    {
        auto xi = x.convert_to_int32();
        return xi.convert_low_to_double();
    }
    static simd_type eval_int64(const simd_type& x)
    {
        auto xi = x.convert_to_int64();
        return xi.convert_to_double();
    }
};

template<class Tag>
struct test_cast<ms::simd<double, 256, Tag>>
{
    using simd_type = ms::simd<double, 256, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xf = x.convert_to_float();

        return simd_type(xf.convert_to_double());
    };

    static simd_type eval_int32(const simd_type& x)
    {
        auto xi = x.convert_to_int32();
        return simd_type(xi.convert_low_to_double(), xi.convert_high_to_double());
    }
    static simd_type eval_int64(const simd_type& x)
    {
        auto xi = x.convert_to_int64();
        return xi.convert_to_double();
    }
};

template<int Bits, class Tag>
struct test_cast<ms::simd_compl<float, Bits, Tag>>
{
    using simd_type = ms::simd_compl<float, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xd_l = x.convert_low_to_double();
        auto xd_h = x.convert_high_to_double();

        return simd_type(xd_l.convert_to_float(), xd_h.convert_to_float());
    }
};

template<class Tag>
struct test_cast<ms::simd_compl<double, 128, Tag>>
{
    using simd_type = ms::simd_compl<double, 128, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_float = typename simd_type::simd_float;
        simd_float xf = x.convert_to_float();

        return simd_type(xf.convert_low_to_double());
    };
};

template<class Tag>
struct test_cast<ms::simd_compl<double, 256, Tag>>
{
    using simd_type = ms::simd_compl<double, 256, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_float = typename simd_type::simd_float;
        simd_float xf = x.convert_to_float();

        return simd_type(xf.convert_to_double());
    };
};

//
template<int Bits, class Tag>
struct test_cast<ms::simd<int32_t, Bits, Tag>>
{
    using simd_type = ms::simd<int32_t, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xd_l = x.convert_low_to_int64();
        auto xd_h = x.convert_high_to_int64();

        return simd_type(xd_l.convert_to_int32(), xd_h.convert_to_int32());
    }

    static simd_type eval_float(const simd_type& x)
    {
        auto xi = x.convert_to_float();
        return xi.convert_to_int32();
    }

    static simd_type eval_double(const simd_type& x)
    {
        auto xd_l = x.convert_low_to_double();
        auto xd_h = x.convert_high_to_double();

        return simd_type(xd_l.convert_to_int32(), xd_h.convert_to_int32());
    }
};

template<int Bits, class Tag>
struct test_scal_cast<ms::simd<int32_t, Bits, Tag>>
{
    using simd_type = ms::simd<int32_t, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xd_l = x.convert_to_int64();
        return xd_l.convert_to_int32();
    }

    static simd_type eval_float(const simd_type& x)
    {
        auto xi = x.convert_to_float();
        return xi.convert_to_int32();
    }

    static simd_type eval_double(const simd_type& x)
    {
        auto xd_l = x.convert_to_double();
        return xd_l.convert_to_int32();
    }
};

template<int Bits, class Tag>
struct test_scal_cast<ms::simd<int64_t, Bits, Tag>>
{
    using simd_type = ms::simd<int64_t, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_int32 = typename simd_type::simd_int32;
        simd_int32 xf = x.convert_to_int32();

        return xf.convert_to_int64();
    };

    static simd_type eval_float(const simd_type& x)
    {
        auto xf = x.convert_to_float();
        return xf.convert_to_int64();
    }

    static simd_type eval_double(const simd_type& x)
    {
        auto xd = x.convert_to_double();
        return xd.convert_to_int64();
    }
};

template<int Bits>
struct test_scal_cast<ms::simd<int64_t, Bits, ms::sse_tag>>
{
    using simd_type = ms::simd<int64_t, Bits, ms::sse_tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_int32 = typename simd_type::simd_int32;
        simd_int32 xf = x.convert_to_int32();

        return xf.convert_low_to_int64();
    };

    static simd_type eval_float(const simd_type& x)
    {
        auto xf = x.convert_to_float();
        return xf.convert_low_to_int64();
    }

    static simd_type eval_double(const simd_type& x)
    {
        auto xd = x.convert_to_double();
        return xd.convert_to_int64();
    }
};

template<class Tag>
struct test_cast<ms::simd<int64_t, 128, Tag>>
{
    using simd_type = ms::simd<int64_t, 128, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xf = x.convert_to_int32();

        return simd_type(xf.convert_low_to_int64());
    };

    static simd_type eval_float(const simd_type& x)
    {
        auto xf = x.convert_to_float();
        return simd_type(xf.convert_low_to_int64());
    }

    static simd_type eval_double(const simd_type& x)
    {
        auto xd = x.convert_to_double();
        return xd.convert_to_int64();
    }
};

template<class Tag>
struct test_cast<ms::simd<int64_t, 256, Tag>>
{
    using simd_type = ms::simd<int64_t, 256, Tag>;

    static simd_type eval(const simd_type& x)
    {
        auto xf = x.convert_to_int32();

        return simd_type(xf.convert_to_int64());
    };

    static simd_type eval_float(const simd_type& x)
    {
        auto xf = x.convert_to_float();
        return simd_type(xf.convert_to_int64());
    }

    static simd_type eval_double(const simd_type& x)
    {
        auto xd = x.convert_to_double();
        return xd.convert_to_int64();
    }
};

struct Func_reverse
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return reverse(x); 
    }

    static std::string name()
    { 
        return "reverse"; 
    };
};

struct Func_cast
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_cast<T>::eval(x); 
    }

    static std::string name()
    { 
        return "cast"; 
    };
};

struct Func_scal_cast
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_scal_cast<T>::eval(x); 
    }

    static std::string name()
    { 
        return "scal_cast"; 
    };
};

struct Func_cast_float
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_cast<T>::eval_float(x); 
    }

    static std::string name()
    { 
        return "cast_float"; 
    };
};

struct Func_cast_double
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_cast<T>::eval_double(x); 
    }

    static std::string name()
    { 
        return "cast_double"; 
    };
};

struct Func_scal_cast_float
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_scal_cast<T>::eval_float(x); 
    }

    static std::string name()
    { 
        return "cast_float"; 
    };
};

struct Func_scal_cast_double
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_scal_cast<T>::eval_double(x); 
    }

    static std::string name()
    { 
        return "cast_double"; 
    };
};

struct Func_cast_int32
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_cast<T>::eval_int32(x); 
    }

    static std::string name()
    { 
        return "cast_int32"; 
    };
};

struct Func_cast_int64
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_cast<T>::eval_int64(x); 
    }

    static std::string name()
    { 
        return "cast_int64"; 
    };
};

struct Func_scal_cast_int32
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_scal_cast<T>::eval_int32(x); 
    }

    static std::string name()
    { 
        return "scal_cast_int32"; 
    };
};

struct Func_scal_cast_int64
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return test_scal_cast<T>::eval_int64(x); 
    }

    static std::string name()
    { 
        return "scal_cast_int64"; 
    };
};

struct Func_extract_low
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return T(x.extract_low()); 
    }

    static std::string name()
    { 
        return "extr_low"; 
    };
};

struct Func_extract_high
{    
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return T(x.extract_high()); 
    }

    static std::string name()
    { 
        return "extr_high"; 
    };
};

struct Func_uminus
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return -x; 
    }

    static std::string name()
    { 
        return "uminus"; 
    };
};

struct Func_any_nan
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return T(any_nan(x)); 
    }

    static std::string name()
    { 
        return "any_nan"; 
    };
};

struct Func_is_nan
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return is_nan(x); 
    }

    static std::string name()
    { 
        return "is_nan"; 
    };
};

struct Func_is_finite
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return is_finite(x); 
    }

    static std::string name()
    { 
        return "is_finite"; 
    };
};

struct Func_bit_not
{
    template<class T>    
    force_inline static T eval(const T& x1)
    { 
        return bitwise_not(x1); 
    }
    
    static std::string name()
    { 
        return "bit_not"; 
    };
};

struct Func_any
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return T(any(lt(x, make_zero<T>::eval()))); 
    }

    static std::string name()
    { 
        return "any"; 
    };
};

struct Func_all
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return T(all(lt(x, make_zero<T>::eval()))); 
    }

    static std::string name()
    { 
        return "all"; 
    };
};

struct Func_shift_left
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return shift_left(x, 10); 
    }

    static std::string name()
    { 
        return "shl"; 
    };
};

struct Func_shift_left2
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return shift_left(x, 40); 
    }

    static std::string name()
    { 
        return "shl2"; 
    };
};

struct Func_shift_right
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return shift_right(x, 10); 
    }

    static std::string name()
    { 
        return "shr"; 
    };
};

struct Func_shift_right2
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return shift_right(x, 40); 
    }

    static std::string name()
    { 
        return "shr2"; 
    };
};

struct Func_shift_right_a
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return shift_right_arithmetic(x, 10); 
    }

    static std::string name()
    { 
        return "shra"; 
    };
};

struct Func_shift_right_a2
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return shift_right_arithmetic(x, 40); 
    }

    static std::string name()
    { 
        return "shra2"; 
    };
};

template<class Val>
struct Func_gather_int
{
    const Val* m_ptr;

    Func_gather_int(const Val* ptr)
        :m_ptr(ptr)
    {};

    template<int Bits, class Tag>    
    force_inline matcl::simd::simd<Val, Bits, Tag> 
    eval(const matcl::simd::simd<int32_t, Bits, Tag>& x) const
    { 
        return matcl::simd::simd<Val, Bits, Tag>::gather(m_ptr, x); 
    }

    static std::string name()
    { 
        return "gtr_int"; 
    };
};

template<class Val>
struct Func_gather_int32
{
    const Val* m_ptr;

    Func_gather_int32(const Val* ptr)
        :m_ptr(ptr)
    {};

    template<int Bits, class Tag>    
    force_inline matcl::simd::simd<Val, Bits, Tag> 
    eval(const matcl::simd::simd<int32_t, Bits, Tag>& x) const
    { 
        return matcl::simd::simd<Val, Bits, Tag>::gather(m_ptr, x); 
    }

    static std::string name()
    { 
        return "gtr_int32"; 
    };
};

template<class Val>
struct Func_gather_int64
{
    const Val* m_ptr;

    Func_gather_int64(const Val* ptr)
        :m_ptr(ptr)
    {};

    template<int Bits, class Tag>    
    force_inline matcl::simd::simd<Val, Bits, Tag> 
    eval(const matcl::simd::simd<int64_t, Bits, Tag>& x) const
    { 
        return matcl::simd::simd<Val, Bits, Tag>::gather(m_ptr, x); 
    }

    static std::string name()
    { 
        return "gtr_int64"; 
    };
};

template<class Val>
struct Func_pow2ki_int32
{};

template<>
struct Func_pow2ki_int32<float>
{
    template<int Bits, class Tag>    
    force_inline matcl::simd::simd<float, Bits, Tag> 
    eval(const matcl::simd::simd<int32_t, Bits, Tag>& x) const
    { 
        return pow2ki(x); 
    }

    static std::string name()
    { 
        return "pow2ki_int32"; 
    };
};

template<>
struct Func_pow2ki_int32<double>
{
    template<int Bits, class Tag>    
    force_inline matcl::simd::simd<double, Bits, Tag> 
    eval(const matcl::simd::simd<int32_t, Bits, Tag>& x) const
    { 
        (void)x;
        return matcl::simd::simd<double, Bits, Tag>::zero();
    }

    static std::string name()
    { 
        return "pow2ki_int32"; 
    };
};

template<class Val>
struct Func_pow2ki_int64
{};

template<>
struct Func_pow2ki_int64<double>
{
    template<int Bits, class Tag>    
    force_inline matcl::simd::simd<double, Bits, Tag> 
    eval(const matcl::simd::simd<int64_t, Bits, Tag>& x) const
    { 
        return pow2ki(x); 
    }

    static std::string name()
    { 
        return "pow2ki_int64"; 
    };
};

template<>
struct Func_pow2ki_int64<float>
{
    template<int Bits, class Tag>    
    force_inline matcl::simd::simd<float, Bits, Tag> 
    eval(const matcl::simd::simd<int64_t, Bits, Tag>& x) const
    { 
        (void)x;
        return matcl::simd::simd<float, Bits, Tag>::zero();
    }

    static std::string name()
    { 
        return "pow2ki_int64"; 
    };
};

struct Func_pow2k
{  
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return pow2k(x); 
    }

    static std::string name()
    { 
        return "pow2k"; 
    };
};

struct Func_conj
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return conj(x); 
    }

    static std::string name()
    { 
        return "conj"; 
    };
};

struct Func_abs
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return abs(x); 
    }

    static std::string name()
    { 
        return "abs"; 
    };
};

struct Func_signbit_base
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return signbit_base(x); 
    }

    static std::string name()
    { 
        return "sign_base"; 
    };
};

struct Func_hor_sum
{
    template<class T>    
    force_inline static T eval(const T& x) 
    { 
        return T(horizontal_sum(x)); 
    }

    static std::string name()
    { 
        return "hor_sum"; 
    };
};

struct Func_hor_min
{
    template<class T>    
    force_inline static T eval(const T& x) 
    { 
        return T(horizontal_min(x)); 
    }

    static std::string name()
    { 
        return "hor_min"; 
    };
};

struct Func_hor_max
{
    template<class T>    
    force_inline static T eval(const T& x) 
    { 
        return T(horizontal_max(x)); 
    }

    static std::string name()
    { 
        return "hor_max"; 
    };
};

struct Func_sqrt
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return sqrt(x); 
    }

    static std::string name()
    { 
        return "sqrt"; 
    };
};

struct Func_exp
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return exp(x); 
    }

    template<class T>    
    force_inline static T eval_base(const T& x)
    { 
        return std::exp(x); 
    }

    static matcl::mp_float eval_mp(const matcl::mp_float& x)
    { 
        return exp(x); 
    }

    static std::string name()
    { 
        return "exp"; 
    };
};

struct Func_sin
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return sin(x); 
    }

    template<class T>    
    force_inline static T eval_base(const T& x)
    { 
        return std::sin(x); 
    }

    static matcl::mp_float eval_mp(const matcl::mp_float& x)
    { 
        return sin(x); 
    }

    static std::string name()
    { 
        return "sin"; 
    };
};

struct Func_cos
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return cos(x); 
    }

    template<class T>    
    force_inline static T eval_base(const T& x)
    { 
        return std::cos(x); 
    }

    static matcl::mp_float eval_mp(const matcl::mp_float& x)
    { 
        return cos(x); 
    }

    static std::string name()
    { 
        return "cos"; 
    };
};

struct Func_tan
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return tan(x); 
    }

    template<class T>    
    force_inline static T eval_base(const T& x)
    { 
        return std::tan(x); 
    }

    static matcl::mp_float eval_mp(const matcl::mp_float& x)
    { 
        return tan(x); 
    }

    static std::string name()
    { 
        return "tan"; 
    };
};

struct Func_cot
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return cot(x); 
    }

    template<class T>    
    force_inline static T eval_base(const T& x)
    { 
        return T(1) / std::tan(x); 
    }

    static matcl::mp_float eval_mp(const matcl::mp_float& x)
    { 
        return cot(x); 
    }

    static std::string name()
    { 
        return "cot"; 
    };
};

struct Func_log
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return log(x); 
    }

    template<class T>    
    force_inline static T eval_base(const T& x)
    { 
        return std::log(x); 
    }

    static matcl::mp_float eval_mp(const matcl::mp_float& x)
    { 
        return log(x); 
    }

    static std::string name()
    { 
        return "log"; 
    };
};

struct Func_fraction
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return fraction(x); 
    }

    static std::string name()
    { 
        return "fraction"; 
    };
};

struct Func_iexponent
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return cast<T>::eval(iexponent(x)); 
    }

    static std::string name()
    { 
        return "iexponent"; 
    };
};

struct Func_exponent
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return exponent(x); 
    }

    static std::string name()
    { 
        return "exponent"; 
    };
};

struct Func_round
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return round(x); 
    }

    static std::string name()
    { 
        return "round"; 
    };
};

struct Func_floor
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return floor(x); 
    }
    
    static std::string name()
    { 
        return "floor"; 
    };
};

struct Func_ceil
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return ceil(x); 
    }
    
    static std::string name() 
    { 
        return "ceil"; 
    };
};

struct Func_trunc
{
    template<class T>    
    force_inline static T eval(const T& x)
    { 
        return trunc(x); 
    }
    
    static std::string name()
    { 
        return "trunc"; 
    };
};

struct Func_mult
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 * x2; 
    }

    static std::string name()
    { 
        return "mult"; 
    };
};

struct Func_mult_RC
{
    template<class T1, class T2>    
    force_inline static T2 eval(const T1& x1, const T2& x2)
    { 
        return x1 * x2; 
    }

    static std::string name()
    { 
        return "mult rc"; 
    };
};

struct Func_mult_CR
{
    template<class T1, class T2>    
    force_inline static T1 eval(const T1& x1, const T2& x2)
    { 
        return x1 * x2; 
    }

    static std::string name()
    { 
        return "mult cr"; 
    };
};

struct Func_div
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 / x2; 
    }

    static std::string name()
    { 
        return "div"; 
    };
};

struct Func_div_RC
{
    template<class T1, class T2>    
    force_inline static T2 eval(const T1& x1, const T2& x2)
    { 
        return x1 / x2; 
    }

    static std::string name()
    { 
        return "div rc"; 
    };
};

struct Func_div_CR
{
    template<class T1, class T2>    
    force_inline static T1 eval(const T1& x1, const T2& x2)
    { 
        return x1 / x2; 
    }

    static std::string name()
    { 
        return "div cr"; 
    };
};

struct Func_plus
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 + x2; 
    }
    
    static std::string name()
    { 
        return "plus"; 
    };
};

struct Func_minus
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return x1 - x2; 
    }

    static std::string name()
    { 
        return "minus"; 
    };
};

struct Func_sub_add
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return sub_add(x1, x2); 
    }
    
    static std::string name()
    { 
        return "sub_add"; 
    };
};

struct Func_max
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return max(x1, x2); 
    }
    
    static std::string name()
    { 
        return "max"; 
    };
};

struct Func_min
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return min(x1, x2); 
    }

    static std::string name()
    { 
        return "min"; 
    };
};

struct Func_eeq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return eeq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "eeq"; 
    };
};

struct Func_neq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return neq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "neq"; 
    };
};

struct Func_leq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return leq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "leq"; 
    };
};

struct Func_geq
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return geq(x1, x2); 
    }
    
    static std::string name()
    { 
        return "geq"; 
    };
};

struct Func_lt
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return lt(x1, x2); 
    }

    static std::string name()
    { 
        return "lt"; 
    };
};

struct Func_gt
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return gt(x1, x2); 
    }

    static std::string name()
    { 
        return "gt"; 
    };
};

struct Func_fma_f
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fma_f(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fma_f"; 
    };
};

struct Func_fma_a
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fma_a(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fma_a"; 
    };
};

struct Func_if_then_else
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        T cond  = gt(x1, T::zero());
        return if_then_else(cond, x2, x3); 
    }

    static std::string name()
    { 
        return "if_then_else"; 
    };
};

struct Func_fnma_f
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fnma_f(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fnma_f"; 
    };
};

struct Func_fnma_a
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fnma_a(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fnma_a"; 
    };
};

struct Func_fms_f
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fms_f(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fms_f"; 
    };
};

struct Func_fms_a
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fms_a(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fms_a"; 
    };
};

struct Func_fnms_f
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fnms_f(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fnms_f"; 
    };
};

struct Func_fnms_a
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2, const T& x3)
    { 
        return fnms_a(x1, x2, x3); 
    }

    static std::string name()
    { 
        return "fnms_a"; 
    };
};

struct Func_bit_or
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return bitwise_or(x1, x2); 
    }
    
    static std::string name()
    { 
        return "bit_or"; 
    };
};

struct Func_bit_xor
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return bitwise_xor(x1, x2); 
    }
    
    static std::string name()
    { 
        return "bit_xor"; 
    };
};

struct Func_bit_and
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return bitwise_and(x1, x2); 
    }
    
    static std::string name()
    { 
        return "bit_and"; 
    };
};

struct Func_bit_andnot
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        return bitwise_andnot(x1, x2); 
    }
    
    static std::string name()
    { 
        return "bit_andnot"; 
    };
};

struct Func_if_zero_else
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        T cond  = gt(x1, T::zero());
        return if_zero_else(cond, x2); 
    }
    
    static std::string name()
    { 
        return "if_zero_else"; 
    };
};

struct Func_if_nan_else
{
    template<class T>    
    force_inline static T eval(const T& x1, const T& x2)
    { 
        T cond  = gt(x1, T::zero());
        return if_nan_else(cond, x2); 
    }
    
    static std::string name()
    { 
        return "if_nan_else"; 
    };
};

struct Func_select_1
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x)
    { 
        return x.select<0,1,2,3>(); 
    }
    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x)
    { 
        return x.select<0,1,2,3>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x)
    { 
        return x.select<0,1>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x)
    { 
        return x.select<0,1>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<float, 256, Tag> eval(const ms::simd<float, 256, Tag>& x)
    { 
        return x.select<0,1,2,3,4,5,6,7>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 256, Tag> eval(const ms::simd<int32_t, 256, Tag>& x)
    { 
        return x.select<0,1,2,3,4,5,6,7>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 256, Tag> eval(const ms::simd<double, 256, Tag>& x)
    { 
        return x.select<0,1,2,3>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 256, Tag> eval(const ms::simd<int64_t, 256, Tag>& x)
    { 
        return x.select<0,1,2,3>(); 
    }

    static std::string name()
    { 
        return "select 1"; 
    };
};

struct Func_select_2
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x)
    { 
        return x.select<3,2,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x)
    { 
        return x.select<3,2,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x)
    { 
        return x.select<1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x)
    { 
        return x.select<1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<float, 256, Tag> eval(const ms::simd<float, 256, Tag>& x)
    { 
        return x.select<7,6,5,4,3,2,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 256, Tag> eval(const ms::simd<int32_t, 256, Tag>& x)
    { 
        return x.select<7,6,5,4,3,2,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 256, Tag> eval(const ms::simd<double, 256, Tag>& x)
    { 
        return x.select<3,2,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 256, Tag> eval(const ms::simd<int64_t, 256, Tag>& x)
    { 
        return x.select<3,2,1,0>(); 
    }

    static std::string name()
    { 
        return "select 2"; 
    };
};

struct Func_select_3
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x)
    { 
        return x.select<2,3,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x)
    { 
        return x.select<2,3,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x)
    { 
        return x.select<1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x)
    { 
        return x.select<1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<float, 256, Tag> eval(const ms::simd<float, 256, Tag>& x)
    { 
        return x.select<4,5,6,7,0,1,2,3>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 256, Tag> eval(const ms::simd<int32_t, 256, Tag>& x)
    { 
        return x.select<4,5,6,7,0,1,2,3>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 256, Tag> eval(const ms::simd<double, 256, Tag>& x)
    { 
        return x.select<2,3,0,1>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 256, Tag> eval(const ms::simd<int64_t, 256, Tag>& x)
    { 
        return x.select<2,3,0,1>(); 
    }

    static std::string name()
    { 
        return "select 3"; 
    };
};

struct Func_select_4
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x)
    { 
        return x.select<0,1,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x)
    { 
        return x.select<0,1,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x)
    { 
        return x.select<0,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x)
    { 
        return x.select<0,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<float, 256, Tag> eval(const ms::simd<float, 256, Tag>& x)
    { 
        return x.select<0,1,2,3,3,2,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 256, Tag> eval(const ms::simd<int32_t, 256, Tag>& x)
    { 
        return x.select<0,1,2,3,3,2,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 256, Tag> eval(const ms::simd<double, 256, Tag>& x)
    { 
        return x.select<0,1,1,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 256, Tag> eval(const ms::simd<int64_t, 256, Tag>& x)
    { 
        return x.select<0,1,1,0>(); 
    }

    static std::string name()
    { 
        return "select 4"; 
    };
};

struct Func_select_5
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x)
    { 
        return x.select<2,3,3,2>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x)
    { 
        return x.select<2,3,3,2>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x)
    { 
        return x.select<1,1>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x)
    { 
        return x.select<1,1>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<float, 256, Tag> eval(const ms::simd<float, 256, Tag>& x)
    { 
        return x.select<4,5,6,7,7,6,5,4>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 256, Tag> eval(const ms::simd<int32_t, 256, Tag>& x)
    { 
        return x.select<4,5,6,7,7,6,5,4>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 256, Tag> eval(const ms::simd<double, 256, Tag>& x)
    { 
        return x.select<2,3,3,2>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 256, Tag> eval(const ms::simd<int64_t, 256, Tag>& x)
    { 
        return x.select<2,3,3,2>(); 
    }

    static std::string name()
    { 
        return "select 5"; 
    };
};

struct Func_select_6
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x)
    { 
        return x.select<1,0,2,3>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x)
    { 
        return x.select<1,0,2,3>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x)
    { 
        return x.select<1,1>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x)
    { 
        return x.select<1,1>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<float, 256, Tag> eval(const ms::simd<float, 256, Tag>& x)
    { 
        return x.select<3,2,1,0,4,5,6,7>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 256, Tag> eval(const ms::simd<int32_t, 256, Tag>& x)
    { 
        return x.select<3,2,1,0,4,5,6,7>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 256, Tag> eval(const ms::simd<double, 256, Tag>& x)
    { 
        return x.select<0,1,3,2>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 256, Tag> eval(const ms::simd<int64_t, 256, Tag>& x)
    { 
        return x.select<0,1,3,2>(); 
    }

    static std::string name()
    { 
        return "select 6"; 
    };
};

struct Func_select_7
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x)
    { 
        return x.select<1,1,2,2>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x)
    { 
        return x.select<1,1,2,2>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x)
    { 
        return x.select<0,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x)
    { 
        return x.select<0,0>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<float, 256, Tag> eval(const ms::simd<float, 256, Tag>& x)
    { 
        return x.select<2,2,3,3,6,6,7,7>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 256, Tag> eval(const ms::simd<int32_t, 256, Tag>& x)
    { 
        return x.select<2,2,3,3,6,6,7,7>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 256, Tag> eval(const ms::simd<double, 256, Tag>& x)
    { 
        return x.select<1,1,3,3>(); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 256, Tag> eval(const ms::simd<int64_t, 256, Tag>& x)
    { 
        return x.select<1,1,3,3>(); 
    }

    static std::string name()
    { 
        return "select 7"; 
    };
};

struct Func_combine_1
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<0,1,0,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<0,1,0,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<0,0>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<0,0>(x, y); 
    }

    static std::string name()
    { 
        return "combine 1"; 
    };
};

struct Func_combine_2
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<5,6,5,7>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<5,6,5,7>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<0,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<0,1>(x, y); 
    }

    static std::string name()
    { 
        return "combine 2"; 
    };
};

struct Func_combine_3
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<0,1,7,6>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<0,1,7,6>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<0,2>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<0,2>(x, y); 
    }

    static std::string name()
    { 
        return "combine 3"; 
    };
};

struct Func_combine_4
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<7,6,0,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<7,6,0,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<1,0>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<1,0>(x, y); 
    }

    static std::string name()
    { 
        return "combine 4"; 
    };
};

struct Func_combine_5
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<0,4,1,5>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<0,4,1,5>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<1,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<1,1>(x, y); 
    }

    static std::string name()
    { 
        return "combine 5"; 
    };
};

struct Func_combine_6
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<4,0,5,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<4,0,5,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<1,2>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<1,2>(x, y); 
    }

    static std::string name()
    { 
        return "combine 6"; 
    };
};

struct Func_combine_7
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<2,6,3,7>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<2,6,3,7>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<1,3>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<1,3>(x, y); 
    }

    static std::string name()
    { 
        return "combine 7"; 
    };
};

struct Func_combine_8
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<6,2,7,3>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<6,2,7,3>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<0,3>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<0,3>(x, y); 
    }

    static std::string name()
    { 
        return "combine 8"; 
    };
};

struct Func_combine_9
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<3,2,7,6>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<3,2,7,6>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<2,0>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<2,0>(x, y); 
    }

    static std::string name()
    { 
        return "combine 9"; 
    };
};

struct Func_combine_10
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<7,6,3,2>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<7,6,3,2>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<2,1>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<2,1>(x, y); 
    }

    static std::string name()
    { 
        return "combine 10"; 
    };
};

struct Func_combine_11
{
    template<class Tag>    
    force_inline static 
    ms::simd<float, 128, Tag> eval(const ms::simd<float, 128, Tag>& x, 
                                  const ms::simd<float, 128, Tag>& y)
    { 
        return ms::simd<float, 128, Tag>::combine<0,4,7,3>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int32_t, 128, Tag> eval(const ms::simd<int32_t, 128, Tag>& x, 
                                  const ms::simd<int32_t, 128, Tag>& y)
    { 
        return ms::simd<int32_t, 128, Tag>::combine<0,4,7,3>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<double, 128, Tag> eval(const ms::simd<double, 128, Tag>& x, 
                                  const ms::simd<double, 128, Tag>& y)
    { 
        return ms::simd<double, 128, Tag>::combine<2,3>(x, y); 
    }

    template<class Tag>    
    force_inline static 
    ms::simd<int64_t, 128, Tag> eval(const ms::simd<int64_t, 128, Tag>& x, 
                                  const ms::simd<int64_t, 128, Tag>& y)
    { 
        return ms::simd<int64_t, 128, Tag>::combine<2,3>(x, y); 
    }

    static std::string name()
    { 
        return "combine 11"; 
    };
};

}
