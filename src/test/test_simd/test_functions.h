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

namespace test_functions
{

namespace ms = matcl::simd;

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

template<int Bits, class Tag>
struct test_cast<ms::simd<float, Bits, Tag>>
{
    using simd_type = ms::simd<float, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_double = typename simd_type::simd_double;
        simd_double xd_l = x.convert_low_to_double();
        simd_double xd_h = x.convert_high_to_double();

        return simd_type(xd_l.convert_to_float(), xd_h.convert_to_float());
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
};

template<class Tag>
struct test_cast<ms::simd<double, 128, Tag>>
{
    using simd_type = ms::simd<double, 128, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_float = typename simd_type::simd_float;
        simd_float xf = x.convert_to_float();

        return simd_type(xf.convert_low_to_double());
    };
};

template<class Tag>
struct test_cast<ms::simd<double, 256, Tag>>
{
    using simd_type = ms::simd<double, 256, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_float = typename simd_type::simd_float;
        simd_float xf = x.convert_to_float();

        return simd_type(xf.convert_to_double());
    };
};

template<int Bits, class Tag>
struct test_cast<ms::simd_compl<float, Bits, Tag>>
{
    using simd_type = ms::simd_compl<float, Bits, Tag>;

    static simd_type eval(const simd_type& x)
    {
        using simd_double = typename simd_type::simd_double;
        simd_double xd_l = x.convert_low_to_double();
        simd_double xd_h = x.convert_high_to_double();

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
        return bitwise_or(signbit_base(x), T::one()); 
    }

    static std::string name()
    { 
        return "sign_base"; 
    };
};

struct Func_sum_all
{
    template<class T>    
    force_inline static T eval(const T& x) 
    { 
        return T(sum_all(x)); 
    }

    static std::string name()
    { 
        return "sum_all"; 
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

}
