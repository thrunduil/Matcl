/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-simd/basic_functions.h"
#include "matcl-simd/details/func/simd_func_def.h"
#include "matcl-simd/details/arch/simd_func_impl.h"
#include "matcl-core/float/float_binary_rep.h"
#include "matcl-core/IO/scalar_io.h"

namespace matcl { namespace simd
{

namespace ms = matcl::simd;

template<>
force_inline
float true_value<float>::get()
{
    return matcl::hex_float(0xFFFFFFFF);
};

template<>
force_inline
int32_t true_value<int32_t>::get()
{
    return 0xFFFFFFFF;
};

template<>
force_inline
double true_value<double>::get()
{
    return matcl::hex_double(0xFFFFFFFFFFFFFFFF);
};

template<>
force_inline
int64_t true_value<int64_t>::get()
{
    return 0xFFFFFFFFFFFFFFFF;
};

template<class T>
force_inline
T false_value<T>::get()
{
    return T();
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::signbit_base(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_signbit_base<Val, Bits, Simd_tag>::eval(x);    
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::bitwise_or(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_bitwise_or<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::bitwise_and(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_bitwise_and<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::bitwise_xor(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_bitwise_xor<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::bitwise_andnot(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_bitwise_andnot<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::bitwise_not(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_bitwise_not<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::shift_left(const simd<Val, Bits, Simd_tag>& x, unsigned int count)
{
    return details::simd_shift_left<Val, Bits, Simd_tag>::eval(x, count);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::shift_right(const simd<Val, Bits, Simd_tag>& x, unsigned int count)
{
    return details::simd_shift_right<Val, Bits, Simd_tag>::eval(x, count);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::shift_right_arithmetic(const simd<Val, Bits, Simd_tag>& x, unsigned int count)
{
    return details::simd_shift_right_arithmetic<Val, Bits, Simd_tag>::eval(x, count);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator*(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_mult<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator/(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_div<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator+(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_plus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator-(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_minus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator-(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_uminus<Val, Bits, Simd_tag>::eval(x);
};

//
template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::mult(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_mult<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::div(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_div<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::plus(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_plus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::minus(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_minus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::uminus(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_uminus<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::abs(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_abs<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::sub_add(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_sub_add<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fma_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fma_f<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fms_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fms_f<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fnma_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fnma_f<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fnms_f(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fnms_f<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fma_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fma_a<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fms_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fms_a<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fnma_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fnma_a<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fnms_a(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return details::simd_fnms_a<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
Val ms::horizontal_sum(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_horizontal_sum<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
Val ms::horizontal_min(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_horizontal_min<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
Val ms::horizontal_max(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_horizontal_max<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::max(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_max<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::min(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_min<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::eeq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_eeq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::neq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_neq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::lt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_lt<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::gt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_gt<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::leq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_leq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::geq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_geq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::is_nan(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_is_nan<Val, Bits, Simd_tag>::eval(x);
}

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::is_finite(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_is_finite<Val, Bits, Simd_tag>::eval(x);
}

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::sqrt(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_sqrt<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::round(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_round<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::floor(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_floor<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::ceil(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_ceil<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::trunc(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_trunc<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline bool
ms::any_nan(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_any_nan<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline bool
ms::all(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_all<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline bool
ms::any(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_any<Val, Bits, Simd_tag>::eval(x);
};

//-----------------------------------------------------------------------
//                   CONDITIONAL FUNCTIONS
//-----------------------------------------------------------------------
template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::if_nan_else(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& val_false)
{
    return bitwise_or(val_false, test);
}

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::if_zero_else(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& val_false)
{
    return bitwise_andnot(test, val_false);
}

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::if_then_else(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& val_true,
             const simd<Val, Bits, Simd_tag>& val_false)
{
    return details::simd_if_then_else<Val, Bits, Simd_tag>::eval(test, val_true, val_false);
}

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::if_add(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& x,
             const simd<Val, Bits, Simd_tag>& y)
{
    return x + bitwise_and(test, y);
}

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::if_sub(const simd<Val, Bits, Simd_tag>& test, const simd<Val, Bits, Simd_tag>& x,
             const simd<Val, Bits, Simd_tag>& y)
{
    return x - bitwise_and(test, y);
}

//-----------------------------------------------------------------------
//                   LOGICAL FUNCTIONS
//-----------------------------------------------------------------------

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::operator &&(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return bitwise_and(x, y);
}

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::operator ||(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return bitwise_or(x, y);
}

template<class Val, int Bits, class Simd_tag>
force_inline simd<Val, Bits, Simd_tag> 
ms::operator !(const simd<Val, Bits, Simd_tag>& x)
{
    return bitwise_not(x);
}

//-----------------------------------------------------------------------
//                   MISCELLANEOUS FUNCTIONS
//-----------------------------------------------------------------------

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::reverse(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_reverse<Val, Bits, Simd_tag>::eval(x);
};

template<class Val_ret, class Val, int Bits, class Simd_tag>
force_inline simd<Val_ret, Bits, Simd_tag>
reinterpret_as(const simd<Val, Bits, Simd_tag>& x)
{
    return details::simd_reinterpret_as<Val_ret, Val, Bits, Simd_tag>::eval(x);
}

//-----------------------------------------------------------------------
//                   IO FUNCTIONS
//-----------------------------------------------------------------------
template<class Val, int Bits, class Simd_tag>
std::ostream& ms::operator<<(std::ostream& os, const simd<Val, Bits, Simd_tag>& x)
{
    int vec_size    = simd<Val, Bits, Simd_tag>::vector_size;
    const Val* ptr  = x.get_raw_ptr();

    os << "{";
    
    matcl::details::saveload_scalar_helper::eval_print(os, ptr[0]);

    for (int i = 1; i < vec_size; ++i)
    {
        os << ", ";
        matcl::details::saveload_scalar_helper::eval_print(os, ptr[i]);
    };

    os << "}";

    return os;        
};

template<class Val, int Bits, class Simd_tag>
std::istream& ms::operator>>(std::istream& is, simd<Val, Bits, Simd_tag>& v)
{
    using simd_type = simd<Val, Bits, Simd_tag>;
    int vec_size    = simd_type::vector_size;    

    char c  = 0;

    // consume whitespaces
    while (is)
    {
        is.get(c);

        if (c != ' ' && c != '\t'  && c != '\n')
            break;
    }
    
    Val nan         = std::numeric_limits<Val>::quiet_NaN();    
    v               = simd_type(nan);

    if (is.good() == false)
        return is;

    Val* ptr        = v.get_raw_ptr();

    if (c != '{')
    {
        // this is an error
        is.setstate(std::ios::failbit);
        return is;
    }

    char sep    = ' ';
    char fin    = ' ';

    bool ok     = true;
    
    ok &= matcl::details::saveload_scalar_helper::eval_load(is, ptr[0]);

    for (int i = 0; i < vec_size; ++i)
    {
        is >> sep;

        if (sep != ',')
        {
            // this is an error
            is.setstate(std::ios::failbit);
            return is;
        };

        ok &= matcl::details::saveload_scalar_helper::eval_load(is, ptr[i]);
    };

    is >> fin;

    if (fin != '}')
    {
        // this is an error
        is.setstate(std::ios::failbit);
        return is;
    }

    if (ok == false)
    {
        is.setstate(std::ios::failbit);
        return is;
    };

    return is;
};

}}
