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

#include "matcl-simd/basic_complex_functions.h"
#include "matcl-simd/details/func/simd_func_complex_def.h"
#include "matcl-simd/details/arch/simd_func_complex_impl.h"

namespace matcl { namespace simd
{

namespace ms = matcl::simd;

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::reverse(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return details::simd_compl_reverse<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::conj(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return details::simd_compl_conj<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::operator*(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_mult<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator*(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_mult<Val, Bits, Simd_tag>::eval_rc(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator*(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_mult<Val, Bits, Simd_tag>::eval_cr(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::operator/(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_div<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator/(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_div<Val, Bits, Simd_tag>::eval_rc(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator/(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_div<Val, Bits, Simd_tag>::eval_cr(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator+(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_plus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator-(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return details::simd_compl_minus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::operator-(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return details::simd_compl_uminus<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
typename simd_compl<Val, Bits, Simd_tag>::value_type
ms::horizontal_sum(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return details::simd_compl_horizontal_sum<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
std::ostream& ms::operator<<(std::ostream& os, const simd_compl<Val, Bits, Simd_tag>& x)
{
    int vec_size    = simd_compl<Val, Bits, Simd_tag>::vector_size;
    const Val* ptr  = x.get_raw_ptr();

    os << "{" << ptr[0];

    for (int i = 1; i < vec_size; ++i)
        os << ", " << ptr[i];

    os << "}";

    return os;
};

// return true if at least element in the vector x is NAN
template<class Val, int Bits, class Simd_tag>
force_inline bool
any_nan(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return any_nan(x.data);
}

}}
