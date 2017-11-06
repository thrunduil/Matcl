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

#include "matcl-simd/func/simd_func.h"
#include "matcl-simd/func/simd_func_def.h"
#include "matcl-simd/arch/simd_func_impl.h"

namespace matcl { namespace simd
{

namespace ms = matcl::simd;

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>
ms::reverse(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_reverse<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator*(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_mult<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator/(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_div<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator+(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_plus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator-(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_minus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::operator-(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_uminus<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::abs(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_abs<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::sub_add(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_sub_add<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fma(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return simd_fma<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::fms(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd<Val, Bits, Simd_tag>& z)
{
    return simd_fms<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
Val ms::sum_all(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_sum_all<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::max(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_max<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::min(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_min<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::eeq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_eeq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::neq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_neq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::lt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_lt<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::gt(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_gt<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag> 
ms::leq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_leq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::geq(const simd<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_geq<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::sqrt(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_sqrt<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::round(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_round<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::floor(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_floor<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::ceil(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_ceil<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd<Val, Bits, Simd_tag>  
ms::trunc(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_trunc<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline bool
any_inf(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_any_inf<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline bool
any_nan(const simd<Val, Bits, Simd_tag>& x)
{
    return simd_any_nan<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
std::ostream& ms::operator<<(std::ostream& os, const simd<Val, Bits, Simd_tag>& x)
{
    int vec_size    = simd<Val, Bits, Simd_tag>::vector_size;

    os << "{" << x.get(0);

    for (int i = 1; i < vec_size; ++i)
        os << ", " << x.get(i);

    os << "}";

    return os;
};

}}
