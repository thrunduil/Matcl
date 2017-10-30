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

#include "matcl-simd/func/simd_func_complex.h"
#include "matcl-simd/func/simd_func_complex_def.h"
#include "matcl-simd/arch/simd_func_complex_impl.h"

namespace matcl { namespace simd
{

namespace ms = matcl::simd;

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::reverse(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return simd_compl_reverse<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::operator*(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return simd_compl_mult<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator*(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return simd_compl_mult<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator*(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y)
{
    return simd_compl_mult<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::operator/(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return simd_compl_div<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator+(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return simd_compl_plus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::operator-(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y)
{
    return simd_compl_minus<Val, Bits, Simd_tag>::eval(x, y);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag> 
ms::operator-(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return simd_compl_uminus<Val, Bits, Simd_tag>::eval(x);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::fma(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return simd_compl_fma<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::fma(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return simd_compl_fma<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::fms(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return simd_compl_fms<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
ms::fms(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return simd_compl_fms<Val, Bits, Simd_tag>::eval(x, y, z);
};

template<class Val, int Bits, class Simd_tag>
force_inline
typename simd_compl<Val, Bits, Simd_tag>::value_type
ms::sum_all(const simd_compl<Val, Bits, Simd_tag>& x)
{
    return simd_compl_sum_all<Val, Bits, Simd_tag>::eval(x);
};

}}
