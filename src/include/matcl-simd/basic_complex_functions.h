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

#include "matcl-simd/simd_general.h"
#include "matcl-simd/details/complex/simd_complex_impl.h"

namespace matcl { namespace simd
{

// vector of elements in reverse order
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag> 
reverse(const simd_compl<Val, Bits, Simd_tag>& x);

// complex conjugate
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag> 
conj(const simd_compl<Val, Bits, Simd_tag>& x);

// vector multiply x * y
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag> 
operator*(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y);

// vector multiply x * y, where x is a real vector
// for x = [x1, x2, ...] and y = [re1, im1, ...] returns [x1*re1, x2*im1, ...]
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag>
operator*(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y);

// vector multiply x * y, where y is a real vector
// for x = [re1, im1, ...] and y = [y1, y2, ...] returns [y1*re1, y2*im1, ...]
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag>
operator*(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector division x / y
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag> 
operator/(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y);

// vector division x / y, where x is a real vector
// for x = [x_1, d_1, ...] and y = [re_1, im_1, ...] returns [c1_r, c1_i, ...], 
// where cj = complex(cj_r, cj_i) and cj = complex(x_i, 0) / complex(re_i, im_i)
// (elements in the x vector on even positions are ignored)
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag> 
operator/(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y);

// vector division x / y, where y is a real vector
// for x = [re1, im1, ...] and y = [y1, y2, ...] returns [re1/y1, im1/y2, ...]
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag> 
operator/(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y);

// vector add x + y
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag>
operator+(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y);

// vector subtract x - y
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag>
operator-(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y);

// vector unary minus x
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag> 
operator-(const simd_compl<Val, Bits, Simd_tag>& x);

// sum of all elements stored in the vector x
template<class Val, int Bits, class Simd_tag>
typename simd_compl<Val, Bits, Simd_tag>::value_type
horizontal_sum(const simd_compl<Val, Bits, Simd_tag>& x);

// evaluate x*y + z
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fma_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return x * y + z;
};

// evaluate x * y + z, where x is a real vector; see operator* for definition
// of x * y, where x is a real vector and y is a complex vector
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag>
fma_f(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return x * y + z;
};

// evaluate x * y + z, where y is a real vector; see operator* for definition
// of x * y, where y is a real vector and x is a complex vector
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fma_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return x * y + z;
};

// evaluate x * y - z
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fms_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return x * y - z;
};

// evaluate x * y - z, where x is a real vector; see operator* for definition
// of x * y, where x is a real vector and y is a complex vector
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fms_f(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return x * y - z;
};

// evaluate x * y - z, where y is a real vector; see operator* for definition
// of x * y, where y is a real vector and x is a complex vector
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fms_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return x * y - z;
};

// evaluate -x*y + z
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fnma_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return z - x * y;
};

// evaluate -x * y + z, where x is a real vector; see operator* for definition
// of x * y, where x is a real vector and y is a complex vector
template<class Val, int Bits, class Simd_tag>
simd_compl<Val, Bits, Simd_tag>
fnma_f(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return z - x * y;
};

// evaluate -x * y + z, where y is a real vector; see operator* for definition
// of x * y, where y is a real vector and x is a complex vector
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fnma_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return z - x * y;
};

// evaluate -x * y - z
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fnms_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return -(x * y + z);
};

// evaluate -x * y - z, where x is a real vector; see operator* for definition
// of x * y, where x is a real vector and y is a complex vector
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fnms_f(const simd<Val, Bits, Simd_tag>& x, const simd_compl<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return -(x * y + z);
};

// evaluate -x * y - z, where y is a real vector; see operator* for definition
// of x * y, where y is a real vector and x is a complex vector
template<class Val, int Bits, class Simd_tag>
force_inline
simd_compl<Val, Bits, Simd_tag>
fnms_f(const simd_compl<Val, Bits, Simd_tag>& x, const simd<Val, Bits, Simd_tag>& y, 
                       const simd_compl<Val, Bits, Simd_tag>& z)
{
    return -(x * y + z);
};

// return true if at least element in the vector x is +INF
template<class Val, int Bits, class Simd_tag>
bool
any_inf(const simd_compl<Val, Bits, Simd_tag>& x);

// return true if at least element in the vector x is NAN
template<class Val, int Bits, class Simd_tag>
bool
any_nan(const simd_compl<Val, Bits, Simd_tag>& x);

// print content of a vector to a stream
template<class Val, int Bits, class Simd_tag>
std::ostream& operator<<(std::ostream& os, const simd_compl<Val, Bits, Simd_tag>& x);

}}
