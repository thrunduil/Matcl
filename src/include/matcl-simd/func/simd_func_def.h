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

#include "matcl-simd/simd_general.h"

namespace matcl { namespace simd
{

namespace ms = matcl::simd;

template<class Val, int Bits, class Simd_tag>
struct simd_reverse{};

template<class Val, int Bits, class Simd_tag>
struct simd_mult{};

template<class Val, int Bits, class Simd_tag>
struct simd_div{};

template<class Val, int Bits, class Simd_tag>
struct simd_plus{};

template<class Val, int Bits, class Simd_tag>
struct simd_minus{};

template<class Val, int Bits, class Simd_tag>
struct simd_uminus{};

template<class Val, int Bits, class Simd_tag>
struct simd_abs{};

template<class Val, int Bits, class Simd_tag>
struct simd_sub_add{};

template<class Val, int Bits, class Simd_tag>
struct simd_fma{};

template<class Val, int Bits, class Simd_tag>
struct simd_fms{};

template<class Val, int Bits, class Simd_tag>
struct simd_sum_all{};

template<class Val, int Bits, class Simd_tag>
struct simd_max{};

template<class Val, int Bits, class Simd_tag>
struct simd_min{};

template<class Val, int Bits, class Simd_tag>
struct simd_eeq{};

template<class Val, int Bits, class Simd_tag>
struct simd_neq{};

template<class Val, int Bits, class Simd_tag>
struct simd_lt{};

template<class Val, int Bits, class Simd_tag>
struct simd_gt{};

template<class Val, int Bits, class Simd_tag>
struct simd_leq{};

template<class Val, int Bits, class Simd_tag>
struct simd_geq{};

template<class Val, int Bits, class Simd_tag>
struct simd_sqrt{};

template<class Val, int Bits, class Simd_tag>
struct simd_round{};

template<class Val, int Bits, class Simd_tag>
struct simd_floor{};

template<class Val, int Bits, class Simd_tag>
struct simd_ceil{};

template<class Val, int Bits, class Simd_tag>
struct simd_trunc{};

template<class Val, int Bits, class Simd_tag>
struct simd_any_nan{};

template<class Val, int Bits, class Simd_tag>
struct simd_any{};

template<class Val, int Bits, class Simd_tag>
struct simd_all{};

}}
