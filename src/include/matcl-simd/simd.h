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

#include "matcl-simd/arch/reg_128/simd_128.h"
#include "matcl-simd/arch/reg_128/simd_128_compl.h"

#include "matcl-simd/arch/reg_256/simd_256.h"
#include "matcl-simd/arch/reg_256/simd_256_compl.h"

#include "matcl-simd/simd/simd_scalar_func.h"
#include "matcl-simd/simd/simd_128_func.h"
#include "matcl-simd/simd/simd_256_func.h"
#include "matcl-simd/simd/simd_128_compl_func.h"
#include "matcl-simd/simd/simd_256_compl_func.h"

#include "matcl-simd/impl/simd_128_compl_func.inl"
#include "matcl-simd/impl/simd_256_func.inl"
#include "matcl-simd/impl/simd_256_compl_func.inl"
#include "matcl-simd/impl/simd_scal_func.inl"

#include "matcl-simd/arch/reg_128/simd_128.inl"
#include "matcl-simd/arch/reg_128/simd_128_compl.inl"
#include "matcl-simd/arch/reg_256/simd_256.inl"
#include "matcl-simd/arch/reg_256/simd_256_compl.inl"

#include "matcl-simd/simd_utils.h"