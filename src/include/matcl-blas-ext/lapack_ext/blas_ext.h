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

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-ext/config_blas_ext.h"

namespace matcl { namespace lapack
{

// additional functions with blas-like interface

// BLAS Level1
//--------------------------------------------------------------------
//                      SET
//--------------------------------------------------------------------

// Compute Y(1:n) = alpha
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
set(i_type n, V alpha, V* y, i_type incy);

//--------------------------------------------------------------------
//                      YAC
//--------------------------------------------------------------------
// Compute Y = alpha * X
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
yax(i_type n, V alpha, const V* x, i_type incx, V* y, i_type incy);

#if 0

//=======================   YAX_CONJ    =====================================================

// Compute Y = alpha * conj(X)
template<class V> BLAS_EXT_EXPORT
typename details::enable_if_valid<void,V>::type
yax_conj(i_type n, V alpha, const V* x,i_type incx, V* y, i_type incy);

BLAS_EXT_EXPORT void syax_conj(const i_type* n, const s_type* alpha,const s_type* x,const i_type* incx, 
                             s_type* y, const i_type* incy );
BLAS_EXT_EXPORT void dyax_conj(const i_type* n, const d_type* alpha,const d_type* x,const i_type* incx, 
                             d_type* y, const i_type* incy );
BLAS_EXT_EXPORT void cyax_conj(const i_type* n, const c_type* alpha,const c_type* x,const i_type* incx, 
                             c_type* y, const i_type* incy );
BLAS_EXT_EXPORT void zyax_conj(const i_type* n, const z_type* alpha,const z_type* x,const i_type* incx, 
                             z_type* y, const i_type* incy );

#endif

};};
