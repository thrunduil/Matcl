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

#include "blas_ext.inl"

namespace matcl { namespace lapack
{

//--------------------------------------------------------------------
//                      SET
//--------------------------------------------------------------------

template BLAS_EXT_EXPORT void
lapack::set<d_type>(i_type n, d_type alpha, d_type* y, i_type incy);

template BLAS_EXT_EXPORT void
lapack::set<s_type>(i_type n, s_type alpha, s_type* y, i_type incy);

template BLAS_EXT_EXPORT void
lapack::set<z_type>(i_type n, z_type alpha, z_type* y, i_type incy);

template BLAS_EXT_EXPORT void
lapack::set<c_type>(i_type n, c_type alpha, c_type* y, i_type incy);

//--------------------------------------------------------------------
//                      YAX
//--------------------------------------------------------------------
template BLAS_EXT_EXPORT void
lapack::yax<d_type>(i_type n, d_type alpha, const d_type* x, i_type incx, d_type* y, i_type incy);

template BLAS_EXT_EXPORT void
lapack::yax<s_type>(i_type n, s_type alpha, const s_type* x, i_type incx, s_type* y, i_type incy);

template BLAS_EXT_EXPORT void
lapack::yax<z_type>(i_type n, z_type alpha, const z_type* x, i_type incx, z_type* y, i_type incy);

template BLAS_EXT_EXPORT void
lapack::yax<c_type>(i_type n, c_type alpha, const c_type* x, i_type incx, c_type* y, i_type incy);

};};
