/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2017 - 2018
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
#include <complex>

#ifdef BLAS_EXT_EXPORTS
    #define BLAS_EXT_EXPORT  __declspec(dllexport)
#else
    #define BLAS_EXT_EXPORT  __declspec(dllimport)
#endif

#ifdef __unix__
 	#undef BLAS_EXT_EXPORT
    #define BLAS_EXT_EXPORT
#endif

namespace matcl { namespace lapack
{

using s_type    = float;
using d_type    = double;
using c_type    = std::complex<float>;
using z_type    = std::complex<double>;
using i_type    = int;
using l_type    = char;

using sel_fun   = char (*)(...);

};};