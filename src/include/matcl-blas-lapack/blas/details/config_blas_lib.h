/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#pragma once

#include "matcl-blas-lapack/blas_loader/clapack_types.h"
#include "matcl-blas-lapack/lapack/details/clapack_lapack_declarations.h"
#include "matcl-blas-lapack/blas_loader/blas_loader.h"

namespace matcl { namespace lapack
{

using s_type_bl = s_type_wr;
using d_type_bl = d_type_wr;
using c_type_bl = c_type_wr;
using z_type_bl = z_type_wr;
using i_type_bl = i_type_wr;
using sel_fun_bl= sel_fun_wr;
using l_type_bl = l_type_wr;

inline s_type_bl* _c(s_type* x)                 { return x; };
inline d_type_bl* _c(d_type* x)                 { return x; };
inline c_type_bl* _c(c_type* x)                 { return reinterpret_cast<c_type_bl*>(x); };
inline z_type_bl* _c(z_type* x)                 { return reinterpret_cast<z_type_bl*>(x); };
inline i_type_bl* _c(i_type* x)                 { return reinterpret_cast<i_type_bl*>(x); };
inline l_type_bl* _cl(l_type* x)                { return reinterpret_cast<l_type_bl*>(x); };
inline sel_fun_bl _c(sel_fun x)                 { return reinterpret_cast<sel_fun_bl>(x); };
inline char*      _c(char* x)                   { return x; };

inline const s_type_bl* _c(const s_type* x)     { return x; };
inline const d_type_bl* _c(const d_type* x)     { return x; };
inline const c_type_bl* _c(const c_type* x)     { return reinterpret_cast<const c_type_bl*>(x); };
inline const z_type_bl* _c(const z_type* x)     { return reinterpret_cast<const z_type_bl*>(x); };
inline const i_type_bl* _c(const i_type* x)     { return reinterpret_cast<const i_type_bl*>(x); };
inline const l_type_bl* _cl(const l_type* x)    { return reinterpret_cast<const l_type_bl*>(x); };
inline const char*      _c(const char* x)       { return x; };

inline s_type_bl* _rc(s_type* x)                { return x; };
inline d_type_bl* _rc(d_type* x)                { return x; };
inline c_type_bl* _rc(c_type* x)                { return reinterpret_cast<c_type_bl*>(x); };
inline z_type_bl* _rc(z_type* x)                { return reinterpret_cast<z_type_bl*>(x); };
inline i_type_bl* _rc(i_type* x)                { return reinterpret_cast<i_type_bl*>(x); };
inline l_type_bl* _rcl(l_type* x)               { return reinterpret_cast<l_type_bl*>(x); };
inline sel_fun_bl _rc(sel_fun x)                { return reinterpret_cast<sel_fun_bl>(x); };
inline char*       _rc(char* x)                 { return x; };

inline s_type_bl* _rc(const s_type* x)          { return _c(const_cast<s_type*>(x)); };
inline d_type_bl* _rc(const d_type* x)          { return _c(const_cast<d_type*>(x)); };
inline c_type_bl* _rc(const c_type* x)          { return _c(const_cast<c_type*>(x)); };
inline z_type_bl* _rc(const z_type* x)          { return _c(const_cast<z_type*>(x)); };
inline i_type_bl* _rc(const i_type* x)          { return _c(const_cast<i_type*>(x)); };
inline l_type_bl* _rcl(const l_type* x)         { return _cl(const_cast<l_type*>(x)); };
inline char*      _rc(const char* x)            { return _c(const_cast<char*>(x)); };

#define BLAS_NAME(x) raw_blas_lapack::x##_
#define LAPACK_NAME(x) raw_blas_lapack::x##_

};};
