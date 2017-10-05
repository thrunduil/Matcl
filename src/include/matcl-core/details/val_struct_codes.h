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

#include "matcl-core/config.h"
#include "matcl-core/general/fwd_decls.h"
#include "matcl-core/memory/global_objects.h"
#include "matcl-core/matrix/enums.h"
#include "matcl-core/matrix/scalar_types.h"

namespace matcl { namespace details
{

template<class T>   
            struct value_to_code                    {};
template<>  struct value_to_code<Integer>           {static const value_code value = value_code::v_integer;};
template<>  struct value_to_code<Float>             {static const value_code value = value_code::v_float;};
template<>  struct value_to_code<Real>              {static const value_code value = value_code::v_real;};
template<>  struct value_to_code<Float_complex>     {static const value_code value = value_code::v_float_complex;};
template<>  struct value_to_code<Complex>           {static const value_code value = value_code::v_complex;};
template<>  struct value_to_code<Object>            {static const value_code value = value_code::v_object;};

template<class T> 
            struct struct_to_code                   {};
template<>  struct struct_to_code<struct_dense>     {static const struct_code value = struct_code::struct_dense;};
template<>  struct struct_to_code<struct_sparse>    {static const struct_code value = struct_code::struct_sparse;};
template<>  struct struct_to_code<struct_banded>    {static const struct_code value = struct_code::struct_banded;};
template<>  struct struct_to_code<struct_scalar>    {static const struct_code value = struct_code::struct_scalar;};

template<matcl::value_code v>    
            struct code_to_value                                {};
template<>  struct code_to_value<value_code::v_integer>         {using type = Integer;};
template<>  struct code_to_value<value_code::v_float>           {using type = Float;};
template<>  struct code_to_value<value_code::v_real>            {using type = Real;};
template<>  struct code_to_value<value_code::v_float_complex>   {using type = Float_complex;};
template<>  struct code_to_value<value_code::v_complex>         {using type = Complex;};
template<>  struct code_to_value<value_code::v_object>          {using type = Object;};

template<matcl::struct_code v> 
            struct code_to_struct                               {};
template<>  struct code_to_struct<struct_code::struct_dense>    {using type = struct_dense;};
template<>  struct code_to_struct<struct_code::struct_sparse>   {using type = struct_sparse;};
template<>  struct code_to_struct<struct_code::struct_banded>   {using type = struct_banded;};
template<>  struct code_to_struct<struct_code::struct_scalar>   {using type = struct_scalar;};

};};