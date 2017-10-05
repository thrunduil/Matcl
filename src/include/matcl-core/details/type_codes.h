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

#include "matcl-core/memory/global_objects.h"
#include "matcl-core/details/val_struct_codes.h"

namespace matcl { namespace details
{

namespace utils
{
    using integer_scalar        = matcl::raw::Matrix<Integer,struct_scalar>;
    using real_scalar           = matcl::raw::Matrix<Real,struct_scalar>;
    using float_scalar          = matcl::raw::Matrix<Float,struct_scalar>;
    using complex_scalar        = matcl::raw::Matrix<Complex,struct_scalar>;
    using float_complex_scalar  = matcl::raw::Matrix<Float_complex,struct_scalar>;
    using object_scalar         = matcl::raw::Matrix<Object,struct_scalar>;
};

template<class T>
            struct type_to_code                             {};
template<>  struct type_to_code<raw::integer_dense>         {static const mat_code value = mat_code::integer_dense;};
template<>  struct type_to_code<raw::real_dense>            {static const mat_code value = mat_code::real_dense;};
template<>  struct type_to_code<raw::float_dense>           {static const mat_code value = mat_code::float_dense;};
template<>  struct type_to_code<raw::complex_dense>         {static const mat_code value = mat_code::complex_dense;};
template<>  struct type_to_code<raw::float_complex_dense>   {static const mat_code value = mat_code::float_complex_dense;};
template<>  struct type_to_code<raw::object_dense>          {static const mat_code value = mat_code::object_dense;};

template<>  struct type_to_code<raw::integer_sparse>        {static const mat_code value = mat_code::integer_sparse;};
template<>  struct type_to_code<raw::real_sparse>           {static const mat_code value = mat_code::real_sparse;};
template<>  struct type_to_code<raw::float_sparse>          {static const mat_code value = mat_code::float_sparse;};
template<>  struct type_to_code<raw::complex_sparse>        {static const mat_code value = mat_code::complex_sparse;};
template<>  struct type_to_code<raw::float_complex_sparse>  {static const mat_code value = mat_code::float_complex_sparse;};
template<>  struct type_to_code<raw::object_sparse>         {static const mat_code value = mat_code::object_sparse;};

template<>  struct type_to_code<raw::integer_band>          {static const mat_code value = mat_code::integer_band;};
template<>  struct type_to_code<raw::real_band>             {static const mat_code value = mat_code::real_band;};
template<>  struct type_to_code<raw::float_band>            {static const mat_code value = mat_code::float_band;};
template<>  struct type_to_code<raw::complex_band>          {static const mat_code value = mat_code::complex_band;};
template<>  struct type_to_code<raw::float_complex_band>    {static const mat_code value = mat_code::float_complex_band;};
template<>  struct type_to_code<raw::object_band>           {static const mat_code value = mat_code::object_band;};

template<>  struct type_to_code<utils::integer_scalar>      {static const mat_code value = mat_code::integer_scalar;};
template<>  struct type_to_code<utils::float_scalar>        {static const mat_code value = mat_code::float_scalar;};
template<>  struct type_to_code<utils::real_scalar>         {static const mat_code value = mat_code::real_scalar;};
template<>  struct type_to_code<utils::float_complex_scalar>{static const mat_code value = mat_code::float_complex_scalar;};
template<>  struct type_to_code<utils::complex_scalar>      {static const mat_code value = mat_code::complex_scalar;};
template<>  struct type_to_code<utils::object_scalar>       {static const mat_code value = mat_code::object_scalar;};

template<>  struct type_to_code<Integer>                    {static const mat_code value = mat_code::integer_scalar;};
template<>  struct type_to_code<Float>                      {static const mat_code value = mat_code::float_scalar;};
template<>  struct type_to_code<Real>                       {static const mat_code value = mat_code::real_scalar;};
template<>  struct type_to_code<Float_complex>              {static const mat_code value = mat_code::float_complex_scalar;};
template<>  struct type_to_code<Complex>                    {static const mat_code value = mat_code::complex_scalar;};
template<>  struct type_to_code<Object>                     {static const mat_code value = mat_code::object_scalar;};

template<matcl::mat_code v> 
            struct code_to_type                     {};
template<>  struct code_to_type<mat_code::integer_dense>  {using type = raw::integer_dense;};
template<>  struct code_to_type<mat_code::real_dense>     {using type = raw::real_dense;};
template<>  struct code_to_type<mat_code::float_dense>    {using type = raw::float_dense;};
template<>  struct code_to_type<mat_code::complex_dense>  {using type = raw::complex_dense;};
template<>  struct code_to_type<mat_code::float_complex_dense>  
                                                          {using type = raw::float_complex_dense;};
template<>  struct code_to_type<mat_code::object_dense>   {using type = raw::object_dense;};

template<>  struct code_to_type<mat_code::integer_sparse> {using type = raw::integer_sparse;};
template<>  struct code_to_type<mat_code::real_sparse>    {using type = raw::real_sparse;};
template<>  struct code_to_type<mat_code::float_sparse>   {using type = raw::float_sparse;};
template<>  struct code_to_type<mat_code::complex_sparse> {using type = raw::complex_sparse;};
template<>  struct code_to_type<mat_code::float_complex_sparse> 
                                                          {using type = raw::float_complex_sparse;};
template<>  struct code_to_type<mat_code::object_sparse>  {using type = raw::object_sparse;};

template<>  struct code_to_type<mat_code::integer_band>   {using type = raw::integer_band;};
template<>  struct code_to_type<mat_code::real_band>      {using type = raw::real_band;};
template<>  struct code_to_type<mat_code::float_band>     {using type = raw::float_band;};
template<>  struct code_to_type<mat_code::complex_band>   {using type = raw::complex_band;};
template<>  struct code_to_type<mat_code::float_complex_band>   
                                                          {using type = raw::float_complex_band;};
template<>  struct code_to_type<mat_code::object_band>    {using type = raw::object_band;};

template<>  struct code_to_type<mat_code::integer_scalar> {using type = Integer;};
template<>  struct code_to_type<mat_code::float_scalar>   {using type = Float;};
template<>  struct code_to_type<mat_code::real_scalar>    {using type = Real;};
template<>  struct code_to_type<mat_code::float_complex_scalar> 
                                                          {using type = Float_complex;};
template<>  struct code_to_type<mat_code::complex_scalar> {using type = Complex;};
template<>  struct code_to_type<mat_code::object_scalar>  {using type = Object;};

};};

#define macro_first_matrix_type_code (int)mat_code::integer_dense
#define macro_last_matrix_type_code  (int)mat_code::object_band

#define macro_first_scalar_type_code (int)mat_code::integer_scalar
#define macro_last_scalar_type_code  (int)mat_code::object_scalar

#define MACRO_FOREACH_CODE(macro,arg1,arg2)                                                 \
    macro(::matcl::details::code_to_type<matcl::mat_code::integer_dense>::type,arg1,arg2)   \
    macro(::matcl::details::code_to_type<matcl::mat_code::real_dense>::type,arg1,arg2)      \
    macro(::matcl::details::code_to_type<matcl::mat_code::complex_dense>::type,arg1,arg2)   \
    macro(::matcl::details::code_to_type<matcl::mat_code::object_dense>::type,arg1,arg2)    \
    macro(::matcl::details::code_to_type<matcl::mat_code::integer_sparse>::type,arg1,arg2)  \
    macro(::matcl::details::code_to_type<matcl::mat_code::real_sparse>::type,arg1,arg2)     \
    macro(::matcl::details::code_to_type<matcl::mat_code::complex_sparse>::type,arg1,arg2)  \
    macro(::matcl::details::code_to_type<matcl::mat_code::object_sparse>::type,arg1,arg2)   \
    macro(::matcl::details::code_to_type<matcl::mat_code::integer_band>::type,arg1,arg2)    \
    macro(::matcl::details::code_to_type<matcl::mat_code::real_band>::type,arg1,arg2)       \
    macro(::matcl::details::code_to_type<matcl::mat_code::complex_band>::type,arg1,arg2)    \
    macro(::matcl::details::code_to_type<matcl::mat_code::object_band>::type,arg1,arg2)

#define MACRO_FOREACH_CODE_ALL(macro,arg1,arg2)                                             \
    macro(::matcl::details::code_to_type<matcl::mat_code::integer_dense>::type,arg1,arg2)   \
    macro(::matcl::details::code_to_type<matcl::mat_code::float_dense>::type,arg1,arg2)     \
    macro(::matcl::details::code_to_type<matcl::mat_code::real_dense>::type,arg1,arg2)      \
    macro(::matcl::details::code_to_type<matcl::mat_code::float_complex_dense>::type,arg1,arg2)\
    macro(::matcl::details::code_to_type<matcl::mat_code::complex_dense>::type,arg1,arg2)   \
    macro(::matcl::details::code_to_type<matcl::mat_code::object_dense>::type,arg1,arg2)    \
    macro(::matcl::details::code_to_type<matcl::mat_code::integer_sparse>::type,arg1,arg2)  \
    macro(::matcl::details::code_to_type<matcl::mat_code::float_sparse>::type,arg1,arg2)     \
    macro(::matcl::details::code_to_type<matcl::mat_code::real_sparse>::type,arg1,arg2)     \
    macro(::matcl::details::code_to_type<matcl::mat_code::float_complex_sparse>::type,arg1,arg2)\
    macro(::matcl::details::code_to_type<matcl::mat_code::complex_sparse>::type,arg1,arg2)  \
    macro(::matcl::details::code_to_type<matcl::mat_code::object_sparse>::type,arg1,arg2)   \
    macro(::matcl::details::code_to_type<matcl::mat_code::integer_band>::type,arg1,arg2)    \
    macro(::matcl::details::code_to_type<matcl::mat_code::float_band>::type,arg1,arg2)      \
    macro(::matcl::details::code_to_type<matcl::mat_code::real_band>::type,arg1,arg2)       \
    macro(::matcl::details::code_to_type<matcl::mat_code::float_complex_band>::type,arg1,arg2)\
    macro(::matcl::details::code_to_type<matcl::mat_code::complex_band>::type,arg1,arg2)    \
    macro(::matcl::details::code_to_type<matcl::mat_code::object_band>::type,arg1,arg2)
