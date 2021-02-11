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

#include <math.h>

namespace matcl
{

// code of scalar type stored in matrices
enum class value_code : int
{ 
    v_integer = 1,  v_float = 2,    v_real = 3, v_float_complex = 4, 
    v_complex = 5,  v_object = 6
};

// code of internal representation of a matrix
enum class struct_code : int
{ 
    struct_dense = 0, struct_sparse, struct_banded, struct_scalar
};

// file open modes
enum class open_mode : int
{
    readonly, readwrite, readwrite_create
};

// thread modes
enum class thread_mode : int
{
    single_thread, multi_thread, serialized
};

// type of transposition
enum class trans_type : int
{
    no_trans, trans, conj_trans
};

// type of transposition or conjugate
enum class trans_type_ext
{
    no_trans, trans, conj_trans, conj
};

// matrix representation tags
struct struct_dense{};
struct struct_sparse{};
struct struct_banded{};
struct struct_scalar{};

// determine how scalars, dense or sparse matrices are printed
enum class disp_mode : int
{
    default,        // take from global settings
    standard,       // display mathod determinded by matrix struct_type
    all_dense,      // all matrices and scalars are displayed as dense matrix
    matrix_dense,   // all matrices are displayed as dense matrix
    scalar_dense,   // all scalars are displayed as dense matrix
    size            // only for counting number of options
};

// determine how strings are aligned
enum class align_type : int
{
    left,           // justify text to left
    center,         // justify text to center
    right           // justify text to right
};

// code, that describe concrete representation of matrix given by type of 
// scalar stored in matrix and matrix representation.
enum class mat_code : int
{
    integer_scalar, float_scalar,   real_scalar,    float_complex_scalar,   
    complex_scalar, object_scalar,

    integer_dense,  float_dense,    real_dense,     float_complex_dense,
    complex_dense,  object_dense,

    integer_sparse, float_sparse,   real_sparse,    float_complex_sparse,   
    complex_sparse, object_sparse,

    integer_band,   float_band,     real_band,      float_complex_band,
    complex_band,   object_band,

    last_code
};

// classify floating point number
enum class fp_type : int
{
    fp_infinite     = FP_INFINITE,  // positive or negative infinity (overflow)
    fp_nan          = FP_NAN,       // not-a-number
    fp_zero         = FP_ZERO,      // value of zero
    fp_subnormal    = FP_SUBNORMAL, // sub-normal value (underflow)
    fp_normal       = FP_NORMAL,    // normal value (none of the above)
    fp_unknown      = 5,            // unable to determine (function returned unrecognized
                                    // value, or object does not represent a floating point value)
};

};
