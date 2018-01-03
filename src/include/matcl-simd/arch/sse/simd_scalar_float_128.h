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

#include "matcl-simd/config.h"
#include "matcl-simd/simd_general.h"
#include "matcl-simd/default_simd.h"

#include <emmintrin.h>

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE SCALAR FLOAT
//-------------------------------------------------------------------

// vector of one single precision scalar
template<>
class alignas(16) simd<float, 128, scalar_sse_tag>
{
    public:
        // implementation type
        using impl_type     = __m128;

        // type of stored elements
        using value_type    = float;

        // simd tag
        using simd_tag      = scalar_sse_tag;

        // number of bits
        static const int
        number_bits         = 128;

        // type of vector storing half of elements
        using simd_half     = simd<float, 128, scalar_sse_tag>;

        // simd type of the same size storing float values
        using simd_float    = simd<float, 128, scalar_sse_tag>;

        // simd type of the same size storing double values
        using simd_double   = simd<double, 128, scalar_sse_tag>;

        // simd type of the same size storing int32_t values
        using simd_int32    = simd<int32_t, 128, scalar_sse_tag>;

        // simd type of the same size storing int64_t values
        using simd_int64    = simd<int64_t, 128, scalar_sse_tag>;

        // simd type of the same kind storing vector of float values
        using simd_vector    = simd<float, 128, sse_tag>;

    public:
        // number of elements in the vector
        static const int 
        vector_size         = 1;

    public:
        // internal representation
        impl_type           data;

    public:
        // construct uninitialized vector
        simd() = default;

        // construct vector with all elements equal to val
        explicit simd(float val);

        // construct from representation
        simd(const impl_type& v);

        // conversion between simd types
        explicit simd(const simd<float, 128, sse_tag>& s);

        // conversion between simd types
        explicit simd(const simd<float, 128, nosimd_tag>& s);        

        // conversion form simd scalar; set all elements to s.first()
        explicit simd(const simd<float, 128, scalar_nosimd_tag>& s);

        // copy constructor
        simd(const simd<float, 128, scalar_sse_tag>& s) = default;

    public:
        // connstruct vector with all elements set to 0.0
        static simd     zero();

        // connstruct vector with all elements set to -0.0
        static simd     minus_zero();

        // connstruct vector with all elements set to 1.0
        static simd     one();

        // connstruct vector with all elements set to -1.0
        static simd     minus_one();

    public:
        // construct vector with all elements equal to arr[0]
        static simd     broadcast(const float* arr);

        // construct vector with all elements equal to arr
        static simd     broadcast(const float& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd     load(const float* arr, std::true_type aligned);
        static simd     load(const float* arr, std::false_type not_aligned = std::false_type());

        // gather singe-precision (32-bit) floating-point elements from memory using 
        // 32-bit indices, i.e. i-th element of resulting vector is arr[ind[i]]
        static simd     gather(const float* arr, const simd_int32& ind);

        // gather two singe-precision (32-bit) floating-point elements from memory using 
        // 64-bit indices, i.e. i-th element of resulting vector is arr[ind[i]];
        static simd     gather(const float* arr, const simd_int64& ind);

        // set the first element in the vector to v and set 0 to all other elements
        static simd     set_lower(float v);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(float* arr, std::true_type aligned) const;
        void            store(float* arr, std::false_type not_aligned = std::false_type()) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(float* arr) const;

        // get i-th element from the vector; pos is 0-based
        float           get(int pos) const;

        // return the first element in the vector; equivalent to get(0), 
        // but possibly faster
        float           first() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, float val);

        // return pointer to the first element in the vector
        const float*    get_raw_ptr() const;
        float*          get_raw_ptr();

    public:
        // cast all elements to double
        simd_double     convert_to_double() const;

        // convert elements to int32_t, rounding is performed according
        // to current rounding mode (usually round to nearest ties to even)
        simd_int32      convert_to_int32() const;

        // convert elements to int64_t, rounding is performed according
        // to current rounding mode (usually round to nearest ties to even)
        simd_int64      convert_to_int64() const;

        // reinterpret cast to vector of double of the same kind
        simd_double     reinterpret_as_double() const;

        // reinterpret cast to vector of int32 of the same kind
        simd_int32      reinterpret_as_int32() const;

        // reinterpret cast to vector of int64 of the same kind
        simd_int64      reinterpret_as_int64() const;

        // reinterpret cast to simd of the same type storing vectors
        simd_vector     as_vector() const;

    public:
        // plus assign operator
        simd&           operator+=(const simd& x);

        // minus assign operator
        simd&           operator-=(const simd& x);

        // multiply assign operator
        simd&           operator*=(const simd& x);

        // divide assign operator
        simd&           operator/=(const simd& x);
};

}}
