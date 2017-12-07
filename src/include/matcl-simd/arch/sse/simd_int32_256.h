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
#include "matcl-simd/arch/sse/simd_int32_128.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE INT32_T
//-------------------------------------------------------------------

// vector of eight int32_t scalars
template<>
class alignas(32) simd<int32_t, 256, sse_tag>
{
    public:
        // type of vector storing half of elements
        using simd_half     = simd<int32_t, 128, sse_tag>;

        // implementation type
        using impl_type     = simd_half[2];

        // type of stored elements
        using value_type    = int32_t;

        // simd type storing int64_t values of size 256 bits
        using simd_256_int64 = simd<int64_t, 256, sse_tag>;

        // simd type storing 32-bit integers of size 256 bits
        using simd_256_int32 = simd<int32_t, 256, sse_tag>;

        // simd type storing float values of size 256 bits
        using simd_256_float = simd<float, 256, sse_tag>;

        // simd type storing double values of size 256 bits
        using simd_256_double = simd<double, 256, sse_tag>;

    public:
        // number of elements in the vector
        static const int 
        vector_size         = sizeof(impl_type) / sizeof(value_type);    

    public:
        // internal representation
        impl_type           data;

    public:
        // construct uninitialized vector
        simd() = default;

        // construct vector with all elements equal to val
        explicit simd(int32_t val);

        // construct vector with first four elements and last four elements
        // copied from lo_hi
        explicit simd(const simd_half& lo_hi);

        // construct vector with i-th element set to vi
        simd(int32_t v0, int32_t v1, int32_t v2, int32_t v3, int32_t v4, 
             int32_t v5, int32_t v6, int32_t v7);

        // construct from representation
        simd(const impl_type& v);

        // construct vector with first four elements copied from lo
        // and last four elements copied from hi
        simd(const simd_half& lo, const simd_half& hi);

        // conversion between simd types
        explicit simd(const simd<int32_t, 256, nosimd_tag>& s);
        explicit simd(const simd<int32_t, 256, avx_tag>& s);

        // copy constructor
        simd(const simd<int32_t, 256, sse_tag>& s) = default;

    public:
        // connstruct vector with all elements set to 0
        static simd     zero();

        // connstruct vector with all elements set to 1
        static simd     one();

        // connstruct vector with all elements set to -1
        static simd     minus_one();

    public:
        // construct vector with all elements equal to arr[0]
        static simd     broadcast(const int32_t* arr);

        // construct vector with all elements equal to arr
        static simd     broadcast(const int32_t& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd     load(const int32_t* arr, std::true_type aligned);
        static simd     load(const int32_t* arr, std::false_type not_aligned = std::false_type());

        // gather 32-bit integer elements from memory using 32-bit indices, 
        // i.e. i-th element of resulting vector is arr[ind[i]]
        static simd     gather(const int32_t* arr, const simd_256_int32& ind);

        // gather 32-bit floating-point elements from memory using 64-bit indices, 
        // i.e. i-th element of resulting vector is arr[ind[i]];
        // last four elements are the same as the first four elements
        static simd     gather(const int32_t* arr, const simd_256_int64& ind);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(int32_t* arr, std::true_type aligned) const;
        void            store(int32_t* arr, std::false_type not_aligned = std::false_type()) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(int32_t* arr) const;

        // get i-th element from the vector; pos is 0-based
        int32_t         get(int pos) const;

        // return the first element in the vector; equivalent to get(0), 
        // but possibly faster
        int32_t         first() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, int32_t val);

        // return pointer to the first element in the vector
        const int32_t*  get_raw_ptr() const;
        int32_t*        get_raw_ptr();

        // return simd storing first four elements
        simd_half       extract_low() const;

        // return simd storing last four elements
        simd_half       extract_high() const;

    public:
        // convert the first four elements to int64_t
        simd_256_int64  convert_low_to_int64() const;

        // convert the last four elements to int64_t
        simd_256_int64  convert_high_to_int64() const;

        // convert elements to float
        simd_256_float  convert_to_float() const;

        // convert the first four elements to double
        simd_256_double convert_low_to_double() const;

        // convert the last four elements to double
        simd_256_double convert_high_to_double() const;

        // reinterpret cast to vector of double of the same kind
        simd_256_double reinterpret_as_double() const;

        // reinterpret cast to vector of float of the same kind
        simd_256_float  reinterpret_as_float() const;

        // reinterpret cast to vector of int64 of the same kind
        simd_256_int64  reinterpret_as_int64() const;

    public:
        // plus assign operator
        simd&           operator+=(const simd& x);

        // minus assign operator
        simd&           operator-=(const simd& x);

        // multiply assign operator
        simd&           operator*=(const simd& x);
};

}}
