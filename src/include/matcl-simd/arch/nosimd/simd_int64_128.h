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

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          GENERIC INT64_T
//-------------------------------------------------------------------

// vector of two int64_t scalars
template<>
class alignas(16) simd<int64_t, 128, nosimd_tag>
{
    public:
        // implementation type
        using impl_type     = int64_t[2];

        // type of stored elements
        using value_type    = int64_t;

    public:
        // number of elements in the vector
        static const int 
        vector_size         = sizeof(impl_type) / sizeof(value_type);    

        // simd type storing 32-bit integers of size 128 bits
        using simd_128_int32 = simd<int32_t, 128, nosimd_tag>;

        // simd type storing 64-bit integers of size 128 bits
        using simd_128_int64 = simd<int64_t, 128, nosimd_tag>;

        // simd type storing float values of size 128 bits
        using simd_128_float  = simd<float, 128, nosimd_tag>;

        // simd type storing double values of size 128 bits
        using simd_128_double = simd<double, 128, nosimd_tag>;

    public:
        // internal representation
        impl_type           data;

    public:
        // construct uninitialized vector
        simd() = default;

        // construct vector with all elements equal to val
        explicit simd(int32_t val);

        // construct vector with all elements equal to val
        explicit simd(int64_t val);

        // set the first element in the vector to v1 and the second element to v2
        simd(int64_t v1, int64_t v2);

        // construct vector with first element copied from lo
        // and last element copied from hi; only lower part of lo and hi is used
        simd(const simd& lo, const simd& hi);

        // construct from representation
        simd(const impl_type& v);

        // conversion between simd types
        explicit simd(const simd<int64_t, 128, sse_tag>& s);

        // copy constructor
        simd(const simd<int64_t, 128, nosimd_tag>& s) = default;

    public:
        // connstruct vector with all elements set to 0
        static simd     zero();

        // connstruct vector with all elements set to 1
        static simd     one();

        // connstruct vector with all elements set to -1
        static simd     minus_one();

    public:
        // construct vector with all elements equal to arr[0]
        static simd     broadcast(const int64_t* arr);

        // construct vector with all elements equal to arr
        static simd     broadcast(const int64_t& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd     load(const int64_t* arr, std::true_type aligned);
        static simd     load(const int64_t* arr, std::false_type not_aligned = std::false_type());

        // gather 64-bit integer elements from memory using 32-bit indices, 
        // i.e. i-th element of resulting vector is arr[ind[i]]
        static simd     gather(const int64_t* arr, const simd_128_int32& ind);

        // gather 64-bit integer elements from memory using 64-bit indices, 
        // i.e. i-th element of resulting vector is arr[ind[i]]
        static simd     gather(const int64_t* arr, const simd_128_int64& ind);

        // set the first element in the vector to v and set 0 to all other elements
        static simd     set_lower(int64_t v);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(int64_t* arr, std::true_type aligned) const;
        void            store(int64_t* arr, std::false_type not_aligned = std::false_type()) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(int64_t* arr) const;

        // get i-th element from the vector; pos is 0-based
        int64_t         get(int pos) const;

        // get i-th element from the vector; Pos is 0-based
        template<int Pos>
        int64_t         get() const;

        // return the first element in the vector; equivalent to get(0), 
        // but possibly faster
        int64_t         first() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, int64_t val);

        // return pointer to the first element in the vector
        const int64_t*  get_raw_ptr() const;
        int64_t*        get_raw_ptr();

        // set i-th element of the vector; Pos is 0-based
        template<int Pos>
        void            set(int64_t val);

        // return simd storing first element
        simd            extract_low() const;

        // return simd storing last element
        simd            extract_high() const;

    public:
        // convert elements to int32_t and store the result in the lower part
        simd_128_int32  convert_to_int32() const;

        // reinterpret cast to vector of double of the same kind
        simd_128_double reinterpret_as_double() const;

        // reinterpret cast to vector of float of the same kind
        simd_128_float  reinterpret_as_float() const;

        // reinterpret cast to vector of int32 of the same kind
        simd_128_int32  reinterpret_as_int32() const;

    public:
        // plus assign operator
        simd&           operator+=(const simd& x);

        // minus assign operator
        simd&           operator-=(const simd& x);

        // multiply assign operator
        simd&           operator*=(const simd& x);
};

}}
