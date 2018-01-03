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
//                          GENERIC SCALAR int32_t
//-------------------------------------------------------------------

// vector of one int32_t scalar
template<>
class alignas(int32_t) simd<int32_t, 128, scalar_nosimd_tag>
{
    public:
        // implementation type
        using impl_type     = int32_t;

        // type of stored elements
        using value_type    = int32_t;

        // simd tag
        using simd_tag      = scalar_nosimd_tag;

        // number of bits
        static const int
        number_bits         = 128;

        // type of vector storing half of elements
        using simd_half     = simd<int32_t, 128, scalar_nosimd_tag>;

        // simd type of the same size storing float values
        using simd_float    = simd<float, 128, scalar_nosimd_tag>;

        // simd type of the same size storing double values
        using simd_double   = simd<double, 128, scalar_nosimd_tag>;

        // simd type of the same size storing int32_t values
        using simd_int32    = simd<int32_t, 128, scalar_nosimd_tag>;

        // simd type of the same size storing int64_t values
        using simd_int64    = simd<int64_t, 128, scalar_nosimd_tag>;

        // simd type of the same kind storing vector of int32_t integers
        using simd_vector    = simd<int32_t, 128, nosimd_tag>;

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
        explicit simd(int32_t val);

      #if MATCL_ARCHITECTURE_HAS_SSE2
        // conversion between simd types
        explicit simd(const simd<int32_t, 128, sse_tag>& s);
      #endif

      #if MATCL_ARCHITECTURE_HAS_SSE2
        // conversion form simd scalar; set all elements to s.first()
        explicit simd(const simd<int32_t, 128, scalar_sse_tag>& s);
      #endif

        // conversion between simd types
        explicit simd(const simd<int32_t, 128, nosimd_tag>& s);        

        // copy constructor
        simd(const simd<int32_t, 128, scalar_nosimd_tag>& s) = default;

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
        static simd     gather(const int32_t* arr, const simd_int32& ind);

        // gather 32-bit integer elements from memory using 64-bit indices, 
        // i.e. i-th element of resulting vector is arr[ind[i]];
        // last two elements are set to zero
        static simd     gather(const int32_t* arr, const simd_int64& ind);

        // set the first element in the vector to v and set 0 to all other elements
        static simd     set_lower(int32_t v);

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

    public:
        // convert all elements to int64_t
        simd_int64      convert_to_int64() const;

        // convert elements to float
        simd_float      convert_to_float() const;

        // convert all elements to double
        simd_double     convert_to_double() const;

        // reinterpret cast to vector of double of the same kind
        simd_double     reinterpret_as_double() const;

        // reinterpret cast to vector of float of the same kind
        simd_float      reinterpret_as_float() const;

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
};

}}
