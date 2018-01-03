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
//                          GENERIC SCALAR DOUBLE
//-------------------------------------------------------------------

// vector of one double precision scalar
template<>
class alignas(double) simd<double, 128, scalar_nosimd_tag>
{
    public:
        // implementation type
        using impl_type     = double;

        // type of stored elements
        using value_type    = double;

        // simd tag
        using simd_tag      = scalar_nosimd_tag;

        // number of bits
        static const int
        number_bits         = 128;

        // type of vector storing half of elements
        using simd_half     = simd<double, 128, scalar_nosimd_tag>;

        // simd type of the same size storing float values
        using simd_float    = simd<float, 128, scalar_nosimd_tag>;

        // simd type of the same size storing double values
        using simd_double   = simd<double, 128, scalar_nosimd_tag>;

        // simd type of the same size storing int32_t values
        using simd_int32    = simd<int32_t, 128, scalar_nosimd_tag>;

        // simd type of the same size storing int64_t values
        using simd_int64    = simd<int64_t, 128, scalar_nosimd_tag>;

        // simd type of the same kind storing vector of double values
        using simd_vector    = simd<double, 128, nosimd_tag>;

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

        // construct vector with all elements equal to val
        explicit simd(int64_t val);

        // construct vector with all elements equal to val
        explicit simd(float val);

        // construct vector with all elements equal to val
        explicit simd(double val);

      #if MATCL_ARCHITECTURE_HAS_SSE2
        // conversion form simd scalar; set all elements to s.first()
        explicit simd(const simd<double, 128, sse_tag>& s);
     #endif

      #if MATCL_ARCHITECTURE_HAS_SSE2
        // conversion form simd scalar; set all elements to s.first()
        explicit simd(const simd<double, 128, scalar_sse_tag>& s);
      #endif

        // conversion form simd scalar; set all elements to s.first()
        explicit simd(const simd<double, 128, nosimd_tag>& s);

        // copy constructor
        simd(const simd<double, 128, scalar_nosimd_tag>& s) = default;

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
        static simd     broadcast(const double* arr);

        // construct vector with all elements equal to arr
        static simd     broadcast(const double& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd     load(const double* arr, std::true_type aligned);
        static simd     load(const double* arr, std::false_type not_aligned = std::false_type());

        // gather double-precision (64-bit) floating-point elements from memory using 
        // 32-bit indices, i.e. i-th element of resulting vector is arr[ind[i]]
        static simd     gather(const double* arr, const simd_int32& ind);

        // gather double-precision (64-bit) floating-point elements from memory using 
        // 64-bit indices, i.e. i-th element of resulting vector is arr[ind[i]]
        static simd     gather(const double* arr, const simd_int64& ind);

        // set the first element in the vector to v and set 0 to all other elements
        static simd     set_lower(double v);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(double* arr, std::true_type aligned) const;
        void            store(double* arr, std::false_type not_aligned = std::false_type()) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(double* arr) const;

        // get i-th element from the vector; pos is 0-based
        double          get(int pos) const;

        // return the first element in the vector; equivalent to get(0), but possibly
        // faster
        double          first() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, double val);

        // return pointer to the first element in the vector
        const double*   get_raw_ptr() const;
        double*         get_raw_ptr();

    public:
        // convert elements to float and store the result in the lower part
        simd_float      convert_to_float() const;

        // convert elements to int32_t, rounding is performed according
        // to current rounding mode (usually round to nearest ties to even)
        simd_int32      convert_to_int32() const;

        // convert elements to int64_t, rounding is performed according
        // to current rounding mode (usually round to nearest ties to even)
        simd_int64      convert_to_int64() const;

        // reinterpret cast to vector of floats of the same kind
        simd_float      reinterpret_as_float() const;

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
