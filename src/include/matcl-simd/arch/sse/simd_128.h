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

#include <emmintrin.h>

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 DOUBLE
//-------------------------------------------------------------------
// vector of two double precision scalars
template<>
class alignas(16) simd<double, 128, sse_tag>
{
    public:
        // implementation type
        using impl_type     = __m128d;

        // type of stored elements
        using value_type    = double;

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
        explicit simd(Integer val);

        // construct vector with all elements equal to val
        explicit simd(float val);

        // construct vector with all elements equal to val
        explicit simd(double val);

        // set the first element in the vector to v1 and the second element to v2
        simd(double v1, double v2);

        // construct from representation
        simd(const impl_type& v);

    public:
        // connstruct vector with all elements set to zero
        static simd     zero();

        // construct vector with all elements equal to arr[0]
        static simd     broadcast(const double* arr);

        // construct vector with all elements equal to arr
        static simd     broadcast(const double& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd     load(const double* arr, std::true_type aligned);
        static simd     load(const double* arr, std::false_type not_aligned);

        // set the first element in the vector to v and set 0 to all other elements
        static simd     set_lower(double v);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(double* arr, std::true_type aligned) const;
        void            store(double* arr, std::false_type not_aligned) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(double* arr) const;

        // get i-th element from the vector; pos is 0-based
        double          get(int pos) const;

        // get i-th element from the vector; Pos is 0-based
        template<int Pos>
        double          get() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, double val);

        // set i-th element of the vector; Pos is 0-based
        template<int Pos>
        void            set(double val);
};

//-------------------------------------------------------------------
//                          SSE2 FLOAT
//-------------------------------------------------------------------
// vector of four single precision scalars
template<>
class alignas(16) simd<float, 128, sse_tag>
{
    public:
        // implementation type
        using impl_type     = __m128;

        // type of stored elements
        using value_type    = float;

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
        explicit simd(float val);

        // construct vector with i-th element set to vi
        simd(float v0, float v1, float v2, float v3);

        // construct from representation
        simd(const impl_type& v);

    public:
        // connstruct vector with all elements set to zero
        static simd     zero();

        // construct vector with all elements equal to arr[0]
        static simd     broadcast(const float* arr);

        // construct vector with all elements equal to arr
        static simd     broadcast(const float& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd     load(const float* arr, std::true_type aligned);
        static simd     load(const float* arr, std::false_type not_aligned);

        // set the first element in the vector to v and set 0 to all other elements
        static simd     set_lower(float v);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(float* arr, std::true_type aligned) const;
        void            store(float* arr, std::false_type not_aligned) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(float* arr) const;

        // get i-th element from the vector; pos is 0-based
        float           get(int pos) const;

        // get i-th element from the vector; Pos is 0-based
        template<int Pos>
        float           get() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, float val);

        // set i-th element of the vector; Pos is 0-based
        template<int Pos>
        void            set(float val);
};

}}
