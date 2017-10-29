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
#include "matcl-simd/simd/simd_general.h"
#include "matcl-simd/simd_complex.h"

#include <immintrin.h>

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          AVX DOUBLE
//-------------------------------------------------------------------
// vector of four double precision scalars
template<>
class alignas(32) simd<double, reg_256>
{
    public:
        // implementation type
        using impl_type     = __m256d;

        // type of stored elements
        using value_type    = double;

        // type of vector storing half of elements
        using simd_half     = simd<double, reg_128>;

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
        explicit simd(const double& val);

        // construct vector with first two elements and last two elements
        // copied from lo_hi
        explicit simd(const simd_half& lo_hi);

        // construct vector with i-th element set to vi
        simd(double v0, double v1, double v2, double v3);

        // construct vector with first two elements copied from lo
        // and last two elements copied from hi
        simd(simd_half lo, simd_half hi);

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

        // load in reverse order, arr points to the last element to load
        static simd     load_reverse(const double* arr, std::true_type aligned);
        static simd     load_reverse(const double* arr, std::false_type not_aligned);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(double* arr, std::true_type aligned) const;
        void            store(double* arr, std::false_type not_aligned) const;

        // store in reverse order, arr points to the last element of store
        void            store_reverse(double* arr, std::true_type aligned) const;
        void            store_reverse(double* arr, std::false_type not_aligned) const;

        // store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(double* arr) const;

        // get i-th element from the vector; pos is 0-based
        double          get(int pos) const;

        // return simd storing first two elements
        simd_half       extract_low() const;

        // return simd storing last two elements
        simd_half       extract_high() const;
};

//-------------------------------------------------------------------
//                          AVX FLOAT
//-------------------------------------------------------------------
// vector of eight single precision scalars
template<>
class alignas(32) simd<float, reg_256>
{
    public:
        // implementation type
        using impl_type     = __m256;

        // type of stored elements
        using value_type    = float;

        // type of vector storing half of elements
        using simd_half     = simd<float, reg_128>;

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
        explicit simd(const float& val);

        // construct vector with first four elements and last four elements
        // copied from lo_hi
        explicit simd(const simd_half& lo_hi);

        // construct vector with i-th element set to vi
        simd(float v0, float v1, float v2, float v3, float v4, float v5, float v6, float v7);

        // construct from representation
        simd(const impl_type& v);

        // construct vector with first four elements copied from lo
        // and last four elements copied from hi
        simd(simd_half lo, simd_half hi);

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

        //load in reverse order, arr points to the last element to load
        static simd     load_reverse(const float* arr, std::true_type aligned);
        static simd     load_reverse(const float* arr, std::false_type not_aligned);

    public:
        // store elements in arr; arr must have length at least vector_size
        void            store(float* arr, std::true_type aligned) const;
        void            store(float* arr, std::false_type not_aligned) const;

        //store in reverse order, arr points to the last element of store
        void            store_reverse(float* arr, std::true_type aligned) const;
        void            store_reverse(float* arr, std::false_type not_aligned) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void            scatter(float* arr) const;

        // get i-th element from the vector; pos is 0-based
        float           get(int pos) const;

        // return simd storing first four elements
        simd_half       extract_low() const;

        // return simd storing last four elements
        simd_half       extract_high() const;
};

}}
