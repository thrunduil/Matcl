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

#include "matcl-simd/arch/reg_128/simd_128.h"
#include "matcl-simd/simd_complex.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          SSE2 DOUBLE COMPLEX
//-------------------------------------------------------------------
// vector of one double precision complex scalar
template<>
class alignas(16) simd_compl<double, reg_128>
{
    public:
        // implementation type
        using impl_type     = simd<double, reg_128>;

        // type of stored elements
        using value_type    = simd_double_complex;

        // type of real and imaginary part of stored elements
        using real_type     = double;

    public:
        // number of elements in the vector
        static const int 
        vector_size         = sizeof(impl_type) / sizeof(real_type) / 2;    

    public:
        // internal representation
        impl_type           data;

    public:
        // construct uninitialized vector
        simd_compl() = default;

        // construct vector storing complex(re, 0)
        explicit simd_compl(Integer re);

        // construct vector storing complex(re, 0)
        explicit simd_compl(float re);

        // construct vector storing complex(re, 0)
        explicit simd_compl(double re);

        // construct vector storing val
        explicit simd_compl(const simd_double_complex& val);

        // construct vector storing val
        explicit simd_compl(const simd_single_complex& val);

        // construct vector storing complex(re, im)
        simd_compl(double re, double im);

        // construct from representation
        simd_compl(const impl_type& v);

    public:
        // connstruct vector with all elements set to zero
        static simd_compl   zero();

        // construct vector storing arr[0]
        static simd_compl   broadcast(const simd_double_complex* arr);

        // construct vector storing arr
        static simd_compl   broadcast(const simd_double_complex& arr);

        // construct vector storing complex(arr[0], 0)
        static simd_compl   broadcast(const double* arr);

        // construct vector storing complex(arr, 0)
        static simd_compl   broadcast(const double& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd_compl   load(const simd_double_complex* arr, std::true_type aligned);
        static simd_compl   load(const simd_double_complex* arr, std::false_type not_aligned);

        //load in reverse order, arr points to the last element to load
        static simd_compl   load_reverse(const simd_double_complex* arr, std::true_type aligned);
        static simd_compl   load_reverse(const simd_double_complex* arr, std::false_type not_aligned);

    public:
        // store elements in arr; arr must have length at least vector_size
        void                store(simd_double_complex* arr, std::true_type aligned) const;
        void                store(simd_double_complex* arr, std::false_type not_aligned) const;

        //store in reverse order, arr points to the last element of store
        void                store_reverse(simd_double_complex* arr, std::true_type aligned) const;
        void                store_reverse(simd_double_complex* arr, std::false_type not_aligned) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void                scatter(simd_double_complex* arr) const;

        // get i-th element from the vector; pos is 0-based
        simd_double_complex get(int pos) const;
};

//-------------------------------------------------------------------
//                          SSE2 FLOAT COMPLEX
//-------------------------------------------------------------------
// vector of two single precision complex scalars
template<>
class alignas(16) simd_compl<float, reg_128>
{
    public:
        // implementation type
        using impl_type     = simd<float, reg_128>;

        // type of stored elements
        using value_type    = simd_single_complex;

        // type of real and imaginary part of stored elements
        using real_type     = float;

    public:
        // number of elements in the vector
        static const int 
        vector_size         = sizeof(impl_type) / sizeof(real_type) / 2;    

    public:
        // internal representation
        impl_type           data;

    public:
        // construct uninitialized vector 
        simd_compl() = default;

        // construct vector with all elements equal to complex(re,0)
        explicit simd_compl(float re);

        // construct vector with all elements equal to val
        explicit simd_compl(const simd_single_complex& val);

        // construct vector with all elements equal to complex(re, im)
        simd_compl(float re, float im);

        // construct vector with i-th element set to vi
        simd_compl(const simd_single_complex& v0, const simd_single_complex& v1);

        // construct vector with i-th element set to complex(re_i, im_i)
        simd_compl(float re_0, float im_0, float re_1, float im_1);

        // construct from representation
        simd_compl(const impl_type& v);

    public:
        // connstruct vector with all elements set to zero
        static simd_compl   zero();

        // construct vector with all elements equal to arr[0]
        static simd_compl   broadcast(const simd_single_complex* arr);

        // construct vector with all elements equal to arr
        static simd_compl   broadcast(const simd_single_complex& arr);

        // construct vector with all elements equal to complex(arr[0], 0)
        static simd_compl   broadcast(const float* arr);

        // construct vector with all elements equal to complex(arr, 0)
        static simd_compl   broadcast(const float& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd_compl   load(const simd_single_complex* arr, std::true_type aligned);
        static simd_compl   load(const simd_single_complex* arr, std::false_type not_aligned);

        //load in reverse order, arr points to the last element to load
        static simd_compl   load_reverse(const simd_single_complex* arr, std::true_type aligned);
        static simd_compl   load_reverse(const simd_single_complex* arr, std::false_type not_aligned);

    public:
        // store elements in arr; arr must have length at least vector_size
        void                store(simd_single_complex* arr, std::true_type aligned) const;
        void                store(simd_single_complex* arr, std::false_type not_aligned) const;

        //store in reverse order, arr points to the last element of store
        void                store_reverse(simd_single_complex* arr, std::true_type aligned) const;
        void                store_reverse(simd_single_complex* arr, std::false_type not_aligned) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void                scatter(simd_single_complex* arr) const;

        // get i-th element from the vector; pos is 0-based
        simd_single_complex get(int pos) const;
};

}}
