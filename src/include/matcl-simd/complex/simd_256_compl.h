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

#include "matcl-simd/arch/simd_impl.h"
#include "matcl-simd/complex/simd_complex.h"
#include "matcl-simd/details/helpers.h"

namespace matcl { namespace simd
{

//-------------------------------------------------------------------
//                          DOUBLE COMPLEX
//-------------------------------------------------------------------
// vector of two double precision complex scalars
template<class Simd_tag>
class alignas(32) simd_compl<double, 256, Simd_tag>
{
    private:
        using simd_half_tag = typename details::simd_half_tag<Simd_tag>::type;

    public:
        // implementation type
        using impl_type     = simd<double, 256, Simd_tag>;

        // type of stored elements
        using value_type    = simd_double_complex;

        // type of real and imaginary part of stored elements
        using real_type     = double;

        // type of vector storing half of elements
        using simd_half     = simd_compl<double, 128, simd_half_tag>;

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

        // construct vector with all elements equal to complex(re, 0)
        explicit simd_compl(Integer re);

        // construct vector with all elements equal to complex(re, 0)
        explicit simd_compl(float re);

        // construct vector with all elements equal to complex(re, 0)
        explicit simd_compl(double re);

        // construct vector with all elements equal to val
        explicit simd_compl(const simd_double_complex& val);

        // construct vector with all elements equal to val
        explicit simd_compl(const simd_single_complex& val);

        // construct vector with first and last element copied from lo_hi
        explicit simd_compl(const simd_half& lo_hi);

        // construct vector with all elements equal to complex(re, im)
        simd_compl(double re, double im);

        // construct vector with i-th element set to vi
        simd_compl(const simd_double_complex& v0, const simd_double_complex& v1);

        // construct vector with i-th element set to complex(re_i, im_i)
        simd_compl(double re_0, double im_0, double re_1, double im_1);

        // construct vector with first element copied from lo
        // and last element copied from hi
        simd_compl(const simd_half& lo, const simd_half& hi);

        // construct from representation
        simd_compl(const impl_type& v);

    public:
        // connstruct vector with all elements set to zero
        static simd_compl   zero();

        // construct vector with all elements equal to arr[0]
        static simd_compl   broadcast(const simd_double_complex* arr);

        // construct vector with all elements equal to complex(arr[0], 0)
        static simd_compl   broadcast(const double* arr);

        // construct vector with all elements equal to arr
        static simd_compl   broadcast(const simd_double_complex& arr);

        // construct vector with all elements equal to complex(arr, 0)
        static simd_compl   broadcast(const double& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd_compl   load(const simd_double_complex* arr, std::true_type aligned);
        static simd_compl   load(const simd_double_complex* arr, std::false_type not_aligned);

    public:
        // store elements in arr; arr must have length at least vector_size
        void                store(simd_double_complex* arr, std::true_type aligned) const;
        void                store(simd_double_complex* arr, std::false_type not_aligned) const;

        // store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void                scatter(simd_double_complex* arr) const;

        // get i-th element from the vector; pos is 0-based
        simd_double_complex get(int pos) const;

        // get i-th element from the vector; Pos is 0-based
        template<int Pos>
        simd_double_complex get() const;

        // return simd storing first element
        simd_half           extract_low() const;

        // return simd storing last element
        simd_half           extract_high() const;
};

//-------------------------------------------------------------------
//                          FLOAT COMPLEX
//-------------------------------------------------------------------
// vector of four single precision complex scalars
template<class Simd_tag>
class alignas(32) simd_compl<float, 256, Simd_tag>
{
    private:
        using simd_half_tag = typename details::simd_half_tag<Simd_tag>::type;

    public:
        // implementation type
        using impl_type     = simd<float, 256, Simd_tag>;

        // type of stored elements
        using value_type    = simd_single_complex;

        // type of real and imaginary part of stored elements
        using real_type     = float;

        // type of vector storing half of elements
        using simd_half     = simd_compl<float, 128, simd_half_tag>;

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

        // construct vector with all elements equal to complex(re, 0)
        explicit simd_compl(float re);

        // construct vector with all elements equal to val
        explicit simd_compl(const simd_single_complex& val);

        // construct vector with first two and last element two copied from lo_hi
        explicit simd_compl(const simd_half& lo_hi);

        // construct vector with all elements equal to complex(re, im)
        simd_compl(float re, float im);

        // construct vector with i-th element set to vi
        simd_compl(const simd_single_complex& v0, const simd_single_complex& v1,
                   const simd_single_complex& v2, const simd_single_complex& v3);

        // construct vector with i-th element set to complex(re_i, im_i)
        simd_compl(float re_0, float im_0, float re_1, float im_1,
                   float re_2, float im_2, float re_3, float im_3);

        // construct vector with first two elements copied from lo
        // and last two elements copied from hi
        simd_compl(const simd_half& lo, const simd_half& hi);

        // construct from representation
        force_inline
        simd_compl(const impl_type& v);

    public:
        // connstruct vector with all elements set to zero
        static simd_compl   zero();

        // construct vector with all elements equal to arr[0]
        static simd_compl   broadcast(const simd_single_complex* arr);
        static simd_compl   broadcast(const float* arr);

        // construct vector with all elements equal to arr
        static simd_compl   broadcast(const simd_single_complex& arr);
        static simd_compl   broadcast(const float& arr);

        // construct vector with elements copied from arr; arr must have length
        // at least vector_size
        static simd_compl   load(const simd_single_complex* arr, std::true_type aligned);
        static simd_compl   load(const simd_single_complex* arr, std::false_type not_aligned);

    public:
        // store elements in arr; arr must have length at least vector_size
        void                store(simd_single_complex* arr, std::true_type aligned) const;
        void                store(simd_single_complex* arr, std::false_type not_aligned) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void                scatter(simd_single_complex* arr) const;

        // get i-th element from the vector; pos is 0-based
        simd_single_complex get(int pos) const;

        // get i-th element from the vector; Pos is 0-based
        template<int Pos>
        simd_single_complex get() const;

        // return simd storing first two elements
        simd_half           extract_low() const;

        // return simd storing last two elements
        simd_half           extract_high() const;
};

}}
