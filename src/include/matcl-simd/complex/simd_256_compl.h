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

#include "matcl-simd/details/arch/simd_impl.h"
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

        // simd tag
        using simd_tag      = Simd_tag;

        // number of bits
        static const int
        number_bits         = 256;

        // type of real and imaginary part of stored elements
        using real_type     = double;

        // type of vector storing half of elements
        using simd_half     = simd_compl<double, 128, simd_half_tag>;

        // simd type storing float values
        using simd_float   = simd_compl<float, 128, simd_half_tag>;

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
        explicit simd_compl(int32_t re);

        // construct vector with all elements equal to complex(re, 0)
        explicit simd_compl(int64_t re);

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

        // conversion between simd types
        template<class Tag>
        explicit simd_compl(const simd_compl<double, 256, Tag>& s);

        // copy constructor
        simd_compl(const simd_compl<double, 256, Simd_tag>& s) = default;

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
        static simd_compl   load(const simd_double_complex* arr, 
                                std::false_type not_aligned = std::false_type());

    public:
        // store elements in arr; arr must have length at least vector_size
        void                store(simd_double_complex* arr, std::true_type aligned) const;
        void                store(simd_double_complex* arr, 
                                std::false_type not_aligned = std::false_type()) const;

        // store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void                scatter(simd_double_complex* arr) const;

        // get i-th element from the vector; pos is 0-based
        simd_double_complex get(int pos) const;

        // return the first element in the vector; equivalent to get(0), 
        // but possibly faster
        simd_double_complex first() const;

        // set i-th element of the vector; pos is 0-based
        void                set(int pos, const simd_double_complex& val);

        // return pointer to the first element in the vector
        const simd_double_complex*
                            get_raw_ptr() const;
        simd_double_complex*
                            get_raw_ptr();

        // return simd storing first element
        simd_half           extract_low() const;

        // return simd storing last element
        simd_half           extract_high() const;

        // convert elements to float complex
        simd_float          convert_to_float() const;

    public:
        // plus assign operator
        simd_compl&         operator+=(const simd_compl& x);

        // minus assign operator
        simd_compl&         operator-=(const simd_compl& x);

        // multiply assign operator
        simd_compl&         operator*=(const simd_compl& x);

        // divide assign operator
        simd_compl&         operator/=(const simd_compl& x);
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

        // simd type of the same kind storing double complex values
        using simd_double   = simd_compl<double, 256, Simd_tag>;

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
        simd_compl(const impl_type& v);

        // conversion between simd types
        template<class Tag>
        explicit simd_compl(const simd_compl<float, 256, Tag>& s);

        // copy constructor
        simd_compl(const simd_compl<float, 256, Simd_tag>& s) = default;

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
        static simd_compl   load(const simd_single_complex* arr, 
                                std::false_type not_aligned = std::false_type());

    public:
        // store elements in arr; arr must have length at least vector_size
        void                store(simd_single_complex* arr, std::true_type aligned) const;
        void                store(simd_single_complex* arr, 
                                std::false_type not_aligned = std::false_type()) const;

        //store elements with in array with stepping (Step = 1,-1 is not optimized)
        template<int Step>
        void                scatter(simd_single_complex* arr) const;

        // get i-th element from the vector; pos is 0-based
        simd_single_complex get(int pos) const;

        // return the first element in the vector; equivalent to get(0), 
        // but possibly faster
        simd_single_complex first() const;

        // set i-th element of the vector; pos is 0-based
        void                set(int pos, const simd_single_complex& val);

        // return pointer to the first element in the vector
        const simd_single_complex*
                            get_raw_ptr() const;
        simd_single_complex*
                            get_raw_ptr();

        // return simd storing first two elements
        simd_half           extract_low() const;

        // return simd storing last two elements
        simd_half           extract_high() const;

        // convert the first four elements to double
        simd_double         convert_low_to_double() const;

        // convert the last four elements to double
        simd_double         convert_high_to_double() const;

    public:
        // plus assign operator
        simd_compl&         operator+=(const simd_compl& x);

        // minus assign operator
        simd_compl&         operator-=(const simd_compl& x);

        // multiply assign operator
        simd_compl&         operator*=(const simd_compl& x);

        // divide assign operator
        simd_compl&         operator/=(const simd_compl& x);
};

}}
