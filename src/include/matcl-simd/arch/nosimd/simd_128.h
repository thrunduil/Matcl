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
//                          GENERIC DOUBLE
//-------------------------------------------------------------------

// vector of two double precision scalars
template<>
class alignas(16) simd<double, 128, nosimd_tag>
{
    public:
        // implementation type
        using impl_type     = double[2];

        // type of stored elements
        using value_type    = double;

    public:
        // number of elements in the vector
        static const int 
        vector_size         = sizeof(impl_type) / sizeof(value_type);    

        // simd type of the same kind storing float values
        using simd_float   = simd<float, 128, nosimd_tag>;

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

        // construct vector with first element copied from lo
        // and last element copied from hi; only lower part of lo and hi is used
        simd(const simd& lo, const simd& hi);

        // construct from representation
        simd(const impl_type& v);

        // conversion between simd types
        explicit simd(const simd<double, 128, sse_tag>& s);

        // copy constructor
        simd(const simd<double, 128, nosimd_tag>& s) = default;

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

        // return the first element in the vector; equivalent to get(0), 
        // but possibly faster
        double          first() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, double val);

        // return pointer to the first element in the vector
        const double*   get_raw_ptr() const;
        double*         get_raw_ptr();

        // set i-th element of the vector; Pos is 0-based
        template<int Pos>
        void            set(double val);

        // cast elements to float complex and store the result in the lower part
        simd_float      cast_to_float() const;

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

//-------------------------------------------------------------------
//                          GENERIC FLOAT
//-------------------------------------------------------------------

// vector of four single precision scalars
template<>
class alignas(16) simd<float, 128, nosimd_tag>
{
    public:
        // implementation type
        using impl_type     = float[4];

        // type of stored elements
        using value_type    = float;

        // simd type of the same kind storing double values
        using simd_double   = simd<double, 128, nosimd_tag>;

        // simd type storing 4 double values
        using simd_double_2 = simd<double, 256, nosimd_tag>;

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

        // construct vector with first two elements copied from lo
        // and last two elements copied from hi; only lower part of lo and
        // hi is used
        simd(const simd& lo, const simd& hi);

        // construct from representation
        simd(const impl_type& v);

        // conversion between simd types
        explicit simd(const simd<float, 128, sse_tag>& s);

        // copy constructor
        simd(const simd<float, 128, nosimd_tag>& s) = default;

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
        static simd     load(const float* arr, std::false_type not_aligned);

        // set the first element in the vector to v and set 0 to all other elements
        static simd     set_lower(float v);

        // cast the first two elements to double
        simd_double     cast_low_to_double() const;

        // cast the last two elements to double
        simd_double     cast_high_to_double() const;

        // cast all elements to double
        simd_double_2   cast_to_double() const;

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

        // return the first element in the vector; equivalent to get(0), 
        // but possibly faster
        float           first() const;

        // set i-th element of the vector; pos is 0-based
        void            set(int pos, float val);

        // set i-th element of the vector; Pos is 0-based
        template<int Pos>
        void            set(float val);

        // return pointer to the first element in the vector
        const float*    get_raw_ptr() const;
        float*          get_raw_ptr();

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
