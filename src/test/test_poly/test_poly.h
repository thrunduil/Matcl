/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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

#include "matcl-core/config.h"
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-core/float/twofold.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-simd/details/poly/utils.h"

#include <string>

namespace matcl { namespace simd { namespace details
{

template<>
struct broadcast<mp_float>
{
    template<class Coef>
    force_inline
    static mp_float eval(const Coef* arr)
    {
        return mp_float(arr[0]);
    }

    template<class Coef>
    force_inline
    static mp_float eval(const Coef& arr)
    {
        return mp_float(arr);
    }
};

}}}

namespace matcl { namespace test
{

void test_poly();
void test_poly_cond();
void test_poly_dyn(bool pow2);

class test_polynomials
{
    public:
        template<class T>
        void    make(formatted_disp& dm);

        template<class T>
        void    make_cond(formatted_disp& dm);

        template<class T>
        void    make_dynamic(formatted_disp& dm);

        template<class T>
        void    make_dynamic_pow2(formatted_disp& dm);

    private:
        template<class T, class TS, int Size>
        void    make2(formatted_disp& dm);

        template<class T, class TS>
        void    make_cond2(formatted_disp& dm, int size);

        template<class T, class TS>
        void    make_dynamic(formatted_disp& dm, int size);

        template<class T, class Poly>
        bool    test_equal_ulp(int size, const T* x, const T* res1, const mp_float* res2,
                    double max_err, double& max_error);

        template<class T, class Poly>
        bool    test_equal_cond(int size, const T* x, const T* res1, const mp_float* res2,
                    double max_err, double& max_error);

        template<class T, class Poly>
        bool    test_equal_cond_comp(int size, const T* x, const mp_float* res2, 
                    double& max_error);

        template<class T>
        bool    test_equal_cond(T res1, const mp_float& res2, T cond, double max_err, 
                                double& error);
        template<class T>
        bool    test_equal_cond_comp(T res1, const mp_float& res2, T res_err, double& error);

        template<class T>
        bool    test_equal_ulp(T res1, const mp_float& res2, double max_err, double& error);
};

template<class T, int Size>
struct random_polynomial
{
    static const int size   = Size;
    static T polynomial[size];

    static void rand();
    static int get_size();
};

template<class T, int Max_size>
struct random_polynomial_dyn
{    
    static int  m_size;
    static T    polynomial[Max_size];

    static void rand(int size);
    static int  get_size();
};

template<class T>
T rand_scalar();

template<class T, int Size>
T random_polynomial<T,Size>::polynomial[size];

template<class T, int Size>
void random_polynomial<T,Size>::rand()
{
    for (int i = 0; i < size; ++i)
        polynomial[i] = rand_scalar<T>();
}

template<class T, int Size>
int random_polynomial<T,Size>::get_size()
{
    return size;
}

template<class T, int Max_size>
T random_polynomial_dyn<T,Max_size>::polynomial[Max_size];

template<class T, int Max_size>
int random_polynomial_dyn<T,Max_size>::m_size = 0;

template<class T, int Max_size>
void random_polynomial_dyn<T, Max_size>::rand(int size0)
{
    int size    = std::min(size0, Max_size);
    m_size      = size;

    for (int i = 0; i < size; ++i)
        polynomial[i] = rand_scalar<T>();
}

template<class T, int Size>
int random_polynomial_dyn<T,Size>::get_size()
{
    return m_size;
}

template<class T, class Poly>
struct Func_horn_1
{   
    template<class TS>
    force_inline
    static TS eval(const TS& arg)
    {
        TS res   = simd::horner<Poly::size>(arg, Poly::polynomial);
        return res;
    };
};

template<class T, class Poly>
struct Func_horn_2
{   
    template<class TS>
    force_inline
    static TS eval(const TS& arg)
    {
        TS res   = simd::horner(arg, Poly::get_size(), Poly::polynomial);
        return res;
    };
};

template<class T, class Poly>
struct Func_horn_cond
{   
    template<class TS>
    force_inline
    static TS eval(const TS& arg)
    {
        TS err;
        TS res   = simd::horner_and_error(arg, Poly::get_size(), Poly::polynomial, err);
        return res;
    };
};

template<class T, class Poly>
struct Func_horn_comp_cond
{   
    template<class TS>
    force_inline
    static TS eval(const TS& arg)
    {
        TS err;
        bool exact;
        TS res   = simd::compensated_horner_and_error(arg, Poly::get_size(), 
                    Poly::polynomial, err, exact);
        return res;
    };
};

template<class T, class Poly>
struct Func_estrin
{
    template<class TS>
    force_inline
    static TS eval(const TS& arg)
    {
        TS res   = simd::estrin<Poly::size>(arg, Poly::polynomial);
        return res;
    };
};

template<class T, class Poly>
struct Func_estrin2
{
    template<class TS>
    force_inline
    static TS eval(const TS& arg)
    {
        TS res   = simd::estrin(arg, Poly::get_size(), Poly::polynomial);
        return res;
    };
};

template<class T, class Poly>
struct Func_horn_twofold
{   
    template<class TS>
    force_inline
    static TS eval(const TS& arg)
    {
        auto res   = simd::compensated_horner<TS, T>(arg, Poly::size, Poly::polynomial);
        return res.value;
    };
};

template<class Poly>
struct Func_horn_mp
{
    template<class TS>
    static mp_float eval(const TS& arg)
    {
        precision prec  = precision(4*53);
        mp_float x      = mp_float(arg, prec);
        mp_float res    = simd::horner(x, Poly::get_size(), Poly::polynomial);

        return res;
    };
};

}}
