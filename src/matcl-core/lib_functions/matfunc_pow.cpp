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

#include "matcl-core/details/scalfunc_complex.h"

#include <cfloat>

namespace matcl { namespace raw { namespace details
{

namespace mr = matcl :: raw;

inline Complex exp_complex(const Complex& arg)
{
    return scal_func::exp(arg);
};

inline Float_complex exp_complex(const Float_complex& arg)
{
    return scal_func::exp(arg);
};

template<class Ret,class T1,class T2,bool isc_1,bool isc_2>
struct pow_complex_helper{};

template<class Ret,class T1,class T2>
struct pow_complex_helper<Ret,T1,T2,true,true>
{
    static Ret eval(const T1& arg1, const T2& arg2)
    {
        if (imag(arg2) == 0)
        {
            using ST2 = typename T2::real_type;
            return pow_complex_helper<Ret,T1,ST2,true,false>::eval(arg1, real(arg2));
        }
        else if (imag(arg1) == 0)
        {
            using ST1 = typename T1::real_type;
            return pow_complex_helper<Ret,ST1,T2,false,true>::eval(real(arg1), arg2);
        }
        else
        {
            Ret tmp = matcl::details::mul_c(arg2 , scal_func::log(Ret(arg1)));
            return exp_complex(tmp);
        };
    };
};

template<class Ret,class T1,class T2>
struct pow_complex_helper<Ret,T1,T2,true,false>
{
    static Ret eval(const T1& arg1, const T2& arg2)
    {
        using Ret_R = typename md::real_type<Ret>::type;

        //x^0 = 1 for any x including nan
        if (arg2 == T2())
            return Ret(1,0);
        if (arg2 == 1)
            return Ret(arg1);

        if (arg2 == -1)
        {
            using T1P   = typename md::unify_types<Ret_R, T1>::type;
            return matcl::details::inv_c(T1P(arg1));
        };

        if (imag(arg1) == 0)
        {
            if (real(arg1) >= 0)
                return Ret(std::pow(Ret_R(real(arg1)), Ret_R(arg2)));

            if (scal_func::isregular(arg2) == false || scal_func::finite(real(arg1)) == false)
                return Ret(std::pow(Ret_R(real(arg1)), Ret_R(arg2)));
        };
    
        Ret tmp = matcl::details::mul_c(arg2 , scal_func::log(Ret(arg1)));
        return exp_complex(tmp);
    };
};

template<class Ret,class T1>
struct pow_complex_helper<Ret,T1,Integer,true,false>
{
    using real_type = typename md::real_type<Ret>::type;

    static Ret eval(const T1& arg1, Integer arg2)
    {
        using Ret_R             = typename md::real_type<Ret>::type;

        //x^0 = 1 for any x including nan
        if (arg2 == 0)
            return Ret(1.0, 0.0);
        if (arg2 == 1)
            return Ret(arg1);

        if (arg2 == -1)
            return md::inv_c(Ret(arg1));

        if (matcl::details::eeq_c(arg1 , T1(1.,0.)))
            return Ret(1.,0.);

        if (imag(arg1) == 0)
            return std::pow(real(arg1), arg2);

        if (scal_func::isregular(arg1) == false)
        {
            if (matcl::details::eeq_c(arg1 , T1()))
            {
                if (arg2 == 0)
                    return Ret(1.,0.);
                else if (arg2 > 0)
                    return Ret(0.,0.);
                else
                    return Ret(constants::inf<real_type>(),0.);
            }
            
            //algorith based on multiplication given incorrect
            //results in case of not finite numbers
            return pow_complex_helper<Ret, T1, Real, true, false>
                ::eval(arg1, Real(arg2));
        };

        Ret tmp                 = arg1;
        unsigned long count     = arg2;

        if (arg2 < 0)
        {
            count = 0 - count;
        };

        for (Ret value = Ret(1); ; tmp *= tmp)
        {	
            // fold in _Left ^ (2 ^ count) as needed
            if ((count & 1) != 0)
            {
                value *= tmp;
            };
            if ((count >>= 1) == 0)
            {
                if (arg2 < 0)
                {
                    auto tmp2 = matcl::details::div_c(Ret_R(1),value);
                    return Ret(tmp2);
                }
                else
                {
                    return Ret(value);
                }
            };
        };
    };
};

template<class Ret,class T1,class T2>
struct pow_complex_helper<Ret,T1,T2,false,true>
{
    static Ret eval(const T1& arg1, const T2& arg2)
    {
        using Ret_R = typename md::real_type<Ret>::type;

        //x^0 = 1 for any x including nan
        if (md::eeq_c(arg2, 0))
            return Ret(1.0, 0.0);
        if (md::eeq_c(arg2, 1))
            return Ret(arg1);
        else if (md::eeq_c(arg2, -1))
            return md::inv_c(Ret(arg1));

        if (imag(arg2) == 0)
        {
            if (arg1 >= 0)
                return Ret(std::pow(arg1, real(arg2)));
            if (scal_func::isregular(real(arg2)) == false || scal_func::finite(arg1) == false)
                return Ret(std::pow(arg1, real(arg2)));

            Ret tmp     = scal_func::log(Ret(arg1));
            tmp         = matcl::details::mul_c(real(arg2) , tmp);
            return exp_complex(tmp);
        };

        if (arg1 > T1())
        {
            if (scal_func::finite(arg1) == false 
                        || scal_func::finite(arg2) == false)
            {
                Ret tmp = matcl::details::mul_c(arg2, scal_func::log(arg1));
                return exp_complex(tmp);
            };

            Ret_R   x_re    = arg1;
            Ret_R   y_re    = real(arg2);
            Ret_R   y_im    = imag(arg2);

            // a ^ (x + iy) = a ^ x * exp(i * log(a) * y)
            Ret_R exp_re    = scal_func::pow(x_re, y_re);
            Ret_R mult_im   = y_im * scal_func::log(x_re);
            Ret exp_im      = scal_func::expi(mult_im);

            Ret_R r_re      = scal_func::mult_zero(exp_re , real(exp_im));
            Ret_R r_im      = scal_func::mult_zero(exp_re , imag(exp_im));
            Ret ret         = Ret(r_re, r_im);

            return ret;
        }
        else
        {
            Ret tmp     = scal_func::log(Ret(arg1));
            tmp         = matcl::details::mul_c(arg2 , tmp);
            return exp_complex(tmp);
        };
    };
};

template<class Ret, class T1, class T2>
Ret pow_complex<Ret, T1,T2>::eval(const T1& arg1, const T2& arg2)
{
    static const bool isc_1 = md::is_complex<T1>::value;
    static const bool isc_2 = md::is_complex<T2>::value;

    return pow_complex_helper<Ret,T1,T2,isc_1,isc_2>::eval(arg1,arg2);
};

};};};

template struct matcl::raw::details::pow_complex<matcl::Float_complex, matcl::Float_complex,matcl::Float_complex>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Float_complex,matcl::Complex>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Float_complex,matcl::Real>;
template struct matcl::raw::details::pow_complex<matcl::Float_complex, matcl::Float_complex,matcl::Float>;
template struct matcl::raw::details::pow_complex<matcl::Float_complex, matcl::Float_complex,matcl::Integer>;

template struct matcl::raw::details::pow_complex<matcl::Float_complex, matcl::Float_complex,matcl::Integer>;

template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Complex,matcl::Float_complex>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Complex,matcl::Complex>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Complex,matcl::Real>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Complex,matcl::Float>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Complex,matcl::Integer>;

template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Real,matcl::Float_complex>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Real,matcl::Complex>;

template struct matcl::raw::details::pow_complex<matcl::Float_complex, matcl::Float,matcl::Float_complex>;
template struct matcl::raw::details::pow_complex<matcl::Complex, matcl::Float,matcl::Complex>;
