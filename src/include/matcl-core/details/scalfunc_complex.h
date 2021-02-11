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
#include "matcl-core/lib_functions/constants.h"
#include "matcl-core/details/complex_details.h"
#include "matcl-core/details/scalfunc_real.h"
#include "matcl-core/details/utils.h"

#include <cfloat>

#pragma warning(push)
#pragma warning(disable:4127)   //conditional expression is constant

#pragma warning(pop)

namespace matcl { namespace raw { namespace details
{

namespace md    = matcl::details;

template<class Ret, class T1, class T2>
struct MATCL_CORE_EXPORT pow_complex
{
    static Ret eval(const T1& arg1, const T2& arg2);
};

template<>
struct is_zero_helper<Float_complex>
{
    force_inline
    static bool eval(const Float_complex& val)
    { 
        return real(val) == 0.0f && imag(val) == 0.0f; 
    };
};

template<>
struct is_zero_helper<Complex>
{
    force_inline
    static bool eval(const Complex& val)
    { 
        return real(val) == 0.0 && imag(val) == 0.0; 
    };
};

template<>
struct is_one_helper<Float_complex>
{
    force_inline
    static bool eval(const Float_complex& val)
    { 
        return real(val) == 1.0f && imag(val) == 0.0f; 
    };
};

template<>
struct is_one_helper<Complex>
{
    force_inline
    static bool eval(const Complex& val)
    { 
        return real(val) == 1.0 && imag(val) == 0.0; 
    };
};

namespace scal_func
{
    namespace mrd = matcl::raw::details;

    //not available for complex; throw exception
    bool MATCL_CORE_EXPORT           signbit(const Complex& x);
    bool MATCL_CORE_EXPORT           signbit(const Float_complex& x);
    Integer MATCL_CORE_EXPORT        ifloor(const Complex& x);
    Integer MATCL_CORE_EXPORT        ifloor(const Float_complex& x);
    Integer MATCL_CORE_EXPORT        iceil(const Complex& x);
    Integer MATCL_CORE_EXPORT        iceil(const Float_complex& x);
    Integer MATCL_CORE_EXPORT        itrunc(const Complex& x);
    Integer MATCL_CORE_EXPORT        itrunc(const Float_complex& x);
    Integer MATCL_CORE_EXPORT        iround(const Complex& x);
    Integer MATCL_CORE_EXPORT        iround(const Float_complex& x);
    Integer MATCL_CORE_EXPORT        isign(const Complex& x);
    Integer MATCL_CORE_EXPORT        isign(const Float_complex& x);
    Real MATCL_CORE_EXPORT           frexp(const Complex& x, Integer& n);
    Float MATCL_CORE_EXPORT          frexp(const Float_complex& x, Integer& n);
    Real MATCL_CORE_EXPORT           modf(const Complex& x, Real& n);
    Float MATCL_CORE_EXPORT          modf(const Float_complex& x, Float& n);
    Complex MATCL_CORE_EXPORT        cbrt(const Complex& arg);
    Float_complex MATCL_CORE_EXPORT  cbrt(const Float_complex& arg);

    Complex MATCL_CORE_EXPORT        sqrt1pm1(const Complex& arg);
    Float_complex MATCL_CORE_EXPORT  sqrt1pm1(const Float_complex& arg);    

    template<class Real>
    force_inline
    Real mult_zero(Real r, Real arg)
    {
        if (is_zero(arg))
            return arg;
        return r * arg;
    };

    //--------------------------------------------------------------------
    template<class Compl>
    struct complex_sign
    {
        using Real_t    = typename md::real_type<Compl>::type;

        static Compl eval(const Compl& x)
        {
            if (scal_func::isregular(x) == false)
            {
                if (mrd::is_zero(x) == true || scal_func::isnan(x) == true)
                    return x;

                bool inf_re = scal_func::isinf(real(x));
                bool inf_im = scal_func::isinf(imag(x));

                if (inf_re == true && inf_im == false)
                    return Compl(scal_func::sign(real(x)));
                if (inf_re == false && inf_im == true)
                    return Compl(Real_t(0), scal_func::sign(imag(x)));
                else
                    return Compl(constants::nan<Real_t>());
            };

            Real_t xa   = scal_func::abs(x);
            Real_t re   = real(x) / xa;
            Real_t im   = imag(x) / xa;

            return Compl(re, im);
        };
    };

    //--------------------------------------------------------------------
    template<class Compl>
    struct complex_exp
    {
        using Real_t    = typename md::real_type<Compl>::type;

        static Compl eval(const Compl& x)
        {
            Real_t x_re   = real(x);
            Real_t x_im   = imag(x);

            if (is_zero(x_im) == true)
                return scal_func:: exp(x_re);

            Real_t exp_re = scal_func::exp(x_re);

            Real_t sin_im, cos_im;
            
            scal_func::sin_cos(x_im, sin_im, cos_im);

            Real_t r_re     = mult_zero(exp_re , cos_im);
            Real_t r_im     = mult_zero(exp_re , sin_im);

            return Compl(r_re, r_im);
        };

        static Compl eval_m1(const Compl& x)
        {
            Real_t x_re   = real(x);
            Real_t x_im   = imag(x);

            if (is_zero(x_im) == true)
                return scal_func:: expm1(x_re);            

            // re = exp(xre) * cos(xim) - 1
            // im = exp(xre) * sin(xim)

            // error in calculating re is:
            //  e1  = |exp(xre) * cos(xim)|/|exp(xre) * cos(xim) - 1|*1.5 ulp
            //      <= [1 + 1/|exp(xre) * cos(xim) - 1|] * 1.5 ulp

            // re = exp(xre) * cos(xim) - 1 = expm1(xre) * cos(xim) - 2*sin(xim/2)^2
            // error in calculating re is
            //  e2 = c1 * 1.5 ulp + c2 * 1.5 ulp
            //  c1 = |expm1(xre) * cos(xim)|/|re|; c2 = |2*sin(xim/2)^2|/|re|
            //  e2 <= [1 + (1+|cos(xim)|) * |cos(xim) - 1|/|exp(xre)*cos(xim) - 1|] * 1.5 ulp
            // since (1+|cos(xim)|) * |cos(xim) - 1| <= 4, e2 is not much worse than e1
            // but can be much better, when cos(xim) ~ 1

            // select e1, when |xim| <= 1 and |xre| > 1.1, then 
            //  1/|exp(xre) * cos(xim) - 1| < 2 and e1 <= 4.5 ulp
            if (abs(x_re) > 1.1 && abs(x_im) < 1.0)
            {
                Compl ret = eval(x);
                return Compl(real(ret) - 1, imag(ret));
            }

            Real_t exp_re = scal_func::exp(x_re);

            Real_t sin_im, cos_im;
            
            scal_func::sin_cos(x_im, sin_im, cos_im);
            
            Real_t r_im     = mult_zero(exp_re , sin_im);

            Real_t exp1     = scal_func::expm1(x_re);
            Real_t sin2     = Real_t(2) * abs2(sin(x_im/2));
            Real_t r_re     = mult_zero(exp1 , cos_im) - sin2;

            return Compl(r_re, r_im);
        };

    };

    //--------------------------------------------------------------------
    template<class Compl>
    struct complex_log
    {
        using Real_t    = typename md::real_type<Compl>::type;

        // eval log(x)
        static Compl eval(const Compl& x)
        {
            Real_t x_re     = real(x);
            Real_t x_im     = imag(x);

            Real_t abs_xre  = abs(x_re);
            Real_t abs_xim  = abs(x_im);

            Real_t min_x    = Real_t(0.85);
            Real_t max_x    = Real_t(1.40);

            bool select_x   = (abs_xre >= min_x && abs_xre <= max_x && abs_xim < 0.5);
            bool select_y   = (abs_xim >= min_x && abs_xim <= max_x && abs_xre < 0.5);

            Real_t r_im     = scal_func::atan2(x_im, x_re);
            Real_t hyp      = scal_func::abs(x);    //1.0 ulp

            if (scal_func::finite(hyp) == false)
            {
                //log(inf) - inf; log(nan) = nan
                return Complex(hyp, r_im);
            };    

            if (select_x || select_y)
            {       
                // ep   = 2.0 * abs(e / hyp_m_1)
                // err1 = (2.5 + 2.5 * abs(ep)) * a + 0.5;
                // err2 = 1.0 * abs(1/log(hyp)) + 0.5
                // a    = |hyp_m_1/(log(hyp_m_1 + 1)*(hyp_m_1 + 1))|
                //
                // where err1 is error of this version; err2 is error of the second version 
                // measured in ulp;
                // err1/err2 < 1 for 0.85 <= x_re <= 1.4; -0.5 <= x_im <= 0.5

                if (select_x)
                {
                    Real_t e        = 1 - abs_xre;                          //0.5 ulp
                    Real_t hyp_m_1  = scal_func::abs2(x_im) 
                                    + scal_func::abs2(e) - 2 * e;           //2.0 + (2.0 + 0.5)*ep + 0.5
                    Real_t r_re     = scal_func::log1p(hyp_m_1) / Real_t(2);//(2.5 + 2.5*M)*a + 0.5
                    return Compl(r_re, r_im);
                }
                else
                {
                    Real_t e        = 1 - abs_xim;                          //0.5
                    Real_t hyp_m_1  = scal_func::abs2(x_re) 
                                    + scal_func::abs2(e) - 2 * e;           //2.0 + (2.0 + 0.5)*ep + 0.5
                    Real_t r_re     = scal_func::log1p(hyp_m_1) / Real_t(2);//(2.5 + 2.5*M)*a + 0.5
                    return Compl(r_re, r_im);
                }                
            }
            else
            {
                Real_t r_re = scal_func::log(hyp);   //1.0 * abs(1/log(hyp)) + 1
                return Compl(r_re, r_im);
            };
        };

        //eval log(1+x)
        static Compl eval_1p(const Compl& x)
        {
            Real_t x_re     = real(x);
            Real_t x_im     = imag(x);

            Real_t abs_xre  = abs(x_re);
            Real_t abs_xim  = abs(x_im);            

            // ep    = 2.0 * abs(x_re / hyp_m_1)
            // err1  = (1.5 + 1*ep) * a + 0.5
            // err2  = 2.0 * abs(1/log(hyp)) + 0.5
            // a     = |hyp_m_1/(log(hyp_m_1 + 1)*(hyp_m_1 + 1))|
            // 
            // |err1| < |err2| for x_re > -0.35

            bool select_1   = (x_re >= -0.35);
            (void)abs_xre;
            (void)abs_xim;

            if (select_1 == false)
            {
                //the standard version is more accurate
                return eval(x + Real_t(1));
            };

            Real_t r_im     = scal_func::atan2(x_im, Real_t(1) + x_re);

            // ep   = 2.0 * abs(x_re / hyp_m_1)
            // err  = 1.0 + ep * (1.0 + 0.0) + 0.5 = 1.5 + 1*ep
            Real_t hyp_m_1  = scal_func::abs2(x_im) + scal_func::abs2(x_re) 
                            + 2 * x_re;
            // err  = (1.5 + 1*ep) * 1.5 + 0.5 (assumption: hyp_m_1 > -0.5, 
            //                                  then ep[log1p] <= 1.5)
            Real_t r_re     = scal_func::log1p(hyp_m_1) / Real_t(2);
            return Compl(r_re, r_im);
        };
    };

    //--------------------------------------------------------------------
    template<class Compl>
    struct complex_sin
    {
        using Real_t    = typename md::real_type<Compl>::type;

        static Compl eval(const Compl& x)
        {
            //sin(z) = sin(x) * cosh(y) + i cos(x) * sinh(y)
            Real_t x_re   = real(x);
            Real_t x_im   = imag(x);

            if (is_zero(x_im) == true)
                return scal_func::sin(x_re);

            Real_t sin_re, cos_re;            
            scal_func::sin_cos(x_re, sin_re, cos_re);

            Real_t sinh_im, cosh_im;            
            scal_func::sinh_cosh(x_im, sinh_im, cosh_im);

            Real_t r_re     = mult_zero(cosh_im, sin_re);
            Real_t r_im     = mult_zero(sinh_im, cos_re);

            return Compl(r_re, r_im);
        };
    };

    //--------------------------------------------------------------------
    template<class Compl>
    struct complex_cos
    {
        using Real_t    = typename md::real_type<Compl>::type;

        static Compl eval(const Compl& x)
        {
            //cos(z) = cos(x) * cosh(y) - i sin(x) * sinh(y)
            Real_t x_re   = real(x);
            Real_t x_im   = imag(x);

            if (is_zero(x_im) == true)
                return scal_func::cos(x_re);

            Real_t sin_re, cos_re;            
            scal_func::sin_cos(x_re, sin_re, cos_re);

            Real_t sinh_im, cosh_im;            
            scal_func::sinh_cosh(x_im, sinh_im, cosh_im);

            Real_t r_re     = mult_zero(cosh_im, cos_re);
            Real_t r_im     = mult_zero(sinh_im, sin_re);

            return Compl(r_re, -r_im);
        };
    };

    //--------------------------------------------------------------------
    template<class Compl>
    struct complex_sinh
    {
        using Real_t    = typename md::real_type<Compl>::type;

        static Compl eval(const Compl& x)
        {
            //sinh(z) = sinh(x) * cos(y) + i cosh(x) * sin(y)
            Real_t x_re   = real(x);
            Real_t x_im   = imag(x);

            if (is_zero(x_im) == true)
                return scal_func::sinh(x_re);

            Real_t sin_im, cos_im;            
            scal_func::sin_cos(x_im, sin_im, cos_im);

            Real_t sinh_re, cosh_re;            
            scal_func::sinh_cosh(x_re, sinh_re, cosh_re);

            Real_t r_re     = mult_zero(sinh_re, cos_im);
            Real_t r_im     = mult_zero(cosh_re, sin_im);

            return Compl(r_re, r_im);
        };
    };

    //--------------------------------------------------------------------
    template<class Compl>
    struct complex_cosh
    {
        using Real_t    = typename md::real_type<Compl>::type;

        static Compl eval(const Compl& x)
        {
            //cosh(z) = cosh(x) * cos(y) + i sinh(x) * sin(y)
            Real_t x_re   = real(x);
            Real_t x_im   = imag(x);

            if (is_zero(x_im) == true)
                return scal_func::cosh(x_re);

            Real_t sin_im, cos_im;            
            scal_func::sin_cos(x_im, sin_im, cos_im);

            Real_t sinh_re, cosh_re;            
            scal_func::sinh_cosh(x_re, sinh_re, cosh_re);

            Real_t r_re     = mult_zero(cosh_re, cos_im);
            Real_t r_im     = mult_zero(sinh_re, sin_im);

            return Compl(r_re, r_im);
        };
    };

    //--------------------------------------------------------------------
    force_inline bool isnan(const Float_complex &x)
    {
        return scal_func::isnan(real(x)) || scal_func::isnan(imag(x));
    }
    force_inline bool isnan(const Complex &x)
    {
        return scal_func::isnan(real(x)) || scal_func::isnan(imag(x));
    }

    //--------------------------------------------------------------------
    force_inline bool finite(const Complex &x)
    {
        return scal_func::finite(real(x)) && scal_func::finite(imag(x));
    };
    force_inline bool finite(const Float_complex &x)
    {
        return scal_func::finite(real(x)) && scal_func::finite(imag(x));
    };

    //--------------------------------------------------------------------
    force_inline bool isinf(const Complex& x)
    {
        return (scal_func::isinf(real(x)) || scal_func::isinf(imag(x))) ? 1 : 0;
    };
    force_inline bool isinf(const Float_complex& x)
    {
        return (scal_func::isinf(real(x)) || scal_func::isinf(imag(x))) ? 1 : 0;
    };

    //--------------------------------------------------------------------
    force_inline bool isregular(const Float_complex &x)
    {
        return scal_func::finite(real(x)) && scal_func::finite(imag(x)) 
                    && mrd::is_zero(x) == false;
    }
    force_inline bool isregular(const Complex &x)
    {
        return scal_func::finite(real(x)) && scal_func::finite(imag(x)) 
                    && mrd::is_zero(x) == false;
    }

    //--------------------------------------------------------------------
    force_inline bool is_int(const Float_complex &x)
    {
        return imag(x) == 0 && scal_func::is_int(real(x));
    }
    force_inline bool is_int(const Complex &x)
    {
        return imag(x) == 0 && scal_func::is_int(real(x));
    }

    //--------------------------------------------------------------------
    force_inline bool is_real(const Float_complex &x)
    {
        return imag(x) == 0;
    }
    force_inline bool is_real(const Complex &x)
    {
        return imag(x) == 0;
    }

    //--------------------------------------------------------------------
    force_inline double nextabove(const Complex& x)
    {
        return scal_func::nextabove(real(x));
    }
    force_inline float nextabove(const Float_complex& x)
    {
        return scal_func::nextabove(real(x));
    }

    //--------------------------------------------------------------------
    force_inline double nextbelow(const Complex& x)
    {
        double re = scal_func::nextbelow(real(x));
        return re;
    }
    force_inline float nextbelow(const Float_complex& x)
    {
        float re = scal_func::nextbelow(real(x));
        return re;
    }

    //--------------------------------------------------------------------
    force_inline double abs(const Complex& arg)		
    {	
        //in VS hypot is slightly faster than std::abs but less accurate
        //(around 1.5ulp to 1.0 ulp)
        return std::abs(arg.value);
    };
    force_inline float abs(const Float_complex& arg)		
    {	
        return std::abs(arg.value);
    };

    //--------------------------------------------------------------------
    force_inline double eps(const Complex &a) 
    { 
        return scal_func::eps(scal_func::abs(a));
    };
    force_inline float eps(const Float_complex &a) 
    { 
        return scal_func::eps(scal_func::abs(a));
    };

    //--------------------------------------------------------------------
    force_inline Complex sqrt(const Complex& arg)		
    {	
        return Complex(std::sqrt(arg.value));
    };
    force_inline Float_complex sqrt(const Float_complex& arg)		
    {	
        return Float_complex(std::sqrt(arg.value));
    };

    //--------------------------------------------------------------------
    force_inline double abs2(const Complex& arg)		
    {	
        return arg.value.real() * arg.value.real() 
                    + arg.value.imag() * arg.value.imag();	
    };
    force_inline float abs2(const Float_complex& arg)		
    {	
        return arg.value.real() * arg.value.real() 
                    + arg.value.imag() * arg.value.imag();	
    };

    //--------------------------------------------------------------------
    force_inline Complex conj(const Complex& arg)		
    {	
        return Complex(real(arg), -imag(arg));
    };
    force_inline Float_complex conj(const Float_complex& arg)		
    {	
        return Float_complex(real(arg), -imag(arg));
    };

    //--------------------------------------------------------------------
    force_inline double arg(const Complex &a) 
    { 
        return std::arg(a.value); 
    };
    force_inline float arg(const Float_complex &a) 
    { 
        return std::arg(a.value); 
    };

    //--------------------------------------------------------------------
    force_inline Complex log(const Complex& arg)		
    {	
        return complex_log<Complex>::eval(arg);
    };
    force_inline Float_complex log(const Float_complex& arg)		
    {	
        return complex_log<Complex>::eval(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex log10(const Complex& a)		
    {	
        return matcl::details::mul_c(scal_func::log(a), constants::log10e());
    };
    force_inline Float_complex log10(const Float_complex& a)		
    {	
        return matcl::details::mul_c(scal_func::log(a), constants::f_log10e());
    };

    //--------------------------------------------------------------------
    force_inline Complex log2(const Complex &a)
    {
        return matcl::details::mul_c(scal_func::log(a), constants::log2e());
    };
    force_inline Float_complex log2(const Float_complex &a)
    {
        return matcl::details::mul_c(scal_func::log(a), constants::f_log2e());
    };

    //--------------------------------------------------------------------
    force_inline Complex scal_func::log1p(const Complex& arg)		
    {	
        return complex_log<Complex>::eval_1p(arg);
    };
    force_inline Float_complex scal_func::log1p(const Float_complex& arg)		
    {	
        return complex_log<Float_complex>::eval_1p(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex tan(const Complex& arg)		
    {	
        if (is_zero(imag(arg)) == true)
            return scal_func::tan(real(arg));

        return Complex(std::tan(arg.value));
    };
    force_inline Float_complex tan(const Float_complex& arg)		
    {	
        if (is_zero(imag(arg)) == true)
            return scal_func::tan(real(arg));

        return Float_complex(std::tan(arg.value));
    };

    //--------------------------------------------------------------------
    force_inline Complex cos(const Complex& arg)		
    {	
        return complex_cos<Complex>::eval(arg);
    };
    force_inline Float_complex cos(const Float_complex& arg)		
    {	
        return complex_cos<Float_complex>::eval(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex sec(const Complex &x) 
    { 
        return matcl::details::inv_c(scal_func::cos(x));
    };
    force_inline Float_complex sec(const Float_complex &x) 
    { 
        return matcl::details::inv_c(scal_func::cos(x));
    };

    //--------------------------------------------------------------------
    force_inline Complex sin(const Complex& arg)		
    {	
        return complex_sin<Complex>::eval(arg);
    };
    force_inline Float_complex sin(const Float_complex& arg)		
    {	
        return complex_sin<Float_complex>::eval(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex cot(const Complex &x) 
    { 
        return matcl::details::inv_c(scal_func::tan(x)); 
    };
    force_inline Float_complex cot(const Float_complex &x) 
    { 
        return matcl::details::inv_c(scal_func::tan(x)); 
    };

    //--------------------------------------------------------------------
    force_inline Complex exp(const Complex& arg)		
    {
        return complex_exp<Complex>::eval(arg);
    };
    force_inline Float_complex exp(const Float_complex& arg)		
    {
        return complex_exp<Float_complex>::eval(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex expm1(const Complex& arg)		
    {
        return complex_exp<Complex>::eval_m1(arg);
    };
    force_inline Float_complex expm1(const Float_complex& arg)		
    {
        return complex_exp<Float_complex>::eval_m1(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex expi(const Complex& arg)		
    {
        return scal_func::exp(Complex(-imag(arg), real(arg)));
    };
    force_inline Float_complex expi(const Float_complex& arg)		
    {
        return scal_func::exp(Float_complex(-imag(arg), real(arg)));
    };

    //--------------------------------------------------------------------
    force_inline Complex exp2(const Complex& arg)		
    {
        double y_re    = real(arg);
        double y_im    = imag(arg);

        // a ^ (x + iy) = a ^ x * exp(i * log(a) * y)
        double v_log2   = constants::ln2();
        double exp_re   = scal_func::exp2(y_re);
        double mult_im  = y_im * v_log2;
        Complex exp_im  = scal_func::expi(mult_im);

        double r_re     = mult_zero(exp_re , real(exp_im));
        double r_im     = mult_zero(exp_re , imag(exp_im));
        Complex ret     = Complex(r_re, r_im);
        return ret;
    };
    force_inline Float_complex exp2(const Float_complex& arg)		
    {
        float y_re      = real(arg);
        float y_im      = imag(arg);

        // a ^ (x + iy) = a ^ x * exp(i * log(a) * y)
        float v_log2    = constants::f_ln2();
        float exp_re    = scal_func::exp2(y_re);
        float mult_im   = y_im * v_log2;

        Float_complex exp_im  = scal_func::expi(mult_im);
        float r_re      = mult_zero(exp_re , real(exp_im));
        float r_im      = mult_zero(exp_re , imag(exp_im));
        Float_complex ret = Float_complex(r_re, r_im);
        return ret;
    };

    //--------------------------------------------------------------------
    force_inline Complex exp10(const Complex& arg)		
    {
        double y_re    = real(arg);
        double y_im    = imag(arg);

        // a ^ (x + iy) = a ^ x * exp(i * log(a) * y)
        double v_log10  = constants::ln10();
        double exp_re   = scal_func::exp10(y_re);
        double mult_im  = y_im * v_log10;
        Complex exp_im  = scal_func::expi(mult_im);
        double r_re     = mult_zero(exp_re , real(exp_im));
        double r_im     = mult_zero(exp_re , imag(exp_im));
        Complex ret     = Complex(r_re, r_im);
        return ret;
    };
    force_inline Float_complex exp10(const Float_complex& arg)		
    {
        float y_re      = real(arg);
        float y_im      = imag(arg);

        // a ^ (x + iy) = a ^ x * exp(i * log(a) * y)
        float v_log10   = constants::f_ln10();
        float exp_re    = scal_func::exp10(y_re);
        float mult_im   = y_im * v_log10;

        Float_complex exp_im  = scal_func::expi(mult_im);
        float r_re          = mult_zero(exp_re , real(exp_im));
        float r_im          = mult_zero(exp_re , imag(exp_im));
        Float_complex ret   = Float_complex(r_re, r_im);
        return ret;
    };

    //--------------------------------------------------------------------
    force_inline Complex sign(const Complex &a) 
    {
        return complex_sign<Complex>::eval(a);
    };
    force_inline Float_complex sign(const Float_complex &a) 
    {
        return complex_sign<Float_complex>::eval(a);
    };

    //--------------------------------------------------------------------
    force_inline Complex floor(const Complex &a)
    {
        return Complex(scal_func::floor(matcl::real(a)), 
                       scal_func::floor(matcl::imag(a)));
    }
    force_inline Float_complex floor(const Float_complex &a)
    {
        return Float_complex(scal_func::floor(matcl::real(a)), 
                             scal_func::floor(matcl::imag(a)));
    }

    //--------------------------------------------------------------------
    force_inline Complex ceil(const Complex &a)
    {
        return Complex(scal_func::ceil(matcl::real(a)),
                       scal_func::ceil(matcl::imag(a)));
    }
    force_inline Float_complex ceil(const Float_complex &a)
    {
        return Float_complex(scal_func::ceil(matcl::real(a)),
                             scal_func::ceil(matcl::imag(a)));
    }

    //--------------------------------------------------------------------
    force_inline Complex trunc(const Complex &a)
    {
        return Complex(scal_func::trunc(matcl::real(a)),
                       scal_func::trunc(matcl::imag(a)));
    }
    force_inline Float_complex trunc(const Float_complex &a)
    {
        return Float_complex(scal_func::trunc(matcl::real(a)),
                             scal_func::trunc(matcl::imag(a)));
    }

    //--------------------------------------------------------------------
    force_inline Complex round(const Complex &a)
    {
        return Complex(scal_func::round(matcl::real(a)), 
                       scal_func::round(matcl::imag(a)));
    }
    force_inline Float_complex round(const Float_complex &a)
    {
        return Float_complex(scal_func::round(matcl::real(a)), 
                             scal_func::round(matcl::imag(a)));
    }

    //--------------------------------------------------------------------
    force_inline Complex csc(const Complex &x) 
    { 
        return matcl::details::inv_c(scal_func::sin(x));
    };
    force_inline Float_complex csc(const Float_complex &x) 
    { 
        return matcl::details::inv_c(scal_func::sin(x));
    };

    //--------------------------------------------------------------------
    force_inline Complex asec(const Complex &x) 
    { 
        return scal_func::acos(matcl::details::inv_c(x));
    };
    force_inline Float_complex asec(const Float_complex &x) 
    { 
        return scal_func::acos(matcl::details::inv_c(x));
    };

    //--------------------------------------------------------------------
    force_inline Complex acsc(const Complex &x) 
    { 
        return scal_func::asin(matcl::details::inv_c(x)); 
    };
    force_inline Float_complex acsc(const Float_complex &x) 
    { 
        return scal_func::asin(matcl::details::inv_c(x)); 
    };

    //--------------------------------------------------------------------
    force_inline Complex tanh(const Complex& arg)		
    {	
        if (is_zero(imag(arg)) == true)
            return scal_func::tanh(real(arg));

        return Complex(std::tanh(arg.value));
    };
    force_inline Float_complex tanh(const Float_complex& arg)		
    {	
        if (is_zero(imag(arg)) == true)
            return scal_func::tanh(real(arg));

        return Float_complex(std::tanh(arg.value));
    };

    //--------------------------------------------------------------------
    force_inline Complex coth(const Complex &x) 
    { 
        return matcl::details::inv_c(scal_func::tanh(x)); 
    };
    force_inline Float_complex coth(const Float_complex &x) 
    { 
        return matcl::details::inv_c(scal_func::tanh(x)); 
    };

    //--------------------------------------------------------------------
    force_inline Complex cosh(const Complex& arg)		
    {	
        return complex_cosh<Complex>::eval(arg);
    };
    force_inline Float_complex cosh(const Float_complex& arg)		
    {	
        return complex_cosh<Float_complex>::eval(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex sech(const Complex &x) 
    { 
        return matcl::details::inv_c(scal_func::cosh(x)); 
    }
    force_inline Float_complex sech(const Float_complex &x) 
    { 
        return matcl::details::inv_c(scal_func::cosh(x)); 
    }

    //--------------------------------------------------------------------
    force_inline Complex sinh(const Complex& arg)		
    {	
        return complex_sinh<Complex>::eval(arg);
    };
    force_inline Float_complex sinh(const Float_complex& arg)		
    {	
        return complex_sinh<Float_complex>::eval(arg);
    };

    //--------------------------------------------------------------------
    force_inline Complex csch(const Complex &x) 
    { 
        return matcl::details::inv_c(scal_func::sinh(x)); 
    }
    force_inline Float_complex csch(const Float_complex &x) 
    { 
        return matcl::details::inv_c(scal_func::sinh(x)); 
    }

    //--------------------------------------------------------------------
    force_inline Complex acot(const Complex &x)
    {
        Complex tmp = matcl::details::div_c(1.,x);
        return scal_func::atan(tmp);
    };
    force_inline Float_complex acot(const Float_complex &x)
    {
        Float_complex tmp = matcl::details::div_c(1.f,x);
        return scal_func::atan(tmp);
    };

    //--------------------------------------------------------------------
    force_inline Complex acoth(const Complex &x)
    {
        Complex tmp1 = matcl::details::div_c(1.,x);
        return scal_func::atanh(tmp1);        
    }
    force_inline Float_complex acoth(const Float_complex &x)
    {
        Float_complex tmp1 = matcl::details::div_c(1.f,x);
        return scal_func::atanh(tmp1);        
    }

    //--------------------------------------------------------------------
    force_inline Complex asech(const Complex &x)
    {
        return scal_func::acosh(matcl::details::inv_c(x));
    }
    force_inline Float_complex asech(const Float_complex &x)
    {
        return scal_func::acosh(matcl::details::inv_c(x));
    }

    //--------------------------------------------------------------------
    force_inline Complex acsch(const Complex &x)
    {
        return scal_func::asinh(matcl::details::inv_c(x));
    }
    force_inline Float_complex acsch(const Float_complex &x)
    {
        return scal_func::asinh(matcl::details::inv_c(x));
    }
    
    //--------------------------------------------------------------------
    force_inline bool isnormal(const Complex& x)
    {
        //for x = Complex(a, 0) we must have: isnormal(x) == isnormal(a)
        bool is_norm_1  = scal_func::isnormal(real(x));
        bool is_norm_2  = scal_func::isnormal(imag(x));
        bool is_zero_1  = mrd::is_zero(real(x));
        bool is_zero_2  = mrd::is_zero(imag(x));

        return is_norm_1 && (is_norm_2 || is_zero_2)
                || (is_norm_1 || is_zero_1) && is_norm_2;
    }
    force_inline bool isnormal(const Float_complex& x)
    {
        //for x = Complex(a, 0) we must have: isnormal(x) == isnormal(a)
        bool is_norm_1  = scal_func::isnormal(real(x));
        bool is_norm_2  = scal_func::isnormal(imag(x));
        bool is_zero_1  = mrd::is_zero(real(x));
        bool is_zero_2  = mrd::is_zero(imag(x));

        return is_norm_1 && (is_norm_2 || is_zero_2)
                || (is_norm_1 || is_zero_1) && is_norm_2;
    }
}

}}}
