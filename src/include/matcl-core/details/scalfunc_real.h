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

#include "matcl-core/config.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-core/general/machine.h"

#if MATCL_ARCHITECTURE_HAS_FMA
    #include "matcl-core/details/simd/simd.h"
#endif

#include <cfloat>
#include <fenv.h> 

namespace matcl { namespace raw { namespace details
{

namespace md    = matcl::details;

template<class T>
struct is_zero_helper;

template<>
struct is_zero_helper<Integer>
{
    force_inline
    static bool eval(Integer val)       { return val == 0; };
};

template<>
struct is_zero_helper<Float>
{
    force_inline
    static bool eval(Float val)         { return val == 0.0f; };
};

template<>
struct is_zero_helper<Real>
{
    force_inline
    static bool eval(Real val)          { return val == 0.0; };
};

template<class T>
struct is_one_helper;

template<>
struct is_one_helper<Integer>
{
    force_inline
    static bool eval(Integer val)       { return val == 1; };
};

template<>
struct is_one_helper<Float>
{
    force_inline
    static bool eval(Float val)         { return val == 1.0f; };
};

template<>
struct is_one_helper<Real>
{
    force_inline
    static bool eval(Real val)          { return val == 1.0; };
};

//--------------------------------------------------------------------
template<class T> 
force_inline bool is_zero(const T& val)
{
    return is_zero_helper<T>::eval(val);
};

template<class T> force_inline
bool is_one(const T& val)
{ 
    return is_one_helper<T>::eval(val);
};

namespace scal_func
{
    namespace mrd = matcl::raw::details;

    double MATCL_CORE_EXPORT         sqrt1pm1(double arg);
    float MATCL_CORE_EXPORT          sqrt1pm1(float arg);
    Complex MATCL_CORE_EXPORT        sqrt1pm1(const Complex& arg);
    Float_complex MATCL_CORE_EXPORT  sqrt1pm1(const Float_complex& arg);

    double MATCL_CORE_EXPORT         asinh(double x);
    float MATCL_CORE_EXPORT          asinh(float x);
    double MATCL_CORE_EXPORT         acosh(double x);
    float MATCL_CORE_EXPORT          acosh(float x);
    double MATCL_CORE_EXPORT         atanh(double x);
    float MATCL_CORE_EXPORT          atanh(float x);
    Float_complex MATCL_CORE_EXPORT  asin(const Float_complex&);
    Complex MATCL_CORE_EXPORT        asin(const Complex&);
    Float_complex MATCL_CORE_EXPORT  acos(const Float_complex&);
    Complex MATCL_CORE_EXPORT        acos(const Complex&);
    Float_complex MATCL_CORE_EXPORT  atan(const Float_complex&);
    Complex MATCL_CORE_EXPORT        atan(const Complex&);
    Float_complex MATCL_CORE_EXPORT  asinh(const Float_complex&);
    Complex MATCL_CORE_EXPORT        asinh(const Complex&);
    Complex MATCL_CORE_EXPORT        atanh(const Complex& x);
    Float_complex MATCL_CORE_EXPORT  atanh(const Float_complex& x);
    Complex MATCL_CORE_EXPORT        acosh(const Complex& x);
    Float_complex MATCL_CORE_EXPORT  acosh(const Float_complex& x);

    float MATCL_CORE_EXPORT          fma_dekker(float x, float y, float z);
    double MATCL_CORE_EXPORT         fma_dekker(double x, double y, double z);

    float MATCL_CORE_EXPORT         float_distance(float x, float y);
    double MATCL_CORE_EXPORT        float_distance(double x, double y);

    Float_complex               log1p(const Float_complex&);
    Complex                     log1p(const Complex&);
    Float_complex               log(const Float_complex&);
    Complex                     log(const Complex&);
    Float_complex               log10(const Float_complex&);
    Complex                     log10(const Complex&);
    Float_complex               log2(const Float_complex&);
    Complex                     log2(const Complex&);
    Float_complex               sqrt(const Float_complex&);
    Complex                     sqrt(const Complex&);
    Float_complex               asec(const Float_complex&);
    Complex                     asec(const Complex&);
    Float_complex               acsc(const Float_complex&);
    Complex                     acsc(const Complex&);
    Float_complex               acoth(const Float_complex&);
    Complex                     acoth(const Complex&);
    Float_complex               asech(const Float_complex&);
    Complex                     asech(const Complex&);

    //--------------------------------------------------------------------
    //                      test range
    //--------------------------------------------------------------------
    template<class T>
    force_inline
    bool test_range_sqrt(const T& A)
    {
        if (A >= 0)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_sqrt1pm1(const T& A)
    {
        if (A >= -1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_log(const T& A)
    {
        if (A >= 0)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_log1p(const T& A)
    {
        if (A >= -1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_log2(const T& A)
    {
        if (A >= 0)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_log10(const T& A)
    {
        if (A >= 0)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_asin(const T& A)
    {
        if (A >= -1 && A <= 1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_acos(const T& A)
    {
        if (A >= -1 && A <= 1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_asec(const T& A)
    {
        if (A <= -1 || A >= 1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_acsc(const T& A)
    {
        if (A <= -1 || A >= 1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_acosh(const T& A)
    {
        if (A >= 1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_atanh(const T& A)
    {
        if (A >= -1 && A <= 1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_acoth(const T& A)
    {
        if (A <= -1 || A >= 1)
            return true;
        else
            return false;
    };

    template<class T>
    force_inline
    bool test_range_asech(const T& A)
    {
        if (A > 0 && A <= 1)
            return true;
        else
            return false;
    };

    //--------------------------------------------------------------------
    force_inline double abs(double x) 
    { 
        return std::abs(x);
    };
    force_inline float abs(float x) 
    { 
        return std::abs(x);
    };
    force_inline Integer abs(Integer x) 
    { 
        return x < 0 ? - x : x; 
    };

    //--------------------------------------------------------------------
    force_inline double abs2(double x) 
    { 
        return x*x; 
    };
    force_inline float abs2(float x) 
    { 
        return x*x; 
    };
    force_inline Integer abs2(Integer x) 
    { 
        return x*x; 
    };

    //--------------------------------------------------------------------
    force_inline bool isnan(Integer)
    {
        return false;
    };
    force_inline bool isnan(double x)
    {
        return x != x;
        //return (std::isnan(x)) ? 1 : 0;
    };
    force_inline bool isnan(float x)
    {
        return x != x;
        //return (std::isnan(x)) ? 1 : 0;
    };

    //--------------------------------------------------------------------
    force_inline bool finite(double x)
    {
        return x-x == 0;
        //std::isfinite is incredibly slow (dll function call)
        //return (std::isfinite(x)) ? 1 : 0;
    };
    force_inline bool finite(float x)
    {
        return x-x == 0;
        //std::isfinite is incredibly slow (dll function call)
        //return (std::isfinite(x)) ? 1 : 0;
    };

    //--------------------------------------------------------------------
    force_inline bool isinf(double x)
    {
        return (scal_func::isnan(x) || scal_func::finite(x)) ? 0 : 1;
    };
    force_inline bool isinf(float x)
    {
        return (scal_func::isnan(x) || scal_func::finite(x)) ? 0 : 1;
    };

    //--------------------------------------------------------------------
    force_inline bool isregular(double x)
    {
        return scal_func::finite(x) == true && x != 0.0;
    };
    force_inline bool isregular(float x)
    {
        return scal_func::finite(x) == true && x != 0.0f;
    };

    //--------------------------------------------------------------------
    force_inline bool is_int(double x)
    {
        return x == (Integer)x;
    };
    force_inline bool is_int(float x)
    {
        return x == (Integer)x;
    };

    //--------------------------------------------------------------------
    force_inline bool is_real(double)
    {
        return true;
    };
    force_inline bool is_real(float)
    {
        return true;
    };

    //--------------------------------------------------------------------
    force_inline double copysign(double x, double y)
    {
        return std::copysign(x,y);
    }
    force_inline double copysign(double x, float y)
    {
        return std::copysign(x,y);
    }
    force_inline float copysign(float x, float y)
    {
        return std::copysign(x,y);
    }
    force_inline double copysign(float x, double y)
    {
        return std::copysign(x,y);
    }

    //--------------------------------------------------------------------
    force_inline double nextabove(double x)
    {
        return std::nextafter(x,constants::inf());
    }
    force_inline float nextabove(float x)
    {
        return std::nextafter(x,constants::f_inf());
    }

    //--------------------------------------------------------------------
    force_inline double nextbelow(double x)
    {
        return std::nextafter(x,-constants::inf());
    }
    force_inline float nextbelow(float x)
    {
        return std::nextafter(x,-constants::f_inf());
    }

    //--------------------------------------------------------------------
    force_inline double eps(double x)
    {
        if (x < 0)
            return scal_func::nextabove(-x) + x;
        else
            return scal_func::nextabove( x) - x;
    }
    force_inline double eps(Integer a)
    {
        return scal_func::eps(Real(a));
    }
    force_inline float eps(float x)
    {
        if (x < 0)
            return scal_func::nextabove(-x) + x;
        else
            return scal_func::nextabove( x) - x;
    }

    //--------------------------------------------------------------------
    force_inline bool signbit(Integer x)
    {
        return (x < 0);
    };
    force_inline bool signbit(const float& x)
    {
        //return std::signbit(x);
        return reinterpret_cast<const double_binary_rep*>(&x)->sign != 0;
    };
    force_inline bool signbit(const double& x)
    {
        //return std::signbit(x);
        return reinterpret_cast<const double_binary_rep*>(&x)->sign != 0;
    };

    //--------------------------------------------------------------------
    force_inline Float sqrt(Float arg)		
    {
        return std::sqrt(arg);
    };
    force_inline Real sqrt(Real arg)		
    {
        return std::sqrt(arg);
    };

    //--------------------------------------------------------------------
    force_inline Float_complex sqrt_c(Float arg)		
    {
        if (test_range_sqrt(arg))
            return Float_complex(scal_func::sqrt(arg));
        else
            return scal_func::sqrt(Float_complex(arg));
    };
    force_inline Complex sqrt_c(Real arg)		
    {
        if (test_range_sqrt(arg))
            return Complex(scal_func::sqrt(arg));
        else
            return scal_func::sqrt(Complex(arg));
    };

    //--------------------------------------------------------------------
    force_inline double cbrt(double arg)
    {
        return std::cbrt(arg);
    }
    force_inline float cbrt(float arg)
    {
        return std::cbrt(arg);
    }

    //--------------------------------------------------------------------
    force_inline Float_complex sqrt1pm1_c(float arg)		
    {
        if (test_range_sqrt1pm1(arg))
            return Float_complex(scal_func::sqrt1pm1(arg));
        else
            return scal_func::sqrt1pm1(Float_complex(arg));
    };
    force_inline Complex sqrt1pm1_c(double arg)		
    {
        if (test_range_sqrt1pm1(arg))
            return Complex(scal_func::sqrt1pm1(arg));
        else
            return scal_func::sqrt1pm1(Complex(arg));
    };

    //--------------------------------------------------------------------
    force_inline double arg(Integer a)
    {
        return (a >= 0) ? 0.0 : constants::pi();
    }
    force_inline double arg(double a)
    {
        if (scal_func::isnan(a))
            return a;

        return scal_func::signbit(a) == false ? 0. : constants::pi();
    }
    force_inline float arg(float a)
    {
        if (scal_func::isnan(a))
            return a;

        return scal_func::signbit(a) == false ? 0.f : constants::f_pi();
    }

    //--------------------------------------------------------------------
    force_inline Float log(Float arg)		
    {
        return std::log(arg);
    };
    force_inline Real log(Real arg)		
    {
        return std::log(arg);
    };

    //--------------------------------------------------------------------
    force_inline Float_complex log_c(Float arg)		
    {
        Float re    = scal_func::log(scal_func::abs(arg));
        Float im    = arg >= 0 ? 0.0f : constants::f_pi();
        return Float_complex(re, im);
    };
    force_inline Complex log_c(Real arg)		
    {
        Real re     = scal_func::log(scal_func::abs(arg));
        Real im     = arg >= 0 ? 0.0f : constants::pi();

        return Complex(re, im);
    };

    //--------------------------------------------------------------------
    force_inline double log1p(double arg)        
    {    
        return std::log1p(arg);
    };
    force_inline float log1p(float arg)        
    {    
        return std::log1p(arg);
    };

    //--------------------------------------------------------------------
    force_inline Float_complex log1p_c(Float arg)		
    {
        if (test_range_log1p(arg))
            return Float_complex(scal_func::log1p(arg));
        else
            return log(Float_complex(arg + 1.f));
    };
    force_inline Complex log1p_c(Real arg)		
    {
        if (test_range_log1p(arg))
            return Complex(scal_func::log1p(arg));
        else
            return log(Complex(arg + 1));
    };

    //--------------------------------------------------------------------
    force_inline Float log10(Float arg)		
    {
        return std::log10(arg);
    };
    force_inline Real log10(Real arg)		
    {
        return std::log10(arg);
    };

    //--------------------------------------------------------------------
    force_inline Float_complex log10_c(Float arg)		
    {
        Float re    = scal_func::log10(scal_func::abs(arg));
        Float im    = arg >= 0 ? 0.0f : constants::f_pi() * constants::f_log10e();

        return Float_complex(re, im);
    };
    force_inline Complex log10_c(Real arg)		
    {
        Real re     = scal_func::log10(scal_func::abs(arg));
        Real im     = arg >= 0 ? 0.0f : constants::pi() * constants::log10e();

        return Complex(re, im);
    };

    //--------------------------------------------------------------------
    force_inline double log2(double a) 
    { 
        return std::log2(a);
    };
    force_inline float log2(float a) 
    { 
        return std::log2(a);
    };

    //--------------------------------------------------------------------
    force_inline Complex log2_c(double arg) 
    { 
        Real re     = scal_func::log2(scal_func::abs(arg));
        Real im     = arg >= 0 ? 0.0f : constants::pi() * constants::log2e();

        return Complex(re, im);
    };
    force_inline Float_complex log2_c(float arg) 
    { 
        Float re    = scal_func::log2(scal_func::abs(arg));
        Float im    = arg >= 0 ? 0.0f : constants::f_pi() * constants::f_log2e();

        return Float_complex(re, im);
    };

    //--------------------------------------------------------------------
    force_inline double tan(double x)
    {
        return std::tan(x);
    };
    force_inline float tan(float x)
    {
        return std::tan(x);
    };

    //--------------------------------------------------------------------
    force_inline double cos(double x)
    {
        return std::cos(x);
    };
    force_inline float cos(float x)
    {
        return std::cos(x);
    };

    //--------------------------------------------------------------------
    force_inline double sec(double x) 
    { 
        return 1. / std::cos(x); 
    };
    force_inline float sec(float x) 
    { 
        return 1.f / std::cos(x); 
    };

    //--------------------------------------------------------------------
    force_inline double sin(double x)
    {
        return std::sin(x);
    };
    force_inline float sin(float x)
    {
        return std::sin(x);
    };

    //--------------------------------------------------------------------
    force_inline void sin_cos(double x, double& v_sin, double& v_cos)
    {
        v_sin = scal_func::sin(x);
        v_cos = scal_func::cos(x);
    };
    force_inline void sin_cos(float x, float& v_sin, float& v_cos)
    {
        v_sin = scal_func::sin(x);
        v_cos = scal_func::cos(x);
    };

    //--------------------------------------------------------------------
    force_inline double cot(double x)
    {
        return 1.0/ scal_func::tan(x);
    };
    force_inline float cot(float x)
    {
        return 1.0f / scal_func::tan(x);
    };

    //--------------------------------------------------------------------
    force_inline Float exp(Float arg)		
    {
        return std::exp(arg);
    };
    force_inline Real exp(Real arg)		
    {
        return std::exp(arg);
    };

    //--------------------------------------------------------------------
    force_inline Float expm1(Float arg)		
    {
        return std::expm1(arg);
    };
    force_inline Real expm1(Real arg)		
    {
        return std::expm1(arg);
    };

    //--------------------------------------------------------------------
    force_inline Float_complex expi(Float arg)		
    {
        Float si, co;
        scal_func::sin_cos(arg, si, co);
        return Float_complex(co, si);
    };
    force_inline Complex expi(Real arg)		
    {
        Real si, co;
        scal_func::sin_cos(arg, si, co);
        return Complex(co, si);
    };

    //--------------------------------------------------------------------
    force_inline double exp2(Integer x)
    {
        return std::exp2(x);
    };
    force_inline double exp2(double x)
    {
        return std::exp2(x);
    };
    force_inline float exp2(float x)
    {
        return std::exp2(x);
    };

    //--------------------------------------------------------------------
    force_inline double exp10(Integer x)
    {
        return std::pow(10.0, Real(x));
    };
    force_inline double exp10(double x)
    {
        return std::pow(10.0, x);
    };
    force_inline float exp10(float x)
    {
        return std::pow(10.0f, x);
    };

    //--------------------------------------------------------------------
    force_inline double sign(double a)
    {
        //if (scal_func::isnan(a)) return a;
        if (a > 0)      return 1.;
        else if (a < 0) return -1.;
        else            return a;
    };
    force_inline float sign(float a)
    {
        //if (scal_func::isnan(a)) return a;
        if (a > 0.f)        return 1.f;
        else if (a < 0.f)   return -1.f;
        else                return a;
    };
    force_inline Integer sign(Integer a)
    {
        if (a > 0) return 1;
        if (a < 0) return -1;
        return 0;
    };

    //--------------------------------------------------------------------
    force_inline Integer isign(double x)
    {
        if (is_zero(x))                 return 0;
        else if (scal_func::signbit(x)) return -1;
        else                            return 1;
    };
    force_inline Integer isign(float x)
    {
        if (is_zero(x))                 return 0;
        else if (scal_func::signbit(x)) return -1;
        else                            return 1;
    };
    force_inline Integer isign(Integer a)
    {
        if (a > 0) return 1;
        if (a < 0) return -1;
        return 0;
    };

    //--------------------------------------------------------------------
    force_inline double floor(double a) 
    {
        return std::floor(a);
    }
    force_inline float floor(float a) 
    {
        return std::floor(a);
    }

    //--------------------------------------------------------------------
    force_inline Integer ifloor(Integer a)
    {
        return a;
    };
    force_inline Integer ifloor(double a)
    {
        return (Integer)floor(a);
    };
    force_inline Integer ifloor(float a)
    {
        return (Integer)floor(a);
    };

    //--------------------------------------------------------------------
    force_inline double ceil(double a) 
    {
        return std::ceil(a);
    }
    force_inline float ceil(float a) 
    {
        return std::ceil(a);
    }

    //--------------------------------------------------------------------
    force_inline Integer iceil(Integer a)
    {
        return a;
    };
    force_inline Integer iceil(double a)
    {
        return (Integer)ceil(a);
    };
    force_inline Integer iceil(float a)
    {
        return (Integer)ceil(a);
    };

    //--------------------------------------------------------------------
    force_inline double trunc(double a) 
    {
        return std::trunc(a);
    }
    force_inline float trunc(float a) 
    {
        return std::trunc(a);
    }

    //--------------------------------------------------------------------
    force_inline Integer itrunc(Integer a)
    {
        return a;
    };
    force_inline Integer itrunc(double a)
    {
        return (Integer)trunc(a);
    };
    force_inline Integer itrunc(float a)
    {
        return (Integer)trunc(a);
    };

    //--------------------------------------------------------------------
    force_inline double round(double a) 
    { 
        ::fesetround(FE_TONEAREST);
        return std::nearbyint(a);
    }
    force_inline float round(float a) 
    { 
        ::fesetround(FE_TONEAREST);
        return std::nearbyint(a);
    }

    //--------------------------------------------------------------------
    force_inline Integer iround(Integer a)
    {
        return a;
    };
    force_inline Integer iround(double a)
    {
        return (Integer)round(a);
    };
    force_inline Integer iround(float a)
    {
        return (Integer)round(a);
    };

    //--------------------------------------------------------------------
    force_inline double csc(double x) 
    { 
        return 1. / scal_func::sin(x); 
    };
    force_inline float csc(float x) 
    { 
        return 1.f / scal_func::sin(x); 
    };

    //--------------------------------------------------------------------
    force_inline double asin(double cm)
    {
        return std::asin(cm);
    };
    force_inline float asin(float cm)
    {
        return std::asin(cm);
    };

    //--------------------------------------------------------------------
    force_inline Complex asin_c(double arg)
    {
        if (test_range_asin(arg))
            return Complex(scal_func::asin(arg));
        else
            return scal_func::asin(Complex(arg));
    };
    force_inline Float_complex asin_c(float arg)
    {
        if (test_range_asin(arg))
            return Float_complex(scal_func::asin(arg));
        else
            return scal_func::asin(Float_complex(arg));
    };

    //--------------------------------------------------------------------
    force_inline double acos(double cm)
    {
        return std::acos(cm);
    };
    force_inline float acos(float cm)
    {
        return std::acos(cm);
    };

    //--------------------------------------------------------------------
    force_inline Complex acos_c(double arg)
    {
        if (test_range_acos(arg))
            return Complex(scal_func::acos(arg));
        else
            return scal_func::acos(Complex(arg));
    };
    force_inline Float_complex acos_c(float arg)
    {
        if (test_range_acos(arg))
            return Float_complex(scal_func::acos(arg));
        else
            return scal_func::acos(Float_complex(arg));
    };

    //--------------------------------------------------------------------
    force_inline double atan(double x)
    {
        return std::atan(x);
    };
    force_inline float atan(float x)
    {
        return std::atan(x);
    };

    //--------------------------------------------------------------------
    force_inline double acot(double x) 
    { 
        return scal_func::atan(1./x);
    };
    force_inline float acot(float x) 
    { 
        return scal_func::atan(1.f/x);
    };

    //--------------------------------------------------------------------
    force_inline double asec(double x) 
    { 
        return std::acos(1./x); 
    };
    force_inline float asec(float x) 
    { 
        return std::acos(1.f/x); 
    };

    //--------------------------------------------------------------------
    force_inline Complex asec_c(double arg) 
    { 
        if (test_range_asec(arg))
            return Complex(scal_func::asec(arg));
        else
            return scal_func::asec(Complex(arg));
    };
    force_inline Float_complex asec_c(float arg) 
    { 
        if (test_range_asec(arg))
            return Float_complex(scal_func::asec(arg));
        else
            return scal_func::asec(Float_complex(arg));
    };

    //--------------------------------------------------------------------
    force_inline double acsc(double x) 
    { 
        return std::asin(1./x); 
    };
    force_inline float acsc(float x) 
    { 
        return std::asin(1.f/x); 
    };

    //--------------------------------------------------------------------
    force_inline Complex acsc_c(double arg) 
    { 
        if (test_range_acsc(arg))
            return Complex(scal_func::acsc(arg));
        else
            return scal_func::acsc(Complex(arg));
    };
    force_inline Float_complex acsc_c(float arg) 
    { 
        if (test_range_acsc(arg))
            return Float_complex(scal_func::acsc(arg));
        else
            return scal_func::acsc(Float_complex(arg));
    };

    //--------------------------------------------------------------------
    force_inline Float tanh(Float arg)		
    {
        return std::tanh(arg);
    };
    force_inline Real tanh(Real arg)		
    {
        return std::tanh(arg);
    };

    //--------------------------------------------------------------------
    force_inline double coth(double x) 
    { 
        return 1./std::tanh(x); 
    };
    force_inline float coth(float x) 
    { 
        return 1.f/std::tanh(x); 
    };

    //--------------------------------------------------------------------
    force_inline double cosh(double x)		
    {	
        return std::cosh(x);
    };
    force_inline float cosh(float x)		
    {	
        return std::cosh(x);
    };

    //--------------------------------------------------------------------
    force_inline double sech(double x) 
    { 
        return 1./std::cosh(x); 
    };
    force_inline float sech(float x) 
    { 
        return 1.f/std::cosh(x); 
    };

    //--------------------------------------------------------------------
    force_inline double sinh(double x) 
    { 
        return std::sinh(x); 
    };
    force_inline float sinh(float x) 
    { 
        return std::sinh(x); 
    };

    //--------------------------------------------------------------------
    force_inline void sinh_cosh(double x, double& sinh, double& cosh)
    {
        sinh = scal_func::sinh(x);
        cosh = scal_func::cosh(x);
    }
    force_inline void sinh_cosh(float x, float& sinh, float& cosh)
    {
        sinh = scal_func::sinh(x);
        cosh = scal_func::cosh(x);
    }

    //--------------------------------------------------------------------
    force_inline double csch(double x) 
    { 
        return 1./std::sinh(x); 
    }
    force_inline float csch(float x) 
    { 
        return 1.f/std::sinh(x); 
    }

    //--------------------------------------------------------------------
    force_inline double acoth(double x)
    {
        return scal_func::atanh(1./x);
    }
    force_inline float acoth(float x)
    {
        return scal_func::atanh(1.f/x);
    }

    //--------------------------------------------------------------------
    force_inline Complex acoth_c(double arg)
    {
        if (test_range_acoth(arg))
            return Complex(scal_func::acoth(arg));
        else
            return scal_func::acoth(Complex(arg));
    }
    force_inline Float_complex acoth_c(float arg)
    {
        if (test_range_acoth(arg))
            return Float_complex(scal_func::acoth(arg));
        else
            return scal_func::acoth(Float_complex(arg));
    }

    //--------------------------------------------------------------------
    force_inline double asech(double x) 
    { 
        return scal_func::acosh(1./x); 
    }
    force_inline float asech(float x) 
    { 
        return scal_func::acosh(1.f/x); 
    }

    //--------------------------------------------------------------------
    force_inline Complex acosh_c(double arg) 
    { 
        if (test_range_acosh(arg))
            return Complex(scal_func::acosh(arg));
        else
            return scal_func::acosh(Complex(arg));
    }
    force_inline Float_complex acosh_c(float arg) 
    { 
        if (test_range_acosh(arg))
            return Float_complex(scal_func::acosh(arg));
        else
            return scal_func::acosh(Float_complex(arg));
    }

    //--------------------------------------------------------------------
    force_inline Complex atanh_c(double arg) 
    { 
        if (test_range_atanh(arg))
            return Complex(scal_func::atanh(arg));
        else
            return scal_func::atanh(Complex(arg));
    }
    force_inline Float_complex atanh_c(float arg) 
    { 
        if (test_range_atanh(arg))
            return Float_complex(scal_func::atanh(arg));
        else
            return scal_func::atanh(Float_complex(arg));
    }

    //--------------------------------------------------------------------
    force_inline Complex asech_c(double arg) 
    { 
        if (test_range_asech(arg))
            return Complex(scal_func::asech(arg));
        else
            return scal_func::asech(Complex(arg));
    }
    force_inline Float_complex asech_c(float arg) 
    { 
        if (test_range_asech(arg))
            return Float_complex(scal_func::asech(arg));
        else
            return scal_func::asech(Float_complex(arg));
    }

    //--------------------------------------------------------------------
    force_inline double acsch(double x) 
    { 
        return scal_func::asinh(1./x); 
    }
    force_inline float acsch(float x) 
    { 
        return scal_func::asinh(1.f/x); 
    }
    
    //--------------------------------------------------------------------
    force_inline bool isnormal(double x)
    {
        return std::isnormal(x);
    }
    force_inline bool isnormal(float x)
    {
        return std::isnormal(x);
    }
    force_inline bool isnormal(Integer x)
    {
        return std::isnormal(Real(x));
    };

    //--------------------------------------------------------------------
    force_inline Integer ilogb(Integer x)
    {
        return std::ilogb(x);
    };
    force_inline Integer ilogb(double x)
    {
        return std::ilogb(x);
    };
    force_inline Integer ilogb(float x)
    {
        return std::ilogb(x);
    };

    //--------------------------------------------------------------------
    force_inline double logb(Integer x)
    {
        return std::logb(x);
    };
    force_inline double logb(double x)
    {
        return std::logb(x);
    };
    force_inline float logb(float x)
    {
        return std::logb(x);
    };

    //--------------------------------------------------------------------
    force_inline double hypot(double x, double y)
    {
        return std::hypot(x,y);
    };
    force_inline float hypot(float x, float y)
    {
        return std::hypot(x,y);
    };

    //--------------------------------------------------------------------
    force_inline double atan2(double x, double y)
    {
        return std::atan2(x,y);
    };
    force_inline float atan2(float x, float y)
    {
        return std::atan2(x,y);
    };

    //--------------------------------------------------------------------
    force_inline double modf(double x, double& y)
    {
        return std::modf(x,&y);
    };
    force_inline float modf(float x, float& y)
    {
        return std::modf(x,&y);
    };

    //--------------------------------------------------------------------
    force_inline double frexp(double x, Integer& y)
    {
        return std::frexp(x,&y);
    };
    force_inline float frexp(float x, Integer& y)
    {
        return std::frexp(x,&y);
    };

    //--------------------------------------------------------------------
    force_inline double fma_f(double x, double y, double z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fma_f(x, y, z);
        #else
            return x * y + z;
        #endif
    };

    force_inline float fma_f(float x, float y, float z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fma_f(x, y, z);
        #else
            return x * y + z;
        #endif
    };

    //--------------------------------------------------------------------
    force_inline double fms_f(double x, double y, double z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fms_f(x, y, z);
        #else
            return x * y - z;
        #endif
    };
    force_inline float fms_f(float x, float y, float z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fms_f(x, y, z);
        #else
            return x * y - z;
        #endif
    };

    //--------------------------------------------------------------------
    force_inline double fma_a(double x, double y, double z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fma_f(x, y, z);
        #else
            return fma_dekker(x, y, z);
        #endif
    };

    force_inline float fma_a(float x, float y, float z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fma_f(x, y, z);
        #else
            return fma_dekker(x, y, z);
        #endif
    };

    //--------------------------------------------------------------------
    force_inline double fms_a(double x, double y, double z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fms_f(x, y, z);
        #else
            return fma_dekker(x, y, -z);
        #endif
    };
    force_inline float fms_a(float x, float y, float z)
    {
        // do not use std::fma, this function is incredibly slow on VS

        #if MATCL_ARCHITECTURE_HAS_FMA
            return simd::fms_f(x, y, z);
        #else
            return fma_dekker(x, y, -z);
        #endif
    };

    //--------------------------------------------------------------------
    template<class T>
    struct dot2_a_impl
    {
        static T eval(T a, T b, T c, T d)
        {
            //Kahan's algorithm

            //error is 2 ulp as shown in
            //C.-P. Jeannerod, N. Louvet, and J.-M. Muller, “Further analysis of
            //Kahan’s algorithm for the accurate computation of 2 × 2 determinants”

            T w = c * d;
            T e = scal_func::fms_a(c, d, w);  //e = c * d - w
            T f = scal_func::fma_a(a, b, w);  //f = a * b + w 
            T g = f + e;
            return g;
        };
    };

    force_inline double dot2_a(double a, double b, double c, double d)
    {
        return dot2_a_impl<double>::eval(a,b,c,d);
    };
    
    force_inline float dot2_a(float a, float b, float c, float d)
    {
        return dot2_a_impl<float>::eval(a,b,c,d);
    };

    //--------------------------------------------------------------------
    force_inline float pow(float x, float y)
    {
        return std::pow(x, y);
    }
    force_inline double pow(double x, double y)
    {
        return std::pow(x, y);
    }

    force_inline float pow(float x, int y)
    {
        return std::pow(x, y);
    }
    force_inline double pow(double x, int y)
    {
        return std::pow(x, y);
    }

    //--------------------------------------------------------------------
    force_inline float ldexp(float x, Integer y)
    {
        return std::ldexp(x, y);
    }
    force_inline double ldexp(double x, Integer y)
    {
        return std::ldexp(x, y);
    }

    //--------------------------------------------------------------------
    force_inline float scalbn(float x, Integer y)
    {
        return std::scalbn(x, y);
    }
    force_inline double scalbn(double x, Integer y)
    {
        return std::scalbn(x, y);
    }

    //--------------------------------------------------------------------
    inline fp_type int_to_fptype(int code)
    {
        switch(code)
        {
            case FP_INFINITE:   return fp_type::fp_infinite;
            case FP_NAN:        return fp_type::fp_nan;
            case FP_ZERO:       return fp_type::fp_zero;
            case FP_SUBNORMAL:  return fp_type::fp_subnormal;
            case FP_NORMAL:     return fp_type::fp_normal;
            default:
                return fp_type::fp_unknown;
        }
    };

    force_inline fp_type fpclassify(float x)
    {
        int ret = std::fpclassify(x);
        return int_to_fptype(ret);
    }
    force_inline fp_type fpclassify(double x)
    {
        int ret = std::fpclassify(x);
        return int_to_fptype(ret);
    }
}

}}}
