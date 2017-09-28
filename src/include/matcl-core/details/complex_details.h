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

#include "matcl-core/matrix/complex_type.h"
#include "matcl-core/details/scalfunc_real.h"

namespace matcl { namespace details
{

namespace mrd = matcl :: raw :: details;

template<class T>   
            struct base_scalar                      { };
template<>  struct base_scalar<bool>                { };
template<>  struct base_scalar<signed char>         { using type = double;   };
template<>  struct base_scalar<unsigned char>       { using type = double;   };
template<>  struct base_scalar<signed short>        { using type = double;   };
template<>  struct base_scalar<unsigned short>      { using type = double;   };
template<>  struct base_scalar<signed int>          { using type = double;   };
template<>  struct base_scalar<unsigned int>        { using type = double;   };
template<>  struct base_scalar<signed long>         { using type = double;   };
template<>  struct base_scalar<unsigned long>       { using type = double;   };
template<>  struct base_scalar<signed long long>    { using type = double;   };
template<>  struct base_scalar<unsigned long long>  { using type = double;   };
template<>  struct base_scalar<float>               { using type = float ;   };
template<>  struct base_scalar<double>              { using type = double;   };
template<>  struct base_scalar<long double>         { using type = double;   };

template<class T>
            struct base_scalar_code                 {};
template<>  struct base_scalar_code<float>          { static const int code = 1;};
template<>  struct base_scalar_code<double>         { static const int code = 2;};

template<int code>  
            struct code_base_scalar                 {};
template<>  struct code_base_scalar<1>              { using type = float ;   };
template<>  struct code_base_scalar<2>              { using type = double;   };

template<class T, class S>
struct max_type_base 
{
    static const int 
    T_code      = base_scalar_code<typename base_scalar<T>::type>::code;

    static const int 
    S_code      = base_scalar_code<typename base_scalar<S>::type>::code;
    
    static const int 
    M_code      = (T_code > S_code) ? T_code : S_code;

    using type  = typename code_base_scalar<M_code>::type;
};

template<class T, class S>
struct convert_complex
{};

template<class T, class S>
struct convert_complex<complex<T>, complex<S>>
{
    using type_1 = complex<T>;
    using type_2 = complex<S>;

    force_inline
    static type_1 eval(const type_2& arg)
    {
        return type_1(real(arg),imag(arg));
    };
};

template<class T>
struct convert_complex<complex<T>, complex<T>>
{
    using type_1 = complex<T>;
    using type_2 = complex<T>;

    force_inline
    static type_1 eval(const type_2& arg)
    {
        return arg;
    };
};

// implement unary minus for complex values
template<class T> struct uminus_helper{};

template<class T>
struct uminus_helper<complex<T>>
{
    using return_type = complex<T>;

    force_inline
    static return_type eval(const complex<T>& a)
    {
        T a_re  = real(a);
        T a_im  = imag(a);

        return return_type(-a_re, -a_im);
    };
};

// implement binary plus for complex values
template<class T, class S> 
struct plus_helper {};

template<class T, class S>
struct plus_helper<complex<T>,complex<S>>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const complex<T>& a, const complex<S>& b)
    {
        real_type a_re  = real_type(real(a));
        real_type a_im  = real_type(imag(a));
        real_type b_re  = real_type(real(b));
        real_type b_im  = real_type(imag(b));

        real_type r_re  = a_re + b_re;
        real_type r_im  = a_im + b_im;

        return return_type(r_re, r_im);
    };
};

template<class T, class S>
struct plus_helper<complex<T>,S>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const complex<T>& a, const S& b)
    {
        real_type a_re  = real_type(real(a));
        real_type a_im  = real_type(imag(a));
        real_type b_re  = real_type(b);

        real_type r_re  = a_re + b_re;
        real_type r_im  = a_im;

        return return_type(r_re, r_im);
    };
};

template<class T, class S>
struct plus_helper<T,complex<S>>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const T& a, const complex<S>& b)
    {
        real_type a_re  = real_type(a);
        real_type b_re  = real_type(real(b));
        real_type b_im  = real_type(imag(b));

        real_type r_re  = a_re + b_re;
        real_type r_im  = b_im;

        return return_type(r_re, r_im);
    };
};

// implement binary minus for complex values
template<class T, class S> 
struct minus_helper {};

template<class T, class S>
struct minus_helper<complex<T>,complex<S>>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const complex<T>& a, const complex<S>& b)
    {
        real_type a_re  = real_type(real(a));
        real_type a_im  = real_type(imag(a));
        real_type b_re  = real_type(real(b));
        real_type b_im  = real_type(imag(b));

        real_type r_re  = a_re - b_re;
        real_type r_im  = a_im - b_im;

        return return_type(r_re, r_im);
    };
};

template<class T, class S>
struct minus_helper<complex<T>,S>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const complex<T>& a, const S& b)
    {
        real_type a_re  = real_type(real(a));
        real_type a_im  = real_type(imag(a));
        real_type b_re  = real_type(b);

        real_type r_re  = a_re - b_re;
        real_type r_im  = a_im;

        return return_type(r_re, r_im);
    };
};

template<class T, class S>
struct minus_helper<T,complex<S>>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const T& a, const complex<S>& b)
    {
        real_type a_re  = real_type(a);
        real_type b_re  = real_type(real(b));
        real_type b_im  = real_type(imag(b));

        real_type r_re  = a_re - b_re;
        real_type r_im  = -b_im;

        return return_type(r_re, r_im);
    };
};

// implement multiplication for complex values
template<class T, class S> 
struct mul_helper {};

template<class T, class S>
struct mul_helper<complex<T>,complex<S>>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const complex<T>& a, const complex<S>& b)
    {
        real_type a_re  = real_type(real(a));        
        real_type b_re  = real_type(real(b));        
        real_type a_im  = real_type(imag(a));
        real_type b_im  = real_type(imag(b));

        namespace ms = mrd::scal_func;

        real_type r_re  = ms::fms(a_re, b_re, a_im * b_im);
        real_type r_im  = ms::fma(a_re, b_im, a_im * b_re);        

        if (ms::isnan(r_re) == true && ms::isnan(r_im) == true)
            return recover_nan(a, b, r_re, r_im);
        else
            return return_type(r_re, r_im);
    };

    static return_type recover_nan(const complex<T>& a, const complex<S>& b,
                                   real_type r_re, real_type r_im)
    {
        //code taken from C99 standard; ISO/IEC 9899:TC2; WG14/N1124
        //notice that correct inf/nan handling is not required neither
        //by C not by C++ standard

        bool recalc = false;

        real_type a_re  = real_type(real(a));        
        real_type b_re  = real_type(real(b));        
        real_type a_im  = real_type(imag(a));
        real_type b_im  = real_type(imag(b));

        namespace ms = mrd::scal_func;

        real_type one   = real_type(1.0);
        real_type zero  = real_type(0.0);

        if (ms::isinf(a_re) || ms::isinf(b_re) ) 
        { 
            // z is infinite
            a_re        = ms::copysign(ms::isinf(a_re) ? one : zero, a_re);
            a_im        = ms::copysign(ms::isinf(a_im) ? one : zero, a_im);

            if (ms::isnan(b_re)) 
                b_re    = ms::copysign(zero, b_re);
            if (ms::isnan(b_im)) 
                b_im    = ms::copysign(zero, b_im);

            recalc = true;
        }

        if (ms::isinf(b_re) || ms::isinf(b_im) ) 
        { 
            // w is infinite
            b_re        = ms::copysign(ms::isinf(b_re) ? one : zero, b_re);
            b_im        = ms::copysign(ms::isinf(b_im) ? one : zero, b_im);

            if (isnan(a_re)) 
                a_re    = ms::copysign(zero, a_re);
            if (isnan(a_im)) 
                a_im    = ms::copysign(zero, a_im);

            recalc = true;
        }

        if (recalc == false)
        {
            if (ms::isinf(a_re * b_re) || ms::isinf(a_im * b_im) 
                || ms::isinf(a_re * b_im) || ms::isinf(a_im * b_re))
            {
                // recover infinities from overflow by changing NaNs to 0
                if (ms::isnan(a_re)) 
                    a_re    = ms::copysign(zero, a_re);
                if (ms::isnan(a_im)) 
                    a_im    = ms::copysign(zero, a_im);
                if (ms::isnan(b_re))
                    b_re    = ms::copysign(zero, b_re);
                if (ms::isnan(b_im))
                    b_im    = ms::copysign(zero, b_im);

                recalc = true;
            }
        }

        if (recalc) 
        {
            r_re = constants::inf<real_type>() * ms::fms(a_re, b_re, - a_im * b_im );
            r_im = constants::inf<real_type>() * ms::fma(a_re, b_im,   a_im * b_re );
        }

        return return_type(r_re, r_im);
    };
};

template<class T, class S>
struct mul_helper<complex<T>,S>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const complex<T>& a, const S& b)
    {
        real_type a_re  = real_type(real(a));
        real_type a_im  = real_type(imag(a));
        real_type b_re  = real_type(b);

        real_type r_re  = a_re * b_re;
        real_type r_im  = a_im * b_re;

        return return_type(r_re, r_im);
    };
};

template<class T, class S>
struct mul_helper<T,complex<S>>
{
    using real_type     = typename max_type_base<T,S>::type;
    using return_type   = complex<real_type>;

    force_inline
    static return_type eval(const T& a, const complex<S>& b)
    {
        real_type a_re  = real_type(a);
        real_type b_re  = real_type(real(b));
        real_type b_im  = real_type(imag(b));

        real_type r_re  = a_re * b_re;
        real_type r_im  = a_re * b_im;

        return return_type(r_re, r_im);
    };
};

// implement division for complex values
template<class T, class S> 
struct div_helper {};

template<class T, class S>
struct div_helper<complex<T>,complex<S>>
{
    using TS            = typename max_type_base<T,S>::type;
    using base_complex  = std::complex<TS>;
    using return_type   = complex<TS>;

    static return_type recover_nan(const complex<T>& a, const complex<S>& b, 
                                   TS ret_re, TS ret_im)
    {
	    TS a_re         = real(a);
        TS a_im         = imag(a);
	    TS b_re         = real(b);
        TS b_im         = imag(b);

        namespace ms = mrd::scal_func;

        if ((b_re == 0.0 && b_im == 0) 
            && (ms::isnan(a_re) == false || ms::isnan(a_im) == false)) 
        {
            TS inf      = constants::inf<TS>();
            ret_re      = ms::copysign(inf, b_re) * a_re;
            ret_im      = ms::copysign(inf, b_re) * a_im;
        }
        else if ((ms::isinf(a_re) == true || ms::isinf(a_im) == true) 
                 && ms::finite(b_re) && ms::finite(b_im)) 
        {
            a_re        = ms::copysign(ms::isinf(a_re) ? TS(1.0) : TS(0.0), a_re);
            a_im        = ms::copysign(ms::isinf(a_im) ? TS(1.0) : TS(0.0), a_im);
            TS inf      = constants::inf<TS>();

            ret_re      = inf * ( a_re * b_re + a_im * b_im );
            ret_im      = inf * ( a_im * b_re - a_re * b_im );
        }
        else if ((ms::isinf(b_re) == true || ms::isinf(b_im) == true) 
                 && ms::finite(a_re) == true && ms::finite(a_im) == true 
                 && ms::isnan(b_re) == false && ms::isnan(b_im) == false) 
        {
            b_re        = ms::copysign(ms::isinf(b_re) ? TS(1.0) : TS(0.0), b_re);
            b_im        = ms::copysign(ms::isinf(b_im) ? TS(1.0) : TS(0.0), b_im);
            ret_re      = TS(0.0) * ( a_re * b_re + a_im * b_im );
            ret_im      = TS(0.0) * ( a_im * b_re - a_re * b_im );
        }

        return return_type(ret_re, ret_im);
    };

    force_inline
    static void scaled_div(TS a_re, TS a_im, TS b_re, TS b_im, 
                           TS& ret_re, TS& ret_im)
    {
        namespace ms = mrd::scal_func;

		TS r        = b_im / b_re;		
        //TS t      = TS(1)/ (b_re + r * b_im);
        TS t        = TS(1) / ms::fma(r, b_im, b_re);        

        //underflow protection; see A Robust Complex Division in Scilab
        //Michael Baudin, Robert L. Smith;

        //use fma instead of sequence of *,+ in order to gain speed and
        //accuracy
        if (r != TS(0))
        {            
            //ret_re    = (a_re + a_im * r) * t;
            ret_re      = ms::fma(a_im, r, a_re) * t;

            //ret_im    = (a_im - a_re * r) * t;
            ret_im      = ms::fma(-a_re, r, a_im) * t;		    
        }
        else
        {
            //ret_re    = (a_re + b_im * (a_im/b_re)) * t;
            ret_re      = ms::fma(b_im, a_im/b_re, a_re) * t;

            //ret_im    = (a_im - b_im * (a_re/b_re)) * t;
            ret_im      = ms::fma(-b_im, a_re/b_re, a_im) * t;		    
        }

        return;
    };

    force_inline
    static return_type eval(const complex<T>& a, const complex<S>& b)
    {
	    TS a_re         = real(a);
        TS a_im         = imag(a);
	    TS b_re         = real(b);
        TS b_im         = imag(b);

        namespace ms = mrd::scal_func;

        TS ret_re, ret_im;

	    if (ms::abs(b_im) < ms::abs(b_re))
        {
            scaled_div(a_re, a_im, b_re, b_im, ret_re, ret_im);
        }
        else
        {
            scaled_div(a_im, a_re, b_im, b_re, ret_re, ret_im);
            ret_im = - ret_im;
        }

        if (ms::isnan(ret_re) == true && ms::isnan(ret_im) == true) 
            return recover_nan(a, b, ret_re, ret_im);
        else
            return return_type(ret_re, ret_im);
    };

    force_inline
    static return_type eval_0(const complex<T>& a, const complex<S>& b)
    {
        if (real(a) == 0 && imag(a) == 0 && real(b) == 0 && imag(b) == 0)
            return return_type(0,0);
        else
            return eval(a,b);
    };

    force_inline
    static return_type eval_1(const complex<T>& a, const complex<S>& b)
    {
        if (real(a) == 0 && imag(a) == 0 && real(b) == 0 && imag(b) == 0)
            return return_type(1,0);
        else
            return eval(a,b);
    };
};

template<class T, class S>
struct div_helper<complex<T>,S>
{
    using TS            = typename max_type_base<T,S>::type;
    using return_type   = complex<TS>;

    force_inline
    static return_type eval(const complex<T>& a, const S& b)
    {
        TS a_re         = TS(real(a));
        TS a_im         = TS(imag(a));
	    TS b_rei        = TS(1) / TS(b);

        TS r_re         = a_re * b_rei;
        TS r_im         = a_im * b_rei;

        return return_type(r_re, r_im);
    };

    force_inline
    static return_type eval_0(const complex<T>& a, const S& b)
    {
        if (real(a) == 0 && imag(a) == 0 && b == 0)
            return return_type(0,0);
        else
            return eval(a,b);
    };

    force_inline
    static return_type eval_1(const complex<T>& a, const S& b)
    {
        if (real(a) == 0 && imag(a) == 0 && b == 0)
            return return_type(1,0);
        else
            return eval(a,b);
    };
};

template<class T, class S>
struct div_helper<T,complex<S>>
{
    using TS            = typename max_type_base<T,S>::type;
    using return_type   = complex<TS>;

    static return_type recover_nan(const T& a, const complex<S>& b, TS ret_re, TS ret_im)
    {
	    TS a_re         = a;
	    TS b_re         = real(b);
        TS b_im         = imag(b);

        namespace ms = mrd::scal_func;

        if ((b_re == 0.0 && b_im == 0) && ms::isnan(a_re) == false)
        {
            TS inf      = constants::inf<TS>();
            ret_re      = ms::copysign(inf, b_re) * a_re;
            ret_im      = TS(0.0);
        }
        else if (ms::isinf(a_re) == true && ms::finite(b_re) && ms::finite(b_im)) 
        {
            a_re        = ms::copysign(ms::isinf(a_re) ? TS(1.0) : TS(0.0), a_re);
            TS inf      = constants::inf<TS>();

            ret_re      = inf * ( a_re * b_re);
            ret_im      = inf * (- a_re * b_im );
        }
        else if ((ms::isinf(b_re) == true || ms::isinf(b_im) == true) 
                 && ms::finite(a_re) == true
                 && ms::isnan(b_re) == false && ms::isnan(b_im) == false) 
        {
            b_re        = ms::copysign(ms::isinf(b_re) ? TS(1.0) : TS(0.0), b_re);
            b_im        = ms::copysign(ms::isinf(b_im) ? TS(1.0) : TS(0.0), b_im);
            ret_re      = TS(0.0) * ( a_re * b_re );
            ret_im      = TS(0.0) * ( - a_re * b_im );
        }

        return return_type(ret_re, ret_im);
    };

    force_inline
    static void scaled_div(TS a_re, TS b_re, TS b_im, TS& ret_re, TS& ret_im)
    {
        namespace ms = mrd::scal_func;

		TS r        = b_im / b_re;
		//TS d      = b_re + r * b_im;
        TS d        = ms::fma(r, b_im, b_re);

        if (r != TS(0))
        {            
		    ret_re  = a_re / d;
		    ret_im  = -ret_re * r;
        }
        else
        {
		    ret_re  = (a_re) / d;
		    ret_im  = (- b_im * (a_re/b_re)) / d;
        }

        return;
    };

    force_inline
    static void scaled_div2(TS a_re, TS b_re, TS b_im, TS& ret_re, TS& ret_im)
    {
        namespace ms = mrd::scal_func;

		TS r        = b_re / b_im;
		//TS d      = b_im + r * b_re;
        TS d        = ms::fma(r, b_re, b_im);

        if (r != TS(0))
        {            
            TS tmp  = a_re / d;
		    ret_re  = tmp * r;
		    ret_im  = -tmp;
        }
        else
        {
		    ret_re  = (b_re * (a_re/ b_im)) / d;
		    ret_im  = -a_re / d;
        }
        return;
    };

    force_inline
    static return_type eval(const T& a, const complex<S>& b)
    {
	    TS a_re         = a;
	    TS b_re         = real(b);
        TS b_im         = imag(b);

        namespace ms = mrd::scal_func;

        TS ret_re, ret_im;

	    if (ms::abs(b_im) < ms::abs(b_re))
            scaled_div(a_re, b_re, b_im, ret_re, ret_im);
        else
            scaled_div2(a_re, b_re, b_im, ret_re, ret_im);

        if (ms::isnan(ret_re) == true && ms::isnan(ret_im) == true) 
            return recover_nan(a, b, ret_re, ret_im);
        else
            return return_type(ret_re, ret_im);
    };

    force_inline
    static return_type eval_0(const T& a, const complex<S>& b)
    {
        if (a == 0 && real(b) == 0 && imag(b) == 0)
            return return_type(0,0);
        else
            return eval(a,b);
    };

    force_inline
    static return_type eval_1(const T& a, const complex<S>& b)
    {
        if (a == 0 && real(b) == 0 && imag(b) == 0)
            return return_type(1,0);
        else
            return eval(a,b);
    };
};

// implement inversion for complex values
template<class T> 
struct inv_helper {};

template<class T>
struct inv_helper<complex<T>>
{
    using return_type = complex<T>;

    force_inline
    static return_type eval(const complex<T>& a)
    {
        return div_helper<T, complex<T>>::eval(T(1.), a);
    };
};

// implement operator< for complex values, lexicographic ordering
template<class T, class S> 
struct lt_helper {};

template<class T, class S>
struct lt_helper<complex<T>,complex<S>>
{
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        if (matcl::real(a) < matcl::real(b))
            return true;
        else if (matcl::real(a) == matcl::real(b))
            return matcl::imag(a) < matcl::imag(b);
        else
            return false;
    };
};

template<class T, class S>
struct lt_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        if (matcl::real(a) < b)
            return true;
        else if (matcl::real(a) == b)
            return matcl::imag(a) < T();
        else
            return false;
    };
};

template<class T, class S>
struct lt_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        if (a < matcl::real(b))
            return true;
        else if (a == matcl::real(b))
            return S() < matcl::imag(b);
        else
            return false;
    };
};

// implement operator> for complex values, lexicographic ordering
template<class T, class S> 
struct gt_helper {};

template<class T, class S>
struct gt_helper<complex<T>,complex<S>>
{
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        if (matcl::real(a) > matcl::real(b))
            return true;
        else if (matcl::real(a) == matcl::real(b))
            return matcl::imag(a) > matcl::imag(b);
        else
            return false;
    };
};

template<class T, class S>
struct gt_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        if (matcl::real(a) > b)
            return true;
        else if (matcl::real(a) == b)
            return matcl::imag(a) > T();
        else
            return false;
    };
};

template<class T, class S>
struct gt_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        if (a > matcl::real(b))
            return true;
        else if (a == matcl::real(b))
            return S() > matcl::imag(b);
        else
            return false;
    };
};

// implement operator<= for complex values, lexicographic ordering
template<class T, class S> 
struct leq_helper {};

template<class T, class S>
struct leq_helper<complex<T>,complex<S>>
{   
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        if (matcl::real(a) < matcl::real(b))
            return true;
        else if (matcl::real(a) == matcl::real(b))
            return matcl::imag(a) <= matcl::imag(b);
        else
            return false;
    };
};

template<class T, class S>
struct leq_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        if (matcl::real(a) < b)
            return true;
        else if (matcl::real(a) == b)
            return matcl::imag(a) <= T();
        else
            return false;
    };
};

template<class T, class S>
struct leq_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        if (a < matcl::real(b))
            return true;
        else if (a == matcl::real(b))
            return S() <= matcl::imag(b);
        else
            return false;
    };
};

// implement operator>= for complex values, lexicographic ordering
template<class T, class S> 
struct geq_helper {};

template<class T, class S>
struct geq_helper<complex<T>,complex<S>>
{
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        if (matcl::real(a) > matcl::real(b))
            return true;
        else if (matcl::real(a) == matcl::real(b))
            return matcl::imag(a) >= matcl::imag(b);
        else
            return false;
    };
};

template<class T, class S>
struct geq_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        if (matcl::real(a) > b)
            return true;
        else if (matcl::real(a) == b)
            return matcl::imag(a) >= T();
        else
            return false;
    };
};

template<class T, class S>
struct geq_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        if (a > matcl::real(b))
            return true;
        else if (a == matcl::real(b))
            return S() >= matcl::imag(b);
        else
            return false;
    };
};

// implement operator== for complex values
template<class T, class S> 
struct eeq_helper {};

template<class T, class S>
struct eeq_helper<complex<T>,complex<S>>
{
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        return (matcl::real(a) == matcl::real(b)) 
            && (matcl::imag(a) == matcl::imag(b));
    };
};

template<class T, class S>
struct eeq_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        return (matcl::real(a) == b) && (matcl::imag(a) == T());
    };
};

template<class T, class S>
struct eeq_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        return (a == matcl::real(b)) && (S() == matcl::imag(b));
    };
};

template<class T, class S>
force_inline
bool eeq_nan_r(T a, S b)
{
    if (mrd::scal_func::isnan(a) == true)
    {
        if (mrd::scal_func::isnan(b) == true)
            return true;
        else
            return false;
    }
    else
    {
        return (a == b);
    }
};

// implement x == y for complex values; nan values are equal
template<class T, class S> 
struct eeq_nan_helper {};

template<class T, class S>
struct eeq_nan_helper<complex<T>,complex<S>>
{
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        return eeq_nan_r(matcl::real(a), matcl::real(b)) 
            && eeq_nan_r(matcl::imag(a), matcl::imag(b));
    };
};

template<class T, class S>
struct eeq_nan_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        return eeq_nan_r(matcl::real(a), b) && (matcl::imag(a) == T());
    };
};

template<class T, class S>
struct eeq_nan_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        return eeq_nan_r(a, matcl::real(b)) && (S() == matcl::imag(b));
    };
};

// implement operator!= for complex values
template<class T, class S> 
struct neq_helper {};

template<class T, class S>
struct neq_helper<complex<T>,complex<S>>
{
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        return (matcl::real(a) != matcl::real(b)) 
                    || (matcl::imag(a) != matcl::imag(b));
    };
};

template<class T, class S>
struct neq_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        return (matcl::real(a) != b) || (matcl::imag(a) != T());
    };
};

template<class T, class S>
struct neq_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        return (a != matcl::real(b)) || (S() != matcl::imag(b));
    };
};

// implement operator!= for complex values; nan values are equal
template<class T, class S> 
struct neq_nan_helper {};

template<class T, class S>
struct neq_nan_helper<complex<T>,complex<S>>
{
    force_inline
    static bool eval(const complex<T>& a, const complex<S>& b)
    {
        return eeq_nan_r(matcl::real(a), matcl::real(b)) == false 
                    || eeq_nan_r(matcl::imag(a), matcl::imag(b)) == false;
    };
};

template<class T, class S>
struct neq_nan_helper<complex<T>,S>
{
    force_inline
    static bool eval(const complex<T>& a, const S& b)
    {
        return eeq_nan_r(matcl::real(a), b) == false || (matcl::imag(a) != T());
    };
};

template<class T, class S>
struct neq_nan_helper<T,complex<S>>
{
    force_inline
    static bool eval(const T& a, const complex<S>& b)
    {
        return eeq_nan_r(a, matcl::real(b)) == false || (S() != matcl::imag(b));
    };
};

template<class T> force_inline
typename details::uminus_helper<T>::return_type
uminus_c(const T& a)
{
    return details::uminus_helper<T>::eval(a);
};

template<class T> force_inline
typename details::inv_helper<T>::return_type
inv_c(const T& a)
{
    return details::inv_helper<T>::eval(a);
};

template<class T, class S> force_inline
typename details::plus_helper<T,S>::return_type
plus_c(const T& a, const S& b)
{
    return details::plus_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
typename details::minus_helper<T,S>::return_type
minus_c(const T& a, const S& b)
{
    return details::minus_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
typename details::mul_helper<T,S>::return_type
mul_c(const T& a, const S& b)
{
    return details::mul_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
typename details::div_helper<T,S>::return_type
div_c(const T& a, const S& b)
{
    return details::div_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
typename details::div_helper<T,S>::return_type
div_0_c(const T& a, const S& b)
{
    return details::div_helper<T,S>::eval_0(a,b);
};

template<class T, class S> force_inline
typename details::div_helper<T,S>::return_type
div_1_c(const T& a, const S& b)
{
    return details::div_helper<T,S>::eval_1(a,b);
};

template<class T, class S> force_inline
bool lt_c(const T& a, const S& b)
{
    return details::lt_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
bool gt_c(const T& a, const S& b)
{
    return details::gt_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
bool leq_c(const T& a, const S& b)
{
    return details::leq_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
bool geq_c(const T& a, const S& b)
{
    return details::geq_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
bool eeq_c(const T& a, const S& b)
{
    return details::eeq_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
bool neq_c(const T& a, const S& b)
{
    return details::neq_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
bool eeq_nan_c(const T& a, const S& b)
{
    return details::eeq_nan_helper<T,S>::eval(a,b);
};

template<class T, class S> force_inline
bool neq_nan_c(const T& a, const S& b)
{
    return details::neq_nan_helper<T,S>::eval(a,b);
};

};};

template<class T> force_inline
matcl::complex<T>& matcl::operator+=(complex<T>& arg, const complex<T>& val)
{
    arg = details::plus_helper<complex<T>,complex<T>>::eval(arg,val);
    return arg;
};

template<class T> force_inline
matcl::complex<T>& matcl::operator-=(complex<T>& arg, const complex<T>& val)
{
    arg = details::minus_helper<complex<T>,complex<T>>::eval(arg,val);
    return arg;
};

template<class T> force_inline
matcl::complex<T>& matcl::operator*=(complex<T>& arg, const complex<T>& val)
{
    arg = details::mul_helper<complex<T>,complex<T>>::eval(arg,val);
    return arg;
};

template<class T> force_inline
matcl::complex<T>& matcl::operator/=(complex<T>& arg, const complex<T>& val)
{
    arg = details::div_helper<complex<T>,complex<T>>::eval(arg,val);
    return arg;
};

template<class T> force_inline
bool matcl::operator==(const complex<T>& arg, const complex<T>& val)
{
    using type = complex<T>;
    return details::eeq_c(arg,val);
};

template<class T> force_inline
bool matcl::operator!=(const complex<T>& arg, const complex<T>& val)
{
    using type = complex<T>;
    return details::neq_c(arg,val);
}

template<class T> force_inline
bool matcl::operator>=(const complex<T>& arg, const complex<T>& val)
{
    using type = complex<T>;
    return details::geq_c(arg,val);
}

template<class T> force_inline
bool matcl::operator<=(const complex<T>& arg, const complex<T>& val)
{
    using type = complex<T>;
    return details::leq_c(arg,val);
}

template<class T> force_inline
bool matcl::operator<(const complex<T>& arg, const complex<T>& val)
{
    using type = complex<T>;
    return details::lt_c(arg,val);
}

template<class T> force_inline
bool matcl::operator>(const complex<T>& arg, const complex<T>& val)
{
    using type = complex<T>;
    return details::gt_c(arg,val);
}

template<class T> force_inline
matcl::complex<T> matcl::operator+(const complex<T>& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::plus_c<type,type>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator+(const complex<T>& arg1, const T& arg2)
{
    using type = complex<T>;
    return details::plus_c<type,T>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator+(const T& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::plus_c<T,type>(arg1,arg2);
}

template<class T> force_inline
matcl::complex<T> matcl::operator-(const complex<T>& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::minus_c<type,type>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator-(const complex<T>& arg1, const T& arg2)
{
    using type = complex<T>;
    return details::minus_c<type,T>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator-(const T& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::minus_c<T,type>(arg1,arg2);
}

template<class T> force_inline
matcl::complex<T> matcl::operator*(const complex<T>& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::mul_c<type,type>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator*(const complex<T>& arg1, const T& arg2)
{
    using type = complex<T>;
    return details::mul_c<type,T>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator*(const T& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::mul_c<T,type>(arg1,arg2);
}

template<class T> force_inline
matcl::complex<T> matcl::operator/(const complex<T>& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::div_c<type,type>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator/(const complex<T>& arg1, const T& arg2)
{
    using type = complex<T>;
    return details::div_c<type,T>(arg1,arg2);
};

template<class T> force_inline
matcl::complex<T> matcl::operator/(const T& arg1, const complex<T>& arg2)
{
    using type = complex<T>;
    return details::div_c<T,type>(arg1,arg2);
}

template<class T> force_inline
matcl::complex<T> matcl::operator-(const complex<T>& arg1)
{
    using type = complex<T>;
    return details::uminus_c<type>(arg1);
}
