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

#include "matcl-mp/func_unary.h"
#include "utils/impl_types.h"
#include "utils/utils.h"
#include "utils/extend_precision.h"
#include "matcl-mp/constants.h"
#include "matcl-mp/func_binary.h"
#include "matcl-core/details/scalfunc_real.h"
#include "error.h"
#include "impl_functions.h"

namespace matcl
{

namespace mmd = matcl::mp::details;
namespace mrd = matcl::raw::details;

static mp_float mult_zero(const mp_float& r, const mp_float& arg, precision req_prec)
{
    req_prec = mmd::result_prec(req_prec, r.get_precision());

    if (is_zero(arg))
        return mp_float(arg, req_prec);

    return mul(r, arg, req_prec);
}

//------------------------------------------------------
//form exp(a * x) where a is real
//template argument Scalar is responsible for calculating a 
//with required precision up to 1 ulp error
template<class Scalar>
struct exp_mult
{
    static mp_complex eval(const Scalar& sc, const mp_complex& x, precision req_prec)
    {
        if (is_zero(imag(x)) == true)
            return sc.eval_real(real(x), req_prec);

        //for INF/NAN call exp to handle these cases
        if (is_regular(x) == false)
        {
            mp_complex ax = mul(x, sc.get(req_prec), req_prec);
            return exp(ax, req_prec);
        }

        precision xprec     = x.get_precision();
        req_prec            = mmd::result_prec(req_prec, xprec);

        // error propagation:
        //  re(exp(x+i*y))  : ep_exp_r  = |x|*a + |y*tan(y)| * b
        //  im(exp(x+i*y))  : ep_exp_i  = |x|*a + |y*cot(y)| * b
        // tan(y) + cot(y)  = 2/sin(2y)

        const mp_float& re  = real(x);
        const mp_float& im  = imag(x);    

        //calculate precision
        Real scal           = sc.get_real();
        Real real_xre       = std::abs(scal * re.cast_float());
        Real mult_re        = real_xre;

        //precision must be large enough to compute sin/cos accurately
        Integer exp_im      = ilogb(im) + 1;
        size_t inc_prec     = 10;   //assumption
        precision req_prec2 = precision((size_t)std::max(exp_im, 0) + req_prec);
        precision min_prec  = precision(req_prec2 + inc_prec);    
        precision int_prec  = min_prec;

        //calculate precision of imaginary part

        for (int iter = 0; ; ++iter)
        {
            mp_float scal_p  = sc.get(min_prec);
            mp_float xim2    = mul(im, mul(2, scal_p, min_prec), min_prec);
            mp_float s       = sin(xim2);

            if (is_zero(s))
            {
                //sinus cannot be exactly zero
                min_prec    = precision(min_prec + min_prec + 4);
                continue;
            };

            mp_float mult_im = abs(xim2 / s);
            Real mult_im_r   = mult_im.cast_float();
                         
            precision np     = mmd::extend_prec_exp_mult(req_prec2, mult_im_r);

            if (np <= min_prec)
            {
                //calculate total precision
                np          = mmd::extend_prec_exp_mult(req_prec2, mult_re + mult_im_r);
                int_prec    = precision(np + 1);
                break;
            };

            if (iter > mmd::get_max_precision_iters())
            {
                mp::error::warn_possibly_inaccurate_result("exp");
                break;
            };

            //restrict precision growth
            if (np > 2 * min_prec || iter > 4)
                np          = precision(2 * min_prec);

            min_prec        = precision(np + 4);
        };

        //do computations
        mp_float scal_p     = sc.get(int_prec);             //1 ulp
        mp_complex x2       = mul(x, scal_p, int_prec);     //1.5 ulp
        return exp(x2, req_prec);
    };
};

//------------------------------------------------------
//form exp(x); special cases are not handled
//template argument Scalar is responsible for calculating x
//with required precision up to ulp_error ulp error
//req_prec cannot be zero
//sin(2 * im(x)) cannot be exactly zero
template<class Scalar>
struct exp_x
{
    static mp_complex eval(const Scalar& sc, precision req_prec)
    {
        // error propagation:
        //  re(exp(x+i*y))  : ep_exp_r  = |x|*a + |y*tan(y)| * b
        //  im(exp(x+i*y))  : ep_exp_i  = |x|*a + |y*cot(y)| * b
        // tan(y) + cot(y)  = 2/sin(2y)

        int iter = 0;
        mp_complex xp;
        Real ulp_error_re;
        Real ulp_error_im;
        precision int_prec  = mmd::extend_prec_dynamic(req_prec, 10.0);

        for(;; ++iter)
        {
            xp              = sc.get(int_prec, ulp_error_re, ulp_error_im);

            if (mrd::scal_func::finite(ulp_error_re + ulp_error_im) == false)
                int_prec    = int_prec * 2;
            else
                break;

            if (iter > mmd::get_max_precision_iters())
            {
                mp::error::warn_possibly_inaccurate_result("exp");
                break;
            };
        }
        
        const mp_float& xre = real(xp);
        const mp_float& xim = imag(xp);    

        //calculate precision
        Real real_xre       = std::abs(xre.cast_float());
        Real mult_re        = real_xre * ulp_error_re;

        //precision must be large enough to compute sin/cos accurately
        Integer exp_im      = ilogb(xim) + 1;
        size_t inc_prec     = 10;   //assumption
        precision req_prec2 = precision((size_t)std::max(exp_im, 0) + int_prec);
        precision min_prec  = precision(req_prec2 + inc_prec);    
        int_prec            = min_prec;

        //calculate precision of imaginary part

        for (; ; ++iter)
        {
            if (iter > mmd::get_max_precision_iters())
            {
                mp::error::warn_possibly_inaccurate_result("exp");
                break;
            };

            if (xp.get_precision() < min_prec)
                xp          = sc.get(min_prec, ulp_error_re, ulp_error_im);

            if (mrd::scal_func::finite(ulp_error_re + ulp_error_im) == false)
            {
                min_prec    = min_prec * 2;
                req_prec2   = min_prec;
                continue;
            }

            mp_float xim2   = ldexp(imag(xp), 1);
            mp_float s      = sin(xim2);

            if (is_zero(s))
            {
                //sinus cannot be exactly zero
                min_prec    = precision(min_prec + min_prec + 4);
                continue;
            };

            mp_float mult_im = abs(xim2 / s);
            Real mult_im_r   = mult_im.cast_float() * ulp_error_im;
                         
            precision np     = mmd::extend_prec_exp_mult(req_prec2, mult_im_r);

            if (np <= min_prec)
            {
                //calculate total precision
                np          = mmd::extend_prec_exp_mult(req_prec2, mult_re + mult_im_r);
                int_prec    = precision(np + 1);
                break;
            };

            //restrict precision growth
            if (np > 2 * min_prec || iter > 4)
                np          = precision(2 * min_prec);

            min_prec        = precision(np + 4);
        };

        //do computations
        if (xp.get_precision() < int_prec)
            xp              = sc.get(int_prec, ulp_error_re, ulp_error_im);

        return exp(xp, req_prec);
    };
};

mp_complex matcl::exp(const mp_complex& x, precision req_prec)
{
    const mp_float& xre = real(x);
    const mp_float& xim = imag(x);

    if (is_zero(xim) == true)
        return exp(xre);

    req_prec            = mmd::result_prec(req_prec, x.get_precision());
    precision int_prec  = mmd::extend_prec_exp(req_prec);

    mp_float exp_re     = exp(xre, int_prec);

    mp_float sin_im, cos_im;
    sin_cos(xim, sin_im, cos_im, int_prec);

    mp_float r_re       = mult_zero(exp_re, cos_im, req_prec);
    mp_float r_im       = mult_zero(exp_re, sin_im, req_prec);

    return mp_complex(r_re, r_im);
}

mp_complex matcl::expm1(const mp_complex& x, precision req_prec)
{
    const mp_float& xre = real(x);
    const mp_float& xim = imag(x);

    if (is_zero(xim) == true)
        return expm1(xre);

    req_prec                = mmd::result_prec(req_prec, x.get_precision());
    Real err0               = 5.0;  //initial assumption, at least 2.5
    precision int_prec      = mmd::extend_prec_dynamic(req_prec, 2*err0);
    precision small_prec    = precision(10);    

    Real tol_small          = 0.05; // not checked if this is an optimal value

    if (abs(xre) < tol_small && abs(xim) < tol_small)
    {
        // re = exp(xre) * cos(xim) - 1 = (exp(xre) -1) * cos(xim) - 2*sin(xim/2)^2

        mp_float r_re, r_im, r_re_p1;

        mp_float exp_re     = exp(xre, int_prec);   //0.5 ulp

        mp_float sin_im, cos_im;
        sin_cos(xim, sin_im, cos_im, int_prec);     //0.5 ulp

        r_im                = mult_zero(exp_re, sin_im, req_prec);  //1.5 ulp

        for (int iter = 0; ; ++iter)
        {
            mp_float exp1   = expm1(xre, int_prec); // 0.5 ulp
            r_re_p1         = mult_zero(exp1, cos_im, int_prec);  //1.5 ulp
            mp_float xim2   = ldexp(xim, -1);       // 0 ulp
            mp_float sin_x  = sin(xim2, int_prec);  // 0.5 ulp
            sin_x           = abs2(sin_x, int_prec);// 1.5 ulp = 0.5 + 2*0.5
            sin_x           = ldexp(sin_x, 1);      // 1.5 ulp

            r_re            = minus(r_re_p1, sin_x, int_prec);  //0.5 + (1+2*M)*1.5

            //calculate required precision
            precision np;

            if (is_zero(r_re) == true)
            {
                np              = precision(2 * int_prec);
            }
            else
            {
                Real min_prop   = div(r_re_p1, r_re, small_prec).cast_float();
                min_prop        = std::abs(min_prop);
                Real err        = 0.5 + (1 + 2.0 * min_prop) * 1.5 + 0.5;
                np              = mmd::extend_prec_dynamic(req_prec, 2.0 * err);
            };

            if (np <= int_prec)
                break;

            if (iter > mmd::get_max_precision_iters())
            {
                mp::error::warn_possibly_inaccurate_result("expm1");
                break;
            };

            //limit precision growth
            if (np > 2 * int_prec || iter > 4)
                np          = precision(2 * int_prec);

            int_prec        = precision(np + 4);
            cos_im          = cos(xim, int_prec);
        };
    
        r_re.set_precision(req_prec);   //0.5 + (1+2*M)*1.5 + 0.5

        return mp_complex(r_re, r_im);
    }
    else
    {
        mp_float r_re, r_im, r_re_p1;

        mp_float exp_re     = exp(xre, int_prec);   //0.5 ulp

        mp_float sin_im, cos_im;
        sin_cos(xim, sin_im, cos_im, int_prec);     //0.5 ulp

        r_im                = mult_zero(exp_re, sin_im, req_prec);  //1.5 ulp

        for (int iter = 0; ; ++iter)
        {
            r_re_p1         = mult_zero(exp_re, cos_im, int_prec);  //1.5 ulp            
            r_re            = minus(r_re_p1, 1, int_prec); //0.5 ulp + M * 1.5 ulp

            //calculate required precision
            precision np;

            if (is_zero(r_re) == true)
            {
                np              = precision(2 * int_prec);
            }
            else
            {
                Real min_prop   = div(r_re_p1, r_re, small_prec).cast_float();
                min_prop        = std::abs(min_prop);
                Real err        = 0.5 + min_prop * 1.5 + 0.5;
                np              = mmd::extend_prec_dynamic(req_prec, 2.0 * err);
            };

            if (np <= int_prec)
                break;

            if (iter > mmd::get_max_precision_iters())
            {
                mp::error::warn_possibly_inaccurate_result("expm1");
                break;
            };

            //limit precision growth
            if (np > 2 * int_prec || iter > 4)
                np          = precision(2 * int_prec);

            int_prec        = precision(np + 4);
            exp_re          = exp(xre, int_prec);
            cos_im          = cos(xim, int_prec);
        };
    
        r_re.set_precision(req_prec);   //0.5 ulp + M * 1.5 ulp + 0.5 ulp

        return mp_complex(r_re, r_im);
    };
}

struct Ln2_scalar
{
    mp_complex eval_real(const mp_float& re, precision req_prec) const
    {
        return mp_complex(exp2(re, req_prec));
    };

    Real get_real() const
    {
        return constants::ln2();
    };

    mp_float get(precision p) const
    {
        return constants::mp_ln2(p);
    }
};

struct Ln10_scalar
{
    mp_complex eval_real(const mp_float& re, precision req_prec) const
    {
        return mp_complex(exp10(re, req_prec));
    };

    Real get_real() const
    {
        return constants::ln10();
    };

    mp_float get(precision p) const
    {
        return constants::mp_ln10(p);
    }
};

mp_complex matcl::exp2(const mp_complex& x, precision req_prec)
{
    return exp_mult<Ln2_scalar>::eval(Ln2_scalar(), x, req_prec);
}

mp_complex matcl::exp10(const mp_complex& x, precision req_prec)
{
    return exp_mult<Ln10_scalar>::eval(Ln10_scalar(), x, req_prec);
}

struct Log_a
{
    const mp_float& m_a;

    Log_a(const mp_float& a)    :m_a(a) {};

    mp_complex eval_real(const mp_float& x, precision req_prec) const
    {
        //eval exp(log(a) * x) = a^x
        return pow(m_a, x, req_prec);
    }

    mp_float get(precision req_prec) const
    {
        return log(m_a, req_prec);
    }

    Real get_real() const
    {
        return log(m_a, precision::precision_double()).cast_float();
    };
};

struct Log_ax
{
    const mp_complex&   m_a;
    const mp_float&     m_b;

    Log_ax(const mp_complex& a, const mp_float& b)
        :m_a(a), m_b(b)
    {};

    mp_complex get(precision p, Real& ulp_error_re, Real& ulp_error_im) const
    {
        ulp_error_re = 2.0;
        ulp_error_im = 2.0;
        return mul(log(m_a, p), m_b, p);
    }
};

//a is not regular or b is not finite; treat special cases as exp(log(a) * b)
mp_complex pow_compl_real_nreg(const mp_complex& a, const mp_float& b, precision p)
{
    return exp(mul(log(a, p), b, p));
}
mp_complex pow_compl_compl_nreg(const mp_complex& a, const mp_complex& b, precision p)
{
    return exp(mul(log(a, p), b, p));
}
mp_complex pow_real_compl_nreg(const mp_float& a, const mp_complex& b, precision p)
{
    return exp(mul(log(a, p), b, p));
}

//im(b) != 0
mp_complex pow_real_compl_reg(const mp_float& a, const mp_complex& b, precision p)
{
    //a and b are regular and a > 0

    if (is_one(a))
        return mp_complex(1.0, p);

    //xi = im(log(a) * b) = log(a) * im(b) != pi/2 * k
    //since log(a) != 0, im(b) != 0 and xi cannot be pi/2 * k exactly

    //exp(log(a) * b)
    return exp_mult<Log_a>::eval(Log_a(a), b, p);   //1.0 ulp
};

//b is not zero, a is not positive real
mp_complex pow_compl_real_reg(const mp_complex& a, const mp_float& b, precision p)
{
    //check when xi = im(log(a) * b) = pi/2 * k:
    //Log(a) = log(|a|) + i Arg(a)
    //xi = Arg(a) * b
    //xi = pi/2 *k for: 1. im(a) = 0, re(a) < 0 2*b integer
    //                  2. im(a) = 0, re(a) > 0 (this case is already handled)
    //                  3. re(a) = 0 b integer
    //                  4. |re(a)| = |im(a)| and b is even integer

    if (is_int(b) == true)
    {
        if (is_zero(imag(a)) == true)
        {
            mp_float res = mp::details::pow_impl(real(a), b, p);    //0.5 ulp
            return mp_complex(res);
        }
        else if (is_zero(real(a)) == true)
        {
            Integer mult = mod(b, 4.0).cast_int();

            mp_float res = mp::details::pow_impl(imag(a), b, p);    //0.5 ulp

            switch(mult)
            {
                case 0: return mp_complex(res, p);                  //0.5 ulp
                case 1: return mp_complex(0.0, res, p);             //0.5 ulp
                case 2: return mp_complex(-res, p);                 //0.5 ulp
                case 3: return mp_complex(0.0, -res, p);            //0.5 ulp
                default:
                    throw std::runtime_error("invalid case");
            };
        }
        else if (abs(real(a)) == abs(imag(a)))
        {
            if (rem(b, 2) == 0)
            {
                //Log(a) = log(sqrt(2)|re|) + i Arg(a)
                mp_float res    = mp::details::pow_impl(abs(real(a)), b, p); //0.5 ulp

                if (b.can_cast_int() == false)
                {
                    //overlow will occur
                    return mp_complex(constants::mp_nan(p), constants::mp_nan(p), p);
                }

                Integer ex      = b.cast_int();
                ex              = ex/2;
                res             = ldexp(res, ex);   //0.5 ulp

                bool sign_eq    = signbit(real(a)) == signbit(imag(a));

                switch (ex % 4)
                {
                    //0.5 ulp
                    case 0: return mp_complex(res, 0, p);
                    case 1: return mp_complex(0, sign_eq ? res : -res, p);
                    case 2: return mp_complex(-res, 0, p);
                    case 3: return mp_complex(0, sign_eq ? -res : res, p);
                    case -1: return mp_complex(0, sign_eq ? -res : res, p);
                    case -2: return mp_complex(-res, 0, p);
                    case -3: return mp_complex(0, sign_eq? res : -res, p);
                    default:
                        throw std::runtime_error("invalid case");
                }
            };
        };
    }
    else if (is_zero(imag(a)) == true)
    {
        //real(a) < 0
        mp_float ex = ldexp(b, 1);

        if (is_int(ex) == true)
        {
            Integer mult    = mod(ex, 4.0).cast_int();
            mp_float res    = mp::details::pow_impl(abs(real(a)), b, p);    //0.5 ulp

            switch(mult)
            {
                //0.5 ulp
                case 0: return mp_complex(res, 0.0, p);
                case 1: return mp_complex(0.0, res, p);
                case 2: return mp_complex(-res,0.0, p);
                case 3: return mp_complex(0.0, -res, p);
                default:
                    throw std::runtime_error("invalid case");
            };
        };
    }

    //general case; sin(2xi) is never exact zero

    //log(a) * b    => 1.5 ulp + 0.5 (so)
    return exp_x<Log_ax>::eval(Log_ax(a, b), p);    //1.0 ulp
}

struct Log_ax_compl
{
    const mp_complex&   m_a;
    const mp_complex&   m_b;

    Log_ax_compl(const mp_complex& a, const mp_complex& b)
        :m_a(a), m_b(b)
    {};

    mp_complex get(precision p, Real& ulp_error_re, Real& ulp_error_im) const
    {
        precision small_prec    = precision(10);

        // re[(x+iy)*(u+iw)]: ep:       = |(u*x)/(u*x - w*y)| * ex + |(w*y)/(u*x - w*y)| * ey
        // im[(x+iy)*(u+iw)]: ep:       = |(w*x)/(u*y + w*x)| * ex + |(u*y)/(u*y + w*x)| * ey

        mp_complex l            = log(m_a, p);      //1 ulp
        mp_complex lb           = mul(l, m_b, p);   // re: 1.0 ulp + M_rr * 1.0 + M_ri * 1.0 + 0.5 (so)
                                                    // im: 1.0 ulp + M_ir * 1.0 + M_ii * 1.0 + 0.5 (so)
        const mp_float& l_re    = real(lb);
        const mp_float& l_im    = imag(lb);

        if (is_zero(l_re) || is_zero(l_im))
            return constants::inf();

        Real M_rr               = div(mul(real(m_b), real(l), small_prec), l_re, small_prec).cast_float();
        Real M_ri               = div(mul(imag(m_b), imag(l), small_prec), l_re, small_prec).cast_float();
        Real M_ir               = div(mul(imag(m_b), real(l), small_prec), l_im, small_prec).cast_float();
        Real M_ii               = div(mul(real(m_b), imag(l), small_prec), l_im, small_prec).cast_float();

        ulp_error_re            = 1.5 + std::abs(M_rr) + std::abs(M_ri);
        ulp_error_im            = 1.5 + std::abs(M_ir) + std::abs(M_ii);

        return lb;
    };
};

mp_complex pow_compl_compl_reg(const mp_complex& a, const mp_complex& b, precision p)
{
    //special cases should aready be handled => b is not real a is not positive real

    if (is_zero(imag(a)) == true && real(a) == -1)
    {
        //(-1)^(x + yi) = exp(-pi*y) * (-1)^x

        Real err        = 10.0; //assumption
        precision ip    = mmd::extend_prec_dynamic(p, 2.0 * err);
        mp_float piy;
        mp_complex m1x;

        for (int iter = 0; iter < 2; ++iter)
        {
            m1x             = pow_compl_real_reg(mp_complex(-1.0, ip), real(b), ip);//1 ulp
            piy             = mul(-constants::mp_pi(ip), imag(b), ip);  //1.5 ulp = 0.5 + 0.5 (prop) + 0.5 (so)
            err             = 2.0 + (0.5 + 1.5 * abs(piy).cast_float());

            precision np    = mmd::extend_prec_dynamic(p, 2.0 * err);

            if (np <= ip)
                break;

            ip              = np + 4;
        };

        piy             = exp(piy, ip);     //0.5 ulp + 1.5 * |piy|
        return mul(piy, m1x, p);            //0.5 + 1 + (0.5 + 1.5 * |piy|) + 0.5(so)
    }
    else if (is_zero(real(a)) == true)
    {
        if (is_one(imag(a)) == true)
        {
            //i^(x + yi)    = exp(-pi/2*y) * i^x

            Real err        = 10.0; //assumption
            precision ip    = mmd::extend_prec_dynamic(p, 2.0 * err);
            mp_float piy;
            mp_complex ix;

            for (int iter = 0; iter < 2; ++iter)
            {
                ix          = pow_compl_real_reg(mp_complex(0.0, 1.0, ip), real(b), ip);  //1 ulp
                piy         = mul(-constants::mp_pi(ip)/2, imag(b), ip);    //1.5 ulp = 0.5 + 0.5 (prop) + 0.5 (so)
                err         = 2.0 + (0.5 + 1.5 * abs(piy).cast_float());

                precision np    = mmd::extend_prec_dynamic(p, 2.0 * err);

                if (np <= ip)
                    break;

                ip              = np + 4;
            };

            piy             = exp(piy, ip);     //0.5 ulp + 1.5 * |piy|
            return mul(piy, ix, p);             //0.5 + 1 + (0.5 + 1.5 * |piy|) + 0.5(so)
        }
        else if (imag(a) == -1)
        {
            //(-i)^(x + yi) = exp(pi/2*y) * (-i)^x
            Real err        = 10.0; //assumption
            precision ip    = mmd::extend_prec_dynamic(p, 2.0 * err);
            mp_float piy;
            mp_complex ix;

            for (int iter = 0; iter < 2; ++iter)
            {
                ix          = pow_compl_real_reg(mp_complex(0.0, -1.0, ip), real(b), ip);  //1 ulp
                piy         = mul(constants::mp_pi(ip)/2, imag(b), ip);     //1.5 ulp = 0.5 + 0.5 (prop) + 0.5 (so)
                err         = 2.0 + (0.5 + 1.5 * abs(piy).cast_float());

                precision np    = mmd::extend_prec_dynamic(p, 2.0 * err);

                if (np <= ip)
                    break;

                ip              = np + 4;
            };

            piy             = exp(piy, ip);     //0.5 ulp + 1.5 * |piy|
            return mul(piy, ix, p);             //0.5 + 1 + (0.5 + 1.5 * |piy|) + 0.5(so)
        };
    };

    //check when im(log(a) * b) = pi/2 * k
    //Log(a) = log(|a|) + i Arg(a)
    //Log(a) * (br + ibi) = [log(|a|)*br - Arg(a) * bi] + [log(|a|)*bi + Arg(a) * br]i
    //[log(|a|)*bi + Arg(a) * br] = pi/2 * k 
    //      => |a| = 1, bi any, => a = 1, -1, i, -i (handled)

    return exp_x<Log_ax_compl>::eval(Log_ax_compl(a,b), p);
}

mp_complex mp::details::pow_real_compl(const mp_float& a, const mp_complex& b, precision p)
{
    p = mmd::result_prec(p, get_precision(a), get_precision(b));

    //x^0 = 1 for any x also NaN
    if (is_zero(b) == true)
        return mp_complex(1, p);
    if (b == 1)
        return mp_complex(a, p);
    if (b == -1)
        return inv(a, p);

    if (is_zero(imag(b)) == true)
    {
        if (a >= 0)
            return pow_impl(a, real(b), p);
        if (is_finite(real(b)) == false || is_finite(a) == false)
            return pow_impl(a, real(b), p);

        //a < 0 and a is regular
        //real(b) is regular
        return pow_compl_real_reg(mp_complex(a), real(b), p);        
    };

    if (a < 0)
    {
        if (is_regular(a) == true && is_regular(b) == true)
            return pow_compl_compl_reg(mp_complex(a), b, p);
        else
            return pow_compl_compl_nreg(mp_complex(a), b, p);
    };

    if (is_regular(a) == true && is_regular(b) == true)
        return pow_real_compl_reg(a, b, p);
    else
        return pow_real_compl_nreg(a, b, p);
}

mp_complex mp::details::pow_compl_real(const mp_complex& a, const mp_float& b, precision p)
{
    p = mmd::result_prec(p, a.get_precision(), b.get_precision());

    //x^0 = 1 for any x also NaN
    if (is_zero(b) == true)
        return mp_complex(1, p);
    if (b == 1)
        return mp_complex(a, p);
    if (b == -1)
        return inv(a, p);

    if (is_zero(imag(a)) == true)
    {
        if (real(a) >= 0)
            return pow_impl(real(a), b, p);
        if (is_regular(b) == false || is_finite(real(a)) == false)
            return pow_impl(real(a), b, p);
    };

    if (is_regular(a) == true && is_regular(b) == true)
        return pow_compl_real_reg(a, b, p);
    else
        return pow_compl_real_nreg(a, b, p);
};

mp_complex mp::details::pow_compl_compl(const mp_complex& a, const mp_complex& b, precision p)
{
    p = mmd::result_prec(p, a.get_precision(), b.get_precision());

    if (is_zero(imag(a)) == true)
        return pow_real_compl(real(a), b, p);

    if (is_zero(imag(b)) == true)
        return pow_compl_real(a, real(b), p);

    if (is_finite(a) == true && is_finite(b) == true)
        return pow_compl_compl_reg(a, b, p);
    else
        return pow_compl_compl_nreg(a, b, p);
};

};
