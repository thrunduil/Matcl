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

namespace matcl
{

namespace mmd = matcl::mp::details;
namespace mrd = matcl::raw::details;

class eval_plus_1
{
    private:
        mp_complex  m_x;

    public:
        eval_plus_1(const mp_complex& x)  : m_x(x) {};

        precision   get_precision() const           { return m_x.get_precision(); };
        bool        is_zero_im() const              { return is_zero(imag(m_x)); };
        bool        is_one_im() const               { return is_one(imag(m_x)); };
        bool        is_mone_im() const              { return imag(m_x) == -1; };
        bool        is_one_re() const               { return is_zero(real(m_x)); };
        bool        is_mone_re() const              { return real(m_x) == -2; };
        bool        is_zero_re() const              { return real(m_x) == -1; };

        mp_complex  eval_prec(precision p) const    { return plus(1, m_x, p); };
};

template<class Eval_x>
mp_complex log_general(const Eval_x& x, Real err_re, precision req_prec)
{
    req_prec                = mmd::result_prec(req_prec, x.get_precision());
    precision min_prec      = mmd::extend_prec_dynamic(req_prec, 2.0 + 2*err_re);
    Real lost_prec          = 2.0 * 20.0;       //initial assumption
    precision int_prec      = mmd::extend_prec_dynamic(req_prec, lost_prec);
    int_prec                = precision(std::max<size_t>(int_prec, min_prec));
    precision small_prec    = precision(10);    // 3 correct digits

    if (x.is_zero_im())
    {
        if (x.is_one_re())
        {
            //log(1) = 0
            return mp_complex(req_prec);
        }
        if (x.is_mone_re())
        {
            //log(-1) = pi*i
            return mp_complex(mp_float(req_prec), constants::mp_pi(req_prec));
        }
    }
    else if (x.is_zero_re())
    {
        if (x.is_one_im())
        {
            //log(1i) = pi/2 * 1i
            return mp_complex(mp_float(req_prec), constants::mp_pi_2(req_prec));
        }
        else if (x.is_mone_im())
        {
            //log(-1i) = -pi/2 * 1i
            return mp_complex(mp_float(req_prec), -constants::mp_pi_2(req_prec));
        }
    };    

    mp_float r_re, x_re, x_im; 

    for(int iter = 0; ; ++iter)
    { 
        mp_complex x_prec       = x.eval_prec(int_prec);
        x_re                    = real(x_prec);
        x_im                    = imag(x_prec);

        mp_float abs_xre        = abs(x_re);
        mp_float abs_xim        = abs(x_im);
    
        mp_float hyp            = hypot(x_re, x_im, int_prec);  // 0.5 ulp + err_re

        if (is_finite(hyp) == false)
        {
            //log(inf) - inf; log(nan) = nan
            mp_float r_im       = atan2(x_im, x_re, req_prec);  // 0.5 ulp + err_re
            return mp_complex(hyp, r_im, req_prec);
        };    

        //r_re is nonzero (precision cannot be too low)
        r_re                = log(hyp, req_prec);   //0.5 + inv_log * (0.5 + err_re)

        Real inv_log        = div(1, r_re, small_prec).cast_float();
        inv_log             = mrd::scal_func::abs(inv_log) + 1;
        lost_prec           = 2.0 * (inv_log * (0.5 + err_re) + 0.5);
        size_t new_prec     = mmd::extend_prec_dynamic(req_prec, lost_prec);

        if (new_prec <= int_prec)
            break;

        if (iter > mmd::get_max_precision_iters())
        {
            mp::error::warn_possibly_inaccurate_result("log");
            break;
        };

        //limit precision growth
        if (new_prec > 2 * int_prec || iter > 4)
            new_prec        = precision(2 * int_prec);

        int_prec            = precision(new_prec + 4);
        hyp                 = hypot(x_re, x_im, int_prec);
    };

    mp_float r_im           = atan2(x_im, x_re, req_prec);  // 0.5 ulp + err_re + 0.5 (so)

    return mp_complex(r_re, r_im);
}

mp_complex matcl::log(const mp_complex& x, precision req_prec)
{
    const mp_float& x_re    = real(x);
    const mp_float& x_im    = imag(x);

    req_prec                = mmd::result_prec(req_prec, x.get_precision());
    Real lost_prec          = 2.0 * 20.0;       //initial assumption
    precision int_prec      = mmd::extend_prec_dynamic(req_prec, lost_prec);
    precision small_prec    = precision(10);    // 3 correct digits

    mp_float abs_xre        = abs(x_re);
    mp_float abs_xim        = abs(x_im);

    mp_float r_im           = atan2(x_im, x_re, req_prec);  // 0.5 ulp

    mp_float hyp, r_re;

    if (is_one(abs_xre) && is_zero(x_im) || is_zero(x_re) && is_one(abs_xim))
    {
        //hypot(x_re, x_im) = 1 and log(hyp) = 0
        return mp_complex(mp_float(req_prec), r_im);
    };

    hyp                     = hypot(x_re, x_im, int_prec);  // 0.5 ulp

    if (is_finite(hyp) == false)
    {
        //log(inf) - inf; log(nan) = nan
        return mp_complex(hyp, r_im, req_prec);
    };    
    
    Real eps_hyp    = 0.05;
    Real eps_x      = 0.05;

    //this is optimal range of the first method
    //bool select_1   = (hyp > 1 - eps_hyp && hyp < 1 + eps_hyp)
    //                && ((abs_xre > 1 - eps_x && abs_xre < 1 + eps_x)
    //                || (abs_xim > 1 - eps_x && abs_xim < 1 + eps_x));

    //however precision expansion is implemented only for the first method
    //so select this method even for xre and xim far from 1; in this case
    //the first method is not much more inaccurate (at most 15x less accurate)
    bool select_1   = (hyp > 1 - eps_hyp && hyp < 1 + eps_hyp);

    (void) eps_x;

    if (select_1)
    {       
        // this range was identified by comparing precision lost for
        // two methods

        bool greater            = mmd::compare_abs(x_re, x_im) == mmd::cmp_type::greater;

        mp_float e, hyp_m_1;

        for(int iter = 0; ; ++iter)
        {
            if (greater)
            {
                //|x_re| > |x_im|
                e               = minus(1, abs_xre, int_prec);//0.5 ulp
                hyp_m_1         = mul(x_im, x_im, int_prec) + mul(e, e, int_prec)
                                - mul(e, 2, int_prec);          //2.5 + 3*M ulp
                r_re            = log1p(hyp_m_1, req_prec);     //(2.5 + 3*M)*a + 0.5 ulp
                r_re            = div(r_re, 2, req_prec);       //(2.5 + 3*M)*a + 0.5 ulp
            }
            else
            {
                //|x_re| <= |x_im|
                e               = minus(1, abs_xim, int_prec);
                hyp_m_1         = mul(x_re, x_re, int_prec) + mul(e, e, int_prec)
                                - mul(e, 2, int_prec);
                r_re            = log1p(hyp_m_1, req_prec);
                r_re            = div(r_re, 2, req_prec);
            }

            size_t new_prec;
            Real error_prop;

            if (is_zero(hyp_m_1))
            {
                //in this case error propagation is huge; just double the precision
                new_prec        = precision(2 * int_prec);
            }
            else
            {
                error_prop      = 2.0 * div(e, hyp_m_1, small_prec).cast_float();

                //ulp error of error_prop
                Real ep_error   = 2.5 + 3*mrd::scal_func::abs(error_prop) + 1.0;
                Real bits_lost  = mrd::scal_func::ceil(mrd::scal_func::log2(ep_error));

                if (bits_lost + small_prec > int_prec)
                {
                    //error_prop is highly inaccurate; double precision
                    new_prec    = int_prec*2;
                }
                else
                {
                    Real ep_logp1   = 1.05;
                    lost_prec       = 2.0 * ((2.5 + 3 * mrd::scal_func::abs(error_prop)) * ep_logp1 + 0.5);
                    new_prec        = mmd::extend_prec_dynamic(req_prec, lost_prec);
                };
            };

            if (new_prec <= int_prec)
                break;

            if (iter > mmd::get_max_precision_iters())
            {
                mp::error::warn_possibly_inaccurate_result("log");
                break;
            };

            //avoid too large precision
            if (new_prec >= 2 * int_prec || iter > 4)
                new_prec        = 2 * int_prec;

            int_prec            = precision(new_prec + 4);
            hyp                 = hypot(x_re, x_im, int_prec);  // 0.5 ulp
        };

        return mp_complex(r_re, r_im, req_prec);
    };

    for(int iter = 0; ; ++iter)
    { 
        //r_re is nonzero (precision cannot be too low)
        r_re                = log(hyp, req_prec);

        Real inv_log        = div(1, r_re, small_prec).cast_float();
        inv_log             = mrd::scal_func::abs(inv_log) + 1;
        lost_prec           = 2.0 * (0.5*inv_log + 0.5);
        size_t new_prec     = mmd::extend_prec_dynamic(req_prec, lost_prec);

        if (new_prec <= int_prec)
            break;

        if (iter > mmd::get_max_precision_iters())
        {
            mp::error::warn_possibly_inaccurate_result("log");
            break;
        };

        //limit precision growth
        if (new_prec > 2 * int_prec || iter > 4)
            new_prec        = precision(2 * int_prec);

        int_prec            = precision(new_prec + 4);
        hyp                 = hypot(x_re, x_im, int_prec);
    };

    return mp_complex(r_re, r_im, req_prec);
}

mp_complex matcl::log1p(const mp_complex& x, precision req_prec)
{
    const mp_float& x_re    = real(x);
    const mp_float& x_im    = imag(x);

    mp_float abs_xre    = abs(x_re);
    mp_float abs_xim    = abs(x_im);

    Real eps_x          = 0.05;

    if (is_zero(abs_xre) && is_zero(x_im))
    {
        //log(1) = 0
        return mp_complex(req_prec);
    };

    // see also eval
    bool select_1       = (abs_xre < eps_x && abs_xim < eps_x);

    if (select_1 == false)
    {
        //err (1+x) = 0.5
        return log_general(eval_plus_1(x), 0.5, req_prec);
    }

    req_prec                = mmd::result_prec(req_prec, x.get_precision());
    Real lost_prec          = 2.0 * 20.0;       //initial assumption
    precision int_prec      = mmd::extend_prec_dynamic(req_prec, lost_prec);
    precision small_prec    = precision(10);    // 3 correct digits
    precision min_prec      = precision(req_prec + 1);

    mp_float r_im           = atan2(x_im, plus(1, x_re, min_prec), req_prec);  // 1.0 ulp + 0.5 (so)
    mp_float r_re;

    for(int iter = 0; ; ++iter)
    {
        //|x_re| > |x_im|
        mp_float hyp_m_1= mul(x_im, x_im, int_prec) + mul(x_re, x_re, int_prec)
                        + mul(x_re, 2, int_prec);       //1.5 + 1*M ulp
        r_re            = log1p(hyp_m_1, req_prec);     //(1.5 + 1*M)*a + 0.5 ulp
        r_re            = div(r_re, 2, req_prec);       //(1.5 + 1*M)*a + 0.5 ulp

        size_t new_prec;
        Real error_prop;

        if (is_zero(hyp_m_1))
        {
            //in this case error propagation is huge; just double the precision
            new_prec        = precision(2 * int_prec);
        }
        else
        {
            error_prop      = 2.0 * div(x_re, hyp_m_1, small_prec).cast_float();

            //ulp error of error_prop
            Real ep_error   = 1.5 + 1*mrd::scal_func::abs(error_prop) + 0.5;
            Real bits_lost  = mrd::scal_func::ceil(mrd::scal_func::log2(ep_error));

            if (bits_lost + small_prec > int_prec)
            {
                //error_prop is highly inaccurate; double precision
                new_prec    = int_prec*2;
            }
            else
            {
                Real ep_logp1   = 1.05;
                lost_prec       = 2.0 * ((1.5 + 1 * mrd::scal_func::abs(error_prop)) * ep_logp1 + 0.5);
                new_prec        = mmd::extend_prec_dynamic(req_prec, lost_prec);
            };
        };

        if (new_prec <= int_prec)
            break;

        if (iter > mmd::get_max_precision_iters())
        {
            mp::error::warn_possibly_inaccurate_result("log1p");
            break;
        };

        //avoid too large precision
        if (new_prec >= 2 * int_prec || iter > 4)
            new_prec        = 2 * int_prec;

        int_prec            = precision(new_prec + 4);
    };

    return mp_complex(r_re, r_im);
}

};
