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

mp_complex matcl::atanh(const mp_complex& z, precision req_prec)
{
    //adapted from boost::atanh

    const mp_float& z_re    = real(z);
    const mp_float& z_im    = imag(z);

    req_prec                = mmd::result_prec(req_prec, z.get_precision());
    Real err_0              = 10.0; //initial assumption, at least 5
    precision int_prec      = mmd::extend_prec_dynamic(req_prec, 2.0 * err_0);
    precision small_prec    = precision(10);
   
    mp_float half_pi    = constants::mp_pi_2(int_prec);
    mp_float one        = mp_float(1, int_prec);

    mp_float x          = abs(z_re);    //0 ulp
    mp_float y          = abs(z_im);    //0 ulp

    mp_float re, im;

    // real(atanh(z)) = log1p(4*x / ((x-1)*(x-1) + y^2))
    // imag(atanh(z)) = atan(2y, (1 - x)(1 + x) - y^2)/2

    // Begin by handling the special cases specified in C99:
    if(is_nan(x))
    {
        if(is_nan(y) == true)
            return mp_complex(x, x, req_prec);
        else if(is_inf(y) == true)
            return mp_complex(0, (signbit(z_im) ? -half_pi : half_pi), req_prec);
        else
            return mp_complex(x, x, req_prec);
    }
    else if (is_nan(y) == true)
    {
        if (is_zero(x) == true)
            return mp_complex(x, y, req_prec);
        if (is_inf(x) == true)
            return mp_complex(0, y, req_prec);
        else
            return mp_complex(y, y, req_prec);
    }
    else if (is_inf(x) == true || is_inf(y) == true)
    {
        return mp_complex(0, (signbit(z_im) ? -half_pi : half_pi), req_prec);
    };

    if (is_one(x) && is_zero(y))
    {
        if (signbit(z_re) == true)
            return mp_complex(-constants::mp_inf(req_prec), mp_float(req_prec));
        else
            return mp_complex(constants::mp_inf(req_prec), mp_float(req_prec));        
    };
    if (is_one(y) && is_zero(x))
    {
        if (signbit(z_im) == true)
            return mp_complex(mp_float(req_prec), -constants::mp_pi_4(req_prec));
        else
            return mp_complex(mp_float(req_prec), constants::mp_pi_4(req_prec));        
    };

    mp_float yy     = abs2(y, int_prec);                //0.5 ulp
    mp_float mxm1   = minus(one, x, int_prec);          //0.5 ulp

    // The real part is given by:    

    mp_float den    = mul(mxm1, mxm1, int_prec);        //2.0 = 1.5 ulp + 0.5 so !
    den             = plus(den, yy, int_prec);          //3.0 = 2.5 ulp + 0.5 so !
    re              = div(ldexp(x,2), den, int_prec);   //4.0 = 3.5 ulp + 0.5 so !
    re              = log1p(re, int_prec);              //4.5 = 0.5 ulp (log1p) + 4.0 ulp (log1p prop) !
    re              = ldexp(re, -2);                    //4.5 ulp

    if(signbit(z_re) == true)
        re  = -re;

    re.set_precision(req_prec);                         //5.0 ulp

    mp_float tmp;

    //im = atan(2y, (1 - x)(1 + x) - y^2)/2

    if (x == 1 || x == -1)
    {
        tmp         = -yy;
        goto lab_im;
    };

    for (int iter = 0; ; ++iter)
    {
        im              = plus(one, x, int_prec);           //0.5 ulp
        tmp             = mul(mxm1, im, int_prec);          //1.5 ulp = 0.5 + 0.5 (prop) + 0.5 (so)
        tmp             = minus(tmp, yy, int_prec);         //err = 0.5 + min_prop * 1.5

        precision np;

        if (is_zero(tmp) == true)
        {
            //exact zero is not possible; increase precision
            np              = precision(2 * int_prec);
        }
        else
        {
            Real min_prop   = div(yy, tmp, small_prec).cast_float();
            min_prop        = std::abs(min_prop);
            Real err        = 0.5 + (0.5 + min_prop * 1.5) + 0.5 + 0.5;
            np              = mmd::extend_prec_dynamic(req_prec, 2.0 * err);
        };

        if (np <= int_prec)
            break;

        if (iter > mmd::get_max_precision_iters())
        {
            mp::error::warn_possibly_inaccurate_result("atanh");
            break;
        };

        //limit precision growth
        if (np > 2 * int_prec || iter > 4)
            np          = precision(2 * int_prec);

        int_prec        = precision(np + 4);

        yy              = abs2(y, int_prec);                //0.5 ulp
        mxm1            = minus(one, x, int_prec);          //0.5 ulp
    };

  lab_im:
    im                  = atan2(ldexp(y, 1), tmp, int_prec);    //0.5 + (0.5 + min_prop * 1.5) + 0.5 (so)
    im                  = ldexp(im,-1);

    if(signbit(z_im) == true)
        im  = -im;

    im.set_precision(req_prec); //0.5 + (0.5 + min_prop * 1.5) + 0.5 + 0.5

    return mp_complex(re, im);
}

};
