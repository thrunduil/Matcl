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

mp_complex matcl::asin(const mp_complex& z, precision req_prec)
{
    // adapted from boost::asin, which is an implementation of
    // "Implementing the Complex Arcsine and Arccosine Functions using
    // Exception Handling." T E Hull, Thomas F Fairgrieve and Ping
    // Tak Peter Tang. ACM Transactions on Mathematical Software, Vol
    // 23, No 3, Sept 1997.
    //exceptional cases ignored

    const mp_float& z_re    = real(z);
    const mp_float& z_im    = imag(z);

    req_prec                = mmd::result_prec(req_prec, z.get_precision());
    Real err_0              = 6.0; //calculations
    precision int_prec      = mmd::extend_prec_dynamic(req_prec, 2.0 * err_0);
    precision small_prec    = precision(10);

    (void)err_0;

    //parameters of the algorithm
    const Real a_crossover  = 10;
    const Real b_crossover  = 0.6417;

    mp_float x          = abs(z_re);    //0 ulp
    mp_float y          = abs(z_im);    //0 ulp

    mp_float re, im;

    // Begin by handling the special cases specified in C99:
    if (is_nan(x) == true)
    {
        if (is_nan(y) == true)
            return mp_complex(x, x, req_prec);

        if (is_inf(y) == true)
        {
            re          = x;
            im          = constants::mp_inf(req_prec);
            goto lab_ret;
        }
        else
        {
            return mp_complex(x, x, req_prec);
        };
    }
    else if (is_nan(y) == true)
    {
        if (is_zero(x) == true)
        {
            re          = mp_float(req_prec);
            im          = y;
            goto lab_ret;
        }
        else if(is_inf(x) == true)
        {
            re          = y;
            im          = constants::mp_inf(req_prec);
            goto lab_ret;
        }
        else
        {
            return mp_complex(y, y, req_prec);
        }
    }
    else if(is_inf(x) == true)
    {
        if (is_inf(y) == true)
        {
            mp_float pi = constants::mp_pi(req_prec);
            mp_float pi4= ldexp(pi, -2);
            re          = pi4;
            im          = constants::mp_inf(req_prec);
            goto lab_ret;
        }
        else
        {
            mp_float pi = constants::mp_pi(req_prec);
            mp_float pi2= ldexp(pi, -1);
            re          = pi2;
            im          = constants::mp_inf(req_prec);
            goto lab_ret;
        }
    }
    else if (is_inf(y) == true)
    {
        re              = mp_float(req_prec);
        im              = constants::mp_inf(req_prec);
        goto lab_ret;
    };

    if (is_zero(y) && (x <= 1))
        return mp_complex(asin(z_re, req_prec), z_im, req_prec);

    {
        mp_float xp1        = plus(1, x, int_prec);     //0.5 ulp
        mp_float xm1        = minus(x, 1, int_prec);    //0.5 ulp

        mp_float yy         = abs2(y, int_prec);        //0.5 ulp

        mp_float r          = abs2(xp1) + yy;           //2.0 ulp = 1.5(abs) + 0.5(add)
        r                   = sqrt(r,int_prec);         //1.5 ulp = 0.5 + 1/2*2.0

        mp_float s          = abs2(xm1) + yy;           //2.0 ulp = 1.5(abs) + 0.5(add)
        s                   = sqrt(s);                  //1.5 ulp = 0.5 + 1/2*2.0

        mp_float a          = ldexp(r + s, -1);         //2.0 ulp
        mp_float b          = div(x, a, int_prec);      //3.0 ulp = 0.5(div) + 2.0(div prop) + 0.5(so)!        

        if (b <= b_crossover)
        {
            // in this region asin error prop <= 1.3
            re              = asin(b, int_prec);        //4.5ulp = 1.3*3 + 0.5
        }
        else
        {
            mp_float apx    = plus(a, x, int_prec);     //2.5 ulp

            if (x <= 1)
            {
                re          = yy /(r + xp1);            //3.5 ulp ([0.5] / [2.5])
                re          = re + (s-xm1);             //4.0 ulp ([3.5] + [2.0]), no cancelation
                re          = ldexp(apx * re, -1);      //7.0 ulp = [2.5] * [4.0]
                re          = sqrt(re);                 //4.0 ulp = 0.5 + 1/2*7
                re          = div(x, re, int_prec);     //4.5 ulp = [0]/[4.0]                
            }
            else
            {
                re          = apx/(r + xp1);            //5.0 ulp = [2.5]/[2.0]
                mp_float tmp= apx/(s+xm1);              //5.0 ulp = [2.5]/[2.0], no cancelation
                re          = ldexp(re + tmp, -1);      //5.5 ulp = [5] + [5]
                re          = sqrt(re, int_prec);       //3.5 ulp = 0.5 + 5.5/2;
                re          = mul(y, re, int_prec);     //4.0 ulp = [0]* [3.5];
                re          = div(x, re, int_prec);     //4.5 ulp = [0]/[4.0]
            }

            re              = atan(re);                 //5.5 ulp = 0.5 + 4.5(func prop) + 0.5(so)!
        }
        
        if(a <= a_crossover)
        {
            mp_float am1    = yy/(r + xp1);             //3.5 ulp = [0.5] / [2.5]
            if (x < 1)
                am1         = am1 + yy/(s - xm1);       //4.0 ulp = [3.5] + [3.5]                
            else
                am1         = am1 + (s + xm1);          //4.0 ulp = [3.5] + [2.5]
            am1             = ldexp(am1, -1);           //4.0 ulp

            im              = plus(a, 1, int_prec);     //2.5 ulp = [2.0] + [0]
            im              = sqrt(am1 * im);           //4.0 ulp = sqrt([7.0])
            im              = am1 + im;                 //4.5 ulp = [4.0] + [4.0];
            im              = log1p(im);                //5.0 ulp = log1p([4.5]) (im > 0)
        }
        else
        {
            im              = minus(abs2(a), 1, int_prec);  //5.0 ulp = [4.5]-[0], no cancel
            im              = a + sqrt(im);             //3.5 ulp = [2.0] + [3.0]

            //in this range log(im) >= 2.9
            im              = log(im);                  //2.0 ulp = 0.5 + [3.5]/2.9
        }

        re.set_precision(req_prec);     //6.0 ulp
        im.set_precision(req_prec);     //5.5 ulp
    };

  lab_ret:
    // set signs
    if(signbit(z_re) == true)
        re  = -re;
    if(signbit(z_im) == true)
        im  = -im;

    return mp_complex(re, im, req_prec);
};

};
