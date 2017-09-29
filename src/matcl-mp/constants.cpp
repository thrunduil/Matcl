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

#include "matcl-mp/constants.h"
#include "matcl-mp/func_unary.h"
#include "matcl-mp/func_binary.h"
#include "matcl-core/lib_functions/constants.h"
#include "utils/impl_types.h"
#include "utils/utils.h"
#include "matcl-mp/cache.h"
#include "utils/extend_precision.h"

namespace matcl { namespace constants
{

namespace mmd = matcl::mp::details;

//-----------------------------------------------------------------
//                      float constants
//-----------------------------------------------------------------
mp_float constants::mp_eps(precision prec)
{
    return eps(mp_float(1.0, prec));
};

mp_float constants::mp_inf(precision prec)
{
    return mp_float(constants::inf(), prec);
};

mp_float constants::mp_nan(precision prec)
{
    return mp_float(constants::nan(), prec);
};

mp_float constants::mp_pi(precision req_prec)
{
    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());
    mp_float ret(req_prec);

    using impl_type = mmd::impl_float;
    impl_type& val  = *mmd::get_ptr<mp_float>::eval(ret.m_data);

    mpfr_const_pi(val, MPFR_RNDN);
    ret.update_debug();
    return ret;
};

mp_float constants::mp_ln2(precision req_prec)
{
    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());
    mp_float ret(req_prec);

    using impl_type = mmd::impl_float;
    impl_type& val  = *mmd::get_ptr<mp_float>::eval(ret.m_data);

    mpfr_const_log2(val, MPFR_RNDN);
    ret.update_debug();
    return ret;
};

mp_float constants::mp_e(precision req_prec)
{
    static matcl::unique_id id("e");

    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());

    if (matcl::cache::has(id, req_prec) == true)
        return matcl::cache::get(id, req_prec);
       
    mp_float val    = exp(mp_float(1, req_prec), req_prec);

    matcl::cache::set(id, req_prec, val);

    return val;
};

mp_float constants::mp_pi_2(precision req_prec)
{
    static matcl::unique_id id("pi_2");

    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());

    if (matcl::cache::has(id, req_prec) == true)
        return matcl::cache::get(id, req_prec);
       
    mp_float val    = div(mp_pi(req_prec), 2, req_prec);

    matcl::cache::set(id, req_prec, val);

    return val;
};

mp_float constants::mp_pi_4(precision req_prec)
{
    static matcl::unique_id id("pi_4");

    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());

    if (matcl::cache::has(id, req_prec) == true)
        return matcl::cache::get(id, req_prec);
       
    mp_float val    = div(mp_pi(req_prec), 4, req_prec);

    matcl::cache::set(id, req_prec, val);

    return val;
};

mp_float constants::mp_ln10(precision req_prec)
{
    static matcl::unique_id id("ln10");

    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());

    if (matcl::cache::has(id, req_prec) == true)
        return matcl::cache::get(id, req_prec);
       
    mp_float val    = log(mp_float(10, req_prec), req_prec);

    matcl::cache::set(id, req_prec, val);

    return val;
};

mp_float constants::mp_log2e(precision req_prec)
{
    static matcl::unique_id id("log2e");

    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());

    if (matcl::cache::has(id, req_prec) == true)
        return matcl::cache::get(id, req_prec);

    precision prec2 = precision(req_prec + mmd::extra_bits_constants());
    mp_float val    = inv(constants::mp_ln2(prec2), req_prec);

    matcl::cache::set(id, req_prec, val);

    return val;
}

mp_float constants::mp_log10e(precision req_prec)
{
    static matcl::unique_id id("log10e");

    req_prec        = mmd::result_prec(req_prec, mp_float::get_default_precision());

    if (matcl::cache::has(id, req_prec) == true)
        return matcl::cache::get(id, req_prec);

    precision prec2 = precision(req_prec + mmd::extra_bits_constants());
    mp_float val    = log10(constants::mp_e(prec2), req_prec);
    matcl::cache::set(id, req_prec, val);

    return val;
};

//-----------------------------------------------------------------
//                      complex constants
//-----------------------------------------------------------------
mp_complex constants::mp_i(precision prec)
{
    return mp_complex(0.0, 1.0, prec);
};

};};