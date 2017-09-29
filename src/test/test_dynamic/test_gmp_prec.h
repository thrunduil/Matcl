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

#include "matcl-dynamic/matcl_dynamic.h"
#include "matcl-mp/matcl_mp.h"
#include "test_gmp.h"

#include <vector>

namespace matcl { namespace test
{

void test_gmp_prec(std::ostream& os);

class gmp_tester_prec
{
    private:
        using prec_vec  = std::vector<precision>;
        using max_vec   = std::vector<double>;
        using dtup_2    = std::tuple<double, double>;

    public:
        void    make(std::ostream& os);

        static mp_float    rand_scalar(precision prec, const max_vec& max);
        static mp_float    rand_scalar(const prec_vec& prec, const max_vec& max);
        static mp_float    rand_scalar(precision prec, double max);
        static mp_complex  rand_scalar_c(precision prec, double max_re, double max_im);
        static mp_complex  rand_scalar_c_log(precision prec, double max_v);
        static mp_complex  rand_scalar_c_exp2(precision prec, double max_re, double max_im);
        static mp_complex  rand_scalar_c_exp10(precision prec, double max_re, double max_im);
        static mp_complex  rand_scalar_c_log1p(precision prec, double max_v);
        static mp_complex  rand_scalar_c_expm1(precision prec, double max_v);
        static mp_complex  rand_scalar_c(const prec_vec& prec, const max_vec& max);
        static void        rand_scalar_c_bin(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                                mp_complex& s2, mp_float& s3);
        static void        rand_scalar_c_pow_cc(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                                mp_complex& s2, mp_float& s3);
        static void        rand_scalar_c_pow_cr(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                                mp_complex& s2, mp_float& s3);
        static void        rand_scalar_c_pow_rc(const prec_vec& prec, const max_vec& max, mp_complex& s1,
                                mp_complex& s2, mp_float& s3);

        static Real         calc_ulp(const mp_float& res, const mp_float& res_ext);

    private:
        void        test_constants(std::ostream& os);
        void        test_real(std::ostream& os, Integer N, const prec_vec& prec, 
                        const max_vec& max);
        void        test_complex(std::ostream& os, Integer N, const prec_vec& prec, 
                        const max_vec& max);
        void        test_complex_special(std::ostream& os, Integer N, const prec_vec& prec);
        void        test_complex_bin_special(std::ostream& os, Integer N, const prec_vec& prec);

        bool        test_constants(std::ostream& os, precision prec);
        bool        test_real(std::ostream& os, Integer N, precision prec, double max);
        bool        test_complex(std::ostream& os, Integer N, precision prec, 
                        double max_re, double max_im);
        bool        test_complex_special(std::ostream& os, Integer N, precision prec, double max_v);
        bool        test_complex_bin_special(std::ostream& os, Integer N, precision prec, 
                        const prec_vec& prec_v, const max_vec& max);

        void        test_bin_real(std::ostream& os, Integer N, const prec_vec& prec,
                                  const max_vec& max);
        bool        test_bin_real(std::ostream& os, Integer N, precision prec, const max_vec& max);
        void        test_bin_complex(std::ostream& os, Integer N, const prec_vec& prec, const max_vec& max);
        void        test_bin_complex(std::ostream& os, Integer N, precision p, const prec_vec& prec,
                                  const max_vec& max);

        template<class Func>
        double      test_scalar(const mp_float& s, precision p, Integer code);

        template<class Func>
        double      test_constant(precision p);

        template<class Func>
        double      test_bin(const mp_float& s1, const mp_float& s2, precision p, Integer code);

        template<class Func>
        dtup_2      test_scalar_c(const mp_complex& s, precision p, Integer code);

        template<class Func, class T1, class T2>
        dtup_2      test_bin_c(const T1& s1, const T2& s2, precision p, Integer code);

        template<class Func>
        bool        test_scalar_func(std::ostream& os, Integer N, precision prec, double max,
                        double accuracy);

        template<class Func>
        bool        test_constant(std::ostream& os, precision p, double accuracy);

        template<class Func>
        bool        test_bin_func(std::ostream& os, Integer N, precision prec, const max_vec& max,
                        double accuracy);

        template<class Func, class Rand>
        bool        test_bin_func_c(std::ostream& os, Integer N, precision p, const prec_vec& prec, 
                        const max_vec& max, Real accuracy);

        template<class Func, class Rand>
        bool        test_scalar_func_c(std::ostream& os, Integer N, precision prec, 
                        double max_re, double max_im, double accuracy);        
};

}};
