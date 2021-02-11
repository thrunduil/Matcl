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

#include "matcl-mp/matcl_mp.h"
#include "matcl-scalar/matcl_scalar.h"
//#include "test_gmp.h"

//#include <vector>

namespace matcl { namespace test
{

void test_scal_accuracy(std::ostream& os);

class scal_accuracy_tester
{
    private:
        using max_vec   = std::vector<double>;

    public:
        void    make(std::ostream& os);

    private:
        template<class Type>
        void    test_unary(std::ostream& os, Integer N, const max_vec& max_v);

        template<class Type>
        void    test_unary_compl(std::ostream& os, Integer N, const max_vec& max_v);

        void    test_binary(std::ostream& os);

        template<class Type, class Func>
        void    test_unary_func(formatted_disp& os, Integer N, const max_vec& max_v);

        template<class Type, class Func>
        void    test_unary_compl_func(formatted_disp& os, Integer N, const max_vec& max_v);

        template<class Type, class Func>
        void    test_unary_func(formatted_disp& os, Integer N, Real max_v, Real& matcl_ulp,
                    Real& std_ulp);

        template<class Type, class Func>
        void    test_unary_compl_func(formatted_disp& os, Integer N, Real max_re, Real max_im,
                    Real& matcl_ulp, Real& std_ulp);

        template<class Type, class Func>
        void    test_scalar(Type s, Integer code, Real& ulp_matcl, Real& ulp_std);

        template<class Type, class Func>
        void    test_scalar_complex(const Type& s, Integer code, Real& ulp_matcl, Real& ulp_std);

    public:
        template<class Type>
        static Type         rand_scalar(double max);

        template<class Type>
        static Real         calc_ulp(Type res, const mp_float& res_ext);

        template<class Type>
        static Real         calc_ulp_complex(const Type& res, const mp_complex& res_ext);

        template<class Type>
        static precision    prec_type();
};

}};
