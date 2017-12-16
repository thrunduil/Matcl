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

#include "test_simd_config.h"
#include "matcl-scalar/IO/formatted_disp.h"
#include "matcl-mp/matcl_mp.h"
#include "accuracy/test_accuracy.h"

namespace matcl { namespace test
{

void    test_math_accuracy();
void    set_rand_denormals(bool val);

class simd_accuracy_tester
{
    private:
        template<class Val>
        using Func      = matcl::test_accuracy_function<Val>;

    public:
        void    make();

        static 
        bool    rand_denormals();

    private:
        template<class Type>
        void    test_unary(int N);

        template<class Type>
        void    test_unary_func(formatted_disp& os, const Func<Type>& f, int N);

        template<class Type>
        void    test_unary_func(formatted_disp& os, const Func<Type>& f, int N, 
                    double& ulp_base, double& ulp_ref);

        template<class Type>
        void    test_scalar(const Func<Type>& f, Type s, int code, double& ulp_base, double& ulp_ref);

    public:
        template<class Type>
        static double   calc_ulp(Type res, const mp_float& res_ext);

        template<class Type>
        static int      precision_lost_subnormal(Type val);
};

}};
