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
#include "scalar_ext.h"
#include <set>

namespace matcl { namespace test
{

namespace mdy = matcl::dynamic;

void test_gmp_object();

class gmp_object_tester
{
    private:
        using type2     = std::pair<mdy::Type, mdy::Type>;
        using type2_set = std::set<type2>;

    public:
        void    make();

    private:
        bool    test_cons_impl();
        double  test_cons_1(mdy::Type t1, mdy::Type t2, const type2_set& missing_cons, 
                    const type2_set& missing_assign);
        bool    test_cons_val_impl();
        double  test_cons_val_1(const Scalar_ext& s1, const Scalar_ext& s2, Integer code);
        bool    test_object_func_impl();
        double  test_func_1(mdy::Type t1);
        bool    test_object_func_val_impl();
        double  test_func_1_val(const Scalar_ext& s1, Integer code);
        bool    test_func_compare_impl();
        double  test_compare_1(mdy::Type t1, mdy::Type t2, const type2_set& missing);
        bool    test_func_uminus_impl();
        double  test_func_uminus_1(mdy::Type t1);
        bool    test_func_uminus_val_impl();
        double  test_func_uminus_val(const Scalar_ext& s1, Integer code);
        bool    test_func_reim_impl();
        double  test_func_reim_1(mdy::Type t1);
        bool    test_func_reim_val_impl();
        double  test_func_reim_val(const Scalar_ext& s1, Integer code);
        bool    test_func_is_impl();
        double  test_func_is(const Scalar_ext& s1, Integer code);
        bool    test_func_next_impl();
        double  test_func_next(const Scalar_ext& s1, Integer code);
        bool    test_func_sign_impl();
        double  test_func_sign(const Scalar_ext& s1, Integer code);
        bool    test_func_eps_impl();
        double  test_func_eps(const Scalar_ext& s1, Integer code);
        bool    test_func_unify_impl();
        double  test_unify_1(mdy::Type t1, mdy::Type t2, Integer code);

        void    test_cons();
        void    test_cons_val();
        void    test_object_func();
        void    test_object_func_val();
        void    test_func_compare();
        void    test_func_uminus();
        void    test_func_uminus_val();
        void    test_func_reim();
        void    test_func_reim_val();
        void    test_func_is();
        void    test_func_next();
        void    test_func_sign();
        void    test_func_eps();

        template<class Func, bool With_mp = true>
        void    test_scalar_func();

        template<class Func, bool With_mp = true>
        double  test_scalar(const Scalar_ext& s, Integer code);

        void    test_sqrt(); 
        void    test_cbrt(); 
        void    test_sqrt_c(); 
        void    test_exp(); 
        void    test_expm1(); 
        void    test_expi(); 
        void    test_exp2(); 
        void    test_exp10(); 
        void    test_log(); 
        void    test_log1p(); 
        void    test_log2(); 
        void    test_log10(); 
        void    test_log_c(); 
        void    test_log1p_c(); 
        void    test_log2_c(); 
        void    test_log10_c(); 
        void    test_sin();
        void    test_cos();
        void    test_tan();    
        void    test_cot();
        void    test_sec();
        void    test_csc();
        void    test_sinh();
        void    test_cosh();
        void    test_tanh();
        void    test_coth();
        void    test_sech();
        void    test_csch();
        void    test_asin();
        void    test_asin_c();
        void    test_acos();
        void    test_acos_c();
        void    test_atan();
        void    test_acot();
        void    test_asec();
        void    test_asec_c();
        void    test_acsc();
        void    test_acsc_c();
        void    test_asinh();
        void    test_acosh();
        void    test_acosh_c();
        void    test_atanh();
        void    test_atanh_c();
        void    test_acoth();
        void    test_acoth_c();
        void    test_asech();
        void    test_asech_c();
        void    test_acsch();
        void    test_inv();
        void    test_invs();
        void    test_floor();
        void    test_ceil();
        void    test_round();
        void    test_trunc();
        void    test_ifloor();
        void    test_iceil();
        void    test_iround();
        void    test_itrunc();
        void    test_sign();
        void    test_op_neg();
        void    test_op_true();

        void    test_ldexp();
        void    test_scalbn();
        void    test_frexp();
        void    test_modf_frac();
        void    test_modf_int();
        void    test_logb();
        void    test_ilogb();

        void    test_func_unify();
        void    test_combinatorics();
};

}};
