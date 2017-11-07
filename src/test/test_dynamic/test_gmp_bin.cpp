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

#include "test_gmp.h"
#include "matcl-mp/matcl_mp.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-scalar/objects/object_functions.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

#include "matcl-scalar/lib_functions/func_binary.h"
#include "matcl-scalar/lib_functions/func_matrix.h"
#include "matcl-scalar/lib_functions/func_forwarding.h"
#include "matcl-scalar/IO/scalar_io.h"

#include "rand_scalars.h"
#include "eval_cons.h"
#include "utils.h"
#include <iostream>

namespace matcl { namespace test
{

namespace mdy = matcl::dynamic;

void test_gmp_bin()
{
    gmp_tester_bin test;
    test.make();
};

void gmp_tester_bin::make()
{
    test_min_obj();

    test_copysign();
    test_nextafter();    

    test_hypot();
    test_atan2();    
    test_mod();
    test_rem();
    test_fdim();

    test_copysign_obj();
    test_nextafter_obj();    
    test_float_distance_obj();    
    test_hypot_obj();
    test_atan2_obj();    
    test_mod_obj();
    test_rem_obj();
    test_fdim_obj();

    test_add();
    test_sub();
    test_mul();
    test_mmul();
    test_div();
    test_div2();
    test_div_0();
    test_div_1();
    test_idiv();
    test_pow();
    test_pow_c();
    test_min();
    test_max();
    test_op_and();
    test_op_or();
    test_op_xor();
    test_elem_and();
    test_elem_or();
    test_elem_xor();

    test_gt();
    test_lt();    
    test_leq();
    test_geq();
    test_neq();
    test_eeq();
    test_neq_nan();
    test_eeq_nan();

    test_add_obj();
    test_sub_obj();
    test_mul_obj();
    test_mmul_obj();
    test_div_obj();
    test_div2_obj();
    test_div_0_obj();
    test_div_1_obj();
    test_idiv_obj();

    test_eeq_obj(); 
    test_gt_obj();
    test_lt_obj();    
    test_leq_obj();
    test_geq_obj();
    test_neq_obj();
    test_neq_nan_obj();
    test_eeq_nan_obj();
    test_pow_obj();
    test_pow_c_obj();

    test_min_obj();
    test_max_obj();
    test_op_and_obj();
    test_op_or_obj();
    test_op_xor_obj();
    test_elem_and_obj();
    test_elem_or_obj();
    test_elem_xor_obj();

    //TODO
    //test_float_distance();
};

void gmp_tester_bin::test_add()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_add(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "add: ok" << "\n";
    else
        out_stream << "add" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_sub()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_sub(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "sub: ok" << "\n";
    else
        out_stream << "sub" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_mul()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_mul(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "mul: ok" << "\n";
    else
        out_stream << "mul" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_mmul()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_mmul(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "mmul: ok" << "\n";
    else
        out_stream << "mmul" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_div()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div: ok" << "\n";
    else
        out_stream << "div" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_div2()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div2(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div2: ok" << "\n";
    else
        out_stream << "div2" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_div_0()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div_0(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div_0: ok" << "\n";
    else
        out_stream << "div_0" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_div_1()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div_1(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div_1: ok" << "\n";
    else
        out_stream << "div_1" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_idiv()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_idiv(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "idiv: ok" << "\n";
    else
        out_stream << "idiv" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_pow()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_pow(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "pow: ok" << "\n";
    else
        out_stream << "pow" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_pow_c()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_pow_c(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "pow_c: ok" << "\n";
    else
        out_stream << "pow_c" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_atan2()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_atan2(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "atan2: ok" << "\n";
    else
        out_stream << "atan2" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_mod()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_mod(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "mod: ok" << "\n";
    else
        out_stream << "mod" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_rem()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_rem(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "rem: ok" << "\n";
    else
        out_stream << "rem" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_min()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_min(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "min: ok" << "\n";
    else
        out_stream << "min" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_max()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_max(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "max: ok" << "\n";
    else
        out_stream << "max" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_op_and()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_op_and(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "op_and: ok" << "\n";
    else
        out_stream << "op_and" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_op_or()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_op_or(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "op_or: ok" << "\n";
    else
        out_stream << "op_or" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_op_xor()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_op_xor(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "op_xor: ok" << "\n";
    else
        out_stream << "op_xor" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_elem_and()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_elem_and(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "elem_and: ok" << "\n";
    else
        out_stream << "elem_and" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_elem_or()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_elem_or(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "elem_or: ok" << "\n";
    else
        out_stream << "elem_or" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_elem_xor()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_elem_xor(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "elem_xor: ok" << "\n";
    else
        out_stream << "elem_xor" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_fdim()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_fdim(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "fdim: ok" << "\n";
    else
        out_stream << "fdim" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_hypot()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_hypot(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "hypot: ok" << "\n";
    else
        out_stream << "hypot" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_nextafter()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_nextafter(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "nextafter: ok" << "\n";
    else
        out_stream << "nextafter" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_float_distance()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_float_distance(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "float_distance: ok" << "\n";
    else
        out_stream << "float_distance" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_copysign()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_copysign(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "copysign: ok" << "\n";
    else
        out_stream << "copysign" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_eeq()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_eeq(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "eeq: ok" << "\n";
    else
        out_stream << "eeq" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_eeq_nan()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_eeq_nan(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "eeq_nan: ok" << "\n";
    else
        out_stream << "eeq_nan" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_copysign_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_copysign_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "copysign obj: ok" << "\n";
    else
        out_stream << "copysign obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_nextafter_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_nextafter_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "nextafter obj: ok" << "\n";
    else
        out_stream << "nextafter obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_float_distance_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_float_distance_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "float_distance obj: ok" << "\n";
    else
        out_stream << "float_distance obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_atan2_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_atan2_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "atan2 obj: ok" << "\n";
    else
        out_stream << "atan2 obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_mod_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_mod_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "mod obj: ok" << "\n";
    else
        out_stream << "mod obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_rem_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_rem_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "rem obj: ok" << "\n";
    else
        out_stream << "rem obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_max_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_max_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "max obj: ok" << "\n";
    else
        out_stream << "max obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_min_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_min_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "min obj: ok" << "\n";
    else
        out_stream << "min obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_op_or_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_op_or_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "op_or obj: ok" << "\n";
    else
        out_stream << "op_or obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_op_xor_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_op_xor_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "op_xor obj: ok" << "\n";
    else
        out_stream << "op_xor obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_op_and_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_op_and_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "op_and obj: ok" << "\n";
    else
        out_stream << "op_and obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_elem_or_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_elem_or_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "elem_or obj: ok" << "\n";
    else
        out_stream << "elem_or obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_elem_xor_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_elem_xor_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "elem_xor obj: ok" << "\n";
    else
        out_stream << "elem_xor obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_elem_and_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_elem_and_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "elem_and obj: ok" << "\n";
    else
        out_stream << "elem_and obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_fdim_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_fdim_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "fdim obj: ok" << "\n";
    else
        out_stream << "fdim obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_hypot_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_hypot_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "hypot obj: ok" << "\n";
    else
        out_stream << "hypot obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_add_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_add_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "add obj: ok" << "\n";
    else
        out_stream << "add obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_sub_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_sub_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "sub obj: ok" << "\n";
    else
        out_stream << "sub obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_mul_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_mul_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "mul obj: ok" << "\n";
    else
        out_stream << "mul obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_mmul_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_mmul_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "mmul obj: ok" << "\n";
    else
        out_stream << "mmul obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_div_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div obj: ok" << "\n";
    else
        out_stream << "div obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_div2_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div2_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div2 obj: ok" << "\n";
    else
        out_stream << "div2 obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_div_0_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div_0_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div_0 obj: ok" << "\n";
    else
        out_stream << "div_0 obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_div_1_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_div_1_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "div_1 obj: ok" << "\n";
    else
        out_stream << "div_1 obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_idiv_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_idiv_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "idiv obj: ok" << "\n";
    else
        out_stream << "idiv obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_pow_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_pow_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "pow obj: ok" << "\n";
    else
        out_stream << "pow obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_pow_c_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_pow_c_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "pow_c obj: ok" << "\n";
    else
        out_stream << "pow_c obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_eeq_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_eeq_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "eeq obj: ok" << "\n";
    else
        out_stream << "eeq obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_eeq_nan_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_eeq_nan_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "eeq_nan obj: ok" << "\n";
    else
        out_stream << "eeq_nan obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_neq()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_neq(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "neq: ok" << "\n";
    else
        out_stream << "neq" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_neq_nan()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_neq_nan(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "neq_nan: ok" << "\n";
    else
        out_stream << "neq_nan" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_neq_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_neq_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "neq obj: ok" << "\n";
    else
        out_stream << "neq obj" << ": FAILED" << "\n";
};
void gmp_tester_bin::test_neq_nan_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_neq_nan_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "neq_nan obj: ok" << "\n";
    else
        out_stream << "neq_nan obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_leq()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_leq(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "leq: ok" << "\n";
    else
        out_stream << "leq" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_leq_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_leq_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "leq obj: ok" << "\n";
    else
        out_stream << "leq obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_geq()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_geq(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "geq: ok" << "\n";
    else
        out_stream << "geq" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_geq_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_geq_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "geq obj: ok" << "\n";
    else
        out_stream << "geq obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_lt()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_lt(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "lt: ok" << "\n";
    else
        out_stream << "lt" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_lt_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_lt_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "lt obj: ok" << "\n";
    else
        out_stream << "lt obj" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_gt()
{
    std::vector<Scalar> scalars1;
    std::vector<Scalar> scalars2;

    Integer N   = 1000;
    rand_scalars::make(scalars1, N, false);
    rand_scalars::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_gt(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "gt: ok" << "\n";
    else
        out_stream << "gt" << ": FAILED" << "\n";
};

void gmp_tester_bin::test_gt_obj()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, false);
    rand_scalars_ext::make(scalars2, N, false);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_gt_obj(scalars1[i],scalars2[i], i);

    if (res == 0.0)
        out_stream << "gt obj: ok" << "\n";
    else
        out_stream << "gt obj" << ": FAILED" << "\n";
};

template<class Derived>
struct eval_operator : eval_scalars<eval_operator<Derived>>
{
    Integer code;

    eval_operator(Integer c) : code(c){};

    template<class T1, class T2>
    bool eval_op(const T1& a1, const T2& a2)
    {
        return Derived::eval(a1,a2);
    };

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        //out_stream << code << "\n";
        if (code == -1)
            disp("break");

        bool res        = eval_op(s1,s2);
        bool is_real_1  = matcl::imag(s1) == 0;
        bool is_real_2  = matcl::imag(s2) == 0;
        bool is_fin_1   = matcl::is_finite(s1) && is_real_1;
        bool is_fin_2   = matcl::is_finite(s2) && is_real_2;

        mp_int      v1_1    = convert_scalar<mp_int>(s1);
        mp_float    v1_2    = convert_scalar<mp_float>(s1);
        mp_rational v1_3    = convert_scalar<mp_rational>(s1);
        mp_complex  v1_4    = convert_scalar<mp_complex>(s1);

        mp_int v2_1         = convert_scalar<mp_int>(s2);
        mp_float v2_2       = convert_scalar<mp_float>(s2);
        mp_rational v2_3    = convert_scalar<mp_rational>(s2);
        mp_complex v2_4     = convert_scalar<mp_complex>(s2);

        double out = 0;

        if (eval_op(v1_1, s2) != res && is_fin_1 == true)
            out     += 1;
        if (eval_op(v1_2, s2) != res && is_real_1 == true)
            out     += 1;
        if (eval_op(v1_3, s2) != res && is_fin_1 == true)
            out     += 1;
        if (eval_op(v1_4, s2) != res)
            out     += 1;

        if (eval_op(s1, v2_1) != res && is_fin_2 == true)
            out     += 1;
        if (eval_op(s1, v2_2) != res && is_real_2 == true)
            out     += 1;
        if (eval_op(s1, v2_3) != res && is_fin_2 == true)
            out     += 1;
        if (eval_op(s1, v2_4) != res)
            out     += 1;

        if (eval_op(v1_1, v2_1) != res && is_fin_1 == true && is_fin_2 == true)
            out     += 1;
        if (eval_op(v1_2, v2_1) != res && is_real_1 == true && is_fin_2 == true)
            out     += 1;
        if (eval_op(v1_3, v2_1) != res && is_fin_1 == true && is_fin_2 == true)
            out     += 1;
        if (eval_op(v1_4, v2_1) != res && is_fin_2 == true)
            out     += 1;

        if (eval_op(v1_1, v2_2) != res && is_fin_1 == true && is_real_2 == true)
            out     += 1;
        if (eval_op(v1_2, v2_2) != res && is_real_1 == true && is_real_2 == true)
            out     += 1;
        if (eval_op(v1_3, v2_2) != res && is_fin_1 == true && is_real_2 == true)
            out     += 1;
        if (eval_op(v1_4, v2_2) != res && is_real_2 == true)
            out     += 1;

        if (eval_op(v1_1, v2_3) != res && is_fin_1 == true && is_fin_2 == true)
            out     += 1;
        if (eval_op(v1_2, v2_3) != res && is_real_1 == true && is_fin_2 == true)
            out     += 1;
        if (eval_op(v1_3, v2_3) != res && is_fin_1 == true && is_fin_2 == true)
            out     += 1;
        if (eval_op(v1_4, v2_3) != res && is_fin_2 == true)
            out     += 1;

        if (eval_op(v1_1, v2_4) != res && is_fin_1 == true)
            out     += 1;
        if (eval_op(v1_2, v2_4) != res && is_real_1 == true)
            out     += 1;
        if (eval_op(v1_3, v2_4) != res && is_fin_1 == true)
            out     += 1;
        if (eval_op(v1_4, v2_4) != res)
            out     += 1;

        if (out != 0)
            out_stream << code << " " << s1 << " " << s2 << "\n";

        return out;
    };
};

struct eval_eeq : eval_operator<eval_eeq>
{
    eval_eeq(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return a1 == a2;
    };
};
struct eval_eeq_nan : eval_operator<eval_eeq_nan>
{
    eval_eeq_nan(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return eeq_nan(a1, a2);
    };
};

struct eval_neq : eval_operator<eval_neq>
{
    eval_neq(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return a1 != a2;
    };
};
struct eval_neq_nan : eval_operator<eval_neq_nan>
{
    eval_neq_nan(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return neq_nan(a1, a2);
    };
};
struct eval_leq : eval_operator<eval_leq>
{
    eval_leq(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return a1 <= a2;
    };
};

struct eval_geq : eval_operator<eval_geq>
{
    eval_geq(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return a1 >= a2;
    };
};

struct eval_lt : eval_operator<eval_lt>
{
    eval_lt(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return a1 < a2;
    };
};

struct eval_gt : eval_operator<eval_gt>
{
    eval_gt(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return a1 > a2;
    };
};
double gmp_tester_bin::test_eeq(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_eeq test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_eeq_nan(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_eeq_nan test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_neq(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_neq test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_neq_nan(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_neq_nan test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_leq(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_leq test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_geq(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_geq test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_lt(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_lt test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_gt(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_gt test(code);
    return test.make(s1, s2);
};

template<class Derived>
struct eval_operator_obj : eval_scalars<eval_operator_obj<Derived>>
{
    Integer code;

    eval_operator_obj(Integer c) : code(c){};

    template<class T1, class T2>
    bool eval_op(const T1& a1, const T2& a2)
    {
        return Derived::eval(a1,a2);
    };

    template<class T1, class T2>
    double eval_scal_func(const T1& s1, const T2& s2)
    {
        //out_stream << code << "\n";
        //if (code == 0)
        //    disp("break");

        bool res_1      = eval_op(s1,s2);

        using type_1    = object_type<T1>;
        using type_2    = object_type<T2>;

        type_1 o1(s1);
        type_2 o2(s2);

        bool res_2      = eval_op(o1, o2);
        bool res_21     = eval_op(s1, o2);
        bool res_22     = eval_op(o1, s2);

        matcl::dynamic::object oo1      = matcl::dynamic::object(o1);
        matcl::dynamic::object oo2      = matcl::dynamic::object(o2);

        bool res_3      = eval_op(oo1, oo2);

        double out      = 0;
        if (res_1 != res_2)
            out         += 1;
        if (res_1 != res_21)
            out         += 1;
        if (res_1 != res_22)
            out         += 1;
        if (res_1 != res_3)
            out         += 1;

        if (out != 0.0)
        {
            out_stream << s1 << " " << s2 << "\n";
            out_stream << code << "\n";
        }
        return out;
    };
};

struct eval_eeq_obj : eval_operator_obj<eval_eeq_obj>
{
    eval_eeq_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)(a1 == a2);
    };
};
struct eval_neq_obj : eval_operator_obj<eval_neq_obj>
{
    eval_neq_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)(a1 != a2);
    };
};
struct eval_eeq_nan_obj : eval_operator_obj<eval_eeq_nan_obj>
{
    eval_eeq_nan_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)eeq_nan(a1, a2);
    };
};
struct eval_neq_nan_obj : eval_operator_obj<eval_neq_nan_obj>
{
    eval_neq_nan_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)neq_nan(a1, a2);
    };
};

struct eval_leq_obj : eval_operator_obj<eval_leq_obj>
{
    eval_leq_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)(a1 <= a2);
    };
};
struct eval_geq_obj : eval_operator_obj<eval_geq_obj>
{
    eval_geq_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)(a1 >= a2);
    };
};
struct eval_lt_obj : eval_operator_obj<eval_lt_obj>
{
    eval_lt_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)(a1 < a2);
    };
};
struct eval_gt_obj : eval_operator_obj<eval_gt_obj>
{
    eval_gt_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return (bool)(a1 > a2);
    };
};
double gmp_tester_bin::test_eeq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_eeq_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_neq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_neq_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_eeq_nan_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_eeq_nan_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_neq_nan_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_neq_nan_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_leq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_leq_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_geq_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_geq_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_lt_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_lt_obj test(code);
    return test.make(s1, s2);
};

double gmp_tester_bin::test_gt_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_gt_obj test(code);
    return test.make(s1, s2);
};

struct eval_op_or : eval_operator<eval_op_or>
{
    eval_op_or(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return op_or(a1, a2);
    };
};
struct eval_op_xor : eval_operator<eval_op_xor>
{
    eval_op_xor(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return op_xor(a1, a2);
    };
};

struct eval_op_and : eval_operator<eval_op_and>
{
    eval_op_and(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return op_and(a1, a2);
    };
};

struct eval_elem_or : eval_operator<eval_elem_or>
{
    eval_elem_or(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return elem_or(a1, a2);
    };
};

struct eval_elem_xor : eval_operator<eval_elem_xor>
{
    eval_elem_xor(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return elem_xor(a1, a2);
    };
};

struct eval_elem_and : eval_operator<eval_elem_and>
{
    eval_elem_and(Integer c) : eval_operator(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return elem_and(a1, a2);
    };
};

double gmp_tester_bin::test_op_or(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_op_or test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_op_xor(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_op_xor test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_op_and(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_op_and test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_elem_or(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_elem_or test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_elem_xor(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_elem_xor test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_elem_and(const Scalar& s1, const Scalar& s2, Integer code)
{
    eval_elem_and test(code);
    return test.make(s1, s2);
};

struct eval_op_or_obj : eval_operator_obj<eval_op_or_obj>
{
    eval_op_or_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return op_or(a1, a2);
    };
};

struct eval_op_xor_obj : eval_operator_obj<eval_op_xor_obj>
{
    eval_op_xor_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return op_xor(a1, a2);
    };
};

struct eval_op_and_obj : eval_operator_obj<eval_op_and_obj>
{
    eval_op_and_obj(Integer c) : eval_operator_obj(c){};

    template<class T1, class T2>
    static bool eval(const T1& a1, const T2& a2)
    {
        return op_and(a1, a2);
    };
};

double gmp_tester_bin::test_op_or_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_op_or_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_op_xor_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_op_xor_obj test(code);
    return test.make(s1, s2);
};
double gmp_tester_bin::test_op_and_obj(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_op_and_obj test(code);
    return test.make(s1, s2);
};


}};
