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

#pragma warning(push)
#pragma warning(disable: 4244)  //conversion from 'const int' to 'matcl::Float', possible loss of data
#pragma warning(disable:4800)//forcing value to bool 'true' or 'false' (performance warning)

#include "test_gmp_object.h"
#include "test_gmp.h"
#include "matcl-mp-obj/mp_object.h"
#include "matcl-dynamic/type.h"
#include "rand_scalars.h"
#include "eval_cons.h"
#include "matcl-scalar/IO/scalar_io.h"

#pragma warning(pop)

#include <iostream>
#include <set>

namespace matcl { namespace test
{

void test_gmp_object()
{
    gmp_object_tester test;
    test.make();
};

void gmp_object_tester::make()
{
    test_cons();
    test_cons_val();
    test_object_func();
    test_object_func_val();
    test_func_compare();
    test_func_uminus();
    test_func_uminus_val();
    test_func_reim();
    test_func_reim_val();
    test_func_is();
    test_func_next();
    test_func_sign();
    test_func_eps();

    test_sqrt();
    test_cbrt();
    test_sqrt_c();
    test_exp();
    test_expm1();
    test_expi();
    test_exp2();
    test_exp10();
    test_log(); 
    test_log2(); 
    test_log10(); 
    test_log1p(); 
    test_log_c(); 
    test_log2_c(); 
    test_log10_c(); 
    test_log1p_c(); 
    test_sin();
    test_cos();
    test_tan();    
    test_cot();
    test_sec();
    test_csc();
    test_sinh();
    test_cosh();
    test_tanh();
    test_coth();
    test_sech();
    test_csch();  
    test_asin();
    test_asin_c();
    test_acos();
    test_acos_c();
    test_atan();
    test_acot();
    test_asec();
    test_asec_c();
    test_acsc();
    test_acsc_c();
    test_asinh();
    test_acosh();
    test_acosh_c();    
    test_atanh();
    test_atanh_c();  
    test_acoth();
    test_acoth_c();
    test_asech();
    test_asech_c();
    test_acsch();
    test_inv();
    test_invs();
    test_floor();
    test_ceil();
    test_round();
    test_trunc();
    test_ifloor();
    test_iceil();
    test_iround();
    test_itrunc();
    test_sign();
    test_op_neg();
    test_op_true();

    test_ldexp();
    test_scalbn();
    test_frexp();
    test_modf_frac();
    test_modf_int();
    test_logb();
    test_ilogb();

    test_func_unify();
    test_combinatorics();
};

void gmp_object_tester::test_cons()
{
    bool res;

    try
    {
        res = test_cons_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_cons EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "constructors: " << "ok" << "\n";
    else
        out_stream << "constructors: " << "FAILED" << "\n";
};
void gmp_object_tester::test_cons_val()
{
    bool res;

    try
    {
        res = test_cons_val_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_cons_val EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "constructors val: " << "ok" << "\n";
    else
        out_stream << "constructors val: " << "FAILED" << "\n";
};

void gmp_object_tester::test_object_func()
{
    bool res;

    try
    {
        res = test_object_func_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_object_func EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "object_func: " << "ok" << "\n";
    else
        out_stream << "object_func: " << "FAILED" << "\n";
};

void gmp_object_tester::test_object_func_val()
{
    bool res;

    try
    {
        res = test_object_func_val_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_object_func_val EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "object_func val: " << "ok" << "\n";
    else
        out_stream << "object_func val: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_compare()
{
    bool res;

    try
    {
        res = test_func_compare_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_compare EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "compare_func: " << "ok" << "\n";
    else
        out_stream << "compare_func: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_unify()
{
    bool res;

    try
    {
        res = test_func_unify_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_unify EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "unify mp: " << "ok" << "\n";
    else
        out_stream << "unify mp: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_uminus()
{
    bool res;

    try
    {
        res = test_func_uminus_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_uminus EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "func_uminus: " << "ok" << "\n";
    else
        out_stream << "func_uminus: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_uminus_val()
{
    bool res;

    try
    {
        res = test_func_uminus_val_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_uminus_val EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "func_uminus val: " << "ok" << "\n";
    else
        out_stream << "func_uminus val: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_reim()
{
    bool res;

    try
    {
        res = test_func_reim_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_reim EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "reim: " << "ok" << "\n";
    else
        out_stream << "reim: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_is()
{
    bool res;

    try
    {
        res = test_func_is_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_is EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "is obj: " << "ok" << "\n";
    else
        out_stream << "is obj: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_next()
{
    bool res;

    try
    {
        res = test_func_next_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_next EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "next obj: " << "ok" << "\n";
    else
        out_stream << "next obj: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_sign()
{
    bool res;

    try
    {
        res = test_func_sign_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_sign EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "sign obj: " << "ok" << "\n";
    else
        out_stream << "sign obj: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_eps()
{
    bool res;

    try
    {
        res = test_func_eps_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_eps EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "eps obj: " << "ok" << "\n";
    else
        out_stream << "eps obj: " << "FAILED" << "\n";
};

void gmp_object_tester::test_func_reim_val()
{
    bool res;

    try
    {
        res = test_func_reim_val_impl();
    }
    catch(std::exception& ex)
    {
        out_stream << "test_func_reim_val EXCEPTION: " << ex.what() << "\n";
    };

    if (res == true)
        out_stream << "reim val: " << "ok" << "\n";
    else
        out_stream << "reim val: " << "FAILED" << "\n";
};

template<class Func, bool With_mp>
double gmp_object_tester::test_scalar(const Scalar_ext& s, Integer code)
{
    eval_scalar_func_templ_obj<Func,With_mp> test1(code);
    double res1 = test1.make(s);
    return res1;
};

template<class Func, bool With_mp>
void gmp_object_tester::test_scalar_func()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_scalar<Func, With_mp>(scalars1[i], i);

    if (res != 0.0)
        out_stream << Func::name() << " obj: FAILED" << "\n";
    else
        out_stream << Func::name() << " obj: ok" << "\n";        
};

void gmp_object_tester::test_ldexp()
{
    return test_scalar_func<Ldexp3p_func_t>();
};
void gmp_object_tester::test_scalbn()
{
    return test_scalar_func<Scalbn3p_func_t>();
};
void gmp_object_tester::test_frexp()
{
    return test_scalar_func<Frexp_func_t>();
};
void gmp_object_tester::test_modf_frac()
{
    return test_scalar_func<Modf_frac_func_t>();
};
void gmp_object_tester::test_modf_int()
{
    return test_scalar_func<Modf_int_func_t>();
};

void gmp_object_tester::test_logb()
{
    return test_scalar_func<Logb_func>();
};
void gmp_object_tester::test_ilogb()
{
    return test_scalar_func<Ilogb_func>();
};

void gmp_object_tester::test_sqrt()
{
    return test_scalar_func<Sqrt_func>();
};
void gmp_object_tester::test_cbrt()
{
    return test_scalar_func<Cbrt_func>();
};

void gmp_object_tester::test_sqrt_c()
{
    return test_scalar_func<Sqrt_c_func>();
};

void gmp_object_tester::test_exp()
{
    return test_scalar_func<Exp_func>();
};
void gmp_object_tester::test_expm1()
{
    return test_scalar_func<Expm1_func>();
};
void gmp_object_tester::test_expi()
{
    return test_scalar_func<Expi_func>();
};

void gmp_object_tester::test_exp2()
{
    return test_scalar_func<Exp2_func>();
};
void gmp_object_tester::test_exp10()
{
    return test_scalar_func<Exp10_func>();
};

void gmp_object_tester::test_log()
{
    return test_scalar_func<Log_func>();
};
void gmp_object_tester::test_log1p()
{
    return test_scalar_func<Log1p_func>();
};
void gmp_object_tester::test_log2()
{
    return test_scalar_func<Log2_func>();
};
void gmp_object_tester::test_log10()
{
    return test_scalar_func<Log10_func>();
};
void gmp_object_tester::test_log_c()
{
    return test_scalar_func<Log_c_func>();
};
void gmp_object_tester::test_log1p_c()
{
    return test_scalar_func<Log1p_c_func>();
};
void gmp_object_tester::test_log2_c()
{
    return test_scalar_func<Log2_c_func>();
};
void gmp_object_tester::test_log10_c()
{
    return test_scalar_func<Log10_c_func>();
};
void gmp_object_tester::test_sin()
{
    return test_scalar_func<Sin_func>();
};
void gmp_object_tester::test_cos()
{
    return test_scalar_func<Cos_func>();
};
void gmp_object_tester::test_tan()
{
    return test_scalar_func<Tan_func>();
};
void gmp_object_tester::test_cot()
{
    return test_scalar_func<Cot_func>();
};
void gmp_object_tester::test_sec()
{
    return test_scalar_func<Sec_func>();
};
void gmp_object_tester::test_csc()
{
    return test_scalar_func<Csc_func>();
};
void gmp_object_tester::test_sinh()
{
    return test_scalar_func<Sinh_func>();
};
void gmp_object_tester::test_cosh()
{
    return test_scalar_func<Cosh_func>();
};
void gmp_object_tester::test_tanh()
{
    return test_scalar_func<Tanh_func>();
};
void gmp_object_tester::test_coth()
{
    return test_scalar_func<Coth_func>();
};
void gmp_object_tester::test_sech()
{
    return test_scalar_func<Sech_func>();
};
void gmp_object_tester::test_csch()
{
    return test_scalar_func<Csch_func>();
};
void gmp_object_tester::test_asin()
{
    return test_scalar_func<Asin_func>();
};
void gmp_object_tester::test_asin_c()
{
    return test_scalar_func<Asin_c_func>();
};
void gmp_object_tester::test_acos()
{
    return test_scalar_func<Acos_func>();
};
void gmp_object_tester::test_acos_c()
{
    return test_scalar_func<Acos_c_func>();
};
void gmp_object_tester::test_atan()
{
    return test_scalar_func<Atan_func>();
};
void gmp_object_tester::test_acot()
{
    return test_scalar_func<Acot_func,false>();
};
void gmp_object_tester::test_asec()
{
    return test_scalar_func<Asec_func,false>();
};
void gmp_object_tester::test_asec_c()
{
    return test_scalar_func<Asec_c_func,false>();
};
void gmp_object_tester::test_acsc()
{
    return test_scalar_func<Acsc_func,false>();
};
void gmp_object_tester::test_acsc_c()
{
    return test_scalar_func<Acsc_c_func,false>();
};
void gmp_object_tester::test_asinh()
{
    return test_scalar_func<Asinh_func>();
};
void gmp_object_tester::test_acosh()
{
    return test_scalar_func<Acosh_func>();
};
void gmp_object_tester::test_acosh_c()
{
    return test_scalar_func<Acosh_c_func>();
};
void gmp_object_tester::test_atanh()
{
    return test_scalar_func<Atanh_func>();
};
void gmp_object_tester::test_atanh_c()
{
    return test_scalar_func<Atanh_c_func>();
};
void gmp_object_tester::test_acoth()
{
    return test_scalar_func<Acoth_func,false>();
};
void gmp_object_tester::test_acoth_c()
{
    return test_scalar_func<Acoth_c_func,false>();
};
void gmp_object_tester::test_asech()
{
    return test_scalar_func<Asech_func,false>();
};
void gmp_object_tester::test_asech_c()
{
    return test_scalar_func<Asech_c_func,false>();
};
void gmp_object_tester::test_acsch()
{
    return test_scalar_func<Acsch_func,false>();
};
void gmp_object_tester::test_inv()
{
    return test_scalar_func<Inv_func>();
};
void gmp_object_tester::test_invs()
{
    return test_scalar_func<Invs_func>();
};

void gmp_object_tester::test_floor()
{
    return test_scalar_func<Floor_func>();
}
void gmp_object_tester::test_combinatorics()
{
    test_scalar_func<Factorial_func>();
    test_scalar_func<Double_factorial_func>();
    test_scalar_func<Binomial_coefficient_func>();
}

void gmp_object_tester::test_ceil()
{
    return test_scalar_func<Ceil_func>();
}
void gmp_object_tester::test_round()
{
    return test_scalar_func<Round_func>();
}
void gmp_object_tester::test_trunc()
{
    return test_scalar_func<Trunc_func>();
}
void gmp_object_tester::test_ifloor()
{
    return test_scalar_func<IFloor_func>();
}

void gmp_object_tester::test_iceil()
{
    return test_scalar_func<ICeil_func>();
}
void gmp_object_tester::test_iround()
{
    return test_scalar_func<IRound_func>();
}
void gmp_object_tester::test_itrunc()
{
    return test_scalar_func<ITrunc_func>();
}
void gmp_object_tester::test_sign()
{
    return test_scalar_func<Sign_func>();
}
void gmp_object_tester::test_op_neg()
{
    return test_scalar_func<Op_neg_func>();
}
void gmp_object_tester::test_op_true()
{
    return test_scalar_func<Op_true_func>();
}

bool gmp_object_tester::test_cons_impl()
{
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();
    mdy::Type mi = MP_int::get_static_type();
    mdy::Type mf = MP_float::get_static_type();
    mdy::Type mc = MP_complex::get_static_type();
    mdy::Type mq = MP_rational::get_static_type();

    std::vector<mdy::Type> types;
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    types.push_back(mi);
    types.push_back(mf);
    types.push_back(mc);
    types.push_back(mq);

    using type2 = std::pair<mdy::Type, mdy::Type>;

    std::set<type2> missing_cons;
    missing_cons.insert(type2(ti,tr));
    missing_cons.insert(type2(ti,tf));
    missing_cons.insert(type2(ti,tc));
    missing_cons.insert(type2(ti,tcf));
    missing_cons.insert(type2(tf,tc));
    missing_cons.insert(type2(tf,tcf));
    missing_cons.insert(type2(tr,tc));
    missing_cons.insert(type2(tr,tcf));

    missing_cons.insert(type2(ti,mf));
    missing_cons.insert(type2(ti,mc));
    missing_cons.insert(type2(ti,mq));
    missing_cons.insert(type2(tr,mc));
    missing_cons.insert(type2(tf,mc));
    missing_cons.insert(type2(mi,tr));
    missing_cons.insert(type2(mi,tf));
    missing_cons.insert(type2(mi,tc));
    missing_cons.insert(type2(mi,tcf));
    missing_cons.insert(type2(mi,mf));
    missing_cons.insert(type2(mi,mq));
    missing_cons.insert(type2(mi,mc));
    missing_cons.insert(type2(mf,tc));
    missing_cons.insert(type2(mf,tcf));
    missing_cons.insert(type2(mf,mc));
    missing_cons.insert(type2(mq,tc));
    missing_cons.insert(type2(mq,tcf));
    missing_cons.insert(type2(mq,mc));

    std::set<type2> missing_assign;

    missing_assign.insert(type2(ti,tr));
    missing_assign.insert(type2(ti,tf));
    missing_assign.insert(type2(ti,tc));
    missing_assign.insert(type2(ti,tcf));
    missing_assign.insert(type2(tf,tc));
    missing_assign.insert(type2(tf,tcf));
    missing_assign.insert(type2(tr,tc));
    missing_assign.insert(type2(tr,tcf));

    missing_assign.insert(type2(ti,mi));
    missing_assign.insert(type2(ti,mf));
    missing_assign.insert(type2(ti,mc));
    missing_assign.insert(type2(ti,mq));
    missing_assign.insert(type2(tr,mi));
    missing_assign.insert(type2(tr,mf));
    missing_assign.insert(type2(tr,mc));
    missing_assign.insert(type2(tr,mq));
    missing_assign.insert(type2(tf,mi));
    missing_assign.insert(type2(tf,mf));
    missing_assign.insert(type2(tf,mc));
    missing_assign.insert(type2(tf,mq));
    missing_assign.insert(type2(tc,mi));
    missing_assign.insert(type2(tc,mf));
    missing_assign.insert(type2(tc,mc));
    missing_assign.insert(type2(tc,mq));
    missing_assign.insert(type2(tcf,mi));
    missing_assign.insert(type2(tcf,mf));
    missing_assign.insert(type2(tcf,mc));
    missing_assign.insert(type2(tcf,mq));

    missing_assign.insert(type2(mi,mf));
    missing_assign.insert(type2(mi,mq));
    missing_assign.insert(type2(mi,mc));
    missing_assign.insert(type2(mf,mc));
    missing_assign.insert(type2(mq,mf));
    missing_assign.insert(type2(mq,mc));

    missing_assign.insert(type2(mi,tf));
    missing_assign.insert(type2(mi,tr));
    missing_assign.insert(type2(mi,tc));
    missing_assign.insert(type2(mi,tcf));

    missing_assign.insert(type2(mf,tc));
    missing_assign.insert(type2(mf,tcf));
    missing_assign.insert(type2(mq,tf));
    missing_assign.insert(type2(mq,tr));
    missing_assign.insert(type2(mq,tc));
    missing_assign.insert(type2(mq,tcf));

    double res = 0;

    for (const auto& pos1 : types)
    for (const auto& pos2 : types)
    {
        res += test_cons_1(pos1, pos2, missing_cons, missing_assign);
    };

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_cons_1(mdy::Type t1, mdy::Type t2, const type2_set& missing_cons,
                                      const type2_set& missing_assign)
{
    mdy::object rhs(t2);
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        mdy::object x = cast(t1, rhs);
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        res += 1;
    }

    try
    {        
        mdy::object t(t1, rhs);

        if (missing_cons.find(type2(t1, t2)) != missing_cons.end())
            res += 1;
    }
    catch(std::exception& ex)
    {
        if (missing_cons.find(type2(t1, t2)) == missing_cons.end())
        {
            disp(ex.what());
            res += 1;
        };
    }

    try
    {        
        mdy::object x = convert(t1, rhs);

        if (missing_cons.find(type2(t1, t2)) != missing_cons.end())
            res += 1;
    }
    catch(std::exception& ex)
    {
        if (missing_cons.find(type2(t1, t2)) == missing_cons.end())
        {
            disp(ex.what());
            res += 1;
        };
    }

    try
    {
        lhs = rhs;

        if (missing_assign.find(type2(t1, t2)) != missing_assign.end())
            res += 1;
    }
    catch(std::exception& ex)
    {
        if (missing_assign.find(type2(t1, t2)) == missing_assign.end())
        {
            disp(ex.what());
            res += 1;
        };
    }

    return res;
};

bool gmp_object_tester::test_cons_val_impl()
{
    std::vector<Scalar_ext> scalars1;
    std::vector<Scalar_ext> scalars2;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);
    rand_scalars_ext::make(scalars2, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_cons_val_1(scalars1[i],scalars2[i], i);

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_cons_val_1(const Scalar_ext& s1, const Scalar_ext& s2, Integer code)
{
    eval_cons_val test1(code);
    double res1 = test1.make(s1, s2);

    eval_assign_val test2(code);
    double res2 = test2.make(s1, s2);

    eval_cast_val test3(code);
    double res3 = test3.make(s1, s2);

    return res1 + res2 + res3; 
};


bool gmp_object_tester::test_object_func_impl()
{
    mdy::Type mi = MP_int::get_static_type();
    mdy::Type mf = MP_float::get_static_type();
    mdy::Type mc = MP_complex::get_static_type();
    mdy::Type mq = MP_rational::get_static_type();

    std::vector<mdy::Type> types;
    types.push_back(mi);
    types.push_back(mf);
    types.push_back(mc);
    types.push_back(mq);

    double res = 0;
    for (const auto& pos1 : types)
    {
        res += test_func_1(pos1);
    };

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_func_1(mdy::Type t1)
{
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        (bool)lhs;
        !lhs;
        lhs.clone();
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        res += 1;
    }

    if (t1 != mdy::Type())
    {
        if (lhs.is_zero() == false)
            res += 1;
        if (lhs.is_null() == true)
            res += 1;
    };

    return res;
};

bool gmp_object_tester::test_object_func_val_impl()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_1_val(scalars1[i], i);

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_func_1_val(const Scalar_ext& s1, Integer code)
{
    eval_object_func_val test1(code);
    double res1 = test1.make(s1);

    return res1;
};

bool gmp_object_tester::test_func_compare_impl()
{
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();

    mdy::Type mi = MP_int::get_static_type();
    mdy::Type mf = MP_float::get_static_type();
    mdy::Type mc = MP_complex::get_static_type();
    mdy::Type mq = MP_rational::get_static_type();

    std::vector<mdy::Type> types;
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    types.push_back(mi);
    types.push_back(mf);
    types.push_back(mc);
    types.push_back(mq);

    using type2 = std::pair<mdy::Type, mdy::Type>;

    std::set<type2> missing;

    double res = 0;
    for (const auto& pos1 : types)
    for (const auto& pos2 : types)
    {
        res += test_compare_1(pos1, pos2, missing);
    };

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_compare_1(mdy::Type t1, mdy::Type t2, const type2_set& missing)
{
    mdy::object rhs(t2);
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        bool res1 = (bool)(lhs == rhs);
        bool res2 = (bool)(lhs != rhs);
        bool res3 = (bool)(lhs < rhs);
        bool res4 = (bool)(lhs > rhs);
        bool res5 = (bool)(lhs <= rhs);
        bool res6 = (bool)(lhs >= rhs);
    }
    catch(std::exception& ex)
    {
        if (missing.find(type2(t1, t2)) == missing.end())
        {
            disp(ex.what());
            res += 1;
        };
    }

    return res;
};

bool gmp_object_tester::test_func_uminus_impl()
{
    mdy::Type mi = MP_int::get_static_type();
    mdy::Type mf = MP_float::get_static_type();
    mdy::Type mc = MP_complex::get_static_type();
    mdy::Type mq = MP_rational::get_static_type();

    std::vector<mdy::Type> types;
    types.push_back(mi);
    types.push_back(mf);
    types.push_back(mc);
    types.push_back(mq);

    double res = 0;
    for (const auto& pos1 : types)
    {
        res += test_func_uminus_1(pos1);
    };

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_func_uminus_1(mdy::Type t1)
{
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        mdy::object r = -lhs;

        if (r.get_type() != lhs.get_type())
            res += 1;
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        res += 1;
    }

    return res;
};

bool gmp_object_tester::test_func_uminus_val_impl()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_uminus_val(scalars1[i], i);

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_func_uminus_val(const Scalar_ext& s1, Integer code)
{
    eval_uminus_val test1(code);
    double res1 = test1.make(s1);

    return res1;
};

bool gmp_object_tester::test_func_reim_impl()
{
    mdy::Type mi = MP_int::get_static_type();
    mdy::Type mf = MP_float::get_static_type();
    mdy::Type mc = MP_complex::get_static_type();
    mdy::Type mq = MP_rational::get_static_type();

    std::vector<mdy::Type> types;
    types.push_back(mi);
    types.push_back(mf);
    types.push_back(mc);
    types.push_back(mq);

    double res = 0;
    for (const auto& pos1 : types)
    {
        res += test_func_reim_1(pos1);
    };

    return (res == 0) ? true : false;
}

static bool is_complex_type(mdy::Type t)
{
    return t == MP_complex::get_static_type();
}

double gmp_object_tester::test_func_reim_1(mdy::Type t1)
{
    mdy::object lhs(t1);

    double res = 0;

    try
    {        
        mdy::object res1 = real(lhs);
        mdy::object res2 = imag(lhs);

        if (res1.get_type() != res2.get_type())
            res += 1;

        if (res1.get_type() != lhs.get_type() && is_complex_type(t1) == false)
            res += 1;
    }
    catch(std::exception& ex)
    {
        disp(ex.what());
        res += 1;
    }

    return res;
};

bool gmp_object_tester::test_func_reim_val_impl()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_reim_val(scalars1[i], i);

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_func_reim_val(const Scalar_ext& s1, Integer code)
{
    eval_reim_val test1(code);
    double res1 = test1.make(s1);

    return res1;
};

bool gmp_object_tester::test_func_unify_impl()
{
    mdy::Type ti = mdy::predefined::type_int();
    mdy::Type tr = mdy::predefined::type_real();
    mdy::Type tf = mdy::predefined::type_float();
    mdy::Type tc = mdy::predefined::type_complex();
    mdy::Type tcf = mdy::predefined::type_float_complex();

    mdy::Type mi = MP_int::get_static_type();
    mdy::Type mf = MP_float::get_static_type();
    mdy::Type mc = MP_complex::get_static_type();
    mdy::Type mq = MP_rational::get_static_type();

    std::vector<mdy::Type> types;
    types.push_back(ti);
    types.push_back(tr);
    types.push_back(tf);
    types.push_back(tc);
    types.push_back(tcf);
    types.push_back(mi);
    types.push_back(mf);
    types.push_back(mc);
    types.push_back(mq);

    double res = 0;
    Integer code = 0;
    for (const auto& pos1 : types)
    for (const auto& pos2 : types)
    {
        res += test_unify_1(pos1, pos2, code);
        ++code;
    };

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_unify_1(mdy::Type t1, mdy::Type t2, Integer code)
{
    //if (code == 45)
    //    out_stream << "break" << "\n";

    mdy::Type tt[]   = {t1, t2};    
    mdy::Type ta     = mdy::operations::return_type(mdy::functions::op_plus::eval(), 2, tt);
    mdy::Type td     = mdy::operations::return_type(mdy::functions::op_div::eval(), 2, tt);
    mdy::Type tr     = MP_rational::get_static_type();

    mdy::Type tu1    = mdy::operations::unify_types(t1, t2);
    mdy::Type tu2    = mdy::operations::unify_types(tu1, tr);

    if (tu1 == ta)
        return 0.0;
    if (tu2 == td)
        return 0.0;

    out_stream << code << "\n";
    out_stream << t1 << " " << t2 << "\n";
    out_stream << tu1 << " " << ta << "\n";
    out_stream << tu2 << " " << td << "\n";

    return 1.0;
}

bool gmp_object_tester::test_func_is_impl()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_is(scalars1[i], i);

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_func_is(const Scalar_ext& s1, Integer code)
{
    eval_is test1(code);
    double res1 = test1.make(s1);

    return res1;
};

bool gmp_object_tester::test_func_next_impl()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_next(scalars1[i], i);

    return (res == 0) ? true : false;
}

bool gmp_object_tester::test_func_sign_impl()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_sign(scalars1[i], i);

    return (res == 0) ? true : false;
}

bool gmp_object_tester::test_func_eps_impl()
{
    std::vector<Scalar_ext> scalars1;

    Integer N   = 1000;
    rand_scalars_ext::make(scalars1, N, true);

    double res = 0.0;

    for (Integer i = 0; i < 1000; ++i)
        res     += test_func_eps(scalars1[i], i);

    return (res == 0) ? true : false;
}

double gmp_object_tester::test_func_next(const Scalar_ext& s1, Integer code)
{
    eval_next test1(code);
    double res1 = test1.make(s1);

    return res1;
};

double gmp_object_tester::test_func_sign(const Scalar_ext& s1, Integer code)
{
    eval_signbit test1(code);
    double res1 = test1.make(s1);

    eval_isign test2(code);
    double res2 = test2.make(s1);

    return res1 + res2;
};

double gmp_object_tester::test_func_eps(const Scalar_ext& s1, Integer code)
{
    test_eval_eps test1(code);
    double res1 = test1.make(s1);

    return res1;
};

}};
