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

#include "test_simd_config.h"
#include "matcl-core/config.h"
#include "matcl-simd/simd_math.h"
#include "matcl-mp/matcl_mp.h"
#include "accuracy/test_accuracy.h"

#include "test_simd_accuracy.h"

namespace test_functions_accuracy
{

template<class Val>
struct Func_exp : matcl::test_accuracy_function<Val>
{
    bool m_rand_denormals;

    Func_exp()
    {
        m_rand_denormals    = matcl::test::simd_accuracy_tester::rand_denormals();

        Val min     = m_rand_denormals ? std::log(std::numeric_limits<Val>::denorm_min())
                        : std::log(std::numeric_limits<Val>::min());

        Val mean    = std::numeric_limits<Val>::epsilon() / Val(256.0);
        Val max     = std::log(std::numeric_limits<Val>::max());

        add_log_range(min, -mean);
        add_log_range(mean, max);        
    };

    virtual std::string name() const override
    { 
        return "exp"; 
    };

    virtual bool rand_denormals() const override
    {
        return m_rand_denormals;
    }
    
    matcl::mp_float eval_mp(const matcl::mp_float& v) const override
    {
        return exp(v);
    }

    Val eval_base(const Val& v) const override
    {
        return matcl::simd::exp(v);
    }

    Val eval_ref(const Val& v) const override
    {
        return std::exp(v);
    }
};

template<class Val>
struct Func_log : matcl::test_accuracy_function<Val>
{
    bool m_rand_denormals;

    Func_log()
    {
        m_rand_denormals    = matcl::test::simd_accuracy_tester::rand_denormals();

        Val min     = m_rand_denormals ? std::numeric_limits<Val>::denorm_min()
                        : std::numeric_limits<Val>::min();

        Val max     = std::numeric_limits<Val>::max();

        add_log_range(min, max);        
    };

    virtual std::string name() const override
    { 
        return "log"; 
    };

    virtual bool rand_denormals() const override
    {
        return m_rand_denormals;
    }
    
    matcl::mp_float eval_mp(const matcl::mp_float& v) const override
    {
        return log(v);
    }

    Val eval_base(const Val& v) const override
    {
        return matcl::simd::log(v);
    }

    Val eval_ref(const Val& v) const override
    {
        return std::log(v);
    }
};

template<class Val>
struct add_pi2_mult
{};

template<>
struct add_pi2_mult<double>
{
    static void eval(matcl::test_accuracy_function<double>& test)
    {
        // worst cases for argument reduction
        test.add_special_value(std::ldexp(6381956970095103.0, 796));

        int max_k       = (1 << 18);

        matcl::precision prec   = matcl::precision(60);
        matcl::mp_float pi2     = matcl::constants::mp_pi(prec);
        
        for (int i = 0; i < max_k; ++i)
        {
            matcl::mp_float xt  = pi2 * matcl::mp_float(i);
            
            double xd           = xt.cast_float();
            test.add_special_value(xd);
        }
    }
};

template<>
struct add_pi2_mult<float>
{
    static void eval(matcl::test_accuracy_function<float>& test)
    {
        // worst cases for argument reduction
        test.add_special_value(std::ldexp(16367173.0f, 72));

        int max_k       = (1 << 18);

        matcl::precision prec   = matcl::precision(60);
        matcl::mp_float pi2     = matcl::constants::mp_pi(prec);
        
        for (int i = 0; i < max_k; ++i)
        {
            matcl::mp_float xt  = pi2 * matcl::mp_float(i);
            
            float xf            = float(xt.cast_float());
            test.add_special_value(xf);
        }
    }
};

template<class Val>
struct Func_sin : matcl::test_accuracy_function<Val>
{
    bool m_rand_denormals;

    Func_sin()
    {
        m_rand_denormals    = matcl::test::simd_accuracy_tester::rand_denormals();

        Val max     = std::numeric_limits<Val>::max();
        Val max_1   = Val(823549.6);

        add_log_range(-max, max);
        add_linear_range(-max_1, max_1);

        add_pi2_mult<Val>::eval(*this);
    };

    virtual std::string name() const override
    { 
        return "sin"; 
    };

    virtual bool rand_denormals() const override
    {
        return m_rand_denormals;
    }
    
    matcl::mp_float eval_mp(const matcl::mp_float& v) const override
    {
        return sin(v);
    }

    Val eval_base(const Val& v) const override
    {
        return matcl::simd::sin(v);
    }

    Val eval_ref(const Val& v) const override
    {
        return std::sin(v);
    }
};

template<class Val>
struct Func_cos : matcl::test_accuracy_function<Val>
{
    bool m_rand_denormals;

    Func_cos()
    {
        m_rand_denormals    = matcl::test::simd_accuracy_tester::rand_denormals();

        Val max     = std::numeric_limits<Val>::max();
        Val max_1   = Val(823549.6);

        add_log_range(-max, max);
        add_linear_range(-max_1, max_1);

        add_pi2_mult<Val>::eval(*this);
    };

    virtual std::string name() const override
    { 
        return "cos"; 
    };

    virtual bool rand_denormals() const override
    {
        return m_rand_denormals;
    }
    
    matcl::mp_float eval_mp(const matcl::mp_float& v) const override
    {
        return cos(v);
    }

    Val eval_base(const Val& v) const override
    {
        return matcl::simd::cos(v);
    }

    Val eval_ref(const Val& v) const override
    {
        return std::cos(v);
    }
};

template<class Val>
struct Func_tan : matcl::test_accuracy_function<Val>
{
    bool m_rand_denormals;

    Func_tan()
    {
        m_rand_denormals    = matcl::test::simd_accuracy_tester::rand_denormals();

        Val max     = std::numeric_limits<Val>::max();
        Val max_1   = Val(823549.6);

        add_log_range(-max, max);
        add_linear_range(-max_1, max_1);

        add_pi2_mult<Val>::eval(*this);
    };

    virtual std::string name() const override
    { 
        return "tan"; 
    };

    virtual bool rand_denormals() const override
    {
        return m_rand_denormals;
    }
    
    matcl::mp_float eval_mp(const matcl::mp_float& v) const override
    {
        return tan(v);
    }

    Val eval_base(const Val& v) const override
    {
        return matcl::simd::tan(v);
    }

    Val eval_ref(const Val& v) const override
    {
        return std::tan(v);
    }
};

template<class Val>
struct Func_cot : matcl::test_accuracy_function<Val>
{
    bool m_rand_denormals;

    Func_cot()
    {
        m_rand_denormals    = matcl::test::simd_accuracy_tester::rand_denormals();

        Val max     = std::numeric_limits<Val>::max();
        Val max_1   = Val(823549.6);

        add_log_range(-max, max);
        add_linear_range(-max_1, max_1);

        add_pi2_mult<Val>::eval(*this);
    };

    virtual std::string name() const override
    { 
        return "cot"; 
    };

    virtual bool rand_denormals() const override
    {
        return m_rand_denormals;
    }
    
    matcl::mp_float eval_mp(const matcl::mp_float& v) const override
    {
        return cot(v);
    }

    Val eval_base(const Val& v) const override
    {
        return matcl::simd::cot(v);
    }

    Val eval_ref(const Val& v) const override
    {
        return Val(1) / std::tan(v);
    }
};

}