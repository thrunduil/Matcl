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

}