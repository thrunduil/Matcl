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

#include "test_accuracy.h"
#include "matcl-scalar/lib_functions/scalar_gen.h"

namespace matcl { namespace details
{

template<class Val>
struct rand_value{};

template<>
struct rand_value<double>
{
    static double eval()
    {
        return matcl::rand();
    }
};

template<>
struct rand_value<float>
{
    static float eval()
    {
        return matcl::frand();
    }
};


}}

namespace matcl
{

template<class Value>
void test_accuracy_function<Value>::add_log_range(Value min, Value max)
{
    if (matcl::is_nan(min) == true || matcl::is_nan(max) == true)
        return;

    if (min > max)
        return add_log_range(max, min);

    min     = std::max(min, -std::numeric_limits<Value>::max());
    max     = std::min(max, std::numeric_limits<Value>::max());

    if (min == max)
    {
        add_special_value(min);
        return;
    }

    Value prod  = min * max;

    if (prod < Value())
    {
        Value min1  = min;
        Value max1  = -min_value();

        Value min2  = min_value();
        Value max2  = max;

        add_log_range(min1, max1);
        add_log_range(min2, max2);
        return;
    }    

    if (prod == 0)
    {
        if (min == Value(0))
        {
            min     = min_value();
            return add_log_range(min, max);
        }

        if (max == Value(0))
        {
            max     = -min_value();
            return add_log_range(min, max);
        }
    }

    add_log_range_split(min, max);
}

template<class Value>
void test_accuracy_function<Value>::add_log_range_split(Value min, Value max)
{
    bool sign   = false;

    if (min < 0)
    {
        sign        = true;
        std::swap(min, max);
        min         = -min;
        max         = -max;
    }

    int n_steps     = 16;
    Value min_l     = std::log(min);
    Value max_l     = std::log(max);
    Value step      = (max_l - min_l) / n_steps;

    Value r1        = min_l;

    for (int i = 0; i < n_steps; ++i)
    {
        Value r2    = r1 + step;
        Value min_s = std::exp(r1);
        Value max_s = std::exp(r2);

        if (i == 0)
            min_s   = min;
        if (i == n_steps - 1)
            max_s   = max;

        if (sign == false)
        {
            m_log_range_vec.push_back(range_type(min_s, max_s));
        }
        else
        {
            m_log_range_vec.push_back(range_type(-max_s, -min_s));
        }

        r1          = r2;
    }
}

template<class Value>
void test_accuracy_function<Value>::add_linear_range(const Value& min, const Value& max)
{
    if (matcl::is_nan(min) == true || matcl::is_nan(max) == true)
        return;

    if (min > max)
        return add_linear_range(max, min);

    min     = std::max(min, -std::numeric_limits<Value>::max());
    max     = std::min(max, std::numeric_limits<Value>::max());

    if (min == max)
    {
        add_special_value(min);
        return;
    }

    m_lin_range_vec.push_back(range_type(min, max));
}

template<class Value>
void test_accuracy_function<Value>::add_special_value(const Value& val)
{
    m_value_vec.push_back(val);
}

template<class Value>
int test_accuracy_function<Value>::add_special_values(std::vector<Value>& vec) const
{
    int N   = (int)vec.size();

    vec.push_back(Value(0));
    vec.push_back(-Value(0));
    vec.push_back(Value(1));
    vec.push_back(Value(-1));
    vec.push_back(Value(2));
    vec.push_back(Value(-2));
    vec.push_back(Value(1.0/2.0));
    vec.push_back(-Value(1.0/2.0));

    vec.push_back(std::numeric_limits<Value>::infinity());
    vec.push_back(-std::numeric_limits<Value>::infinity());
    vec.push_back(std::numeric_limits<Value>::quiet_NaN());
    vec.push_back(std::numeric_limits<Value>::max());
    vec.push_back(-std::numeric_limits<Value>::max());
    vec.push_back(std::numeric_limits<Value>::min());
    vec.push_back(-std::numeric_limits<Value>::min());
    vec.push_back(min_value());
    vec.push_back(-min_value());

    vec.insert(vec.end(), m_value_vec.begin(), m_value_vec.end());

    int n_elems     = (int)vec.size() - N;
    return n_elems;
}

template<class Value>
void test_accuracy_function<Value>::rand_values(std::vector<Value>& vec, int size) const
{
    vec.reserve(vec.size() + size);

    int N           = add_special_values(vec);
    
    int n_log_range = (int)m_log_range_vec.size();
    int n_lin_range = (int)m_lin_range_vec.size();
    int n_range     = n_log_range + n_lin_range;

    if (n_range == 0)
    {
        Value rmin  = -std::numeric_limits<Value>::max();
        Value rmax  = std::numeric_limits<Value>::max();

        for (int i = N; i < size; ++i)
        {
            Value v     = rand_value(rmin, rmax, true);
            vec.push_back(v);
        }
    }
    else
    {
        for (int i = N; i < size; ++i)
        {
            Value v;
            int r_ind   = std::abs(matcl::irand()) % n_range;

            if (r_ind < n_log_range)
            {
                Value rmin  = m_log_range_vec[r_ind].first;
                Value rmax  = m_log_range_vec[r_ind].second;
                v           = rand_value(rmin, rmax, true);
            }
            else
            {
                Value rmin  = m_lin_range_vec[r_ind - n_log_range].first;
                Value rmax  = m_lin_range_vec[r_ind - n_log_range].second;
                v           = rand_value(rmin, rmax, true);
            }

            vec.push_back(v);
        }
    };
}

template<class Value>
Value test_accuracy_function<Value>::rand_value(const Value& min, const Value& max, bool log_range) const
{
    if (log_range == false)
    {
        Value v = min + details::rand_value<Value>::eval() * (max - min);
        return v;
    }

    if (min > Value(0))
    {
        Value l_min     = std::log(min);
        Value l_max     = std::log(max);

        Value v     = l_min + details::rand_value<Value>::eval() * (l_max - l_min);
        v           = std::exp(v);
        return v;
    }
    else
    {
        Value l_min     = std::log(-max);
        Value l_max     = std::log(-min);

        Value v     = l_min + details::rand_value<Value>::eval() * (l_max - l_min);
        v           = std::exp(v);
        return -v;
    }
}

template<class Value>
Value test_accuracy_function<Value>::min_value() const
{
    if (rand_denormals() == true)
        return std::numeric_limits<Value>::denorm_min();
    else
        return std::numeric_limits<Value>::min();
}

template<class Value>
bool test_accuracy_function<Value>::rand_denormals() const
{
    return true;
}

}