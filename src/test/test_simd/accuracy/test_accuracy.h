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

#include "matcl-core/config.h"
#include <vector>

namespace matcl
{

template<class Value>
class test_accuracy_function
{
    private:
        using range_type    = std::pair<Value, Value>;
        using range_vec     = std::vector<range_type>;
        using value_vec     = std::vector<Value>;

    private:
        range_vec           m_log_range_vec;
        range_vec           m_lin_range_vec;
        value_vec           m_value_vec;

    public:
        virtual ~test_accuracy_function(){};

        void                add_log_range(Value min, Value max);
        void                add_linear_range(const Value& min, const Value& max);
        void                add_special_value(const Value& val);

        void                rand_values(std::vector<Value>& vec, int size) const;

        virtual bool        rand_denormals() const;

    public:
        virtual std::string name() const = 0;   

        virtual mp_float    eval_mp(const matcl::mp_float& v) const = 0;

        virtual Value       eval_base(const Value& v) const = 0;

        virtual Value       eval_ref(const Value& v) const = 0;

    private:
        void                add_log_range_split(Value min, Value max);
        int                 add_special_values(std::vector<Value>& vec) const;
        Value               rand_value(const Value& min, const Value& max, bool log_range) const;
        Value               min_value() const;
};

}

#include "test_accuracy.inl"