/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018 - 2021
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

#include <string>
#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/data_struct/scalar.h"

namespace matcl { namespace test
{

class unary_function
{
    protected:
        bool			m_is_error;
        std::string		m_error;

    private:
        virtual Real	eval_mat(const Matrix& mat,bool show_res, int code) = 0;
        virtual Real	eval_scalar(const Scalar& s,bool show_res, int code) = 0;

    public:
        unary_function()							: m_is_error(false) {};

        virtual			~unary_function()			{}; 

        virtual long    n_new_objects()             { return 0; };
        Real			eval(const Matrix& mat,bool show_res, int code);
        Real			eval_scal(const Scalar& s,bool show_res, int code);

        bool			is_error() const			{ return m_is_error;};
        std::string		get_error() const			{ return m_error; };

        template<class Derived>
        Real	        eval_scalar_impl(const Scalar& s1) const;
};

template<class Derived>
Real unary_function::eval_scalar_impl(const Scalar& s) const
{
    switch(s.get_value_code())
    {
        case value_code::v_integer:
        {
            Integer value = s.get_int();
            return dynamic_cast<const Derived*>(this)->eval_scal_func(value);
        }
        case value_code::v_float:
        {
            Float value		= s.get_float();
            return dynamic_cast<const Derived*>(this)->eval_scal_func(value);
        }
        case value_code::v_float_complex:
        {
            Float_complex value = s.get_fcomplex();
            return dynamic_cast<const Derived*>(this)->eval_scal_func(value);
        }
        case value_code::v_real:
        {
            Real value		= s.get_real();
            return dynamic_cast<const Derived*>(this)->eval_scal_func(value);
        }
        case value_code::v_complex:
        {
            Complex value = s.get_complex();
            return dynamic_cast<const Derived*>(this)->eval_scal_func(value);
        }
        default:
        {
            return 0;
        }
    };
};

};};