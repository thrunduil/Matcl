/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2018
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

class bin_function
{
    protected:
        bool			m_is_error;
        std::string		m_error;

    private:
        virtual Real	eval_mat(const Matrix& mat1,const Matrix& mat2, int code) = 0;
        virtual Real	eval_scalar(const Scalar& s1, const Scalar& s2, int code) = 0;

    public:
        bin_function()								: m_is_error(false) {};

        virtual			~bin_function()				{}; 

        Real			eval(const Matrix& mat1,const Matrix& mat2, int code);
        Real			eval_scal(const Scalar& s1, const Scalar& s2, int code);

        bool			is_error() const			{ return m_is_error;};
        std::string		get_error() const			{ return m_error; };

        template<class Derived>
        Real	        eval_scalar_impl(const Scalar& s1, const Scalar& s2);

    protected:
        Real			check_value_type(const Matrix& mat1,const Matrix& mat2,const Matrix& out, 
                                        bool float_return);
        Real			check_value_type(const Matrix& mat1,const Matrix& mat2);
};

template<class Derived>
Real bin_function::eval_scalar_impl(const Scalar& s1, const Scalar& s2)
{
    switch((int)s1.get_value_code()*10+(int)s2.get_value_code())
    {
        case (int)value_code::v_integer*10+(int)value_code::v_integer:
        {
            Integer value1 = s1.get_int();
            Integer value2 = s2.get_int();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_integer*10+(int)value_code::v_float:
        {
            Integer value1	= s1.get_int();
            Float value2	= s2.get_float();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_integer*10+(int)value_code::v_real:
        {
            Integer value1	= s1.get_int();
            Real value2		= s2.get_real();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_integer*10+(int)value_code::v_float_complex:
        {
            Integer value1	        = s1.get_int();
            Float_complex value2	= s2.get_fcomplex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_integer*10+(int)value_code::v_complex:
        {
            Integer value1	= s1.get_int();
            Complex value2	= s2.get_complex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float*10+(int)value_code::v_integer:
        {
            Float value1	= s1.get_float();
            Integer value2	= s2.get_int();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_real*10+(int)value_code::v_integer:
        {
            Real value1		= s1.get_real();
            Integer value2	= s2.get_int();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float*10+(int)value_code::v_float:
        {
            Float value1	= s1.get_float();
            Float value2	= s2.get_float();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_real*10+(int)value_code::v_float:
        {
            Real value1	    = s1.get_real();
            Float value2	= s2.get_float();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float*10+(int)value_code::v_real:
        {
            Float value1	= s1.get_float();
            Real value2	    = s2.get_real();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_real*10+(int)value_code::v_real:
        {
            Real value1		= s1.get_real();
            Real value2		= s2.get_real();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float*10+(int)value_code::v_float_complex:
        {
            Float value1		    = s1.get_float();
            Float_complex value2	= s2.get_fcomplex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_real*10+(int)value_code::v_float_complex:
        {
            Real value1		        = s1.get_real();
            Float_complex value2	= s2.get_fcomplex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float*10+(int)value_code::v_complex:
        {
            Float value1		    = s1.get_float();
            Complex value2	        = s2.get_complex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_real*10+(int)value_code::v_complex:
        {
            Real value1		= s1.get_real();
            Complex value2	= s2.get_complex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float_complex*10+(int)value_code::v_integer:
        {
            Float_complex value1	= s1.get_fcomplex();
            Integer value2	        = s2.get_int();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_complex*10+(int)value_code::v_integer:
        {
            Complex value1	= s1.get_complex();
            Integer value2	= s2.get_int();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float_complex*10+(int)value_code::v_float:
        {
            Float_complex value1	= s1.get_fcomplex();
            Float value2		    = s2.get_float();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_complex*10+(int)value_code::v_float:
        {
            Complex value1	        = s1.get_complex();
            Float value2		    = s2.get_float();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float_complex*10+(int)value_code::v_real:
        {
            Float_complex value1	= s1.get_fcomplex();
            Real value2		        = s2.get_real();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_complex*10+(int)value_code::v_real:
        {
            Complex value1	= s1.get_complex();
            Real value2		= s2.get_real();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float_complex*10+(int)value_code::v_float_complex:
        {
            Float_complex value1	= s1.get_fcomplex();
            Float_complex value2	= s2.get_fcomplex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_complex*10+(int)value_code::v_float_complex:
        {
            Complex value1	        = s1.get_complex();
            Float_complex value2	= s2.get_fcomplex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_float_complex*10+(int)value_code::v_complex:
        {
            Float_complex value1	= s1.get_fcomplex();
            Complex value2	        = s2.get_complex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        case (int)value_code::v_complex*10+(int)value_code::v_complex:
        {
            Complex value1	= s1.get_complex();
            Complex value2	= s2.get_complex();
            return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
        }
        default:
        {
            return 0;
        }
    };
};

};};