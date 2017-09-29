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
#include "scalar.h"
#include "scalar_ext.h"
#include <vector>

namespace matcl { namespace test
{

class rand_scalars
{
    public:
        static void make(std::vector<Scalar>& scal, Integer n, bool fractions);

    public:
        static Scalar          rand_scalar(bool fractions);
        static Integer         rand_int();
        static Float           rand_float(bool fractions);
        static Real            rand_real(bool fractions);
        static Complex         rand_compl(bool fractions);
        static Float_complex   rand_fcompl(bool fractions);
};

class rand_scalars_ext
{
    public:
        static void make(std::vector<Scalar_ext>& scal, Integer n, bool fractions);

    public:
        static Scalar_ext      rand_scalar(bool fractions);
        static mp_int          rand_mpint(bool fractions);
        static mp_float        rand_mpfloat(bool fractions);
        static mp_complex      rand_mpcompl(bool fractions);
        static mp_rational     rand_mprat(bool fractions);
};

template<class Derived>
struct eval_scalars_1
{
    virtual ~eval_scalars_1(){};

    double make(const Scalar& s1)
    {
        switch((int)s1.get_value_code())
        {
            case (int)value_code::v_integer:
            {
                Integer value1 = s1.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_float:
            {
                Float value1    = s1.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_real:
            {
                Real value1    = s1.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_complex:
            {
                Complex value1  = s1.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_float_complex:
            {
                Float_complex value1  = s1.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            default:
            {
                return 1;
            }
        };
    };

    double make(const Scalar_ext& s1)
    {
        using value_code = value_ext_code;

        switch((int)s1.get_value_code())
        {
            case (int)value_code::v_int:
            {
                Integer value1 = s1.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_float:
            {
                Float value1    = s1.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_real:
            {
                Real value1    = s1.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_complex:
            {
                Complex value1  = s1.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_fcomplex:
            {
                Float_complex value1  = s1.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_mpint:
            {
                const mp_int& value1  = s1.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_mpfloat:
            {
                const mp_float& value1  = s1.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_mpcompl:
            {
                const mp_complex& value1  = s1.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            case (int)value_code::v_mprat:
            {
                const mp_rational& value1  = s1.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1);
            }
            default:
            {
                return 1;
            }
        };
    };
};

template<class Derived>
struct eval_scalars
{
    virtual ~eval_scalars(){};

    double make(const Scalar& s1, const Scalar& s2)
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

    double make(const Scalar_ext& s1, const Scalar_ext& s2)
    {
        using value_code = value_ext_code;

        switch((int)s1.get_value_code()*10+(int)s2.get_value_code())
        {
            case (int)value_code::v_int*10+(int)value_code::v_int:
            {
                Integer value1 = s1.get_int();
                Integer value2 = s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_float:
            {
                Integer value1	= s1.get_int();
                Float value2	= s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_real:
            {
                Integer value1	= s1.get_int();
                Real value2		= s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_fcomplex:
            {
                Integer value1	        = s1.get_int();
                Float_complex value2	= s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_complex:
            {
                Integer value1	= s1.get_int();
                Complex value2	= s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_mpint:
            {
                Integer value1	        = s1.get_int();
                const mp_int& value2	= s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_mpfloat:
            {
                Integer value1	        = s1.get_int();
                const mp_float& value2	= s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_mpcompl:
            {
                Integer value1	            = s1.get_int();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_int*10+(int)value_code::v_mprat:
            {
                Integer value1	            = s1.get_int();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            case (int)value_code::v_float*10+(int)value_code::v_int:
            {
                Float value1	= s1.get_float();
                Integer value2	= s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_float:
            {
                Float value1	= s1.get_float();
                Float value2	= s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_real:
            {
                Float value1	= s1.get_float();
                Real value2	    = s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_fcomplex:
            {
                Float value1		    = s1.get_float();
                Float_complex value2	= s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_complex:
            {
                Float value1		    = s1.get_float();
                Complex value2	        = s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_mpint:
            {
                Float value1		    = s1.get_float();
                const mp_int& value2	= s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_mpfloat:
            {
                Float value1		    = s1.get_float();
                const mp_float& value2	= s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_mpcompl:
            {
                Float value1		        = s1.get_float();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_float*10+(int)value_code::v_mprat:
            {
                Float value1		        = s1.get_float();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            case (int)value_code::v_real*10+(int)value_code::v_int:
            {
                Real value1		= s1.get_real();
                Integer value2	= s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_float:
            {
                Real value1	    = s1.get_real();
                Float value2	= s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_real:
            {
                Real value1		= s1.get_real();
                Real value2		= s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_fcomplex:
            {
                Real value1		        = s1.get_real();
                Float_complex value2	= s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_complex:
            {
                Real value1		= s1.get_real();
                Complex value2	= s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_mpint:
            {
                Real value1		        = s1.get_real();
                const mp_int& value2	= s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_mpfloat:
            {
                Real value1		        = s1.get_real();
                const mp_float& value2	= s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_mpcompl:
            {
                Real value1		            = s1.get_real();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_real*10+(int)value_code::v_mprat:
            {
                Real value1		            = s1.get_real();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            
            case (int)value_code::v_fcomplex*10+(int)value_code::v_int:
            {
                Float_complex value1	= s1.get_fcomplex();
                Integer value2	        = s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_float:
            {
                Float_complex value1	= s1.get_fcomplex();
                Float value2		    = s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_real:
            {
                Float_complex value1	= s1.get_fcomplex();
                Real value2		        = s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_fcomplex:
            {
                Float_complex value1	= s1.get_fcomplex();
                Float_complex value2	= s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_complex:
            {
                Float_complex value1	= s1.get_fcomplex();
                Complex value2	        = s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_mpint:
            {
                Float_complex value1	= s1.get_fcomplex();
                const mp_int& value2	= s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_mpfloat:
            {
                Float_complex value1	= s1.get_fcomplex();
                const mp_float& value2	= s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_mpcompl:
            {
                Float_complex value1	    = s1.get_fcomplex();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_fcomplex*10+(int)value_code::v_mprat:
            {
                Float_complex value1	    = s1.get_fcomplex();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            case (int)value_code::v_complex*10+(int)value_code::v_int:
            {
                Complex value1	= s1.get_complex();
                Integer value2	= s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_float:
            {
                Complex value1	        = s1.get_complex();
                Float value2		    = s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_real:
            {
                Complex value1	= s1.get_complex();
                Real value2		= s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_fcomplex:
            {
                Complex value1	        = s1.get_complex();
                Float_complex value2	= s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_complex:
            {
                Complex value1	= s1.get_complex();
                Complex value2	= s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_mpint:
            {
                Complex value1	        = s1.get_complex();
                const mp_int& value2	= s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_mpfloat:
            {
                Complex value1	        = s1.get_complex();
                const mp_float& value2	= s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_mpcompl:
            {
                Complex value1	            = s1.get_complex();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_complex*10+(int)value_code::v_mprat:
            {
                Complex value1	            = s1.get_complex();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            //mpint
            case (int)value_code::v_mpint*10+(int)value_code::v_int:
            {
                const mp_int& value1	= s1.get_mpint();
                Integer value2	        = s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_float:
            {
                const mp_int& value1	= s1.get_mpint();
                Float value2		    = s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_real:
            {
                const mp_int& value1	= s1.get_mpint();
                Real value2		        = s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_fcomplex:
            {
                const mp_int& value1	= s1.get_mpint();
                Float_complex value2	= s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_complex:
            {
                const mp_int& value1	= s1.get_mpint();
                Complex value2	        = s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_mpint:
            {
                const mp_int& value1	= s1.get_mpint();
                const mp_int& value2	= s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_mpfloat:
            {
                const mp_int& value1	= s1.get_mpint();
                const mp_float& value2	= s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_mpcompl:
            {
                const mp_int& value1	= s1.get_mpint();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpint*10+(int)value_code::v_mprat:
            {
                const mp_int& value1	    = s1.get_mpint();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            //mpfloat
            case (int)value_code::v_mpfloat*10+(int)value_code::v_int:
            {
                const mp_float& value1	= s1.get_mpfloat();
                Integer value2	        = s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_float:
            {
                const mp_float& value1	= s1.get_mpfloat();
                Float value2		    = s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_real:
            {
                const mp_float& value1	= s1.get_mpfloat();
                Real value2		        = s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_fcomplex:
            {
                const mp_float& value1	= s1.get_mpfloat();
                Float_complex value2	= s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_complex:
            {
                const mp_float& value1	= s1.get_mpfloat();
                Complex value2	        = s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_mpint:
            {
                const mp_float& value1	= s1.get_mpfloat();
                const mp_int& value2	= s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_mpfloat:
            {
                const mp_float& value1	= s1.get_mpfloat();
                const mp_float& value2	= s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_mpcompl:
            {
                const mp_float& value1	    = s1.get_mpfloat();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpfloat*10+(int)value_code::v_mprat:
            {
                const mp_float& value1	    = s1.get_mpfloat();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            //mpcompl
            case (int)value_code::v_mpcompl*10+(int)value_code::v_int:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                Integer value2	            = s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_float:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                Float value2		        = s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_real:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                Real value2		            = s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_fcomplex:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                Float_complex value2	    = s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_complex:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                Complex value2	            = s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_mpint:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                const mp_int& value2	    = s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_mpfloat:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                const mp_float& value2	    = s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_mpcompl:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mpcompl*10+(int)value_code::v_mprat:
            {
                const mp_complex& value1    = s1.get_mpcomplex();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            //mprat
            case (int)value_code::v_mprat*10+(int)value_code::v_int:
            {
                const mp_rational& value1   = s1.get_mprational();
                Integer value2	            = s2.get_int();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_float:
            {
                const mp_rational& value1   = s1.get_mprational();
                Float value2		        = s2.get_float();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_real:
            {
                const mp_rational& value1   = s1.get_mprational();
                Real value2		            = s2.get_real();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_fcomplex:
            {
                const mp_rational& value1   = s1.get_mprational();
                Float_complex value2	    = s2.get_fcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_complex:
            {
                const mp_rational& value1   = s1.get_mprational();
                Complex value2	            = s2.get_complex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_mpint:
            {
                const mp_rational& value1   = s1.get_mprational();
                const mp_int& value2	    = s2.get_mpint();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_mpfloat:
            {
                const mp_rational& value1   = s1.get_mprational();
                const mp_float& value2	    = s2.get_mpfloat();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_mpcompl:
            {
                const mp_rational& value1   = s1.get_mprational();
                const mp_complex& value2	= s2.get_mpcomplex();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }
            case (int)value_code::v_mprat*10+(int)value_code::v_mprat:
            {
                const mp_rational& value1   = s1.get_mprational();
                const mp_rational& value2	= s2.get_mprational();
                return dynamic_cast<Derived*>(this)->eval_scal_func(value1,value2);
            }

            default:
            {
                return 1;
            }
        };
    };
};

}};
