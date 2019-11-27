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

#include "matcl-matrep/matcl_matrep.h"
#include "test/test_matcl/framework/matrix_set/matrix_set.h"
#include "test/test_matcl/framework/matrix_utils.h"
#include "matcl-linalg/matcl_linalg.h"

namespace matcl { namespace test
{

template<class func_helper>
class scal_func : public unary_function
{
    public:
        virtual Real eval_mat(const Matrix& mat,bool,int code )
        {
            (void)code;
            try
            {		
                Matrix out_full = func_helper::eval(full(mat));
                value_code vc   = matrix_traits::complex_value_type(mat.get_value_code());
                matcl::mat_code nt = matrix_traits::get_matrix_type(vc,mat.get_struct_code());

                Matrix mat_c    = convert(mat,nt);
                Matrix out		= func_helper::eval(mat);
                Matrix out_c	= func_helper::eval(mat_c);
                check_struct(out);
                Real dif        = norm_1(out - out_full)/(norm_1(out) + constants::eps());
                dif             += norm_1(out - out_c)/(norm_1(out) + 1.);

                if (abs(dif) < 1000 * out.length() * constants::eps())
                    dif=0;

                return dif;
            }
            catch(error::scalar_required)
            {
                return mat.cols() == 1 && mat.rows() == 1 ? 1. : 0.;
            }
            catch(const std::exception& ex)
            {
                m_is_error = true;
                m_error = ex.what();
                return 1.;
            }
            catch(...)
            {
                m_is_error = true;
                m_error = "unknown error";
                return 1.;
            };
        };

        virtual Real eval_scalar(const Scalar& mat,bool,int code )
        {
            (void)code;

            Real dif = 0;

            switch (mat.get_value_code())
            {
                case value_code::v_integer:
                {
                    Integer val = mat.get_int();
                    Matrix res1 = func_helper::eval(val);
                    dif         = norm_1(func_helper::eval_scal(val) - res1);

                    if(abs(dif) < 1000.*constants::eps()) 
                        dif = 0.; // Absolute tolerance

                    dif = dif / (norm_1(res1) + constants::eps()*(1. + norm_1(val)));
                    break;
                }
                case value_code::v_float:
                {
                    Float val       = mat.get_float();
                    Matrix res1     = func_helper::eval(val);
                    Float_complex val_c	= val;
                    dif		        = norm_1(func_helper::eval_scal(val) - res1);
                    dif				+=norm_1(func_helper::eval_scal(val) - func_helper::eval_scal(val_c));

                    if(abs(dif) < 1000.*constants::f_eps()) 
                        dif = 0.; // Absolute tolerance

                    dif = dif / (norm_1(res1) + constants::f_eps()*(1. + norm_1(val)));
                    break;
                }
                case value_code::v_real:
                {
                    Real val		= mat.get_real();
                    Matrix res1     = func_helper::eval(val);
                    Complex val_c	= val;
                    dif		        = norm_1(func_helper::eval_scal(val) - res1);
                    dif				+=norm_1(func_helper::eval_scal(val) - func_helper::eval_scal(val_c));

                    if(abs(dif) < 1000.*constants::eps()) 
                        dif = 0.; // Absolute tolerance

                    dif = dif / (norm_1(res1) + constants::eps()*(1. + norm_1(val)));
                    break;
                }
                case value_code::v_float_complex:
                {
                    Float_complex val   = mat.get_fcomplex();
                    Matrix res1         = func_helper::eval(val);
                    dif                 = norm_1(func_helper::eval_scal(val) - res1);

                    if(abs(dif) < 1000.*constants::f_eps()) 
                        dif = 0.; // Absolute tolerance

                    dif = dif / (norm_1(res1) + constants::f_eps()*(1. + norm_1(val)));
                    break;
                }
                case value_code::v_complex:
                {
                    Complex val     = mat.get_complex();
                    Matrix res1     = func_helper::eval(val);
                    dif             = norm_1(func_helper::eval_scal(val) - res1);

                    if(abs(dif) < 1000.*constants::eps())
                        dif = 0.; // Absolute tolerance

                    dif = dif / (norm_1(res1) + constants::eps()*(1. + norm_1(val)));
                    break;
                }
                default:
                {
                    return 0;
                }
            };

            if (dif < 1000 * constants::eps())
                dif = 0;

            return dif;
        };
};

};};
