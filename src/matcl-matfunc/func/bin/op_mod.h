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

#include "matcl-matfunc/func/bin/op_info.h"
#include "matcl-internals/func/raw_func_op_helpers.h"

#pragma warning(push)
#pragma warning(disable: 4127)  // conditional expression is constant

namespace matcl { namespace details
{

namespace mdyf = matcl::dynamic::functions;

template<class Derived, bool Is_mod>
class op_mod_base : public op_info<Derived>
{
    public:
        template<class Val>
        bool is_special_case(const Val& val, int ) const
        { 
            return mrd::is_zero(val); 
        };

        template<class Val>
        void eval_special_case(matcl::Matrix& ret, const Val& v,const Matrix& mat) const
        { 
            ret = zero<Val>(mat,v); 
        };

        template<class Val>
        void eval_special_case(matcl::Matrix& ret ,const Matrix& mat, const Val& v) const
        { 
            ret = zero<Val>(mat,v); 
        };

        template<class Val>
        Matrix zero(const Matrix& mat, const Val& v) const
        {
            Integer m   = mat.rows();
            Integer n   = mat.cols();

            matcl::value_code vt = mat.get_value_code();
        
            if (vt == value_code::v_integer && std::is_same<Val, Integer>::value)
            {
                return ispzeros(m,n);
            }
            else if ((std::is_same<Val,Float>::value || std::is_same<Val,Float_complex>::value)
                     && (vt == value_code::v_float || vt == value_code::v_float_complex) )
            {
                return fspzeros(m,n);
            }
            else if (std::is_same<Val,Object>::value == false && vt != value_code::v_object)
            {
                return spzeros(m,n);
            }

            using ti_type_1 = matcl::ti::ti_object;
            using ti_type_2 = typename ti::get_ti_type<Val>::type;

            ti_type_1 t1    = mat.get_type();
            ti_type_2 t2    = matcl::ti::get_ti(v);

            ti_type_1 ti_r;
            
            if (Is_mod == true)
                ti_r        = matcl::ti::get_return_ti<ti_type_1>(mdyf::mod::eval(),t1,t2);
            else
                ti_r        = matcl::ti::get_return_ti<ti_type_1>(mdyf::rem::eval(),t1,t2);

            return mrd::create_matrix<Object>(ti_r, m, n);
        };

        template<class UT, class TI1, class TI2>
        void eval_special_case_00_impl(matcl::Matrix& ret, Integer m, Integer n, TI1 t1, TI2 t2) const
        {
            using UTR       = md::real_type<UT>::type;
            using ti_type   = matcl::ti::ti_type<UTR>;

            ti_type ti_r;

            if (Is_mod == true)
                ti_r        = matcl::ti::get_return_ti<ti_type>(mdyf::mod::eval(),t1,t2);
            else
                ti_r        = matcl::ti::get_return_ti<ti_type>(mdyf::rem::eval(),t1,t2);

            ret = mrd::create_matrix<UTR>(ti_r, m,n);
        };
};

}}

#pragma warning(pop)