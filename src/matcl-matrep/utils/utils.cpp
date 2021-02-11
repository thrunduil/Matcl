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

#include "matcl-internals/base/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/details/matrix.inl"

namespace matcl { namespace details
{

template<class val_type>
typename unit_matrix<val_type,struct_dense>::matrix_type 
    unit_matrix<val_type,struct_dense>::eval(tinfo_t ret_ti, Integer m)
{
    (void)ret_ti;
    Matrix tmp = matcl::eye(m,matrix_traits::value_code<val_type>::value);
    return tmp.move_impl<matrix_type>().move();
};

unit_matrix<Object,struct_dense>::matrix_type 
    unit_matrix<Object,struct_dense>::eval(tinfo_t ret_ti, Integer m)
{
    Matrix tmp = matcl::eye(ti::ti_object(ret_ti), m);
    return tmp.move_impl<matrix_type>().move();
};

template<class val_type>
typename unit_matrix<val_type,struct_sparse>::matrix_type 
unit_matrix<val_type,struct_sparse>::eval(tinfo_t ret_ti, Integer m)
{
    (void)ret_ti;
    Matrix tmp = matcl::speye(m,matrix_traits::value_code<val_type>::value);
    return tmp.move_impl<matrix_type>().move();
};

unit_matrix<Object,struct_sparse>::matrix_type 
unit_matrix<Object,struct_sparse>::eval(tinfo_t ret_ti, Integer m)
{
    Matrix tmp = matcl::speye(ti::ti_object(ret_ti), m);
    return tmp.move_impl<matrix_type>().move();
};

template<class val_type>
typename unit_matrix<val_type,struct_banded>::matrix_type         
unit_matrix<val_type,struct_banded>::eval(tinfo_t ret_ti, Integer m)
{
    (void)ret_ti;
    Matrix tmp = matcl::beye(m,m,0,0,matrix_traits::value_code<val_type>::value);
    return tmp.move_impl<matrix_type>().move();
};

unit_matrix<Object,struct_banded>::matrix_type 
unit_matrix<Object,struct_banded>::eval(tinfo_t ret_ti, Integer m)
{
    Matrix tmp = matcl::beye(ti::ti_object(ret_ti), m,m,0,0);
    return tmp.move_impl<matrix_type>().move();
};

//----------------------------------------------------------------
//                  trans_manip
//----------------------------------------------------------------

trans_type_ext trans_manip::link_trans(trans_type_ext t1, trans_type_ext t2)
{
    switch (t1)
    {
        case trans_type_ext::no_trans:
        {
            return t2;
        };
        case trans_type_ext::conj:
        {
            switch (t2)
            {
                case trans_type_ext::no_trans:      return trans_type_ext::conj;
                case trans_type_ext::trans:         return trans_type_ext::conj_trans;
                case trans_type_ext::conj_trans:    return trans_type_ext::trans;
                case trans_type_ext::conj:          return trans_type_ext::no_trans;
            }
        };
        case trans_type_ext::trans:
        {
            switch (t2)
            {
                case trans_type_ext::no_trans:      return trans_type_ext::trans;
                case trans_type_ext::trans:         return trans_type_ext::no_trans;
                case trans_type_ext::conj_trans:    return trans_type_ext::conj;
                case trans_type_ext::conj:          return trans_type_ext::conj_trans;
            }
        };
        case trans_type_ext::conj_trans:
        {
            switch (t2)
            {
                case trans_type_ext::no_trans:      return trans_type_ext::conj_trans;
                case trans_type_ext::trans:         return trans_type_ext::conj;
                case trans_type_ext::conj_trans:    return trans_type_ext::no_trans;
                case trans_type_ext::conj:          return trans_type_ext::trans;
            }
        }
        default:
            matcl_assert(0,"unknown case");
            throw;
    };
};

trans_type_ext trans_manip::convert_trans(trans_type t)
{
    switch (t)
    {
        case trans_type::no_trans:      return trans_type_ext::no_trans;
        case trans_type::trans:         return trans_type_ext::trans;
        case trans_type::conj_trans:    return trans_type_ext::conj_trans;

        //this case should never happened
        default:                        return trans_type_ext::no_trans;
    };
};

trans_type trans_manip::convert_trans(trans_type_ext t, bool& conj)
{
    conj = false;
    switch (t)
    {
        case trans_type_ext::no_trans:      return trans_type::no_trans;
        case trans_type_ext::trans:         return trans_type::trans;
        case trans_type_ext::conj_trans:    return trans_type::conj_trans;

        case trans_type_ext::conj:          
        {
            conj = true;
            return trans_type::no_trans;
        }

        //this case should never happened
        default:                        return trans_type::no_trans;
    };
};

//----------------------------------------------------------------
//                  intantiation
//----------------------------------------------------------------
template struct unit_matrix<Integer,struct_dense>;
template struct unit_matrix<Real,struct_dense>;
template struct unit_matrix<Float,struct_dense>;
template struct unit_matrix<Complex,struct_dense>;
template struct unit_matrix<Float_complex,struct_dense>;

template struct unit_matrix<Integer,struct_banded>;
template struct unit_matrix<Real,struct_banded>;
template struct unit_matrix<Float,struct_banded>;
template struct unit_matrix<Complex,struct_banded>;
template struct unit_matrix<Float_complex,struct_banded>;

template struct unit_matrix<Integer,struct_sparse>;
template struct unit_matrix<Real,struct_sparse>;
template struct unit_matrix<Float,struct_sparse>;
template struct unit_matrix<Complex,struct_sparse>;
template struct unit_matrix<Float_complex,struct_sparse>;

};};
