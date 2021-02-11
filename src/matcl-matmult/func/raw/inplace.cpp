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

#include "matcl-internals/func/inplace.h"
#include "matcl-internals/base/utils.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/lib_functions/func_unary.h"
#include "matcl-matmult/func/raw/mmul/mmul_utils.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-internals/func/raw_manip.h"

namespace matcl { namespace raw { namespace inplace
{

namespace mrd = matcl::raw::details;

template<class val_type, class str_type>
struct copy_utril_to_ltril_impl
{};

template<class val_type>
struct copy_utril_to_ltril_impl<val_type, struct_dense>
{
    using matrix_type = Matrix<val_type, struct_dense>;

    static void eval(matrix_type& mat, bool conj_transpose)
    {
        if (conj_transpose == true)
            return eval_conj<true>(mat);
        else
            return eval_conj<false>(mat);
    };

    template<bool Conj>
    static void eval_conj(matrix_type& mat)
    {
        Integer M       = mat.rows();
        Integer N       = mat.cols();
        Integer K       = std::min(M,N);
        Integer mat_ld  = mat.ld();

        //block size
        const Integer NB    = optim_params::block_size_trans;

        //blocked version, reduce cache misses in matrix transpose
        for (Integer ib = 0; ib < K; ib += NB)
        {
            for (Integer jb = 0; jb < K; jb += NB)
            {
                Integer ibl     = std::min(K, ib + NB);
                Integer jbl     = std::min(K, jb + NB);

                val_type* ptr_m = mat.ptr() + ib;
                val_type* ptr_i = mat.ptr() + ib * mat_ld;

                for (Integer i = ib; i < ibl; ++i)
                {                    
                    for (Integer j = jb; j < jbl; ++j)
                    {
                        val_type at = mrd::make_conj<Conj>::eval(ptr_m[j*mat_ld]);
                        mrd::reset_helper(ptr_i[j], at);
                    }

                    ptr_i       += mat_ld;
                    ptr_m       += 1;
                };
            };
        };
    };
};

template<class val_type>
struct copy_utril_to_ltril_impl<val_type, struct_banded>
{
    using matrix_type = Matrix<val_type, struct_banded>;
    static void eval(matrix_type& mat, bool conj_transpose)
    {
        Integer fd      = mat.first_diag();
        Integer ld      = mat.last_diag();
        Integer mat_ld  = mat.ld();

        if (fd == 0)
            return;

        val_type Z = matcl::details::default_value<val_type>(mat.get_type());

        for (Integer d = fd; d < -ld; ++ d)
        {
            Integer s       = mat.diag_length(d);
            val_type* ptr   = mat.rep_ptr() + mat.first_elem_diag(d);

            for (Integer i = 0; i < s; ++i)
            {
                mrd::reset_helper(ptr[0], Z);
                ptr         += mat_ld;
            };
        }
        for (Integer d = -ld; d < 0; ++d)
        {
            Integer su      = mat.diag_length(-d);
            Integer sl      = mat.diag_length(d);
            Integer s       = std::min(su, sl);

            val_type* ptr_l = mat.rep_ptr() + mat.first_elem_diag(d);
            val_type* ptr_u = mat.rep_ptr() + mat.first_elem_diag(-d);

            if (conj_transpose == false)
            {
                for (Integer i = 0; i < s; ++i)
                {
                    mrd::reset_helper(ptr_l[0],ptr_u[0]);
                    ptr_l   += mat_ld;
                    ptr_u   += mat_ld;
                };
            }
            else
            {
                for (Integer i = 0; i < s; ++i)
                {
                    val_type ac    = mrd::conj_helper<val_type>::eval(ptr_u[0]);

                    mrd::reset_helper(ptr_l[0], ac);
                    ptr_l   += mat_ld;
                    ptr_u   += mat_ld;
                };
            };

            for (Integer i = s; i < sl; ++i)
            {
                mrd::reset_helper(ptr_l[0], Z);
                ptr_l       += mat_ld;
            };
        };
    };
};

template<class matrix_type>
void copy_utril_to_ltril<matrix_type>::eval(matrix_type& mat, bool conj_transpose)
{
    using str_type  = typename matrix_type::struct_type;
    using val_type  = typename matrix_type::value_type;

    return copy_utril_to_ltril_impl<val_type, str_type>::eval(mat,conj_transpose);
};

template<class matrix_type>
void make_triu<matrix_type>::eval(matcl::Matrix& ret, matrix_type& mat, Integer d)
{
    mrd::manip_tr_helper<matrix_type>::eval_triu(ret,mat,d,true);
}

template<class matrix_type>
void make_tril<matrix_type>::eval(matcl::Matrix& ret, matrix_type& mat, Integer d)
{
    mrd::manip_tr_helper<matrix_type>::eval_tril(ret,mat,d,true);
};

}}};

MACRO_INSTANTIATE_DENSE_1(matcl::raw::inplace::copy_utril_to_ltril)
MACRO_INSTANTIATE_BAND_1(matcl::raw::inplace::copy_utril_to_ltril)

MACRO_INSTANTIATE_DENSE_1(matcl::raw::inplace::make_tril)
MACRO_INSTANTIATE_BAND_1(matcl::raw::inplace::make_tril)

MACRO_INSTANTIATE_DENSE_1(matcl::raw::inplace::make_triu)
MACRO_INSTANTIATE_BAND_1(matcl::raw::inplace::make_triu)
