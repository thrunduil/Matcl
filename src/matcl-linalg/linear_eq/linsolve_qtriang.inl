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

#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "linsolve_general.inl"
#include "linsolve_triang.inl"

namespace matcl { namespace details
{

//--------------------------------------------------------------------------------
//                  STRUCT TYPES
//--------------------------------------------------------------------------------

template<class S1, class S2, class V1, class V2>
struct linsolve_qtriang_str
{
    using M1    = raw::Matrix<V1,S1>;
    using M2    = raw::Matrix<V2,S2>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, trans_type trans,
                     const matcl::options& opts)
    {
        //it should aready be checked that A, B are real dense

        using V     = typename details::unify_types<V1,V2>::type;
        using VR    = typename details::real_type<V>::type;

        using MD_1  = raw::Matrix<VR,struct_dense>;
        using MD_2  = raw::Matrix<VR,struct_dense>;

        MD_1 Ac     = raw::converter<MD_1, M1>::eval(A);
        MD_2 Bc     = raw::converter<MD_2, M2>::eval(B);

        return linsolve_qtriang_str<struct_dense, struct_dense, VR, VR>::eval(ret, Ac, p,q,Bc, trans,opts);
    };
};

template<class V, bool Is_compl = details::is_complex<V>::value>
struct linsolve_qtriang2_impl
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                     matcl::permvec q, const M2& B, trans_type trans,
                     const matcl::options& opts)
    {
        (void)opts;

        //A and B are converted to real dense matrices

        Integer N       = A.rows();
        Integer Nrhs    = B.cols();

        matcl::lapack::i_type info = 0;

        if (trans != trans_type::no_trans)
            std::swap(p, q);
                
        M2 Bc   = B.make_unique();

        //permute rows of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), 1);
        };
        
        lapack::qtrtrs(get_trans_code(trans), N, Nrhs, lap(A.ptr()), A.ld(), lap(Bc.ptr()), Bc.ld(), 
                        &info);

        if (info)
            throw error::error_singular();

        //permute rows of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, q.length(), lap(q_ptr), -1);
        };		

        if (p.is_id() == false || q.is_id() == false)
        {
            Bc.get_struct().reset();
        }
        else
        {
            if (trans == trans_type::no_trans)
            {
                if (B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::qtriu);
                else
                    Bc.get_struct().reset();
            }
            else
            {
                if (B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::qtril);
                else
                    Bc.get_struct().reset();
            };
        };

        ret = Matrix(Bc,false);
        return;
    };
};

template<class V>
struct linsolve_qtriang2_impl<V,true>
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                     matcl::permvec q, const M2& B, trans_type trans,
                     const matcl::options& opts)
    {
        //it should aready be checked that A, B are real dense
        using VR    = typename details::real_type<V>::type;

        using MD_1  = raw::Matrix<VR,struct_dense>;
        using MD_2  = raw::Matrix<VR,struct_dense>;

        MD_1 Ac     = raw::converter<MD_1, M1>::eval(A);
        MD_2 Bc     = raw::converter<MD_2, M2>::eval(B);

        return linsolve_qtriang2_impl<VR>::eval(ret, Ac, p,q,Bc, trans,opts);
    };
};

template<class V>
struct linsolve_qtriang_str<struct_dense, struct_dense, V, V>
{
    using M1    = raw::Matrix<V,struct_dense>;
    using M2    = raw::Matrix<V,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                     matcl::permvec q, const M2& B, trans_type trans,
                     const matcl::options& opts)
    {
        return linsolve_qtriang2_impl<V>::eval(ret, A, p,q,B,trans,opts);
    };
};

}};