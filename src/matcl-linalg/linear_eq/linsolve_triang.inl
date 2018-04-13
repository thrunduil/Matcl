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
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type2_switch.h"
#include "linsolve_general.inl"
#include "linsolve_sparse.inl"

namespace matcl { namespace details
{

//--------------------------------------------------------------------------------
//                  STRUCT TYPES
//--------------------------------------------------------------------------------
template<class S1, class S2, class V1, class V2>
struct linsolve_triang_str
{};

//--------------------------------------------------------------------------------
//                  DENSE-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_dense,struct_dense,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_dense>;
    using M2    = raw::Matrix<V2,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                     matcl::permvec q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        (void)opts;

        char uplo, diag = 'N';

        if (is_lt)  uplo = 'L';
        else        uplo = 'U';

        Integer N       = A.rows();
        Integer Nrhs    = B.cols();

        matcl::lapack::i_type info = 0;
        
        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_dense>;

        M Ac    = raw::converter<M,M1>::eval(A);
        M Bc    = raw::converter<M,M2>::eval(B.make_unique());

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute rows of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), 1);
        };

        lapack::trtrs(&uplo, get_trans_code(trans), &diag, N, Nrhs, lap(Ac.ptr()), Ac.ld(), 
                      lap(Bc.ptr()), Bc.ld(), &info);

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
                if (is_lt == true && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else if (is_lt == false && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else
                    Bc.get_struct().reset();
            }
            else
            {
                if (is_lt == true && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else if (is_lt == false && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else
                    Bc.get_struct().reset();
            };
        };

        ret = Matrix(Bc,false);
        return;
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, matcl::permvec p, matcl::permvec q, 
                                const M2& B, bool is_lt, trans_type trans, const options& opts)
    {
        (void)opts;

        char uplo, diag = 'N';

        if (is_lt)  uplo = 'L';
        else        uplo = 'U';

        Integer N       = A.rows();
        Integer Mrhs    = B.rows();

        matcl::lapack::i_type info = 0;
        
        using V     = typename unify_types<V1,V2>::type;
        using M     = raw::Matrix<V,struct_dense>;

        M Ac    = raw::converter<M,M1>::eval(A);
        M Bc    = raw::converter<M,M2>::eval(B.make_unique());

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute columns of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, q.length(), lap(q_ptr), 1);
        };

        lapack::trtrs_rev(&uplo, get_trans_code(trans), &diag, N, Mrhs, lap(Ac.ptr()), Ac.ld(), 
                          lap(Bc.ptr()), Bc.ld(), &info);

        if (info)
            throw error::error_singular();

        //permute columns of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), -1);
        };

        if (p.is_id() == false || q.is_id() == false)
        {
            Bc.get_struct().reset();
        }
        else
        {
            if (trans == trans_type::no_trans)
            {
                if (is_lt == true && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else if (is_lt == false && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else
                    Bc.get_struct().reset();
            }
            else
            {
                if (is_lt == true && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else if (is_lt == false && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else
                    Bc.get_struct().reset();
            };
        };

        ret = Matrix(Bc,false);
        return;
    };
};

//--------------------------------------------------------------------------------
//                  DENSE-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_dense,struct_sparse,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_dense>;
    using M2    = raw::Matrix<V2,struct_sparse>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_dense,struct_dense,V1,V2>::eval(ret, A,p,q, Bc,is_lt,trans, opts);
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_dense,struct_dense,V1,V2>::eval_rev(ret,A,p,q,Bc,is_lt,trans,opts);
    };
};

//--------------------------------------------------------------------------------
//                  DENSE-BANDED
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_dense,struct_banded,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_dense>;
    using M2    = raw::Matrix<V2,struct_banded>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_dense,struct_dense,V1,V2>::eval(ret,A,p,q, Bc,is_lt,trans,opts);
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_dense,struct_dense,V1,V2>::eval_rev(ret, A,p,q, Bc,is_lt,trans, opts);
    };
};

//--------------------------------------------------------------------------------
//                  BANDED-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_banded,struct_dense,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_banded>;
    using M2    = raw::Matrix<V2,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                     matcl::permvec q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        (void)opts;

        char uplo;
        if (is_lt)      uplo = 'L';
        else            uplo = 'U';

        Integer N       = A.rows();

        if (A.has_diag(0) == false)
            throw error::error_singular();

        Integer kd      = (is_lt? A.number_subdiagonals() : A.number_superdiagonals());
        Integer Nrhs    = B.cols();

        matcl::lapack::i_type info = 0;

        using V     = typename unify_types<V1,V2>::type;
        using MA    = raw::Matrix<V,struct_banded>;
        using MB    = raw::Matrix<V,struct_dense>;

        MA Ac   = raw::converter<MA,M1>::eval(A);
        MB Bc   = raw::converter<MB,M2>::eval(B.make_unique());

        Integer off = A.first_elem_diag(0);
       
        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute rows of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswp(Bc.cols(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), 1);
        };

        lapack::tbtrs(&uplo, get_trans_code(trans), "N", N, kd, Nrhs, lap(Ac.rep_ptr() + off), Ac.ld(), 
                        lap(Bc.ptr()), Bc.ld(), &info);

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
                if (is_lt == true && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else if (is_lt == false && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else
                    Bc.get_struct().reset();
            }
            else
            {
                if (is_lt == true && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else if (is_lt == false && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else
                    Bc.get_struct().reset();
            };
        };

        ret = Matrix(Bc,true);
        return;
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, matcl::permvec p, matcl::permvec q, 
                                  const M2& B, bool is_lt, trans_type trans, const options& opts)
    {
        (void)opts;

        char uplo;
        if (is_lt)      uplo = 'L';
        else            uplo = 'U';

        if (A.has_diag(0) == false)
            throw error::error_singular();

        Integer N       = A.rows();
        Integer kd      = (is_lt? A.number_subdiagonals() : A.number_superdiagonals());
        Integer Mrhs    = B.rows();

        matcl::lapack::i_type info = 0;

        using V     = typename unify_types<V1,V2>::type;
        using MA    = raw::Matrix<V,struct_banded>;
        using MB    = raw::Matrix<V,struct_dense>;

        MA Ac   = raw::converter<MA,M1>::eval(A);
        MB Bc   = raw::converter<MB,M2>::eval(B.make_unique());

        Integer off = Ac.first_elem_diag(0);
       
        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute columns of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, q.length(), lap(q_ptr), 1);
        };

        lapack::tbtrs_rev(&uplo, get_trans_code(trans), "N", N, kd, Mrhs, lap(Ac.rep_ptr() + off), Ac.ld(), 
                        lap(Bc.ptr()), Bc.ld(), &info);

        if (info)
            throw error::error_singular();

        //permute columns of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswpc(Bc.rows(), lap(Bc.ptr()), Bc.ld(), 1, p.length(), lap(p_ptr), -1);
        };

        if (p.is_id() == false || q.is_id() == false)
        {
            Bc.get_struct().reset();
        }
        else
        {
            if (trans == trans_type::no_trans)
            {
                if (is_lt == true && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else if (is_lt == false && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else
                    Bc.get_struct().reset();
            }
            else
            {
                if (is_lt == true && B.get_struct().is_triu())
                    Bc.get_struct().set(predefined_struct_type::triu);
                else if (is_lt == false && B.get_struct().is_tril())
                    Bc.get_struct().set(predefined_struct_type::tril);
                else
                    Bc.get_struct().reset();
            };
        };

        ret = Matrix(Bc,true);
        return;
    };
};

//--------------------------------------------------------------------------------
//                  BANDED-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_banded,struct_sparse,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_banded>;
    using M2    = raw::Matrix<V2,struct_sparse>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_banded,struct_dense,V1,V2>::eval(ret, A,p,q,Bc,is_lt,trans, opts);
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_banded,struct_dense,V1,V2>::eval_rev(ret, A,p,q,Bc,is_lt,trans, opts);
    };
};

//--------------------------------------------------------------------------------
//                  BANDED-BANDED
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_banded,struct_banded,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_banded>;
    using M2    = raw::Matrix<V2,struct_banded>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_banded,struct_dense,V1,V2>::eval(ret, A,p,q,Bc,is_lt,trans, opts);
    };
    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const options& opts)
    {
        using M2c   = raw::Matrix<V2,struct_dense>;
        M2c Bc      = raw::converter<M2c,M2>::eval(B);

        return linsolve_triang_str<struct_banded,struct_dense,V1,V2>::eval_rev(ret, A,p,q,Bc,is_lt,trans, opts);
    };
};

//--------------------------------------------------------------------------------
//                  SPARSE-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_sparse,struct_dense,V1,V2>
    : public linsolve_triang_sparse<struct_dense, V1, V2>
{};

//--------------------------------------------------------------------------------
//                  SPARSE-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_sparse,struct_sparse,V1,V2>
    : public linsolve_triang_sparse<struct_sparse, V1, V2>
{};

//--------------------------------------------------------------------------------
//                  SPARSE-BAND
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_str<struct_sparse,struct_banded,V1,V2>
    : public linsolve_triang_sparse<struct_banded, V1, V2>
{};

}};