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

#include "matcl-internals/base/utils.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-internals/algs/scatter.h"
#include "matcl-internals/base/sort.h"

namespace matcl { namespace details
{

template<class struct_type, class V1, class V2>
struct linsolve_triang_sparse{};

//--------------------------------------------------------------------------------
//                  TRIANG SPARSE-DENSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_sparse<struct_dense,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_sparse>;
    using M2    = raw::Matrix<V2,struct_dense>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        if (is_lt)      return eval_lt(ret,A,p,q,B,trans,opts);
        else            return eval_ut(ret,A,p,q,B,trans,opts);
    };
    
    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const options& opts)
    {
        if (is_lt)      return eval_lt_rev(ret,A,p,q,B,trans,opts);
        else            return eval_ut_rev(ret,A,p,q,B,trans,opts);
    };
    
    static void eval_lt(matcl::Matrix& ret, const M1& A, matcl::permvec p, matcl::permvec q, 
                                const M2& B, trans_type trans, const options& opts)
    {        
        (void)opts;

        const raw::details::sparse_ccs<V1>& rep = A.rep();
        Integer n               = rep.cols();
        Integer Nrhs            = B.cols();

        const Integer * Ad_c	= rep.ptr_c();
        const Integer * Ad_r	= rep.ptr_r();
        const V1 * Ad_x	        = rep.ptr_x();
      
        for (Integer j = 0; j < n; ++j)
        {
            Integer k = Ad_c[j];

            // First non-zero element may be zero. 
            // Find the first element which is really non-zero
            while(k < Ad_c[j+1] && Ad_x[k] == 0) 
                k++;

            if (k == Ad_c[j+1] || Ad_r[k] != j || Ad_x[k] == 0)
                throw error::error_singular();
        };        

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_dense>;

        M X     = raw::converter<M,M2>::eval(B.make_unique());
        V* x    = X.ptr();

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute rows of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswp(X.cols(), lap(X.ptr()), X.ld(), 1, p.length(), lap(p_ptr), 1);
        };

        if (trans == trans_type::no_trans)
        {
            for (Integer i = 0; i < Nrhs; ++i, x+= X.ld())
            {           
                for (Integer j = 0; j < n; ++j)
                {
                    Integer k = Ad_c[j];
                
                    // avoid 0 elements as above
                    while(Ad_x[k] == 0) 
                        k++;
                
                    V1 val  = Ad_x[k++];
                    x[j]    = x[j]/val;

                    for (; k < Ad_c[j + 1]; ++k)
                        x[Ad_r[k]] -= Ad_x[k] * x[j];
                }
            };

            if (B.get_struct().is_tril())
                X.get_struct().set(predefined_struct_type::tril);
            else
                X.get_struct().reset();
        }
        else if (trans == trans_type::trans)
        {
            for (Integer i = 0; i < Nrhs; ++i, x+= X.ld())
            {           
                for (Integer j = n - 1; j >= 0; --j)
                {
                    Integer k = Ad_c[j];

                    // avoid 0 elements as above
                    while(Ad_x[k] == 0) 
                        k++;

                    V1 val  = Ad_x[k++];
                    
                    for (; k < Ad_c[j + 1]; ++k)
                        x[j] -= Ad_x[k] * x[Ad_r[k]];

                    x[j] = x[j]/val;
                }
            };

            if (X.get_struct().is_triu())
                X.get_struct().set(predefined_struct_type::triu);
            else
                X.get_struct().reset();
        }
        else if (trans == trans_type::conj_trans)
        {
            for (Integer i = 0; i < Nrhs; ++i, x+= X.ld())
            {           
                for (Integer j = n - 1; j >= 0; --j)
                {
                    Integer k = Ad_c[j];

                    // avoid 0 elements as above
                    while(Ad_x[k] == 0) 
                        k++;

                    V1 val  = matcl::conj(Ad_x[k++]);
                    
                    for (; k < Ad_c[j + 1]; ++k)
                        x[j] -= matcl::conj(Ad_x[k]) * x[Ad_r[k]];

                    x[j] = x[j]/val;
                }
            };

            if (X.get_struct().is_triu())
                X.get_struct().set(predefined_struct_type::triu);
            else
                X.get_struct().reset();
        }
        else
        {
            matcl_assert(0,"unknown case");
            throw;
        };

        //permute rows of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswp(X.cols(), lap(X.ptr()), X.ld(), 1, q.length(), lap(q_ptr), -1);
        };	

        if (p.is_id() == false || q.is_id() == false)
            X.get_struct().reset();

        ret = Matrix(X,true);
        return;
    };

    static void eval_ut(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                        matcl::permvec q, const M2& B, trans_type trans,
                        const options& opts)
    {        
        (void)opts;

        const raw::details::sparse_ccs<V1>& rep = A.rep();
        Integer n               = rep.cols();
        Integer Nrhs            = B.cols();

        const Integer * Ad_c	= rep.ptr_c();
        const Integer * Ad_r	= rep.ptr_r();
        const V1 * Ad_x	        = rep.ptr_x();

        for (Integer j = n-1; j >= 0; --j)
        {
            Integer k = Ad_c[j + 1] - 1;

            // Last non-zero element may be zero
            // Find the element which is non-zero
            while(k > Ad_c[j] && Ad_x[k] == 0) 
                k--;

            if (k < Ad_c[j] || Ad_r[k] != j || Ad_x[k] == 0)
                throw error::error_singular();
        };

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_dense>;

        M X     = raw::converter<M,M2>::eval(B.make_unique());
        V* x    = X.ptr();

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute rows of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswp(X.cols(), lap(X.ptr()), X.ld(), 1, p.length(), lap(p_ptr), 1);
        };

        if (trans == trans_type::no_trans)
        {
            for (Integer i = 0; i < Nrhs; ++i, x+= X.ld())
            {           
                for (Integer j = n-1; j >= 0; --j)
                {
                    Integer k = Ad_c[j + 1] - 1;
              
                    // avoid 0 elements as above
                    while(Ad_x[k] == 0) 
                        k--;

                    V1 val  = Ad_x[k--];
                    x[j]    = x[j]/val;

                    for (; k >= Ad_c[j]; --k)
                        x[Ad_r[k]] -= Ad_x[k] * x[j];
                }
            };

            if (X.get_struct().is_triu())
                X.get_struct().set(predefined_struct_type::triu);
            else
                X.get_struct().reset();
        }
        else if (trans == trans_type::trans)
        {
            for (Integer i = 0; i < Nrhs; ++i, x+= X.ld())
            {           
                for (Integer j = 0; j < n; ++j)
                {
                    Integer k = Ad_c[j + 1] - 1;
              
                    // avoid 0 elements as above
                    while(Ad_x[k] == 0) 
                        k--;
                
                    V1 val  = Ad_x[k--];                    

                    for (; k >= Ad_c[j]; --k)
                        x[j] -= Ad_x[k] * x[Ad_r[k]];

                    x[j]    = x[j]/val;
                }
            };

            if (B.get_struct().is_tril())
                X.get_struct().set(predefined_struct_type::tril);
            else
                X.get_struct().reset();
        }
        else if (trans == trans_type::conj_trans)
        {
            for (Integer i = 0; i < Nrhs; ++i, x+= X.ld())
            {           
                for (Integer j = 0; j < n; ++j)
                {
                    Integer k = Ad_c[j + 1] - 1;
              
                    // avoid 0 elements as above
                    while(Ad_x[k] == 0) 
                        k--;
                
                    V1 val  = matcl::conj(Ad_x[k--]);                    

                    for (; k >= Ad_c[j]; --k)
                        x[j] -= matcl::conj(Ad_x[k]) * x[Ad_r[k]];

                    x[j]    = x[j]/val;
                }
            };

            if (B.get_struct().is_tril())
                X.get_struct().set(predefined_struct_type::tril);
            else
                X.get_struct().reset();
        }
        else
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        };

        //permute rows of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswp(X.cols(), lap(X.ptr()), X.ld(), 1, q.length(), lap(q_ptr), -1);
        };		

        if (p.is_id() == false || q.is_id() == false)
            X.get_struct().reset();

        ret = Matrix(X,true);
        return;
    };

    static void eval_lt_rev(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                            matcl::permvec q, const M2& B, trans_type trans,
                            const options& opts)
    {
        (void)opts;

        const raw::details::sparse_ccs<V1>& rep = A.rep();
        Integer n               = rep.cols();
        Integer Mrhs            = B.rows();

        const Integer * Ad_c	= rep.ptr_c();
        const Integer * Ad_r	= rep.ptr_r();
        const V1 * Ad_x	        = rep.ptr_x();
      
        for (Integer j = 0; j < n; ++j)
        {
            Integer k = Ad_c[j];

            // First non-zero element may be zero. 
            // Find the first element which is non-zero
            while(k < Ad_c[j+1] && Ad_x[k] == 0) 
                k++;

            if (k == Ad_c[j+1] || Ad_r[k] != j || Ad_x[k] == 0)
                throw error::error_singular();
        };        

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_dense>;

        M X         = raw::converter<M,M2>::eval(B.make_unique());
        V* x        = X.ptr();

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute columns of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswpc(X.rows(), lap(X.ptr()), X.ld(), 1, q.length(), lap(q_ptr), 1);
        };

        if (trans == trans_type::no_trans)
        {
            V* ptr_col  = x + (n-1) * X.ld();

            for (Integer j = n-1; j >= 0; --j, ptr_col -= X.ld())
            {
                Integer k = Ad_c[j];
                
                // avoid 0 elements as above
                while(Ad_x[k] == 0) 
                    k++;

                V1 val  = Ad_x[k++];

                for (; k < Ad_c[j + 1]; ++k)
                {
                    Integer pos = Ad_r[k];
                    V1 mult     = Ad_x[k];

                    if (mult == 0)
                        continue;

                    V* ptr_pos  = x + pos * X.ld();

                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] -= mult * ptr_pos[i];
                };

                if (val != 1)
                {
                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] = ptr_col[i] / val;
                };
            };

            if (B.get_struct().is_tril())
                X.get_struct().set(predefined_struct_type::tril);
            else
                X.get_struct().reset();            
        }
        else if (trans == trans_type::trans)
        {
            V* ptr_col  = x;

            for (Integer j = 0; j < n; ++j, ptr_col += X.ld())
            {
                Integer k = Ad_c[j];
                
                // avoid 0 elements as above
                while(Ad_x[k] == 0) 
                    k++;

                V1 val  = Ad_x[k++];

                if (val != 1)
                {
                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] = ptr_col[i] / val;
                };

                for (; k < Ad_c[j + 1]; ++k)
                {
                    Integer pos = Ad_r[k];
                    V1 mult     = Ad_x[k];

                    if (mult == 0)
                        continue;

                    V* ptr_pos  = x + pos * X.ld();

                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_pos[i] -= mult * ptr_col[i];
                };
            };

            if (B.get_struct().is_triu())
                X.get_struct().set(predefined_struct_type::triu);
            else
                X.get_struct().reset();
        }
        else if (trans == trans_type::conj_trans)
        {
            V* ptr_col  = x;

            for (Integer j = 0; j < n; ++j, ptr_col += X.ld())
            {
                Integer k = Ad_c[j];
                
                // avoid 0 elements as above
                while(Ad_x[k] == 0) 
                    k++;

                V1 val  = matcl::conj(Ad_x[k++]);

                if (val != 1)
                {
                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] = ptr_col[i] / val;
                };

                for (; k < Ad_c[j + 1]; ++k)
                {
                    Integer pos = Ad_r[k];
                    V1 mult     = matcl::conj(Ad_x[k]);

                    if (mult == 0)
                        continue;

                    V* ptr_pos  = x + pos * X.ld();

                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_pos[i] -= mult * ptr_col[i];
                };
            };

            if (B.get_struct().is_triu())
                X.get_struct().set(predefined_struct_type::triu);
            else
                X.get_struct().reset();
        }
        else
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        };

        //permute columns of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswpc(X.rows(), lap(X.ptr()), X.ld(), 1, p.length(), lap(p_ptr), -1);
        };

        if (p.is_id() == false || q.is_id() == false)
            X.get_struct().reset();

        ret = Matrix(X,true);
        return;
    };

    static void eval_ut_rev(matcl::Matrix& ret, const M1& A, matcl::permvec p, 
                            matcl::permvec q, const M2& B, trans_type trans,
                            const options& opts)
    {
        (void)opts;

        const raw::details::sparse_ccs<V1>& rep = A.rep();
        Integer n               = rep.cols();
        Integer Mrhs            = B.rows();

        const Integer * Ad_c	= rep.ptr_c();
        const Integer * Ad_r	= rep.ptr_r();
        const V1 * Ad_x	        = rep.ptr_x();

        for (Integer j = n-1; j >= 0; --j)
        {
            Integer k = Ad_c[j + 1] - 1;

            // Last non-zero element may be zero
            // Find the element which is really non-zero
            while(k > Ad_c[j] && Ad_x[k] == 0) 
                k--;

            if (k < Ad_c[j] || Ad_r[k] != j || Ad_x[k] == 0)
                throw error::error_singular();
        };

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_dense>;

        M X         = raw::converter<M,M2>::eval(B.make_unique());
        V* x        = X.ptr();        

        if (trans != trans_type::no_trans)
            std::swap(p,q);

        //permute columns of B with q
        if (q.is_id() == false)
        {
            Matrix q_int    = q.to_interchanges_matrix();
            Integer* q_ptr  = q_int.get_array_unique<Integer>();

            lapack::laswpc(X.rows(), lap(X.ptr()), X.ld(), 1, q.length(), lap(q_ptr), 1);
        };

        if (trans == trans_type::no_trans)
        {
            V* ptr_col  = x;

            for (Integer j = 0; j < n; ++j, ptr_col += X.ld())
            {      
                Integer k = Ad_c[j + 1] - 1;

                // avoid 0 elements as above
                while(Ad_x[k] == 0) 
                    k--;

                V1 val  = Ad_x[k--];

                for (; k >= Ad_c[j]; --k)
                {
                    Integer pos = Ad_r[k];
                    V1 mult     = Ad_x[k];

                    if (mult == 0)
                        continue;

                    V* ptr_pos  = x + pos * X.ld();

                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] -= mult * ptr_pos[i];
                };

                if (val != 1)
                {
                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] = ptr_col[i] / val;
                };
            };

            if (X.get_struct().is_triu())
                X.get_struct().set(predefined_struct_type::triu);
            else
                X.get_struct().reset();            
        }
        else if (trans == trans_type::trans)
        {
            V* ptr_col  = x + (n-1)*X.ld();

            for (Integer j = n-1; j >= 0; --j, ptr_col -= X.ld())
            {      
                Integer k = Ad_c[j + 1] - 1;

                // avoid 0 elements as above
                while(Ad_x[k] == 0) 
                    k--;

                V1 val  = Ad_x[k--];

                if (val != 1)
                {
                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] = ptr_col[i] / val;
                };

                for (; k >= Ad_c[j]; --k)
                {
                    Integer pos = Ad_r[k];
                    V1 mult     = Ad_x[k];

                    if (mult == 0)
                        continue;

                    V* ptr_pos  = x + pos * X.ld();

                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_pos[i] -= mult * ptr_col[i];
                };
            };

            if (X.get_struct().is_tril())
                X.get_struct().set(predefined_struct_type::tril);
            else
                X.get_struct().reset();
        }
        else if (trans == trans_type::conj_trans)
        {
            V* ptr_col  = x + (n-1)*X.ld();

            for (Integer j = n-1; j >= 0; --j, ptr_col -= X.ld())
            {      
                Integer k = Ad_c[j + 1] - 1;

                // avoid 0 elements as above
                while(Ad_x[k] == 0) 
                    k--;

                V1 val  = matcl::conj(Ad_x[k--]);

                if (val != 1)
                {
                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_col[i] = ptr_col[i] / val;
                };

                for (; k >= Ad_c[j]; --k)
                {
                    Integer pos = Ad_r[k];
                    V1 mult     = matcl::conj(Ad_x[k]);

                    if (mult == 0)
                        continue;

                    V* ptr_pos  = x + pos * X.ld();

                    for (Integer i = 0; i < Mrhs; ++i)
                        ptr_pos[i] -= mult * ptr_col[i];
                };
            };

            if (X.get_struct().is_tril())
                X.get_struct().set(predefined_struct_type::tril);
            else
                X.get_struct().reset();
        }
        else
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        };

        //permute columns of B with p
        if (p.is_id() == false)
        {
            Matrix p_int    = p.to_interchanges_matrix();
            Integer* p_ptr  = p_int.get_array_unique<Integer>();

            lapack::laswpc(X.rows(), lap(X.ptr()), X.ld(), 1, p.length(), lap(p_ptr), -1);
        };

        if (p.is_id() == false || q.is_id() == false)
            X.get_struct().reset();

        ret = Matrix(X,true);
        return;
    };
};

//--------------------------------------------------------------------------------
//                  TRIANG SPARSE-SPARSE
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_sparse<struct_sparse,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_sparse>;
    using M2    = raw::Matrix<V2,struct_sparse>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        if (is_lt)      return eval_lt(ret,A,p,q,B,trans,opts);
        else            return eval_ut(ret,A,p,q,B,trans,opts);
    };
    
    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const options& opts)
    {
        if (is_lt)      return eval_lt_rev(ret,A,p,q,B,trans,opts);
        else            return eval_ut_rev(ret,A,p,q,B,trans,opts);
    };

    static void eval_lt(matcl::Matrix& ret, const M1& A, matcl::permvec p, matcl::permvec q, 
                        const M2& B, trans_type trans, const options& opts)
    {
        (void)opts;

        Integer n = A.cols();
        Integer m = B.cols();

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_sparse>;

        const raw::details::sparse_ccs<V1>& rep_A = A.rep();
        const raw::details::sparse_ccs<V2>& rep_B = B.rep();

        const Integer * A_c	= rep_A.ptr_c();
        const Integer * A_r	= rep_A.ptr_r();
        const V1 * A_x	    = rep_A.ptr_x();

        for (Integer j = 0; j < n; ++j)
        {
            Integer k = A_c[j];
            
            // first non-zero element may be zero
            // find the first element which is really non-zero
            while(k < A_c[j+1] && A_x[k] == 0) 
                k++;

            if (k == A_c[j+1] || A_r[k] != j || A_x[k] == 0)
                throw error::error_singular();
        };

        M C(A.get_type(),n,m,n);
        raw::details::sparse_ccs<V>& out = C.rep();
        
        using workspace = matcl::details::workspace2<V>;
        using scatter   = matcl::algorithm::scatter;

        workspace           x(C.get_type(), n);
        
        scatter sc_w        = matcl::algorithm::scatter::get(n,m);

        Integer * out_c	    = out.ptr_c();
        Integer * out_r	    = out.ptr_r();
        V * out_x	        = out.ptr_x();

        const Integer * B_c	= rep_B.ptr_c();
        const Integer * B_r	= rep_B.ptr_r();
        const V2 * B_x	    = rep_B.ptr_x();

        if (trans != trans_type::no_trans)
            std::swap(p,q);            

        p                   = invperm(p);

        bool perm_q         = (q.is_id() == false);
        bool perm_p         = (p.is_id() == false);

        const Integer* ptr_q    = perm_q ? q.to_array() : nullptr;
        const Integer* ptr_p    = perm_p ? p.to_array() : nullptr;

        Integer nz          = 0;

        for (Integer j = 0; j < m; ++j)
        {
            out_c[j]        = nz;

            if (B_c[j] == B_c[j+1])
                continue;

            auto mark       = sc_w.next_mark();

            if (nz + n > out.nzmax())
            {
                out.add_memory(out.nzmax() + n);
                out_r	    = out.ptr_r();
                out_x	    = out.ptr_x();
            };

            //get symbolic structure
            Integer nz_old  = nz;
            Integer last    = 0;

            if (perm_p == true)
            {
                for (Integer kb = B_c[j]; kb < B_c[j+1]; ++kb)
                {
                    V val        = B_x[kb];

                    if (val == 0)
                        continue;
                    
                    Integer i   = ptr_p[B_r[kb]] - 1;
                    last        = std::max(last,i);
                    x[i]        = val;
                    sc_w[i]     = mark;       
                    out_r[nz++] = i;
                }
            }
            else
            {
                for (Integer kb = B_c[j]; kb < B_c[j+1]; ++kb)
                {
                    V val        = B_x[kb];

                    if (val == 0)
                        continue;
                    
                    Integer i   = B_r[kb];
                    x[i]        = val;
                    sc_w[i]     = mark;                    
                    out_r[nz++] = i;
                }

                last            = B_r[B_c[j+1]-1];
            };

            if (trans == trans_type::no_trans)
            {
                for (Integer pos = nz_old; pos < nz; ++pos)
                {
                    Integer jl  = out_r[pos];

                    //visit jl column
                    for (Integer kl = A_c[jl]; kl < A_c[jl + 1]; ++kl)
                    {
                        Integer i = A_r[kl];
                        if (sc_w[i] < mark && A_x[kl] != 0)
                        {
                            sc_w[i]     = mark;
                            x[i]        = 0;
                            out_r[nz++] = i;
                        };
                    };
                }

                //sort row indices
                matcl::utils::sort_q(out_r + nz_old, nz - nz_old);

                //numerical phase
                for (Integer pos = nz_old; pos < nz; ++pos)
                {                
                    Integer jl  = out_r[pos];
                    Integer kl  = A_c[jl];

                    // avoid 0 elements as above
                    while(A_x[kl] == 0) 
                        kl++;

                    x[jl]   = x[jl]/A_x[kl];
                    V xj    = x[jl];

                    for (kl++; kl < A_c[jl + 1]; ++kl)
                    {
                        Integer i   = A_r[kl];
                        x[i]        -= xj * A_x[kl];
                    };
                };
            }
            else
            {
                //get symbolic structure                
                for (Integer jl = last; jl >= 0; --jl)
                {       
                    if (sc_w[jl] == mark)
                        continue;

                    for (Integer kl = A_c[jl]; kl < A_c[jl + 1]; ++kl)
                    {
                        Integer pos = A_r[kl];

                        if (sc_w[pos] == mark && A_x[kl] != 0)
                        {
                            sc_w[jl]    = mark;
                            x[jl]       = 0;
                            out_r[nz++] = jl;
                            break;
                        };
                    };
                };

                //sort row indices
                matcl::utils::sort_q(out_r + nz_old, nz - nz_old);

                //numerical phase
                if (trans == trans_type::trans)
                {                
                    for (Integer pos = nz - 1; pos >= nz_old; --pos)
                    {                
                        Integer jl  = out_r[pos];
                        Integer kl  = A_c[jl];

                        // avoid 0 elements
                        while(A_x[kl] == 0) 
                            kl++;

                        V1 val  = A_x[kl++];
                    
                        for (; kl < A_c[jl + 1]; ++kl)
                        {
                            Integer pos2 = A_r[kl];

                            if (sc_w[pos2] == mark)
                                x[jl] -= A_x[kl] * x[pos2];
                        };

                        x[jl]       = x[jl] / val;
                    };
                }
                else
                {
                    for (Integer pos = nz - 1; pos >= nz_old; --pos)
                    {                
                        Integer jl  = out_r[pos];
                        Integer kl  = A_c[jl];

                        // avoid 0 elements
                        while(A_x[kl] == 0) 
                            kl++;

                        V1 val  = matcl::conj(A_x[kl++]);
                    
                        for (; kl < A_c[jl + 1]; ++kl)
                        {
                            Integer pos2 = A_r[kl];

                            if (sc_w[pos2] == mark)
                                x[jl] -= matcl::conj(A_x[kl]) * x[pos2];
                        };

                        x[jl]       = x[jl] / val;
                    };
                };
            };

            for (Integer k = out_c[j]; k < nz; k++)
            {
                out_x[k] = x[out_r[k]];
            };

            //permute rows with q
            if (perm_q == true)
            {
                for (Integer k = out_c[j]; k < nz; k++)
                    out_r[k] = ptr_q[out_r[k]] - 1;

                Integer nz_old2 = out_c[j];
                matcl::utils::sort_q(out_r + nz_old2, out_x + nz_old2, nz - nz_old2);
            };
        };

        out_c[m] = nz;

        out.add_memory(-1);        

        if (perm_q == true || perm_p == true)
        {
            C.get_struct().reset();
        }
        else
        {
            if (trans == trans_type::no_trans)
            {            
                if (B.get_struct().is_tril())
                    C.get_struct().set(predefined_struct_type::tril);
                else
                    C.get_struct().reset();
            }
            else
            {
                if (B.get_struct().is_triu())
                    C.get_struct().set(predefined_struct_type::triu);
                else
                    C.get_struct().reset();
            };
        };

        ret = Matrix(C,true);
        return;
    };

    static void eval_ut(matcl::Matrix& ret, const M1& A, matcl::permvec p, matcl::permvec q, 
                        const M2& B, trans_type trans, const options& opts)
    {       
        (void)opts;

        Integer n = A.rows();
        Integer m = B.cols();

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_sparse>;

        const raw::details::sparse_ccs<V1>& rep_A = A.rep();
        const raw::details::sparse_ccs<V2>& rep_B = B.rep();

        const Integer * A_c	= rep_A.ptr_c();
        const Integer * A_r	= rep_A.ptr_r();
        const V1 * A_x	    = rep_A.ptr_x();

        for (Integer j = n-1; j >= 0; --j)
        {
            Integer k = A_c[j + 1] - 1;

            // last non-zero element may be zero
            // find the element which is really non-zero!
            while(k > A_c[j] && A_x[k] == 0) 
                k--;

            if (k < A_c[j] || A_r[k] != j || A_x[k] == 0)
                throw error::error_singular();
        };                

        M C(A.get_type(),n,m,n);
        raw::details::sparse_ccs<V>& out = C.rep();

        using workspace = matcl::details::workspace2<V>;
        using scatter   = matcl::algorithm::scatter;

        workspace           x(C.get_type(), n);        
        scatter sc_w        = matcl::algorithm::scatter::get(n,m);

        Integer * out_c	    = out.ptr_c();
        Integer * out_r	    = out.ptr_r();
        V * out_x	        = out.ptr_x();

        const Integer * B_c	= rep_B.ptr_c();
        const Integer * B_r	= rep_B.ptr_r();
        const V2 * B_x	    = rep_B.ptr_x();

        if (trans != trans_type::no_trans)
            std::swap(p,q);            

        p                   = invperm(p);

        bool perm_q         = (q.is_id() == false);
        bool perm_p         = (p.is_id() == false);

        const Integer* ptr_q    = perm_q ? q.to_array() : nullptr;
        const Integer* ptr_p    = perm_p ? p.to_array() : nullptr;

        Integer nz          = 0;

        for (Integer j = 0; j < m; ++j)
        {
            out_c[j]        = nz;

            if (B_c[j] == B_c[j+1])
                continue;

            auto mark       = sc_w.next_mark();

            if (nz + n > out.nzmax())
            {
                out.add_memory(out.nzmax() + n);
                out_r	    = out.ptr_r();
                out_x	    = out.ptr_x();
            };

            Integer nz_old  = nz;
            Integer first   = n;

            if (perm_p == true)
            {
                for (Integer kb = B_c[j]; kb < B_c[j+1]; ++kb)
                {
                    V val        = B_x[kb];

                    if (val == 0)
                        continue;
                    
                    Integer i   = ptr_p[B_r[kb]] - 1;
                    first       = std::min(first,i);
                    x[i]        = val;
                    sc_w[i]     = mark;        
                    out_r[nz++] = i;
                }
            }
            else
            {
                for (Integer kb = B_c[j]; kb < B_c[j+1]; ++kb)
                {
                    V val        = B_x[kb];

                    if (val == 0)
                        continue;
                    
                    Integer i   = B_r[kb];
                    x[i]        = val;
                    sc_w[i]     = mark;                    
                    out_r[nz++] = i;
                }

                first           = B_r[B_c[j]];
            };

            if (trans == trans_type::no_trans)
            {
                //get symbolic structure
                for (Integer pos = nz_old; pos < nz; ++pos)
                {
                    Integer jl  = out_r[pos];

                    //visit jl column
                    for (Integer kl = A_c[jl]; kl < A_c[jl + 1]; ++kl)
                    {
                        Integer i = A_r[kl];
                        if (sc_w[i] < mark && A_x[kl] != 0)
                        {
                            sc_w[i]     = mark;
                            x[i]        = 0;
                            out_r[nz++] = i;
                        };
                    };
                }

                //sort row indices
                matcl::utils::sort_q(out_r + nz_old, nz - nz_old);

                //numerical phase
                for (Integer pos = nz - 1; pos >=  nz_old; --pos)
                {                
                    Integer ju  = out_r[pos];
                    Integer ku  = A_c[ju + 1] - 1;

                    // avoid 0 elements as above
                    while(A_x[ku] == 0)
                        ku--;

                    x[ju]   = x[ju]/A_x[ku];
                    V xj    = x[ju];

                    for (--ku; ku >= A_c[ju]; --ku)
                    {
                        Integer i   = A_r[ku];
                        x[i]        -= xj * A_x[ku];
                    };
                };
            }
            else
            {
                //get symbolic structure                
                for (Integer ju = first; ju < n; ++ju)
                {       
                    if (sc_w[ju] == mark)
                        continue;

                    for (Integer ku = A_c[ju+1]-1; ku >= A_c[ju]; --ku)
                    {
                        Integer pos = A_r[ku];

                        if (sc_w[pos] == mark && A_x[ku] != 0)
                        {
                            sc_w[ju]    = mark;
                            x[ju]       = 0;
                            out_r[nz++] = ju;
                            break;
                        };
                    };
                };

                //sort row indices
                matcl::utils::sort_q(out_r + nz_old, nz - nz_old);

                //numerical phase
                if (trans == trans_type::trans)
                {                
                    for (Integer pos = nz_old; pos < nz; ++pos)
                    {                
                        Integer ju  = out_r[pos];
                        Integer ku  = A_c[ju+1] - 1;

                        // avoid 0 elements
                        while(A_x[ku] == 0) 
                            ku--;

                        V1 val  = A_x[ku--];
                    
                        for (; ku >= A_c[ju]; --ku)
                        {
                            Integer pos2 = A_r[ku];

                            if (sc_w[pos2] == mark)
                                x[ju] -= A_x[ku] * x[pos2];
                        };

                        x[ju]       = x[ju] / val;
                    };
                }
                else
                {
                    for (Integer pos = nz_old; pos < nz; ++pos)
                    {                
                        Integer ju  = out_r[pos];
                        Integer ku  = A_c[ju+1] - 1;

                        // avoid 0 elements
                        while(A_x[ku] == 0) 
                            ku--;

                        V1 val  = matcl::conj(A_x[ku--]);
                    
                        for (; ku >= A_c[ju]; --ku)
                        {
                            Integer pos2 = A_r[ku];

                            if (sc_w[pos2] == mark)
                                x[ju] -= matcl::conj(A_x[ku]) * x[pos2];
                        };

                        x[ju]       = x[ju] / val;
                    };
                };
            };

            for (Integer k = out_c[j]; k < nz; ++k)
                out_x[k] = x[out_r[k]];

            //permute rows with q
            if (perm_q == true)
            {
                for (Integer k = out_c[j]; k < nz; k++)
                    out_r[k] = ptr_q[out_r[k]] - 1;

                Integer nz_old2 = out_c[j];
                matcl::utils::sort_q(out_r + nz_old2, out_x + nz_old2, nz - nz_old2);
            };
        }

        out_c[m] = nz;

        out.add_memory(-1);

        if (perm_q == true || perm_p == true)
        {
            C.get_struct().reset();
        }
        else
        {
            if (trans == trans_type::no_trans)
            {
                if (B.get_struct().is_triu())
                    C.get_struct().set(predefined_struct_type::triu);
                else
                    C.get_struct().reset();
            }
            else
            {
                if (B.get_struct().is_tril())
                    C.get_struct().set(predefined_struct_type::tril);
                else
                    C.get_struct().reset();
            };
        };

        ret = Matrix(C,true);
        return;
    };

    static void eval_lt_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                            const matcl::permvec& q, const M2& B, trans_type trans,
                            const options& opts)
    {        
        Integer n = A.cols();
        Integer m = B.rows();

        switch(trans)
        {
            case trans_type::no_trans:
                break;
            case trans_type::trans:
            {
                if (m < n / 2)
                {
                    eval_lt(ret, A, p,q, matcl::trans(matcl::Matrix(B,false)).impl_unique<M2>(), 
                                              trans_type::no_trans, opts);
                    ret = matcl::trans(ret);
                    return;
                }
                else
                {
                    eval_ut_rev(ret, matcl::trans(Matrix(A,false)).impl_unique<M1>(),q,p, B, trans_type::no_trans,
                                opts);
                    return;
                };
            }
            case trans_type::conj_trans:
            {
                if (m < n / 2)
                {
                    eval_lt(ret, A, p,q, matcl::ctrans(matcl::Matrix(B,false)).impl_unique<M2>(), 
                                              trans_type::no_trans, opts);
                    ret = matcl::ctrans(ret);
                    return;
                }
                else
                {
                    eval_ut_rev(ret, matcl::ctrans(Matrix(A,false)).impl_unique<M1>(),q,p, B, trans_type::no_trans,
                                opts);
                    return;
                };
            }
        };
        
        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_sparse>;

        const raw::details::sparse_ccs<V1>& rep_A = A.rep();
        const raw::details::sparse_ccs<V2>& rep_B = B.rep();

        const Integer * A_c	= rep_A.ptr_c();
        const Integer * A_r	= rep_A.ptr_r();
        const V1 * A_x	    = rep_A.ptr_x();

        for (Integer j = 0; j < n; ++j)
        {
            Integer k = A_c[j];
            
            // first non-zero element may be zero
            // find the first element which is really non-zero
            while(k < A_c[j+1] && A_x[k] == 0) 
                k++;

            if (k == A_c[j+1] || A_r[k] != j || A_x[k] == 0)
                throw error::error_singular();
        };

        const Integer * B_c	= rep_B.ptr_c();
        const Integer * B_r	= rep_B.ptr_r();
        const V2 * B_x	    = rep_B.ptr_x();

        bool perm_p         = (p.is_id() == false);
        bool perm_q         = (q.is_id() == false);

        const Integer* q_ptr    = perm_q ? q.to_array() : nullptr;

        if (m == 1)
        {
            M C(A.get_type(), m, n, n);
            raw::details::sparse_ccs<V>& out = C.rep();

            Integer * out_c	    = out.ptr_c();
            Integer * out_r	    = out.ptr_r();
            V * out_x	        = out.ptr_x();

            Integer nz          = 0;

            for (Integer j = 0; j < n; ++j)
            {
                out_c[j]        = nz;
                Integer col_a   = n - 1 - j;
                Integer col_b   = (perm_q? q_ptr[col_a] - 1 : col_a);
                Integer kl      = A_c[col_a];

                // avoid 0 elements as above
                while(A_x[kl] == 0)
                    ++kl;

                V val           = A_x[kl++];
                V x             = 0;

                if (B_c[col_b+1] > B_c[col_b])
                    x           = B_x[B_c[col_b]];

                //numerical phase
                for (Integer ka = kl; ka < A_c[col_a + 1]; ++ka)
                {
                    Integer ja  = n - 1 - A_r[ka];
                    V elem      = A_x[ka];

                    //visit ja column of out matrix
                    for (Integer ko = out_c[ja]; ko < out_c[ja + 1]; ++ko)
                        x       -= out_x[ko] * elem;
                }

                if (x != 0)
                {
                    out_x[nz]   = x / val;
                    out_r[nz]   = 0;
                    ++nz;
                };
            };

            out_c[n] = nz;

            std::reverse(out_x, out_x + nz);
            std::reverse(out_r, out_r + nz);
            std::reverse(out_c, out_c + n + 1);

            for (Integer i = 0; i <= n; ++i)
                out_c[i] = nz - out_c[i];

            out.add_memory(-1);        
       
            if (perm_q == true || perm_p == true)
            {
                C.get_struct().reset();
            }
            else
            {
                if (B.get_struct().is_tril())
                    C.get_struct().set(predefined_struct_type::tril);
                else
                    C.get_struct().reset();
            };

            ret = Matrix(C,true);

            if (perm_p == true)
            {
                ret = ret(matcl::colon(), invperm(p).to_matrix());
            }

            return;
        };

        M C(A.get_type(), m, n, m + n);
        raw::details::sparse_ccs<V>& out = C.rep();

        Integer * out_c	    = out.ptr_c();
        Integer * out_r	    = out.ptr_r();
        V * out_x	        = out.ptr_x();

        using workspace     = matcl::details::workspace2<V>;
        using scatter       = matcl::algorithm::scatter;

        workspace           x(C.get_type(), m);        
        scatter sc_w        = matcl::algorithm::scatter::get(m, n);

        Integer nz          = 0;

        for (Integer j = 0; j < n; ++j)
        {
            out_c[j]        = nz;            

            if (nz + m > out.nzmax())
            {
                out.add_memory(out.nzmax() + m);
                out_r	    = out.ptr_r();
                out_x	    = out.ptr_x();
            };

            Integer col_a   = n - 1 - j;
            Integer col_b   = (perm_q? q_ptr[col_a] - 1 : col_a);
            Integer kl      = A_c[col_a];

            // avoid 0 elements as above
            while(A_x[kl] == 0)
                ++kl;

            V val           = A_x[kl++];

            //get symbolic structure
            Integer nz_old  = nz;
            auto mark       = sc_w.next_mark();

            for (Integer kb = B_c[col_b]; kb < B_c[col_b+1]; ++kb)
            {
                V tmp       = B_x[kb];

                if (tmp == 0)
                    continue;
                    
                Integer i   = B_r[kb];
                x[i]        = tmp;
                sc_w[i]     = mark;                    
                out_r[nz++] = i;
            }

            for (Integer ka = kl; ka < A_c[col_a + 1]; ++ka)
            {
                Integer ja  = n - 1 - A_r[ka];

                if (A_x[ka] == 0)
                    continue;

                //visit ja column
                for (Integer ko = out_c[ja]; ko < out_c[ja + 1]; ++ko)
                {
                    Integer i = out_r[ko];
                    if (sc_w[i] < mark)
                    {
                        sc_w[i]     = mark;
                        x[i]        = 0;
                        out_r[nz++] = i;
                    };
                };
            }

            //sort row indices in reversed order
            matcl::utils::sort_q_rev(out_r + nz_old, nz - nz_old);

            //numerical phase
            for (Integer ka = kl; ka < A_c[col_a + 1]; ++ka)
            {
                Integer ja  = n - 1 - A_r[ka];
                V elem      = A_x[ka];

                //visit ja column of out matrix
                for (Integer ko = out_c[ja]; ko < out_c[ja + 1]; ++ko)
                {
                    Integer i   = out_r[ko];
                    x[i]        -= out_x[ko] * elem;
                };
            }

            for (Integer k = out_c[j]; k < nz; ++k)
                out_x[k] = x[out_r[k]] / val;
        };

        out_c[n] = nz;

        std::reverse(out_x, out_x + nz);
        std::reverse(out_r, out_r + nz);
        std::reverse(out_c, out_c + n + 1);

        for (Integer i = 0; i <= n; ++i)
            out_c[i] = nz - out_c[i];

        out.add_memory(-1);        
       
        if (perm_q == true || perm_p == true)
        {
            C.get_struct().reset();
        }
        else
        {
            if (B.get_struct().is_tril())
                C.get_struct().set(predefined_struct_type::tril);
            else
                C.get_struct().reset();
        };

        ret = Matrix(C,true);

        if (perm_p)
            ret = ret(matcl::colon(), invperm(p).to_matrix());

        return;
    };

    static void eval_ut_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                            const matcl::permvec& q, const M2& B, trans_type trans,
                            const options& opts)
    {        
        Integer n = A.rows();
        Integer m = B.rows();

        switch(trans)
        {
            case trans_type::no_trans:
                break;
            case trans_type::trans:
            {
                if (m < n / 2)
                {
                    eval_ut(ret, A,p,q, matcl::trans(matcl::Matrix(B,false)).impl_unique<M2>(), 
                                              trans_type::no_trans, opts);
                    ret = matcl::trans(ret);
                    return;
                }
                else
                {
                    eval_lt_rev(ret, matcl::trans(Matrix(A,false)).impl_unique<M1>(),q,p, B, trans_type::no_trans,
                                opts);
                    return;
                };
            }
            case trans_type::conj_trans:
            {
                if (m < n / 2)
                {
                    eval_ut(ret, A,p,q, matcl::ctrans(matcl::Matrix(B,false)).impl_unique<M2>(), 
                                              trans_type::no_trans, opts);
                    ret = matcl::ctrans(ret);
                    return;
                }
                else
                {
                    eval_lt_rev(ret, matcl::ctrans(Matrix(A,false)).impl_unique<M1>(),q,p, B, trans_type::no_trans,
                                opts);
                    return;
                };
            }
        };

        using V = typename unify_types<V1,V2>::type;
        using M = raw::Matrix<V,struct_sparse>;

        const raw::details::sparse_ccs<V1>& rep_A = A.rep();
        const raw::details::sparse_ccs<V2>& rep_B = B.rep();

        const Integer * A_c	= rep_A.ptr_c();
        const Integer * A_r	= rep_A.ptr_r();
        const V1 * A_x	    = rep_A.ptr_x();

        for (Integer j = n-1; j >= 0; --j)
        {
            Integer k = A_c[j + 1] - 1;

            // last non-zero element may be zero
            // find the element which is non-zero
            while(k > A_c[j] && A_x[k] == 0) 
                k--;

            if (k < A_c[j] || A_r[k] != j || A_x[k] == 0)
                throw error::error_singular();
        };                

        const Integer * B_c	= rep_B.ptr_c();
        const Integer * B_r	= rep_B.ptr_r();
        const V2 * B_x	    = rep_B.ptr_x();

        bool perm_p         = (p.is_id() == false);
        bool perm_q         = (q.is_id() == false);

        const Integer* q_ptr    = perm_q ? q.to_array() : nullptr;

        if (m == 1)
        {
            M C(A.get_type(), m, n, n);
            raw::details::sparse_ccs<V>& out = C.rep();

            Integer * out_c	    = out.ptr_c();
            Integer * out_r	    = out.ptr_r();
            V * out_x	        = out.ptr_x();

            Integer nz          = 0;

            for (Integer j = 0; j < n; ++j)
            {
                out_c[j]        = nz;
                Integer ku      = A_c[j + 1] - 1;
                Integer jb      = perm_q? q_ptr[j] - 1 : j;

                // avoid 0 elements as above
                while(A_x[ku] == 0)
                    ku--;

                V val           = A_x[ku--];
                V x             = 0;

                if (B_c[jb+1] > B_c[jb])
                    x           = B_x[B_c[jb]];

                //numerical phase
                for (Integer ka = ku; ka >= A_c[j]; --ka)
                {
                    Integer ja  = A_r[ka];
                    V elem      = A_x[ka];

                    //visit ja column of out matrix
                    for (Integer ko = out_c[ja]; ko < out_c[ja + 1]; ++ko)
                        x       -= out_x[ko] * elem;
                }

                if (x != 0)
                {
                    out_x[nz]   = x / val;
                    out_r[nz]   = 0;
                    ++nz;
                };
            }
            out_c[n] = nz;

            out.add_memory(-1);

            if (perm_q == true || perm_p == true)
            {
                C.get_struct().reset();
            }
            else
            {
                if (B.get_struct().is_triu())
                    C.get_struct().set(predefined_struct_type::triu);
                else
                    C.get_struct().reset();
            };

            ret = Matrix(C,true);

            if (perm_p == true)
                ret = ret(colon(), invperm(p).to_matrix());

            return;
        };

        M C(A.get_type(), m, n, m + n);
        raw::details::sparse_ccs<V>& out = C.rep();

        using workspace = matcl::details::workspace2<V>;
        using scatter   = matcl::algorithm::scatter;

        workspace           x(C.get_type(), m);
        
        scatter sc_w        = matcl::algorithm::scatter::get(m, n);

        Integer * out_c	    = out.ptr_c();
        Integer * out_r	    = out.ptr_r();
        V * out_x	        = out.ptr_x();

        Integer nz          = 0;

        for (Integer j = 0; j < n; ++j)
        {
            out_c[j]        = nz;            

            if (nz + m > out.nzmax())
            {
                out.add_memory(out.nzmax() + m);
                out_r	    = out.ptr_r();
                out_x	    = out.ptr_x();
            };

            Integer jb      = (perm_q ? q_ptr[j] - 1: j);
            Integer ku      = A_c[j + 1] - 1;

            // avoid 0 elements as above
            while(A_x[ku] == 0)
                ku--;

            V val           = A_x[ku--];

            //get symbolic structure
            auto mark       = sc_w.next_mark();
            Integer nz_old  = nz;

            for (Integer kb = B_c[jb]; kb < B_c[jb+1]; ++kb)
            {
                V tmp       = B_x[kb];

                if (tmp == 0)
                    continue;
                    
                Integer i   = B_r[kb];
                x[i]        = tmp;
                sc_w[i]     = mark;                    
                out_r[nz++] = i;
            }
        
            for (Integer ka = ku; ka >= A_c[j]; --ka)
            {
                Integer ja  = A_r[ka];

                if (A_x[ka] == 0)
                    continue;

                //visit ja column
                for (Integer ko = out_c[ja]; ko < out_c[ja + 1]; ++ko)
                {
                    Integer i = out_r[ko];
                    if (sc_w[i] < mark)
                    {
                        sc_w[i]     = mark;
                        x[i]        = 0;
                        out_r[nz++] = i;
                    };
                };
            }

            //sort row indices
            matcl::utils::sort_q(out_r + nz_old, nz - nz_old);

            //numerical phase
            for (Integer ka = ku; ka >= A_c[j]; --ka)
            {
                Integer ja  = A_r[ka];
                V elem      = A_x[ka];

                //visit ja column of out matrix
                for (Integer ko = out_c[ja]; ko < out_c[ja + 1]; ++ko)
                {
                    Integer i   = out_r[ko];
                    x[i]        -= out_x[ko] * elem;
                };
            }

            for (Integer k = out_c[j]; k < nz; ++k)
                out_x[k] = x[out_r[k]] / val;
        }

        out_c[n] = nz;

        out.add_memory(-1);

        if (perm_q == true || perm_p == true)
        {
            C.get_struct().reset();
        }
        else
        {
            if (B.get_struct().is_triu())
                C.get_struct().set(predefined_struct_type::triu);
            else
                C.get_struct().reset();
        };

        ret = Matrix(C,true);

        if (perm_p == true)
        {
            ret = ret(colon(), invperm(p).to_matrix());
        };

        return;
    };
};

//--------------------------------------------------------------------------------
//                  SPARSE-BAND
//--------------------------------------------------------------------------------
template<class V1, class V2>
struct linsolve_triang_sparse<struct_banded,V1,V2>
{
    using M1    = raw::Matrix<V1,struct_sparse>;
    using M2    = raw::Matrix<V2,struct_banded>;

    static void eval(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                     const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                     const options& opts)
    {
        Real density = Real(B.nnz())/(B.rows()+1.)/(B.cols()+1.);
     
        if (density < optim_params::max_sparse_density_min)
        {
            using SparseMatrix  = raw::Matrix<V2,struct_sparse>;
            SparseMatrix Bc     = raw::converter<SparseMatrix,M2>::eval(B);

            return linsolve_triang_str<struct_sparse,struct_sparse,V1,V2>::eval(ret, A,p,q,Bc,is_lt,trans, opts);
        }
        else
        {
            using DenseMatrix   = raw::Matrix<V2,struct_dense>;
            DenseMatrix Bc      = raw::converter<DenseMatrix,M2>::eval(B);

            return linsolve_triang_str<struct_sparse,struct_dense,V1,V2>::eval(ret, A,p,q,Bc,is_lt,trans, opts);
        };
    };

    static void eval_rev(matcl::Matrix& ret, const M1& A, const matcl::permvec& p, 
                         const matcl::permvec& q, const M2& B, bool is_lt, trans_type trans,
                         const options& opts)
    {
        Real density = Real(B.nnz())/(B.rows()+1.)/(B.cols()+1.);

        if (density < optim_params::max_sparse_density_min)
        {
            using SparseMatrix  = raw::Matrix<V2,struct_sparse>;
            SparseMatrix Bc     = raw::converter<SparseMatrix,M2>::eval(B);

            return linsolve_triang_str<struct_sparse,struct_sparse,V1,V2>::eval_rev(ret, A,p,q,Bc,is_lt,trans, opts);
        }
        else
        {
            using DenseMatrix   = raw::Matrix<V2,struct_dense>;
            DenseMatrix Bc      = raw::converter<DenseMatrix,M2>::eval(B);

            return linsolve_triang_str<struct_sparse,struct_dense,V1,V2>::eval_rev(ret,A,p,q,Bc,is_lt,trans,opts);
        };
    };
};

};};
