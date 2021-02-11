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

#include "matcl-internals/algs/scatter.h"
#include "matcl-matrep/lib_functions/func_binary.h"
#include "matcl-core/utils/workspace.h"
#include "matcl-internals/base/sort.h"
#include "matcl-linalg/utils/optim_params.h"

namespace matcl { namespace details
{

//--------------------------------------------------------------------------------
//                                  HELPERS
//--------------------------------------------------------------------------------

template<class V>
struct update_helpers
{
    using VR = typename md::real_type<V>::type;

    static bool scaling_chudl(const VR& SIGMA, const VR& dj, const V& pj, 
                                    VR& djbar, VR& C, V& S0, V& S1);
};

template<class V>
bool update_helpers<V>::scaling_chudl(const VR& SIGMA, const VR& dj, const V& pj, 
                                    VR& djbar, VR& C, V& S0, V& S1)
{
    if (pj == V(0.))
    {
        djbar   = dj;
        S0      = V(0.0);
        S1      = V(0.0);
        C       = 1.0;

        return true;
    };

	VR scale    = dj + sqrt(abs(SIGMA)) * abs(pj);
    djbar   	= abs2(dj/scale) + SIGMA * abs2(pj/scale);		
    
    if (djbar <= 0.)
        return false;
    
    djbar  		= scale * sqrt(djbar);
    S0          = pj / djbar;
    S1          = SIGMA * S0;
    C           = dj / djbar;

    return true;
};

template<class V, class S>
struct chol_remove_diag
{
    using Mat           = raw::Matrix<V,S>;
    using DenseMatrix   = raw::Matrix<V,struct_dense>;

    static void eval(matcl::Matrix& ret, const Mat& A, Integer k)
    {
        matcl::Matrix Ad = Matrix(A.get_diag(), false);
        Ad               = delrows(std::move(Ad), k);

        ret = matcl::bdiag(Ad, 0);
    }
    static void eval(matcl::Matrix& ret, const Mat& A, const matcl::Matrix& K)
    {
        matcl::Matrix Ad = Matrix(A.get_diag(), false);
        Ad               = delrows(std::move(Ad), K);

        ret = matcl::bdiag(Ad, 0);
    }
};

template<class V, class S> 
struct chol_update_impl{};

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V> 
struct chol_update_impl<V,struct_dense> 
{
    using Mat = raw::Matrix<V,struct_dense>;

    static void eval_remove(matcl::Matrix& ret, const Mat& A, bool upper, Integer k)
    {
        if (k == A.rows())
        {
            const Mat& tmp = A.resize(A.rows()-1,A.cols()-1);
            ret = matcl::Matrix(tmp,true);
            return;
        };

        if (A.get_struct().is_id())
        {
            const auto& tmp = A.resize(A.rows()-1,A.cols()-1);
            ret = matcl::Matrix(tmp,true);
            return;
        };

        if (is_diag(Matrix(A,false)))
            return chol_remove_diag<V,struct_dense>::eval(ret, A, k);

        Mat A2      = A.make_unique();
        A2.get_struct().reset_value();

        const Mat& Av   = A2.make_view(k+1, A.rows(), k+1, A.cols());
        Matrix v        = upper ? Matrix(A2.make_view(k,k,k+1,A.cols()),false)
                                : Matrix(A2.make_view(k+1,A.rows(),k,k),false);

        {
            Matrix Au;
            eval_update(Au, Av, v, upper, 1., true);
        };

        Matrix Ar   = Matrix(std::move(A2),false);

        ret = delrowscols(std::move(Ar), k,k);
        return;
    };

    static void eval_remove(matcl::Matrix& ret, const Mat& A, bool upper, const matcl::Matrix& K)
    {
        const Integer* K_ptr    = K.get_array<Integer>();
        Integer N               = K.length();
        Integer K0              = K_ptr[0];

        if (K0 == A.cols() - N + 1)
        {
            const Mat& tmp = A.resize(A.rows()-K.length(),A.cols()-K.length());
            ret     = matcl::Matrix(tmp,true);
            return;
        };

        if (A.get_struct().is_id())
        {
            const auto& tmp = A.resize(A.rows()-K.length(),A.cols()-K.length());
            ret     = matcl::Matrix(tmp,true);
            return;
        };
        if (is_diag(Matrix(A,false)) )
            return chol_remove_diag<V,struct_dense>::eval(ret, A, K);

        Mat A2 = A.make_unique();
        A2.get_struct().reset_value();

        for (Integer i = 0; i < N; ++i)
        {
            Integer k       = K_ptr[i];
            const Mat& Av   = A2.make_view(k+1, A.rows(), k+1, A.cols());
            Matrix v    = upper ? Matrix(A2.make_view(k,k,k+1,A.cols()),false)
                                : Matrix(A2.make_view(k+1,A.rows(),k,k),false);

            //update is inplace => A2 is changed
            eval_update(ret, Av, v, upper, 1., true);
        };
        Matrix Au       = delrowscols(matcl::Matrix(A2,false), K, K);
        ret = Au;
        return;
    };

    static void eval_update(matcl::Matrix& ret, const Mat& A, const matcl::Matrix& w0, bool upper, 
                            const Matrix& Sigma, bool inplace)
    {
        using VR        = typename md::real_type<V>::type;

        Integer N       = A.rows();

        //no rvalue overload for w matrix; always make copy
        Matrix w        = convert(w0, Mat::matrix_code);
        Mat& v          = w.get_impl_unique<Mat>();

        V* v_ptr        = v.ptr();
        Integer v_LD    = v.ld();
        Integer K       = v.cols();

        Mat work(A.get_type(), upper ? 2*N : 0, 1);
        V* work_ptr     = work.ptr();
        
        bool inplace2   = A.get_refstr()->is_unique() == true || inplace == true;
        Mat A2          = (inplace2)? Mat(A, Mat::copy_is_safe()) : A.clone();

        VR sig_scal;
        const VR* sig_ptr;
        Integer sig_step;

        if (Sigma.length() == 1)
        {
            sig_scal    = Sigma.get_scalar<VR>();
            sig_ptr     = &sig_scal;
            sig_step    = 0;
        }
        else
        {
            sig_ptr     = Sigma.get_array<VR>();
            sig_step    = 1;
        };

        Integer info;

        V* U_ptr        = A2.ptr();
        Integer LDA     = A2.ld();

        if (upper == true)
        {
            for (Integer i = 0; i < K; ++i)
            {
                lapack::chudu(lap(U_ptr), LDA, N, sig_ptr[0], lap(v_ptr), 1, lap(work_ptr), 
                              details::lap(&info));
                sig_ptr += sig_step;
                v_ptr   += v_LD;

                if (info > 0)
                    throw error::error_invalid_chol_update_nonposdef();
            }
        }
        else
        {
            for (Integer i = 0; i < K; ++i)
            {
                lapack::chudl(lap(U_ptr), LDA, N, sig_ptr[0], lap(v_ptr), 1, details::lap(&info));
                sig_ptr += sig_step;
                v_ptr   += v_LD;

                if (info > 0)
                    throw error::error_invalid_chol_update_nonposdef();
            };
        }

        ret = matcl::Matrix(A2,true);
        return;
    };
};

//--------------------------------------------------------------------------------
//                                  SPARSE
//--------------------------------------------------------------------------------
class rem_heap
{
    private:
        struct item
        {
            Integer         length;
            const Integer*  root;
            size_t          offset;
            bool            const_root;

            item(Integer l, const Integer* r, size_t off, bool is_const)
                :length(l), root(r), offset(off), const_root(is_const)
            {};

            void rebind(const Integer* main_root)
            {
                if (const_root == false)
                    root = main_root + offset;
            };
        };

        using items         = std::vector<item>;
        using map_vec       = std::map<Integer, items>;

    private:
        map_vec             m_vectors;
        const Integer*      m_root_ptr;

    public:
        rem_heap(const Integer* root_ptr);

        void                build_remove(const Matrix& K, const Integer* ptr_c, const Integer* ptr_r);
        template<class VR>
        void                build_update(const VR* sigma_ptr, Integer sigma_step, const Integer* W_c, 
                                    const Integer* W_r, Integer KN);

        void                add_vector(Integer first_pos, const Integer* ptr, Integer len, bool is_const_root);

        Integer             top_index() const;
        Integer             top_length() const;
        void                next_index();

        const Integer*      get_vector(Integer l) const;
        Integer             get_vector_length(Integer l) const;

        Integer             nnz_current_group() const;
        Integer             nnz_all_groups() const;
        
        void                rebind(const Integer* root_r_ptr);

    private:
        Integer             calc_nnz_group(const items& vec) const;
        void                rebind_items(items& vec, const Integer* root_r_ptr);
};

rem_heap::rem_heap(const Integer* root_ptr)
    :m_root_ptr(root_ptr)
{};

void rem_heap::build_remove(const Matrix& K, const Integer* ptr_c, const Integer* ptr_r)
{
    Integer N                   = K.length();
    const Integer* k_ptr        = K.get_array<Integer>();

    for (Integer i = 0; i < N; ++i)
    {
        //indices are 1-based
        Integer k               = k_ptr[i] - 1;
        Integer start_row       = k + 1;

        Integer rs_length       = ptr_c[k+1] - ptr_c[k];
        Integer off             = ptr_c[k];
        const Integer* rs_root  = ptr_r + off;

        add_vector(start_row, rs_root, rs_length, true);
    };
};

template<class VR>
void rem_heap::build_update(const VR* ptr_sig, Integer sig_step, const Integer* W_c, 
                            const Integer* W_r, Integer K)
{
    for (Integer i = 0; i < K; ++i)
    {
        Integer k1              = W_c[i];
        Integer k2              = W_c[i+1];

        Integer rs_length       = k2 - k1;

        if (rs_length == 0 || ptr_sig[i*sig_step] == VR(0.0))
            continue;

        Integer off             = k1;
        const Integer* rs_root  = W_r + off;
        Integer start_row       = rs_root[0];

        add_vector(start_row, rs_root, rs_length, true);
    };
};

void rem_heap::add_vector(Integer start_row, const Integer* rs_root, Integer rs_length, 
                          bool is_const_root)
{
    while (rs_length > 0 && (rs_root[0] < start_row) )
    {
        ++rs_root;
        --rs_length;
    };

    if (rs_length == 0)
        return;

    Integer k_first = rs_root[0];
    size_t offset   = is_const_root? 0 : rs_root - m_root_ptr;

    auto pos = m_vectors.find(k_first);

    if (pos == m_vectors.end())
        m_vectors.insert(pos, map_vec::value_type(k_first, 
                            items{item{rs_length, rs_root, offset,is_const_root}}));
    else
        pos->second.push_back(item{rs_length, rs_root, offset,is_const_root});
};

Integer rem_heap::top_index() const
{
    return m_vectors.begin()->first;
};
Integer rem_heap::top_length() const
{
    if (m_vectors.empty() == true)
        return 0;

    return (int)m_vectors.begin()->second.size();
};

void rem_heap::next_index()
{
    if (m_vectors.size() > 0)
        m_vectors.erase(m_vectors.begin());
};

const Integer* rem_heap::get_vector(Integer l) const
{
    return m_vectors.begin()->second[l].root;
};
Integer rem_heap::get_vector_length(Integer l) const
{
    return m_vectors.begin()->second[l].length;
};

Integer rem_heap::nnz_current_group() const
{
    return calc_nnz_group(m_vectors.begin()->second);
};
Integer rem_heap::nnz_all_groups() const
{
    Integer nz = 0;
    for (const auto& pos : m_vectors)
        nz += calc_nnz_group(pos.second);

    return nz;
};
Integer rem_heap::calc_nnz_group(const items& vec) const
{
    Integer nz = 0;
    for (const auto& pos : vec)
        nz += pos.length;

    return nz;
};

void rem_heap::rebind(const Integer* root_r_ptr)
{
    m_root_ptr = root_r_ptr;

    for (auto& pos : m_vectors)
        rebind_items(pos.second, root_r_ptr);
};

void rem_heap::rebind_items(items& vec, const Integer* root_r_ptr)
{
    for (auto& pos : vec)
        pos.rebind(root_r_ptr);
};

template<class V> 
struct chol_update_impl<V,struct_sparse> 
{
    using Mat       = raw::Matrix<V,struct_sparse>;
    using VP        = pod_type<V>;
    using VR        = typename md::real_type<V>::type;
    using workspace = pod_workspace<VP>;

    static void eval_update(matcl::Matrix& ret, const Mat& A, const matcl::Matrix& w, 
                            bool upper, const Matrix& sigma0, bool inplace)
    {
        using VR    = typename md::real_type<V>::type;

        Integer K   = w.cols();
       
        //special cases
        if (A.rows() == 0 || K == 0)
        {
            ret = matcl::Matrix(A,true);
            return;
        };
         
        if (A.nnz() == 0 && K == 1)
        {
            //simply create new cholesky decomposition
            VR sigma    = sigma0.get_scalar<VR>();

            if (sigma < 0)
                throw error::error_invalid_chol_update_nonposdef();

            if (upper == true)
            {
                ret = sqrt(sigma) * ctrans(vec(w));
                ret.resize(A.rows(), A.cols());
                ret.add_struct(predefined_struct_type::triu);
                return;
            }
            else
            {
                ret = sqrt(sigma) * vec(w);
                ret.resize(A.rows(), A.cols());
                ret.add_struct(predefined_struct_type::tril);
                ret;
            };
        };

        if (upper == true)
            return update_ut(ret, A,w,sigma0,inplace);
        else
            return update_block_lt(ret, A,w,sigma0,inplace);
    };

    static void eval_remove(matcl::Matrix& ret, const Mat& A, bool upper, const matcl::Matrix& K)
    {
        const Integer* K_ptr    = K.get_array<Integer>();
        Integer N               = K.length();
        Integer K0              = K_ptr[0];

        if (K0 == A.cols() - N + 1)
        {
            const Mat& tmp  = A.resize(A.rows()-K.length(),A.cols()-K.length());
            ret             = matcl::Matrix(tmp,true);
            return;
        };

        if (A.get_struct().is_id())
        {
            const Mat& tmp  = A.resize(A.rows()-1,A.cols()-1);
            ret             = matcl::Matrix(tmp,true);
            return;
        };

        if (is_diag(Matrix(A,false)))
            return chol_remove_diag<V,struct_sparse>::eval(ret, A, K);

        if (A.nnz() == 0)
        {
            const Mat& tmp     = A.resize(A.rows()-1,A.cols()-1);
            ret         = matcl::Matrix(tmp,true);
            return;
        };

        if (upper == true)
        {
            //TODO: implement for sparse matrix;
            using MC    = raw::Matrix<V,struct_dense>;
            MC AC       = raw::converter<MC,Mat>::eval(A);
            return chol_update_impl<V,struct_dense>::eval_remove(ret, AC, upper, K);
        };

        bool inplace = false;
        if (A.get_refstr()->is_unique() == true)
            inplace = true;

        Matrix Au;
        eval_remove_block_lt(Au, A, K, inplace);
        ret = delrowscols(std::move(Au),K,K);
        return;
    };

    static void eval_remove_block_lt(matcl::Matrix& ret, const Mat& A, const matcl::Matrix& K, bool inplace)
    {
        const mrd::sparse_ccs<V>& A_rep = A.rep();
        const Integer* A_c              = A_rep.ptr_c();
        const Integer* A_r              = A_rep.ptr_r();

        rem_heap m_rem(A_r);
        m_rem.build_remove(K, A_c, A_r);

        Mat A2 = symbolic_analysis(A, m_rem, inplace);
        A2.get_struct().reset();

        Integer KN              = K.length();
        const Integer* K_ptr    = K.get_array<Integer>();

        for (Integer i = 0; i < KN; ++i)
            update_numeric_remove(A2, K_ptr[i] - 1);

        A2.get_struct().set(predefined_struct_type::tril);
        ret = Matrix(A2,true);
        return;
    };

    static void update_ut(matcl::Matrix& ret, const Mat& A, const matcl::Matrix& w, 
                          const Matrix& sigma, bool)
    {
        //TODO: implement for sparse matrix;
        using MC    = raw::Matrix<V,struct_dense>;
        MC AC       = raw::converter<MC,Mat>::eval(A);
        return chol_update_impl<V,struct_dense>::eval_update(ret, AC, w, true, sigma, true);        
    };

    static void update_block_lt(matcl::Matrix& ret, const Mat& A, const matcl::Matrix& w, 
                          const Matrix& sigma, bool inplace)
    {
        using VR = typename md::real_type<V>::type;

        //at this stage value code of A is unified with value code of w
        Matrix nz_W         = nnz(w, 1);
        Integer max_nnz_row = max_vec(nz_W).get_scalar<Integer>();
        Integer N           = A.rows();
        Real density        = Real(max_nnz_row)/Real(N);

        if (density > linalg_optim_params::max_row_density_chol_update())
        {
            //W is too dense for sparse update
            using Mat_D     = raw::Matrix<V,struct_dense>;
            Mat_D Ac        = raw::converter<Mat_D,Mat>::eval(A);

            bool upper      = false;
            return chol_update_impl<V,struct_dense>::eval_update(ret, Ac, w, upper, sigma, true);
        };

        const Mat& W                    = convert(w, Mat::matrix_code).get_impl<Mat>();
        const mrd::sparse_ccs<V>& W_rep = W.rep();
        const mrd::sparse_ccs<V>& A_rep = A.rep();
        const Integer* W_c              = W_rep.ptr_c();
        const Integer* W_r              = W_rep.ptr_r();
        Integer KN                      = w.cols();

        VR sig_scal;
        const VR* sig_ptr;
        Integer sig_step;

        if (sigma.length() == 1)
        {
            sig_scal    = sigma.get_scalar<VR>();
            sig_ptr     = &sig_scal;
            sig_step    = 0;
        }
        else
        {
            sig_ptr     = sigma.get_array<VR>();
            sig_step    = 1;
        };

        rem_heap m_rem(A_rep.ptr_r());
        m_rem.build_update<VR>(sig_ptr, sig_step, W_c, W_r, KN);

        Mat A2      = symbolic_analysis(A, m_rem, inplace);        
        A2.get_struct().reset();

        for (Integer i = 0; i < KN; ++i)
            update_numeric_prepare(A2, i, W, sig_ptr[i*sig_step]);

        A2.get_struct().set(predefined_struct_type::tril);
        ret = Matrix(A2,true);
        return;
    };

    static Mat symbolic_analysis(const Mat& A, rem_heap& m_rem, bool inplace)
    {
        const mrd::sparse_ccs<V>& A_rep = A.rep();
        const Integer* A_c              = A_rep.ptr_c();
        const Integer* A_r              = A_rep.ptr_r();
        const V* A_x                    = A_rep.ptr_x();

        using scatter = matcl::algorithm::scatter;

        scatter sc = scatter::get(A.rows(), A.cols());
             
        bool fillins                = false;

        if (m_rem.top_length() == 0)
        {
            if (inplace == true || A.get_refstr()->is_unique() == true)
                return Mat(A, Mat::copy_is_safe());
            else 
                return A.copy();
        };

        Integer j                   = m_rem.top_index();

        for(;;)
        {
            Integer pos_A           = A_c[j];
            Integer pos_Al          = A_c[j+1];
            Integer r               = A_r[pos_A];                        

            while (r < j && pos_A < pos_Al)
                r                   = A_r[++pos_A];

            auto mark               = sc.next_mark();

            for (Integer k = pos_A; k < pos_Al; ++k)
            {
                Integer r2          = A_r[k];
                sc[r2]              = mark;
            };   

            //test for fill-ins
            Integer rem_length      = m_rem.top_length();
            for (Integer l = 0; l < rem_length; ++l)
            {
                const Integer* rs   = m_rem.get_vector(l);
                Integer rs_length   = m_rem.get_vector_length(l);

                for (Integer k = 0; k < rs_length; ++k)
                {
                    Integer r2      = rs[k];
                    if (sc[r2] < mark)
                    {
                        fillins         = true;
                        goto lab_fill;
                    };
                };
            };

            //store symbolic structure
            m_rem.add_vector(j + 1, A_r + pos_A + 1, pos_Al - pos_A - 1, false);

            //goto next column
            m_rem.next_index();

            if (m_rem.top_length() == 0)
            {
                if (inplace == true || A.get_refstr()->is_unique() == true)
                    return Mat(A, Mat::copy_is_safe());
                else 
                    return A.copy();
            };
            
            j                       = m_rem.top_index();
        };

      lab_fill:

        //create new matrix, copy, and continue symbolic analysis
        Mat out(A.get_type(), A.rows(), A.cols(), A.nnz() + m_rem.nnz_all_groups());

        mrd::sparse_ccs<V> rep  = out.rep();

        Integer* out_c          = rep.ptr_c();
        Integer* out_r          = rep.ptr_r();
        V* out_x                = rep.ptr_x();

        m_rem.rebind(out_r);

        Integer nz              = 0;
        Integer j_out           = 0;

        for(;;)
        {
            Integer dnz         = A_c[j] - A_c[j_out]
                                + A_c[j + 1] - A_c[j]
                                + m_rem.nnz_current_group();

            auto mark           = sc.next_mark();

            if (nz + dnz > rep.nzmax()) 
            {
                rep.add_memory(rep.nzmax() + dnz);

                out_r           = rep.ptr_r();
                out_x           = rep.ptr_x();

                m_rem.rebind(out_r);
            };

            //copy unmodified elements to out
            for (; j_out < j; ++ j_out)
            {
                out_c[j_out]    = nz;
                for (Integer k = A_c[j_out]; k < A_c[j_out+1]; ++k)
                {
                    out_r[nz]   = A_r[k];
                    out_x[nz++] = A_x[k];
                };
            };

            out_c[j_out]        = nz;
            Integer nz_old      = nz;

            Integer pos_A       = A_c[j];
            Integer pos_Al      = A_c[j+1];
            Integer r           = A_r[pos_A];
            
            while (r < j && pos_A < pos_Al)
                r               = A_r[++pos_A];

            for (Integer k = pos_A; k < pos_Al; ++k)
            {
                Integer r2      = A_r[k];
                sc[r2]          = mark;

                out_r[nz]       = r2;
                out_x[nz++]     = A_x[k];
            };

            bool inserted           = false;
            Integer rem_length      = m_rem.top_length();
            for (Integer l = 0; l < rem_length; ++l)
            {
                const Integer* rs   = m_rem.get_vector(l);
                Integer rs_length   = m_rem.get_vector_length(l);

                for (Integer k = 0; k < rs_length; ++k)
                {
                    Integer r2      = rs[k];
                    if (sc[r2] < mark)
                    {
                        sc[r2]      = mark;

                        out_r[nz]   = r2;
                        out_x[nz++] = 0;

                        inserted    = true;
                    };
                };
            };

            //store results
            if (inserted == true)
                matcl::utils::sort_q(out_r + nz_old, out_x + nz_old,  nz - nz_old);

            //store symbolic structure
            m_rem.add_vector(j + 1, out_r + nz_old + 1, nz - nz_old - 1, false);

            //goto next column
            m_rem.next_index();
            ++ j_out;

            if (m_rem.top_length() == 0)
                break;
            
            j                       = m_rem.top_index();
        };

        //copy unmodified elements to out
        Integer dnz             = A_c[A.cols()] - A_c[j_out];
        if (nz + dnz > rep.nzmax()) 
        {
            rep.add_memory(rep.nzmax() + dnz);

            out_r               = rep.ptr_r();
            out_x               = rep.ptr_x();
        };
        for (; j_out < A.cols(); ++ j_out)
        {
            out_c[j_out]        = nz;
            for (Integer k = A_c[j_out]; k < A_c[j_out+1]; ++k)
            {
                out_r[nz]       = A_r[k];
                out_x[nz++]     = A_x[k];
            };
        };
        out_c[A.cols()]         = nz;

        return out;
    };

    static void update_numeric_remove(const Mat& A, Integer k)
    {
        const mrd::sparse_ccs<V>& rep   = A.rep();
        const Integer* A_c              = rep.ptr_c();
        const Integer* A_r              = rep.ptr_r();
        const V* A_x                    = rep.ptr_x();

        Integer pos                 = A_c[k];

        while ((A_r[pos] <= k || A_x[pos] == 0) && pos < A_c[k+1])
            ++pos;

        Integer length              = A_c[k+1] - pos;
        const Integer* r_ptr        = A_r + pos;
        const V* x_ptr              = A_x + pos;
        Integer start_col           = A_r[pos];

        return update_numeric(A, start_col, length, r_ptr, x_ptr, VR(1.0));
    }

    static void update_numeric_prepare(const Mat& A, Integer k, const Mat& W, const VR& sigma)
    {
        const mrd::sparse_ccs<V>& rep   = W.rep();
        const Integer* W_c              = rep.ptr_c();
        const Integer* W_r              = rep.ptr_r();
        const V* W_x                    = rep.ptr_x();

        Integer pos                 = W_c[k];

        while ((W_x[pos] == 0) && pos < W_c[k+1])
            ++pos;

        Integer length              = W_c[k+1] - pos;
        const Integer* r_ptr        = W_r + pos;
        const V* x_ptr              = W_x + pos;
        Integer start_col           = W_r[pos];

        return update_numeric(A, start_col, length, r_ptr, x_ptr, sigma);
    }

    static void update_numeric(const Mat& A, Integer start_col, Integer length, const Integer* r_ptr, 
                                const V* x_ptr, const VR& sigma)
    {
        if (length == 0 || sigma == VR(0.0))
            return;

        //update_numeric is inplace algorithm
        mrd::sparse_ccs<V> rep      = A.rep();
        Integer* A_c                = rep.ptr_c();
        Integer* A_r                = rep.ptr_r();
        V* A_x                      = rep.ptr_x();

        using VR        = typename md::real_type<V>::type;
        using scatter   = matcl::algorithm::scatter;

        scatter sc      = scatter::get(A.rows(), 1);

        workspace work_w(A.rows());

        auto mark                   = sc.next_mark();
        V* work_ptr                 = reinterpret_cast<V*>(work_w.ptr());

        //copy w to workspace      
        for (Integer i = 0; i < length; ++i)
        {
            Integer r               = r_ptr[i]; 
            work_ptr[r]             = x_ptr[i];
            sc[r]                   = mark;
        };
       
        Integer j                   = start_col;

        for(;;)
        {
            Integer pos_A           = A_c[j];
            Integer pos_Al          = A_c[j+1];
            Integer r               = A_r[pos_A];
            
            while (r < j && pos_A < pos_Al)
                r                   = A_r[++pos_A];

            // diagonal element must always exist
            assert(r == j);

            V pj                    = work_ptr[j];            
            VR dj                   = real(A_x[pos_A]);

            VR djbar, C;
            V S0, S1;

            bool info               = update_helpers<V>::scaling_chudl(sigma, dj, pj, djbar, C, S0, S1);

            if (info == false)
                throw error::error_invalid_chol_update_nonposdef();

            A_x[pos_A]              = djbar;

            for (Integer k2 = pos_A + 1; k2 < pos_Al; ++k2)
            {
                Integer r2          = A_r[k2];
                V Ax                = A_x[k2];
                V tmp, Y;

                if (sc[r2] == mark)
                {
                    Y               = work_ptr[r2];
                    tmp             = C * work_ptr[r2] - S0 * Ax;
                }
                else
                {
                    Y               = V(0.0);
                    tmp             = - S0 * Ax;
                    sc[r2]          = mark;
                };
                
                A_x[k2]             = C * Ax + S1 * Y;                
                work_ptr[r2]        = tmp;
            };

            //goto next column
            if (pos_A + 1 < pos_Al)
                j                   = A_r[pos_A + 1];
            else
                j                   = -1;

            if (j == -1)
                break;
        };
    };
};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class V> 
struct chol_update_impl<V,struct_banded> 
{
    using Mat = raw::Matrix<V,struct_banded>;

    static void eval_remove(matcl::Matrix& ret, const Mat& A, bool upper, Integer k)
    {
        if (k == A.rows())
        {
            const auto& tmp = A.resize(A.rows()-1,A.cols()-1);
            ret = matcl::Matrix(tmp,true);
            return;
        };

        if (A.get_struct().is_id())
        {
            const auto& tmp = A.resize(A.rows()-1,A.cols()-1);
            ret = matcl::Matrix(tmp,true);
            return;
        };

        if (A.get_struct().is_diag() || is_diag(Matrix(A,false)))
        {
            return chol_remove_diag<V,struct_banded>::eval(ret, A, k);
        };

        Real density = Real(A.nnz())/(A.rows()+1.)/(A.cols()+1.);

        if (density < optim_params::max_sparse_density_min)
        {
            using SparseMatrix = raw::Matrix<V,struct_sparse>;
            SparseMatrix smat = raw::converter<SparseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
            return chol_update<V,struct_sparse>::eval_remove(ret, smat, upper, k);
        }
        else
        {
            using DenseMatrix = raw::Matrix<V,struct_dense>;
            DenseMatrix smat = raw::converter<DenseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
            return chol_update<V,struct_dense>::eval_remove(ret, smat, upper, k);
        };
    };

    static void eval_remove(matcl::Matrix& ret, const Mat& A, bool upper, const matcl::Matrix& K)
    {
        const Integer* K_ptr    = K.get_array<Integer>();
        Integer N               = K.length();
        Integer K0              = K_ptr[0];

        if (K0 == A.cols() - N + 1)
        {
            const Mat& tmp = A.resize(A.rows()-K.length(),A.cols()-K.length());
            ret = matcl::Matrix(tmp,true);
            return;
        };

        if (A.get_struct().is_id())
        {
            const auto& tmp = A.resize(A.rows()-K.length(),A.cols()-K.length());
            ret = matcl::Matrix(tmp,true);
            return;
        };

        if (A.get_struct().is_diag() || is_diag(Matrix(A,false)))
        {
            return chol_remove_diag<V,struct_banded>::eval(ret, A, K);
        };

        Real density = Real(A.nnz())/(A.rows()+1.)/(A.cols()+1.);
        if (density < optim_params::max_sparse_density_min)
        {
            using SparseMatrix = raw::Matrix<V,struct_sparse>;
            SparseMatrix smat = raw::converter<SparseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
            return chol_update<V,struct_sparse>::eval_remove(ret, smat, upper, K);
        }
        else
        {
            using DenseMatrix = raw::Matrix<V,struct_dense>;
            DenseMatrix smat = raw::converter<DenseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
            return chol_update<V,struct_dense>::eval_remove(ret, smat, upper, K);
        };
    };

    static void eval_update(matcl::Matrix& ret, const Mat& A, const matcl::Matrix& w, bool upper,
                            const Matrix& sigma, bool)
    {
        Real density = Real(A.nnz())/(A.rows()+1.)/(A.cols()+1.);

        if (density < optim_params::max_sparse_density_min)
        {
            using SparseMatrix = raw::Matrix<V,struct_sparse>;
            SparseMatrix smat = raw::converter<SparseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
            return chol_update_impl<V,struct_sparse>::eval_update(ret, smat,w,upper,sigma,true);
        }
        else
        {
            using DenseMatrix = raw::Matrix<V,struct_dense>;
            DenseMatrix smat = raw::converter<DenseMatrix,raw::Matrix<V,struct_banded>>::eval(A);
            return chol_update_impl<V,struct_dense>::eval_update(ret, smat,w,upper,sigma,true);
        };
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct chol_update
{
    using M = raw::Matrix<V,S>;

    static void eval_remove(matcl::Matrix& ret, const M& A, bool upper, Integer k)
    {
        chol_update_impl<V,S>::eval_remove(ret, A, upper, k);
        ret.add_struct(upper ? predefined_struct_type::triu : predefined_struct_type::tril);
        return;
    };

    static void eval_remove(matcl::Matrix& ret, const M& A, bool upper, const matcl::Matrix& K)
    {
        chol_update_impl<V,S>::eval_remove(ret, A, upper, K);
        ret.add_struct(upper ? predefined_struct_type::triu : predefined_struct_type::tril);
        return;
    };

    static void eval_update(matcl::Matrix& ret, const M& A, const Matrix& w, bool upper, 
                            const Matrix& sigma)
    {
        chol_update_impl<V,S>::eval_update(ret, A, w, upper, sigma, false);
        ret.add_struct(upper ? predefined_struct_type::triu : predefined_struct_type::tril);
        return;
    };
};

template<class S>
struct chol_update<Integer,S>
{
    using M = raw::Matrix<Integer,S>;

    static void eval_remove(matcl::Matrix& ret, const M& A, bool upper, Integer k)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return chol_update_impl<Real,S>::eval_remove(ret, AC, upper, k);
    };

    static void eval_remove(matcl::Matrix& ret, const M& A, bool upper, const matcl::Matrix& K)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return chol_update_impl<Real,S>::eval_remove(ret, AC, upper, K);
    };

    static void eval_update(matcl::Matrix& ret, const M& A, const Matrix& w, bool upper, 
                            const Matrix& sigma)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return chol_update_impl<Real,S>::eval_update(ret, AC, w, upper, sigma, true);
    };
};

};};
