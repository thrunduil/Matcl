/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/func/raw/raw_manip.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-internals/base/utils.h"
#include "matcl-core/details/integer.h"
#include "matcl-matrep/utils/workspace.h"
#include "matcl-scalar/details/scalfunc_helpers.h"
#include "matcl-matrep/details/struct_flag_predefined.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/base/optim_params.h"
#include "matcl-matrep/lib_functions/func_unary.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace raw { namespace details
{

template<class V, class S>
struct all_finite_helper_impl
{};

template<class V>
struct all_finite_helper_impl<V,struct_dense>
{
    using Mat   = raw::Matrix<V,struct_dense>;

    static bool eval(const Mat& A)
    {
        Integer M       = A.rows();
        Integer N       = A.cols();
        Integer A_ld    = A.ld();

        const V* ptr    = A.ptr();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < M; ++i)
            {
                if ((bool)mrd::isfinite_helper<V>::eval(ptr[i]) == false)
                    return false;
            };

            ptr         += A_ld;
        };

        return true;
    };
};

template<class V>
struct all_finite_helper_impl<V,struct_sparse>
{
    using Mat   = raw::Matrix<V,struct_sparse>;

    static bool eval(const Mat& A)
    {
        if (A.nnz() == 0)
            return true;

        const details::sparse_ccs<V>& d = A.rep();

        const Integer* d_c  = d.ptr_c();
        const V* d_x		= d.ptr_x();

        Integer N           = d.cols();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
            {
                if ((bool)mrd::isfinite_helper<V>::eval(d_x[k]) == false)
                    return false;
            };
        };

        return true;
    };
};

template<class V>
struct all_finite_helper_impl<V,struct_banded>
{
    using Mat   = raw::Matrix<V,struct_banded>;

    static bool eval(const Mat& A)
    {
        Integer N       = A.cols();        
        Integer A_ld    = A.ld();
        Integer fd      = A.first_diag();

        if (fd == A.last_diag())
        {
            const V* ptr    = A.rep_ptr() + A.first_elem_diag(fd);
            Integer s       = A.diag_length(fd);

            for (Integer k = 0; k < s; ++k)
            {
                if ((bool)mrd::isfinite_helper<V>::eval(ptr[0]) == false)
                    return false;

                ptr         += A_ld;
            };

            return true;
        };

        const V* ptr    = A.rep_ptr();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer fp  = A.first_elem_pos(j);

            for (Integer k = fr; k <= lr; ++k, ++fp)
            {
                if ((bool)mrd::isfinite_helper<V>::eval(ptr[fp]) == false)
                    return false;
            };

            ptr         += A_ld;
        };

        return true;
    };
};

}}};

namespace matcl { namespace raw
{

namespace gr = matcl::raw;
namespace mrd = matcl::raw::details;
namespace md = matcl::details;

template<class val_type>
Matrix<val_type,struct_dense> full(const Matrix<val_type,struct_sparse>& A)
{
    using SparseMat = Matrix<val_type,struct_sparse>;
    using FullMat   = Matrix<val_type,struct_dense>;
    return converter<FullMat,SparseMat>::eval(A);
};

template<class val_type>
Matrix<val_type,struct_sparse> sparse(const Matrix<val_type,struct_dense>& A)
{
    using SparseMat = Matrix<val_type,struct_sparse>;
    using FullMat   = Matrix<val_type,struct_dense>;
    return converter<SparseMat,FullMat>::eval(A);
};

template<class val_type>
Matrix<val_type,struct_sparse> sparse(const Matrix<val_type,struct_banded>& A)
{
    using SparseMat = Matrix<val_type,struct_sparse>;
    using BandMat   = Matrix<val_type,struct_banded>;
    return converter<SparseMat,BandMat>::eval(A);
};

template<class val_type>
Matrix<val_type,struct_banded> band(const Matrix<val_type,struct_dense>& A)
{
    using BandMat   = Matrix<val_type,struct_banded>;
    using FullMat   = Matrix<val_type,struct_dense>;
    return converter<BandMat,FullMat>::eval(A);
};

template<class val_type>
Matrix<val_type,struct_banded> band(const Matrix<val_type,struct_sparse>& A)
{
    using BandMat   = Matrix<val_type,struct_banded>;
    using SparseMat = Matrix<val_type,struct_sparse>;	
    return converter<BandMat,SparseMat>::eval(A);
};

template<class ret_type,class M,class struct_type>
struct eval_reshape_helper
{};

template<class ret_type,class M>
struct eval_reshape_helper<ret_type,M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& A,Integer r, Integer c)
    {
        ret = matcl::Matrix(A.reshape(r,c),false);
    };
};

template<class ret_type,class M>
struct eval_reshape_helper<ret_type,M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& A,Integer r, Integer c)
    {
        using value_type = typename M::value_type;

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();

        if ((Ar == r) && (Ac == c))
        {
            ret = matcl::Matrix(A,false);
            return;
        };
        if ((double) Ar * (double) Ac != (double) r * (double) c)
            throw error::invalid_reshape(Ar, Ac, r, c);

        if (An == 0)
        {
            auto tmp = matcl::details::zero_matrix<value_type,struct_sparse>
                            ::eval(A.get_type(),r,c);

            ret = matcl::Matrix(tmp,false);
            return;
        };

        ret_type res(A.get_type(),r,c,An);

        details::sparse_ccs<value_type>& d = res.rep();
        const details::sparse_ccs<value_type>& Ad = A.rep();

        const Integer* Ad_r = Ad.ptr_r();
        const Integer* Ad_c = Ad.ptr_c();
        const value_type* Ad_x = Ad.ptr_x();
        Integer off         = Ad.offset();

        Integer* d_r = d.ptr_r();
        Integer* d_c = d.ptr_c();
        value_type* d_x = d.ptr_x();

        for (Integer Aj = 0, j = 1, i = 0, lastAi = 0; Aj < Ac; ++Aj)
        {
            for (Integer k = Ad_c[Aj]; k < Ad_c[Aj + 1]; ++k)
            {
                Integer Ai	= Ad_r[k];
                i			+= Ai - lastAi;
                lastAi		= Ai;
                if (i >= r)
                { 
                    j += i/r;
                    i %= r; 
                };
                mrd::reset_helper(d_x[k-off],Ad_x[k]);
                d_r[k-off] = i;
                ++d_c[j];
            };
            lastAi	-= Ar;
            if (lastAi <= -r) 
            { 
                j -= lastAi/r; 
                lastAi %= r; 
            };
        };

        for (Integer j = 2; j <= c; ++j)
            d_c[j] += d_c[j - 1];

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class ret_type,class M>
struct eval_reshape_helper<ret_type,M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& A,Integer r, Integer c)
    {
        using value_type = typename M::value_type;

        Integer Ar = A.rows(), Ac = A.cols(), An = A.nnz();

        if ((Ar == r) && (Ac == c))
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if ((double) Ar * (double) Ac != (double) r * (double) c)
            throw error::invalid_reshape(Ar, Ac, r, c);

        if (An == 0)
        {
            auto tmp = matcl::details::zero_matrix<value_type,struct_sparse>
                            ::eval(A.get_type(),r,c);

            ret = matcl::Matrix(tmp,false);
            return;
        };

        ret_type res(A.get_type(),r,c,An);

        details::sparse_ccs<value_type>& d = res.rep();

        Integer* d_r		= d.ptr_r();
        Integer* d_c		= d.ptr_c();
        value_type* d_x		= d.ptr_x();
        Integer nz			= 0;

        const value_type* ptr_A = A.rep_ptr();
        Integer A_ld            = A.ld();

        for (Integer Aj = 0, i = 0, j = 1, lastAi = 0; Aj < Ac; ++Aj)
        {
            Integer first_row   = A.first_row(Aj);
            Integer last_row    = A.last_row(Aj);
            Integer pos         = A.first_elem_pos(Aj);

            for (Integer k = first_row; k <= last_row; ++k, ++pos)
            {
                Integer Ai	= k;
                i			+= Ai - lastAi;
                lastAi		= Ai;
                if (i >= r)
                { 
                    j += i/r;
                    i %= r; 
                };

                const value_type& tmp = ptr_A[pos];

                mrd::reset_helper(d_x[nz],tmp);
                d_r[nz] = i;
                ++d_c[j];
                ++nz;
            };
            lastAi	-= Ar;
            if (lastAi <= -r) 
            { 
                j       -= lastAi/r; 
                lastAi  %= r; 
            };

            ptr_A += A_ld;
        };

        for (Integer j = 2; j <= c; ++j)
            d_c[j] += d_c[j - 1];

        d.add_memory(-1);

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class M>
void details::manip_reshape_helper<M>::eval_reshape(matcl::Matrix& ret, const M& A, Integer m, Integer n)
{
    using struct_type = typename M::struct_type;
    return eval_reshape_helper<ret_type_reshape,M,struct_type>::eval(ret, A,m,n);
};

template<class ret_type,class M,class struct_type>
struct eval_vec_helper
{};

template<class ret_type,class M>
struct eval_vec_helper<ret_type,M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& A)
    {
        ret = matcl::Matrix(A.reshape(A.size(),1),false);
    };
};

template<class ret_type,class M>
struct eval_vec_helper<ret_type,M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& A)
    {
        using value_type = typename M::value_type;

        Integer r = A.rows(), c = A.cols(), nnz = A.nnz();
        error::check_size(r,c);

        Integer size = imult(r,c);

        if (nnz == 0 || r == 0 || c == 0)
        {
            auto tmp = matcl::details::zero_matrix<value_type,struct_sparse>
                            ::eval(A.get_type(),size,1);

            ret = matcl::Matrix(tmp,false);
            return;
        };

        const details::sparse_ccs<value_type>& Ad = A.rep();

        ret_type res(A.get_type(),size,1,nnz);
        details::sparse_ccs<value_type>& d = res.rep();

        const Integer* Ad_r = Ad.ptr_r();
        const Integer* Ad_c = Ad.ptr_c();
        const value_type* Ad_x = Ad.ptr_x();
        Integer off         = Ad.offset();

        Integer* d_r = d.ptr_r();
        Integer* d_c = d.ptr_c();
        value_type* d_x = d.ptr_x();

        d_c[0] = 0;
        for (Integer j = 0, jj = 0; j < c; ++j, jj+=r)
        {
            for (Integer i = Ad_c[j]; i < Ad_c[j + 1]; ++i)
            {
                d_r[i-off] = Ad_r[i] + jj;
                mrd::reset_helper(d_x[i-off],Ad_x[i]);
            }
        };
        d_c[1] = nnz;

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_vec_helper<ret_type,M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& A)
    {		
        using VT    = typename M::value_type;

        Integer r   = A.rows();
        Integer c   = A.cols();
        Integer nnz = A.nnz();

        error::check_size(r,c);

        Integer size = imult(r,c);

        if (nnz == 0 || r == 0 || c == 0)
        {
            auto tmp = matcl::details::zero_matrix<VT,struct_sparse>
                            ::eval(A.get_type(),size,1);
            ret = matcl::Matrix(tmp,false);
            return;
        };

        ret_type res(A.get_type(),size,1,nnz);
        details::sparse_ccs<VT>& d = res.rep();

        Integer* d_r    = d.ptr_r();
        Integer* d_c    = d.ptr_c();
        VT* d_x         = d.ptr_x();

        d_c[0]          = 0;
        Integer nz      = 0;
        
        Integer A_ld    = A.ld();
        Integer fd      = A.first_diag();
        Integer ld      = A.last_diag();

        if (fd == ld)
        {
            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(fd);

            Integer s   = A.diag_length(fd);
            Integer r0  = A.first_row_on_diag(fd);

            for (Integer j = 0; j < s; ++j, r0 += r + 1)
            {
                const VT& tmp = ptr_A[0];

                d_r[nz] = r0;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;

                ptr_A += A_ld;
            };
        }
        else
        {
            const VT* ptr_A = A.rep_ptr();

            for (Integer j = 0, jj = 0; j < c; ++j, jj += r)
            {
                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);
                Integer pos_A		= A.first_elem_pos(j);

                for (Integer i = first_row; i <= last_row; ++i, ++pos_A)
                {
                    const VT& tmp   = ptr_A[pos_A];

                    d_r[nz]     = i + jj;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                };

                ptr_A += A_ld;
            };
        };

        d_c[1] = nz;

        d.add_memory(-1);
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class M>
void details::manip_reshape_helper<M>::eval_vec(matcl::Matrix& ret, const M& m)
{
    using struct_type = typename M::struct_type;
    return eval_vec_helper<ret_type_vec,M,struct_type>::eval(ret,m);
};

template<class ret_type,class M,class struct_type>
struct eval_flipud_helper
{};

template<class ret_type,class M,class struct_type>
struct eval_fliplr_helper
{};

template<class ret_type,class M>
struct eval_flipud_helper<ret_type,M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        Integer r = m.rows(), c = m.cols();

        if (r == 1)
        {
            ret = matcl::Matrix(m,false);
            return;
        };

        ret_type res(m.get_type(),r, c);

        if (res.size() == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        const typename M::value_type* ptr_m = m.ptr();
        typename ret_type::value_type* ptr_res = res.ptr();

        Integer res_ld      = res.ld();
        Integer m_ld        = m.ld();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0, k = r-1; i < r; ++i, --k)
            {
                mrd::reset_helper(ptr_res[i],ptr_m[k]);
            };
            ptr_res += res_ld;
            ptr_m   += m_ld;
        }

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_fliplr_helper<ret_type,M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& m)
    {
        Integer r = m.rows(), c = m.cols();        

        if (c == 1)
        {
            ret = matcl::Matrix(m,false);
            return;
        };

        ret_type res(m.get_type(),r, c);
        const typename M::value_type* ptr_m = m.ptr() + imult(c-1,m.ld());
        typename ret_type::value_type* ptr_res = res.ptr();

        if (res.size() == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        Integer res_ld  = res.ld();
        Integer m_ld    = m.ld();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
            {
                mrd::reset_helper(ptr_res[i],ptr_m[i]);
            };
            ptr_res += res_ld;
            ptr_m   -= m_ld;
        }

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_fliplr_helper<ret_type,M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& A)
    {
        using value_type = typename M::value_type;

        Integer r = A.rows(), c = A.cols(), nnz = A.nnz();

        if (nnz == 0 || r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(A.get_type(),r,c,nnz),false);
            return;
        };
        if (c == 1)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        const details::sparse_ccs<value_type>& Ad = A.rep();

        ret_type res(A.get_type(),r,c,nnz);
        details::sparse_ccs<value_type>& d = res.rep();

        const Integer* Ad_r = Ad.ptr_r();
        const Integer* Ad_c = Ad.ptr_c();
        const value_type* Ad_x = Ad.ptr_x();

        Integer* d_r = d.ptr_r();
        Integer* d_c = d.ptr_c();
        value_type* d_x = d.ptr_x();

        d_c[0] = 0;
        for (Integer j = c - 1, pos_r = 0, pos_c = 1; j >= 0; --j, ++pos_c)
        {
            for (Integer i = Ad_c[j]; i < Ad_c[j+1]; ++i)
            {
                d_r[pos_r] = Ad_r[i];
                mrd::reset_helper(d_x[pos_r],Ad_x[i]);
                ++pos_r;
            }
            d_c[pos_c] = pos_r;
        };
        d_c[c] = nnz;

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_fliplr_helper<ret_type,M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& A)
    {
        using VT    = typename M::value_type;

        Integer r   = A.rows();
        Integer c   = A.cols();
        Integer nnz = A.nnz();

        if (nnz == 0 || r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(A.get_type(),r,c),false);
            return;
        };

        ret_type res(A.get_type(),r,c,nnz);
        details::sparse_ccs<VT>& d = res.rep();

        Integer* d_r    = d.ptr_r();
        Integer* d_c    = d.ptr_c();
        VT* d_x         = d.ptr_x();

        Integer nz      = 0;
        
        Integer A_ld    = A.ld();
        Integer fd      = A.first_diag();
        Integer ld      = A.last_diag();

        if (fd == ld)
        {
            Integer pos_c = 1;
            Integer s   = A.diag_length(fd);
            Integer r0  = A.first_row_on_diag(fd);
            Integer c0  = A.first_col_on_diag(fd);

            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(fd);

            for (Integer j = c; j > c0 + s; --j, ++pos_c)
                d_c[pos_c] = 0;

            ptr_A       += (s - 1) * A_ld;
            r0          += s - 1;

            for (Integer j = c0 + s; j > c0; --j, ++pos_c)
            {
                const VT& tmp   = ptr_A[0];

                d_r[nz]     = r0;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;				

                d_c[pos_c]  = nz;
                ptr_A       -= A_ld;
                r0          -= 1;
            };

            for (Integer j = c0; j > 0; --j, ++pos_c)
                d_c[pos_c]  = nz;
        }
        else
        {
            const VT* ptr_A = A.rep_ptr() + imult(c-1, A_ld);

            for (Integer j = c-1, pos_c = 1; j >= 0; --j, ++pos_c)
            {
                Integer first_row   = A.first_row(j);
                Integer last_row    = A.last_row(j);
                Integer pos         = A.first_elem_pos(j);

                for (Integer i = first_row; i <= last_row; ++i, ++pos)
                {
                    const VT& tmp   = ptr_A[pos];
                    d_r[nz]         = i;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;				
                };
                d_c[pos_c]  = nz;
                ptr_A       -= A_ld;
            };
        };

        d_c[c] = nz;

        d.add_memory(-1);
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_flipud_helper<ret_type,M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& A)
    {
        using value_type = typename M::value_type;

        Integer r = A.rows(), c = A.cols(), nnz = A.nnz();

        if (nnz == 0|| r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(A.get_type(), r,c,nnz),false);
            return;
        };
        if (r == 1)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        const details::sparse_ccs<value_type>& Ad = A.rep();

        ret_type res(A.get_type(),r,c,nnz);
        details::sparse_ccs<value_type>& d = res.rep();

        const Integer* Ad_r = Ad.ptr_r();
        const Integer* Ad_c = Ad.ptr_c();
        const value_type* Ad_x = Ad.ptr_x();

        Integer* d_r = d.ptr_r();
        Integer* d_c = d.ptr_c();
        value_type* d_x = d.ptr_x();

        d_c[0] = 0;
        for (Integer j = 0, pos_r = 0, pos_c = 1; j < c; ++j, ++pos_c)
        {
            for (Integer i = Ad_c[j+1]-1; i >= Ad_c[j]; --i)
            {
                d_r[pos_r] = r - 1 - Ad_r[i];
                mrd::reset_helper(d_x[pos_r],Ad_x[i]);
                ++pos_r;
            }
            d_c[pos_c] = pos_r;
        };
        d_c[c] = nnz;

        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_flipud_helper<ret_type,M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& A)
    {
        using VT    = typename M::value_type;

        Integer r   = A.rows();
        Integer c   = A.cols();
        Integer nnz = A.nnz();

        if (nnz == 0|| r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(A.get_type(),r,c,nnz),false);
            return;
        };

        ret_type res(A.get_type(),r,c,nnz);
        details::sparse_ccs<VT>& d = res.rep();

        Integer* d_r    = d.ptr_r();
        Integer* d_c    = d.ptr_c();
        VT* d_x         = d.ptr_x();

        Integer nz	    = 0;
        d_c[0]		    = 0;
        
        Integer A_ld    = A.ld();
        Integer fd      = A.first_diag();
        Integer ld      = A.last_diag();

        if (fd == ld)
        {
            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(fd);
            Integer s   = A.diag_length(fd);
            Integer r0  = A.first_row_on_diag(fd);
            Integer c0  = A.first_col_on_diag(fd);

            for (Integer j = 0; j < c0; ++j)
                d_c[j] = nz;

            for (Integer j = c0; j < c0 + s; ++j, ++r0)
            {
                d_c[j] = nz;

                const VT& tmp = ptr_A[0];

                d_r[nz]		= r - r0 - 1;
                mrd::reset_helper(d_x[nz],tmp);
                ++nz;
                
                ptr_A   += A_ld;
            };

            for (Integer j = c0 + s; j <= c; ++j)
                d_c[j] = nz;
        }
        else
        {
            const VT* ptr_A = A.rep_ptr();

            for (Integer j = 0; j < c; ++j)
            {
                d_c[j] = nz;

                Integer first_row	= A.first_row(j);
                Integer last_row	= A.last_row(j);
                Integer pos			= A.first_elem_pos(j) + last_row - first_row;

                for (Integer i = last_row; i >= first_row; --i, --pos)
                {
                    const VT& tmp   = ptr_A[pos];
                    d_r[nz]		    = r - i - 1;
                    mrd::reset_helper(d_x[nz],tmp);
                    ++nz;
                }				
                ptr_A += A_ld;
            };

            d_c[c] = nz;
        };

        d.add_memory(-1);
        ret = matcl::Matrix(res,false);
        return;
    }; 
};

template<class M>
void details::manip_reshape_helper<M>::eval_flipud(matcl::Matrix& ret, const M& m)
{
    using struct_type = typename M::struct_type;
    return eval_flipud_helper<ret_type_vec,M,struct_type>::eval(ret,m);
};

template<class M>
void details::manip_reshape_helper<M>::eval_fliplr(matcl::Matrix& ret, const M& m)
{
    using struct_type = typename M::struct_type;
    return eval_fliplr_helper<ret_type_vec,M,struct_type>::eval(ret,m);
};

template<class ret_type,class M,class struct_type>
struct eval_repmat_helper
{};

template<class ret_type,class M>
struct eval_repmat_helper<ret_type,M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& mat,Integer m, Integer n)
    {		
        Integer matr = mat.rows(), matc = mat.cols();
        error::check_size(m,matr);
        error::check_size(n,matc);
        Integer resr = imult(m,matr);
        Integer resc = imult(n,matc);
     
        using val_in    = typename M::value_type;
        using val_ret   = typename ret_type::value_type;

        ret_type res(mat.get_type(),resr, resc);
     
        if (res.size() == 0)
        {
            ret = matcl::Matrix(res,false); 
            return;
        };

        const val_in* ptr_mat = mat.ptr();
        val_ret* ptr_res    = res.ptr();
        Integer mat_ld      = mat.ld();
        Integer res_ld      = res.ld();
     
        if (matr == 1 && n == 1) 
        {
            for (Integer j = 0; j < matc; ++j, ptr_mat += mat_ld)
            {
                const val_in& val = ptr_mat[0];
                for (Integer i = 0; i < m; ++i)
                    mrd::reset_helper(ptr_res[i],val);

                ptr_res += res_ld;
            }
    
            ret = matcl::Matrix(res,false);
            return;
        };

        for (Integer j = 0; j < n; ++j)
        {
            ptr_mat = mat.ptr();
            for (Integer matj = 0; matj < matc; ++matj)
            {
                for (Integer i = 0, pos_res = 0; i < m; ++i)
                {
                    for (Integer mati = 0; mati < matr; ++mati, ++pos_res)
                        mrd::reset_helper(ptr_res[pos_res],ptr_mat[mati]);
                };

                ptr_mat += mat_ld;
                ptr_res += res_ld;
            };
        };
     
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_repmat_helper<ret_type,M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& mat,Integer m, Integer n)
    {
        Integer matr = mat.rows(), matc = mat.cols();
        
        error::check_size(m,matr);
        error::check_size(n,matc);

        Integer resr = imult(m,matr), resc = imult(n,matc);
        Integer resnnz = imult(imult(mat.nnz(),m),n);
        
        if (m == 1 && n == 1)
        {
            ret = matcl::Matrix(mat,false);
            return;
        };

        ret_type res(mat.get_type(),resr, resc,resnnz);

        if (resnnz == 0 || resr == 0 || resc == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        using value_type        = typename M::value_type;

        details::sparse_ccs<value_type>& rep	= res.rep();
        const details::sparse_ccs<value_type>& Arep	= mat.rep();

        Integer* d_c			= rep.ptr_c();
        Integer* d_r			= rep.ptr_r();
        value_type* d_x			= rep.ptr_x();

        const Integer* Ad_c		= Arep.ptr_c();
        const Integer* Ad_r		= Arep.ptr_r();
        const value_type* Ad_x	= Arep.ptr_x();

        d_c[0] = 0;
        Integer nz = 0;

        for (Integer j = 0; j < matc; ++j)
        {
            for (Integer k = 1, pos_r = 0; k <= m; ++k, pos_r+= matr)
            {
                for (Integer i = Ad_c[j]; i < Ad_c[j+1]; ++i)
                {
                    Integer p               = Ad_r[i];
                    const value_type& val   = Ad_x[i];

                    d_r[nz] = pos_r + p;
                    mrd::reset_helper(d_x[nz],val);
                    ++nz;
                };

                d_c[j+1] = nz;
            };
        };

        for (Integer k = 2, col = matc + 1; k<= n; ++k)
        {
            for (Integer j = 0; j < matc; ++j, ++col)
            {
                for (Integer i = d_c[j]; i < d_c[j+1]; ++i)
                {
                    Integer p               = d_r[i];
                    const value_type& val   = d_x[i];

                    d_r[nz] = p;
                    mrd::reset_helper(d_x[nz],val);
                    ++nz;
                };

                d_c[col] = nz;
            };
        };

        rep.add_memory(-1);
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class ret_type,class M>
struct eval_repmat_helper<ret_type,M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& mat,Integer m, Integer n)
    {
        if (m == 1 && n == 1)
        {
            ret = matcl::Matrix(mat,false);
            return;
        };

        Integer matr    = mat.rows();
        Integer matc    = mat.cols();

        error::check_size(m,matr);
        error::check_size(n,matc);
        
        Integer resr    = imult(m, matr);
        Integer resc    = imult(n, matc);
        Integer resnnz  = imult(imult(mat.nnz(), m), n); 

        ret_type res(mat.get_type(),resr, resc,resnnz);

        if (resr == 0 || resc == 0 || resnnz == 0)
        {
            ret = matcl::Matrix(res,false);
            return;
        };

        using VT        = typename M::value_type;

        details::sparse_ccs<VT>& rep	= res.rep();

        Integer* d_c	= rep.ptr_c();
        Integer* d_r	= rep.ptr_r();
        VT* d_x         = rep.ptr_x();

        Integer nz      = 0;
        
        Integer mat_ld  = mat.ld();
        Integer fd      = mat.first_diag();
        Integer ld      = mat.last_diag();

        if (fd == ld)
        {
            const VT* ptr_mat   = mat.rep_ptr() + mat.first_elem_diag(fd);

            Integer s       = mat.diag_length(fd);
            Integer r0      = mat.first_row_on_diag(fd);
            Integer c0      = mat.first_col_on_diag(fd);

            for (Integer j = 0; j < c0; ++j)
                d_c[j]      = nz;

            for (Integer j = c0; j < c0 + s; ++j, ++r0)
            {
                d_c[j]          = nz;
                const VT& val   = ptr_mat[0];
                
                for (Integer k = 1, pos_r = r0; k <= m; ++k, pos_r += matr)
                {
                    d_r[nz] = pos_r;
                    mrd::reset_helper(d_x[nz],val);
                    ++nz;					
                };

                ptr_mat     += mat_ld;
            };

            for (Integer j = c0 + s; j <= matc; ++j)
                d_c[j]      = nz;
        }
        else
        {
            const VT* ptr_mat   = mat.rep_ptr();

            for (Integer j = 0; j < matc; ++j)
            {
                d_c[j] = nz;

                for (Integer k = 1, pos_r = 0; k <= m; ++k, pos_r += matr)
                {
                    Integer first_row   = mat.first_row(j);
                    Integer last_row    = mat.last_row(j);
                    Integer pos_A       = mat.first_elem_pos(j);

                    for (Integer i = first_row; i <= last_row; ++i, ++pos_A)
                    {
                        Integer p       = i;
                        const VT& val   = ptr_mat[pos_A];
                        d_r[nz]         = pos_r + p;
                        mrd::reset_helper(d_x[nz],val);
                        ++nz;
                    };
                };		

                ptr_mat += mat_ld;
            };

            d_c[matc]   = nz;
        };

        for (Integer k = 2, col = matc + 1; k<= n; ++k)
        {
            for (Integer j = 0; j < matc; ++j, ++col)
            {
                for (Integer i = d_c[j]; i < d_c[j+1]; ++i)
                {
                    Integer p       = d_r[i];
                    const VT& val   = d_x[i];
                    d_r[nz]         = p;

                    mrd::reset_helper(d_x[nz],val);
                    ++nz;
                };

                d_c[col] = nz;
            };
        };

        rep.add_memory(-1);
        ret = matcl::Matrix(res,false);
        return;
    };
};

template<class M>
void details::manip_reshape_helper<M>::eval_repmat(matcl::Matrix& ret, const M& A, Integer m, Integer n)
{
    using struct_type = typename M::struct_type;
    return eval_repmat_helper<ret_type_repmat,M,struct_type>::eval(ret,A,m,n);
};

template<class ret_type,class M, class struct_type>
struct eval_trans_helper
{};

template<class ret_type,class M>
struct eval_trans_helper<ret_type,M,struct_dense>
{
    using VTR   = typename ret_type::value_type;

    static void eval_mat(VTR* ptr_res, Integer res_ld, const M& m)
    {
        Integer r = m.rows();
        Integer c = m.cols();		

        using VT    = typename M::value_type;        

        Integer m_ld    = m.ld();

        //block size
        const Integer NB    = optim_params::block_size_trans;

        //blocked version, reduce cache misses in matrix transpose
        for (Integer ib = 0; ib < r; ib += NB)
        {
            for (Integer jb = 0; jb < c; jb += NB)
            {
                Integer ibl     = std::min(r, ib + NB);
                Integer jbl     = std::min(c, jb + NB);

                const VT* ptr_m = m.ptr() + ib;
                VTR * ptr_C     = ptr_res + ib * res_ld;

                for (Integer i = ib; i < ibl; ++i)
                {                    
                    for (Integer j = jb; j < jbl; ++j)
                        mrd::reset_helper(ptr_C[j], VTR(ptr_m[j*m_ld]));

                    ptr_C       += res_ld;
                    ptr_m       += 1;
                };
            };
        };

        return;
    };

    static void eval(ret_type& ret, const M& m)
    {
        bool is_A_sq    = m.rows() == m.cols();

        if (m.get_struct().is_symmetric(is_A_sq, is_real_matrix(m)) && is_A_sq)
        {
            ret.assign_to_fresh(raw::converter<ret_type,M>::eval(m));
            return;
        };

        Integer r = m.rows();
        Integer c = m.cols();		

        if (r <= 1 || c <= 1)
        {
            M out = m.reshape(c,r);

            out.set_struct(md::predefined_struct::get_trans(m.get_struct()));
            ret.assign_to_fresh(raw::converter<ret_type,M>::eval(out));

            return;
        }

        ret_type res(m.get_type());
        res.reset_unique(c, r);

        using VT    = typename M::value_type;
        using VTR   = typename ret_type::value_type;

        Integer m_ld    = m.ld();
        Integer res_ld  = res.ld();

        //block size
        const Integer NB    = optim_params::block_size_trans;

        //blocked version, reduce cache misses in matrix transpose
        for (Integer ib = 0; ib < r; ib += NB)
        {
            for (Integer jb = 0; jb < c; jb += NB)
            {
                Integer ibl     = std::min(r, ib + NB);
                Integer jbl     = std::min(c, jb + NB);

                const VT* ptr_m = m.ptr() + ib;
                VTR * ptr_C     = res.ptr() + ib * res_ld;

                for (Integer i = ib; i < ibl; ++i)
                {                    
                    for (Integer j = jb; j < jbl; ++j)
                        mrd::reset_helper(ptr_C[j], VTR(ptr_m[j*m_ld]));

                    ptr_C       += res_ld;
                    ptr_m       += 1;
                };
            };
        };

        res.set_struct(md::predefined_struct::get_trans(m.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };
};

template<class ret_type,class M>
struct eval_trans_helper<ret_type,M,struct_banded>
{
    static void eval(ret_type& ret, const M& A)
    {
        bool is_A_sq    = A.rows() == A.cols();
        if (is_A_sq && A.get_struct().is_symmetric(is_A_sq, is_real_matrix(A)))
        {
            ret.assign_to_fresh(raw::converter<ret_type,M>::eval(A));
            return;
        };

        Integer r       = A.rows();
        Integer c       = A.cols();
        Integer fd      = A.first_diag();
        Integer ld      = A.last_diag();
        Integer A_ld    = A.ld();

        using VT        = typename M::value_type;
        using VTR       = typename ret_type::value_type;        

        if (fd == ld)
        {
            if (r == c && fd == 0)
            {
                ret.assign_to_fresh(raw::converter<ret_type,M>::eval(A));
                return;
            };

            ret_type res(A.get_type(), c, r, -fd, -fd);

            Integer s       = A.diag_length(fd);
            Integer res_ld  = res.ld();

            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(fd);
            VTR* ptr_res    = res.rep_ptr() + res.first_elem_diag(-fd);

            for (Integer i = 0; i < s; ++i)
            {
                mrd::reset_helper<VTR>(ptr_res[0],VTR(ptr_A[0]));
                ptr_A       += A_ld;
                ptr_res     += res_ld;
            };

            res.set_struct(md::predefined_struct::get_trans(A.get_struct()));
            ret.assign_to_fresh(res);

            return;
        };

        ret_type res(A.get_type(),c, r, -ld, -fd);
        
        Integer res_ld      = res.ld();

        for (Integer d = fd; d <= ld; ++d)
        {
            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(d);
            VTR* ptr_res    = res.rep_ptr() + res.first_elem_diag(-d);
            Integer s       = A.diag_length(d);

            for (Integer i = 1; i <= s; ++i)
            {
                mrd::reset_helper(ptr_res[0], VTR(ptr_A[0]));

                ptr_res     += res_ld;
                ptr_A       += A_ld;
            }
        };

        res.set_struct(md::predefined_struct::get_trans(A.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };
};

template<class ret_type,class M>
struct eval_trans_helper<ret_type,M,struct_sparse>
{
    static void eval_reshape(ret_type& ret, const M& A, Integer max_row, Integer max_col,
                             Integer ret_rows, Integer ret_cols)
    {
        using VT        = typename M::value_type;
        using VTR       = typename ret_type::value_type;
        
        Integer r = A.rows();
        Integer c = A.cols();
        Integer n = A.nnz();

        max_row = std::min(max_row,r);
        max_col = std::min(max_col,c);
        max_col = std::min(max_col, ret_rows);
        max_row = std::min(max_row, ret_cols);

        if (n == 0 || r == 0 || c == 0)
        {
            ret_type res(A.get_type(), ret_rows, ret_cols);
            ret.assign_to_fresh(res);
            return;
        };

        details::sparse_ccs<VTR> d(A.get_type(),ret_rows, ret_cols, n);

        const details::sparse_ccs<VT>& Ad = A.rep();
        const Integer * Ad_c    = Ad.ptr_c();
        const Integer * Ad_r    = Ad.ptr_r() + Ad.offset();		
        const VT * Ad_x         = Ad.ptr_x() + Ad.offset();

        Integer * d_r   = d.ptr_r();
        Integer * d_c   = d.ptr_c();
        VTR * d_x       = d.ptr_x();

        matcl::pod_workspace<Integer> p(max_row + 1, 0);

        Integer i, j, k;

        Integer n2 = Ad_c[max_col] - Ad.offset();

        for (i = 0; i < n2; ++i)
        {
            Integer cur_row = Ad_r[i];
            
            if (cur_row >= max_row)
                continue;
            
            ++p[cur_row + 1];
        };

        for (k = 0, i = 1; i <= max_row; ++i)
        {
            p[i]    += k;
            k       = p[i];
            d_c[i]  = p[i];
        }

        Ad_r = Ad.ptr_r();
        Ad_x = Ad.ptr_x();

        for (j = 0; j < max_col; ++j)
        {
            for (i = Ad_c[j]; i < Ad_c[j + 1]; ++i)
            {
                Integer cur_row = Ad_r[i];

                if (cur_row >= max_row)
                    continue;

                k       = p[cur_row]++;
                d_r[k]  = j;

                mrd::reset_helper(d_x[k],VTR(Ad_x[i]));
            }
        };

        Integer nz = d_c[max_row];
        for (j = max_row + 1; j <= ret_cols; ++j)
            d_c[j] = nz;

        ret_type res = sparse_matrix_base<VTR>(d);
        res.set_struct(md::predefined_struct::get_trans(A.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };

    static void eval(ret_type& ret, const M& A)
    {
        bool is_A_sq    = A.rows() == A.cols();

        if (A.get_struct().is_symmetric(is_A_sq, is_real_matrix(A))  && is_A_sq)
        {
            ret.assign_to_fresh(raw::converter<ret_type,M>::eval(A));
            return;
        };

        using VT        = typename M::value_type;
        using VTR       = typename ret_type::value_type;
        
        Integer r = A.rows();
        Integer c = A.cols();
        Integer n = A.nnz();        

        if (n == 0 || r == 0 || c == 0)
        {
            ret_type res(A.get_type(), c, r);
            ret.assign_to_fresh(res);
            return;
        };

        details::sparse_ccs<VTR> d(A.get_type(),c, r, n);

        const details::sparse_ccs<VT>& Ad = A.rep();
        const Integer * Ad_c    = Ad.ptr_c();
        const Integer * Ad_r    = Ad.ptr_r() + Ad.offset();		
        const VT * Ad_x         = Ad.ptr_x() + Ad.offset();

        Integer * d_r   = d.ptr_r();
        Integer * d_c   = d.ptr_c();
        VTR * d_x       = d.ptr_x();

        matcl::pod_workspace<Integer> p(r + 1, 0);

        Integer i, j, k;

        for (i = 0; i < n; ++i)
        {
            ++p[Ad_r[i] + 1];
        };

        for (k = 0, i = 1; i <= r; ++i)
        {
            p[i]    += k;
            k       = p[i];
            d_c[i]  = p[i];
        }

        Ad_r = Ad.ptr_r();
        Ad_x = Ad.ptr_x();

        for (j = 0; j < c; ++j)
        {
            for (i = Ad_c[j]; i < Ad_c[j + 1]; ++i)
            {
                k = p[Ad_r[i]]++;
                d_r[k] = j;
                mrd::reset_helper(d_x[k],VTR(Ad_x[i]));
            }
        };

        ret_type res = sparse_matrix_base<VTR>(d);
        res.set_struct(md::predefined_struct::get_trans(A.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };
};

template<class ret_type,class M, bool is_compl_or_obj, class struct_type>
struct eval_ctrans_helper
{};

template<class ret_type,class M, class struct_type>
struct eval_ctrans_helper<ret_type,M,false,struct_type>
{
    static void eval(ret_type& ret, const M& m)
    {
        return eval_trans_helper<ret_type,M,struct_type>::eval(ret,m);
    };

    static void eval_reshape(ret_type& ret, const M& m, Integer max_row, Integer max_col,
                             Integer ret_rows, Integer ret_cols)
    {
        return eval_trans_helper<ret_type,M,struct_type>
            ::eval_reshape(ret,m,max_row,max_col,ret_rows,ret_cols);
    };
};

template<class ret_type,class M>
struct eval_ctrans_helper<ret_type,M,true,struct_dense>
{
    using VTR   = typename ret_type::value_type;

    static void eval_mat(VTR* ptr_res, Integer res_ld, const M& m)
    {
        Integer r = m.rows();
        Integer c = m.cols();

        using VT    = typename M::value_type;        

        Integer m_ld    = m.ld();

        //block size
        const Integer NB    = optim_params::block_size_trans;

        //blocked version, reduce cache misses in matrix transpose
        for (Integer ib = 0; ib < r; ib += NB)
        {
            for (Integer jb = 0; jb < c; jb += NB)
            {
                Integer ibl     = std::min(r, ib + NB);
                Integer jbl     = std::min(c, jb + NB);

                const VT* ptr_m = m.ptr() + ib;
                VTR * ptr_C     = ptr_res + ib * res_ld;

                for (Integer i = ib; i < ibl; ++i)
                {                    
                    for (Integer j = jb; j < jbl; ++j)
                    {
                        VT ac = mrd::conj_helper<VT>::eval(ptr_m[j*m_ld]);
                        mrd::reset_helper(ptr_C[j], ac);
                    }

                    ptr_C       += res_ld;
                    ptr_m       += 1;
                };
            };
        };

        return;
    };

    static void eval(ret_type& ret, const M& m)
    {
        bool is_A_sq    = m.rows() == m.cols();

        if (m.get_struct().is_hermitian(is_A_sq, is_real_matrix(m)) && is_A_sq)
        {
            ret.assign_to_fresh(raw::converter<ret_type,M>::eval(m));
            return;
        };

        Integer r = m.rows();
        Integer c = m.cols();

        ret_type res(m.get_type(),c, r);

        if (m.size() == 0)
        {
            ret.assign_to_fresh(res);
            return;
        };

        using VT    = typename M::value_type;        

        Integer m_ld    = m.ld();
        Integer res_ld  = res.ld();

        //block size
        const Integer NB    = optim_params::block_size_trans;

        //blocked version, reduce cache misses in matrix transpose
        for (Integer ib = 0; ib < r; ib += NB)
        {
            for (Integer jb = 0; jb < c; jb += NB)
            {
                Integer ibl     = std::min(r, ib + NB);
                Integer jbl     = std::min(c, jb + NB);

                const VT* ptr_m = m.ptr() + ib;
                VTR * ptr_C     = res.ptr() + ib * res_ld;

                for (Integer i = ib; i < ibl; ++i)
                {                    
                    for (Integer j = jb; j < jbl; ++j)
                    {
                        VT ac = mrd::conj_helper<VT>::eval(ptr_m[j*m_ld]);
                        mrd::reset_helper(ptr_C[j], VTR(ac));
                    }

                    ptr_C       += res_ld;
                    ptr_m       += 1;
                };
            };
        };

        res.set_struct(md::predefined_struct::get_ctrans(m.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };
};

template<class ret_type,class M>
struct eval_ctrans_helper<ret_type,M,true,struct_banded>
{
    static void eval(ret_type& ret, const M& A)
    {
        bool is_A_sq    = A.rows() == A.cols();

        if (A.get_struct().is_hermitian(is_A_sq, is_real_matrix(A)) && is_A_sq)
        {
            ret.assign_to_fresh(raw::converter<ret_type,M>::eval(A));
            return;
        };

        Integer r   = A.rows();
        Integer c   = A.cols();
        Integer fd  = A.first_diag();
        Integer ld  = A.last_diag();

        using VT    = typename M::value_type;
        using VTR   = typename ret_type::value_type;
        
        Integer A_ld    = A.ld();

        if (fd == ld && fd == 0)
        {
            static_assert(md::is_complex<VT>::value || std::is_same<VT,Object>::value, "");
            
            Integer s       = A.diag_length(fd);

            ret_type res(A.get_type(), c, r, -fd, -fd);

            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(fd);
            VTR* ptr_res    = res.rep_ptr() + res.first_elem_diag(-fd);
            Integer res_ld  = res.ld();

            for (Integer i = 0; i < s; ++i)
            {
                VT ca       = mrd::conj_helper<VT>::eval(ptr_A[0]);
                mrd::reset_helper<VTR>(ptr_res[0],VTR(ca));
                ptr_A       += A_ld;
                ptr_res     += res_ld;
            };

            res.set_struct(md::predefined_struct::get_ctrans(A.get_struct()));
            ret.assign_to_fresh(res);
            return;
        };
        
        ret_type res(A.get_type(),c, r, -ld, -fd);
        
        Integer res_ld      = res.ld();

        for (Integer d = fd; d <= ld; ++d)
        {
            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(d);
            VTR* ptr_res    = res.rep_ptr() + res.first_elem_diag(-d);
            Integer s       = A.diag_length(d);

            for (Integer i = 1; i <= s; ++i)
            {
                VT ac       = mrd::conj_helper<VT>::eval(ptr_A[0]);
                mrd::reset_helper<VTR>(ptr_res[0], VTR(ac));

                ptr_res     += res_ld;
                ptr_A       += A_ld;
            }
        };

        res.set_struct(md::predefined_struct::get_ctrans(A.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };
};

template<class ret_type,class M>
struct eval_ctrans_helper<ret_type,M,true,struct_sparse>
{
    static void eval_reshape(ret_type& ret, const M& A, Integer max_row, Integer max_col,
                             Integer ret_rows, Integer ret_cols)
    {
        using VT        = typename M::value_type;
        using VTR       = typename ret_type::value_type;
        
        Integer r = A.rows();
        Integer c = A.cols();
        Integer n = A.nnz();

        max_row = std::min(max_row,r);
        max_col = std::min(max_col,c);
        max_col = std::min(max_col, ret_rows);
        max_row = std::min(max_row, ret_cols);

        if (n == 0 || r == 0 || c == 0)
        {
            ret_type res(A.get_type(), ret_rows, ret_cols);
            ret.assign_to_fresh(res);
            return;
        };

        details::sparse_ccs<VTR> d(A.get_type(),ret_rows, ret_cols, n);

        const details::sparse_ccs<VT>& Ad = A.rep();
        const Integer * Ad_c    = Ad.ptr_c();
        const Integer * Ad_r    = Ad.ptr_r() + Ad.offset();		
        const VT * Ad_x         = Ad.ptr_x() + Ad.offset();

        Integer * d_r   = d.ptr_r();
        Integer * d_c   = d.ptr_c();
        VTR * d_x       = d.ptr_x();

        matcl::pod_workspace<Integer> p(max_row + 1, 0);

        Integer i, j, k;

        Integer n2 = Ad_c[max_col] - Ad.offset();

        for (i = 0; i < n2; ++i)
        {
            Integer cur_row = Ad_r[i];
            
            if (cur_row >= max_row)
                continue;
            
            ++p[cur_row + 1];
        };

        for (k = 0, i = 1; i <= max_row; ++i)
        {
            p[i]    += k;
            k       = p[i];
            d_c[i]  = p[i];
        }

        Ad_r = Ad.ptr_r();
        Ad_x = Ad.ptr_x();

        for (j = 0; j < max_col; ++j)
        {
            for (i = Ad_c[j]; i < Ad_c[j + 1]; ++i)
            {
                Integer cur_row = Ad_r[i];

                if (cur_row >= max_row)
                    continue;

                k       = p[cur_row]++;
                d_r[k]  = j;

                VT ac = mrd::conj_helper<VT>::eval(Ad_x[i]);
                mrd::reset_helper(d_x[k],VTR(ac));
            }
        };

        Integer nz = d_c[max_row];
        for (j = max_row + 1; j <= ret_cols; ++j)
            d_c[j] = nz;

        ret_type res = sparse_matrix_base<VTR>(d);
        res.set_struct(md::predefined_struct::get_trans(A.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };

    static void eval(ret_type& ret, const M& A)
    {        
        bool is_A_sq    = A.rows() == A.cols();

        if (A.get_struct().is_hermitian(is_A_sq, is_real_matrix(A)) && is_A_sq)
        {
            ret.assign_to_fresh(raw::converter<ret_type,M>::eval(A));
            return;
        };

        using VTR       = typename ret_type::value_type;
        using VT        = typename M::value_type;
        
        Integer r = A.rows();
        Integer c = A.cols();
        Integer n = A.nnz();

        details::sparse_ccs<VTR> d(A.get_type(),c, r, n);
        
        if (n == 0|| r == 0 || c == 0) 
        {
            ret_type res(A.get_type(), c, r);
            ret.assign_to_fresh(res);
            return;
        };

        const details::sparse_ccs<VT>& Ad = A.rep();

        matcl::pod_workspace<Integer> p(r + 1, 0);

        const Integer* Ad_c     = Ad.ptr_c();
        const Integer* Ad_r     = Ad.ptr_r() + Ad.offset();		
        const VT* Ad_x          = Ad.ptr_x() + Ad.offset();

        Integer* d_r    = d.ptr_r();
        Integer* d_c    = d.ptr_c();
        VTR* d_x        = d.ptr_x();

        Integer i, j, k;

        for (i = 0; i < n; ++i)
        {
            ++p[Ad_r[i] + 1];
        };

        for (k = 0, i = 1; i <= r; ++i)
        {
            p[i]    += k;
            d_c[i]  = p[i];
            k       = p[i];
        };

        Ad_r = Ad.ptr_r();
        Ad_x = Ad.ptr_x();

        for (j = 0; j < c; ++j)
        {
            for (i = Ad_c[j]; i < Ad_c[j + 1]; ++i)
            {
                k       = p[Ad_r[i]]++;
                d_r[k]  = j;

                VT ac = mrd::conj_helper<VT>::eval(Ad_x[i]);
                mrd::reset_helper(d_x[k],VTR(ac));
            }
        };

        ret_type res = sparse_matrix_base<VTR>(d);
        res.set_struct(md::predefined_struct::get_ctrans(A.get_struct()));
        ret.assign_to_fresh(res);
        return;
    };
};

template<class M>
void details::manip_trans_helper<M>::eval_trans(M& ret, const M& m)
{
    using struct_type = typename M::struct_type;
    return eval_trans_helper<ret_type,M,struct_type>::eval(ret, m);
};

template<class M>
void details::manip_trans_helper<M>::eval_ctrans(M& ret, const M& m)
{
    using struct_type = typename M::struct_type;
    static const bool is_compl = md::is_complex<typename M::value_type>::value;
    static const bool is_obj = std::is_same<typename M::value_type,Object>::value;
    return eval_ctrans_helper<ret_type,M,is_compl || is_obj,struct_type>::eval(ret,m);
};

template<class struct_type, class Val_in, class Val_ret>
void details::manip_trans_converter_helper<struct_type, Val_in, Val_ret>::eval_trans(ret_type& ret, const in_type& m)
{
    return eval_trans_helper<ret_type,in_type,struct_type>::eval(ret, m);
};

template<class struct_type, class Val_in, class Val_ret>
void details::manip_trans_converter_helper<struct_type, Val_in, Val_ret>::eval_ctrans(ret_type& ret, const in_type& m)
{
    static const bool is_compl  = md::is_complex<Val_in>::value;
    static const bool is_obj    = std::is_same<Val_in,Object>::value;
    return eval_ctrans_helper<ret_type,in_type,is_compl || is_obj,struct_type>::eval(ret,m);
};

struct estim_bandwidth
{
    template<class Mat>
    static void eval(const Mat& mat, Integer fd0, Integer ld0, Integer& fd, Integer& ld)
    {
        if (fd0 <= 0)
        {
            Integer lde = raw::get_ld(mat,std::max(-fd0,0));
            fd          = std::max(fd0, -lde);
        }
        else if (mat.rows() == 0)
        {
            fd          = 0;
        }
        else
        {
            fd          = std::max(fd0, -mat.rows() + 1);
        }

        if (ld0 >= 0)
        {
            Integer ud  = raw::get_ud(mat,std::max(ld0,0));
            ld          = std::min(ld0, ud);
        }
        else if (mat.cols() == 0)
        {
            ld          = 0;
        }
        else
        {
            ld          = std::min(ld0, mat.cols() - 1);
        }
    };
};

template<class ret_type,class M,class struct_type>
struct eval_tril_helper
{};

template<class ret_type,class M,class struct_type>
struct eval_triu_helper
{};

template<class M,class struct_type>
struct eval_select_band_helper
{};

template<class ret_type,class M>
struct eval_tril_helper<ret_type,M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& A, Integer d, bool rvalue)
    {
        Integer r = A.rows(), c = A.cols();

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (d >= 0 && A.get_struct().is_tril())
        {
            ret = matcl::Matrix(A,false);
            return;            
        };
        if (d >= 1 && A.get_struct().is_hessl())
        {
            ret = matcl::Matrix(A,false);
            return;            
        };

        if (rvalue == true && A.is_unique() == true)
            return eval_inplace(ret, A, d);

        ret_type res(A.get_type(),r, c);

        using value_type    = typename ret_type::value_type;
        using value_type_in = typename M::value_type;

        value_type Z = md::default_value<value_type>(A.get_type());

        value_type* ptr_res         = res.ptr();
        const value_type_in* ptr_A  = A.ptr();
        Integer A_ld                = A.ld();
        Integer res_ld              = res.ld();

        for (Integer j = 0; j < c; ++j)
        {
            Integer mjmdrp1 = (j - d < r) ? j - d : r;
            Integer i;
            for (i = 0; i < mjmdrp1; ++i)
                mrd::reset_helper(ptr_res[i],Z);

            for (; i < r; ++i)
                mrd::reset_helper(ptr_res[i],ptr_A[i]);

            ptr_A   += A_ld;
            ptr_res += res_ld;
        }

        res.set_struct(md::predefined_struct::get_tril(A.get_struct(), d, is_real_matrix(A)));
        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_inplace(matcl::Matrix& ret, const M& A, Integer d)
    {
        Integer r = A.rows(), c = A.cols();

        ret_type res(A);

        using value_type    = typename ret_type::value_type;
        using value_type_in = typename M::value_type;

        value_type Z = md::default_value<value_type>(A.get_type());

        value_type* ptr_res         = res.ptr();
        Integer res_ld              = A.ld();

        for (Integer j = 0; j < c; ++j)
        {
            Integer mjmdrp1 = (j - d < r) ? j - d : r;
            Integer i;
            for (i = 0; i < mjmdrp1; ++i)
                mrd::reset_helper(ptr_res[i],Z);

            ptr_res += res_ld;
        }

        res.set_struct(md::predefined_struct::get_tril(A.get_struct(), d, is_real_matrix(A)));
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class ret_type,class M>
struct eval_triu_helper<ret_type,M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& A, Integer d, bool rvalue)
    {
        Integer r = A.rows(), c = A.cols();		
    
        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (d <= 0 && A.get_struct().is_triu())
        {
            ret = matcl::Matrix(A,false);
            return;
        };        

        if (d <= -1 && A.get_struct().is_hessu())
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (rvalue == true && A.is_unique() == true)
            return eval_inplace(ret, A, d);

        ret_type res(A.get_type(),r, c);

        using value_type    = typename ret_type::value_type;
        using value_type_in = typename M::value_type;
        value_type Z        = md::default_value<value_type>(A.get_type());

        value_type* ptr_res         = res.ptr();
        const value_type_in* ptr_A  = A.ptr();
        Integer A_ld                = A.ld();
        Integer res_ld              = res.ld();

        for (Integer j = 0; j < c; ++j)
        {
            Integer mjmdr = (j - d < r - 1) ? j - d : r - 1;
            Integer i;
            for (i = 0; i <= mjmdr; ++i)
                mrd::reset_helper(ptr_res[i],ptr_A[i]);

            for (; i < r; ++i)
                mrd::reset_helper(ptr_res[i],Z);

            ptr_A   += A_ld;
            ptr_res += res_ld;
        }

        res.set_struct(md::predefined_struct::get_triu(A.get_struct(), d, is_real_matrix(A)));

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_inplace(matcl::Matrix& ret, const M& A, Integer d)
    {
        Integer r = A.rows(), c = A.cols();		

        ret_type res(A);

        using value_type    = typename ret_type::value_type;
        using value_type_in = typename M::value_type;
        value_type Z        = md::default_value<value_type>(A.get_type());

        value_type* ptr_res = res.ptr();
        Integer res_ld      = res.ld();

        for (Integer j = 0; j < c; ++j)
        {
            Integer mjmdr   = (j - d < r - 1) ? j - d : r - 1;
            Integer i       = std::max(mjmdr + 1, 0);

            for (; i < r; ++i)
                mrd::reset_helper(ptr_res[i],Z);

            ptr_res += res_ld;
        }

        res.set_struct(md::predefined_struct::get_triu(A.get_struct(), d, is_real_matrix(A)));

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class M>
struct eval_select_band_helper<M,struct_dense>
{
    static void eval(matcl::Matrix& ret, const M& mat, Integer fd0, Integer ld0)
    {
        Integer r   = mat.rows();
        Integer c   = mat.cols();

        using val_type  = typename M::value_type;
        using ret_type  = raw::Matrix<val_type,struct_banded>;

        if (r == 0 || c == 0)
        {
            ret = matcl::Matrix(ret_type(mat.get_type(),r,c,0,0),false);
            return;
        };

        Integer fd, ld;
        estim_bandwidth::eval(mat, fd0, ld0, fd, ld);

        if (ld < fd)
        {
            using Mat_sparse    = raw::Matrix<val_type, struct_sparse>;
            Mat_sparse res(mat.get_type(), r, c);
            ret = matcl::Matrix(res, false);
            return;
        };

        ret_type res(mat.get_type(),r,c,fd,ld);

        const val_type* ptr = mat.ptr();

        val_type* ptr_ret   = res.rep_ptr();
        Integer mat_ld      = mat.ld();
        Integer res_ld      = res.ld();

        if (fd > 0)
        {
            val_type Z      = md::default_value<val_type>(res.get_type());

            for (Integer d = 0; d < fd; ++d)
            {
                Integer s               = res.diag_length(d);
                val_type* ptr_rd        = ptr_ret + res.first_elem_diag(d);

                for (Integer i = 0; i < s; ++i)
                {
                    mrd::reset_helper(ptr_rd[0], Z);
                    ptr_rd              += res_ld;
                };
            };
        };

        if (ld < 0)
        {
            val_type Z      = md::default_value<val_type>(res.get_type());

            for (Integer d = ld + 1; d <= 0; ++d)
            {
                Integer s               = res.diag_length(d);
                val_type* ptr_rd        = ptr_ret + res.first_elem_diag(d);

                for (Integer i = 0; i < s; ++i)
                {
                    mrd::reset_helper(ptr_rd[0], Z);
                    ptr_rd              += res_ld;
                };
            };
        };

        for (Integer j = 0; j < c; ++j)
        {
            Integer first_row   = res.first_row(j);
            Integer first_pos   = res.first_elem_pos(j) - first_row;
            Integer end         = std::min(j-fd,r-1);

            for (Integer k = std::max(j-ld,0); k <= end; ++k)
                mrd::reset_helper(ptr_ret[first_pos+k],ptr[k]);

            ptr         += mat_ld;
            ptr_ret     += res_ld;
        };

        struct_flag ms  = mat.get_struct();
        struct_flag sf  = ms;
        bool is_mat_sq  = mat.rows() == mat.cols();
        sf              = md::predefined_struct::get_tril(sf, ld, is_real_matrix(mat));
        sf              = md::predefined_struct::get_triu(sf, fd, is_real_matrix(mat));

        if (ms.is_symmetric(is_mat_sq, is_real_matrix(mat)) && fd <= 0 && -fd == ld)
            sf.add(predefined_struct_type::sym);

        if (ms.is_hermitian(is_mat_sq, is_real_matrix(mat)) && fd <= 0 && -fd == ld)
            sf.add(predefined_struct_type::her);

        res.set_struct(sf);
        ret = matcl::Matrix(res,false);
    };
};

template<class ret_type,class M>
struct eval_tril_helper<ret_type,M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& A, Integer D, bool rvalue)
    {
        using value_type = typename M::value_type;

        Integer Ar = A.rows(), Ac = A.cols(), nnz = A.nnz();

        if (D >= 0 && A.get_struct().is_tril())
        {
            ret = matcl::Matrix(A,false);
            return;
        };        

        if (D >= 1 && A.get_struct().is_hessl())
        {
            ret = matcl::Matrix(A,false);
            return;            
        };

        if (nnz == 0 || Ar == 0 || Ac == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (rvalue == true && A.is_unique() == true)
            return eval_inplace(ret, A, D);

        const details::sparse_ccs<value_type>& Ad = A.rep();

        ret_type res(A.get_type(),Ar,Ac,nnz);
        details::sparse_ccs<value_type>& d = res.rep();

        const Integer* Ad_r = Ad.ptr_r();
        const Integer* Ad_c = Ad.ptr_c();
        const value_type* Ad_x = Ad.ptr_x();

        Integer* d_r = d.ptr_r();
        Integer* d_c = d.ptr_c();
        value_type* d_x = d.ptr_x();		

        Integer pos_r = 0;

        for (Integer j = 0; j < Ac; ++j)
        {
            d_c[j] = pos_r;
            Integer row_start = j - D;
            for (Integer i = Ad_c[j]; i < Ad_c[j+1]; ++i)
            {
                Integer p = Ad_r[i];

                if (p < row_start)
                    continue;

                d_r[pos_r] = p;
                mrd::reset_helper(d_x[pos_r],Ad_x[i]);
                ++pos_r;
            }
        };
        d_c[Ac] = pos_r;

        d.add_memory(-1);
        res.set_struct(md::predefined_struct::get_tril(A.get_struct(), D, is_real_matrix(A)));

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_inplace(matcl::Matrix& ret, const M& A, Integer D)
    {
        using value_type = typename M::value_type;

        Integer Ac = A.cols();

        ret_type res(A);
        details::sparse_ccs<value_type>& d = res.rep();

        Integer* d_r    = d.ptr_r();
        Integer* d_c    = d.ptr_c();
        value_type* d_x = d.ptr_x();		

        Integer pos_r   = d.offset();
        Integer j       = 0;

        for (; j < Ac; ++j)
        {
            Integer row_start   = j - D;
            Integer i           = d_c[j];
            Integer i_last      = d_c[j+1];

            for (; i < i_last; ++i)
            {
                Integer p = d_r[i];

                if (p < row_start)
                    goto lab_modify;

                ++pos_r;
            }

            continue;

          lab_modify:

            for (; i < i_last; ++i)
            {
                Integer p = d_r[i];

                if (p < row_start)
                    continue;

                d_r[pos_r] = p;
                mrd::reset_helper(d_x[pos_r],d_x[i]);
                ++pos_r;
            }
            
            ++j;
            break;
        };      

        for (; j < Ac; ++j)
        {            
            Integer row_start   = j - D;
            Integer i_start     = d_c[j];
            Integer i_last      = d_c[j+1];

            d_c[j]              = pos_r;

            for (Integer i = i_start; i < i_last; ++i)
            {
                Integer p = d_r[i];

                if (p < row_start)
                    continue;

                d_r[pos_r] = p;
                mrd::reset_helper(d_x[pos_r],d_x[i]);
                ++pos_r;
            }
        };

        d_c[Ac] = pos_r;

        d.add_memory(-1);
        res.set_struct(md::predefined_struct::get_tril(A.get_struct(), D, is_real_matrix(A)));

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class ret_type,class M>
struct eval_triu_helper<ret_type,M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& A, Integer D, bool rvalue)
    {
        using value_type = typename M::value_type;

        Integer Ar = A.rows(), Ac = A.cols(), nnz = A.nnz();		

        if (D <= 0 && A.get_struct().is_triu())
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (D <= -1 && A.get_struct().is_hessu())
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (nnz == 0 || Ar == 0 || Ac == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (rvalue == true && A.is_unique() == true)
            return eval_inplace(ret, A, D);

        const details::sparse_ccs<value_type>& Ad = A.rep();

        ret_type res(A.get_type(),Ar,Ac,nnz);
        details::sparse_ccs<value_type>& d = res.rep();

        const Integer* Ad_r = Ad.ptr_r();
        const Integer* Ad_c = Ad.ptr_c();
        const value_type* Ad_x = Ad.ptr_x();

        Integer* d_r = d.ptr_r();
        Integer* d_c = d.ptr_c();
        value_type* d_x = d.ptr_x();		

        Integer pos_r = 0;

        for (Integer j = 0; j < Ac; ++j)
        {
            d_c[j] = pos_r;
            Integer row_start = j - D;
            for (Integer i = Ad_c[j]; i < Ad_c[j+1]; ++i)
            {
                Integer p = Ad_r[i];

                if (p > row_start)
                    break;

                d_r[pos_r] = p;
                mrd::reset_helper(d_x[pos_r],Ad_x[i]);
                ++pos_r;
            }
        };
        d_c[Ac] = pos_r;

        d.add_memory(-1);
        res.set_struct(md::predefined_struct::get_triu(A.get_struct(), D, is_real_matrix(A)));

        ret = matcl::Matrix(res,true);
        return;
    };

    static void eval_inplace(matcl::Matrix& ret, const M& A, Integer D)
    {
        using value_type = typename M::value_type;

        Integer Ac = A.cols();
        
        ret_type res(A);
        details::sparse_ccs<value_type>& d = res.rep();

        Integer* d_r    = d.ptr_r();
        Integer* d_c    = d.ptr_c();
        value_type* d_x = d.ptr_x();		

        Integer pos_r   = d.offset();
        Integer j       = 0;

        for (; j < Ac; ++j)
        {
            Integer row_start   = j - D;
            Integer i_start     = d_c[j];
            Integer i_end       = d_c[j+1];
            Integer i           = i_start;

            for (; i < i_end; ++i)
            {
                Integer p = d_r[i];

                if (p > row_start)
                {
                    ++j;
                    goto modif_lab;
                }

                ++pos_r;
            }
        };

      modif_lab:
        for (; j < Ac; ++j)
        {            
            Integer row_start   = j - D;
            Integer i_start     = d_c[j];
            Integer i_end       = d_c[j+1];

            d_c[j]              = pos_r;

            for (Integer i = i_start; i < i_end; ++i)
            {
                Integer p       = d_r[i];

                if (p > row_start)
                    break;

                d_r[pos_r] = p;
                mrd::reset_helper(d_x[pos_r],d_x[i]);
                ++pos_r;
            }
        };
        
        d_c[Ac] = pos_r;

        d.add_memory(-1);
        res.set_struct(md::predefined_struct::get_triu(A.get_struct(), D, is_real_matrix(A)));

        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class M>
struct eval_select_band_helper<M,struct_sparse>
{
    static void eval(matcl::Matrix& ret, const M& mat, Integer fd0, Integer ld0)
    {
        Integer r = mat.rows();
        Integer c = mat.cols();

        using val_type  = typename M::value_type;
        using ret_type  = raw::Matrix<val_type,struct_banded>;

        if (mat.nnz() == 0)
        {
            ret = matcl::Matrix(mat,false);
            return;
        };

        Integer fd, ld;
        estim_bandwidth::eval(mat, fd0, ld0, fd, ld);

        if (ld < fd)
        {
            using Mat_sparse    = raw::Matrix<val_type, struct_sparse>;
            Mat_sparse res(mat.get_type(), r, c);
            ret = matcl::Matrix(res, false);
            return;
        };

        ret_type res(mat.get_type(),r,c,fd,ld);

        val_type* ptr_res   = res.rep_ptr();
        Integer res_ld      = res.ld();
        Integer res_size    = res.impl_size();
        val_type Z          = md::default_value<val_type>(res.get_type());

        for (Integer i = 0; i < res_size; ++i)
            mrd::reset_helper(ptr_res[i],Z);

        const mrd::sparse_ccs<val_type>& Ad = mat.rep();
        const Integer * Ad_c = Ad.ptr_c();
        const Integer * Ad_r = Ad.ptr_r();
        const val_type* Ad_x = Ad.ptr_x();

        for (Integer j = 0; j < c; ++j)
        {
            Integer first_row   = std::max(j-ld,0);
            Integer last_row    = std::min(j-fd,r-1);
            Integer first_pos   = res.first_elem_pos(j) - res.first_row(j);

            for (Integer k = Ad_c[j]; k < Ad_c[j+1]; ++k)
            {
                Integer r2      = Ad_r[k];
                
                if (r2 < first_row || r2 > last_row)
                    continue;

                mrd::reset_helper(ptr_res[first_pos+Ad_r[k]],Ad_x[k]);
            };

            ptr_res += res_ld;
        };

        struct_flag ms  = mat.get_struct();
        struct_flag sf  = ms;
        bool is_mat_sq  = mat.rows() == mat.cols();
        sf              = md::predefined_struct::get_tril(sf, ld, is_real_matrix(mat));
        sf              = md::predefined_struct::get_triu(sf, fd, is_real_matrix(mat));

        if (ms.is_symmetric(is_mat_sq, is_real_matrix(mat)) && fd <= 0 && -fd == ld)
            sf.add(predefined_struct_type::sym);

        if (ms.is_hermitian(is_mat_sq, is_real_matrix(mat)) && fd <= 0 && -fd == ld)
            sf.add(predefined_struct_type::her);

        res.set_struct(sf);
        ret = matcl::Matrix(res,false);
    };
};

template<class ret_type,class M>
struct eval_tril_helper<ret_type,M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& A, Integer D, bool)
    {
        if (A.rows() == 0 || A.cols() == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (D >= 0 && A.get_struct().is_tril())
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (D >= 1 && A.get_struct().is_hessl())
        {
            ret = matcl::Matrix(A,false);
            return;            
        };

        Integer fd  = A.first_diag();
        Integer ld  = A.last_diag();

        if (D >= ld)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        using VT    = typename M::value_type;
        using VTR   = typename ret_type::value_type;

        if (D < fd)
        {
            auto tmp = matcl::details::zero_matrix<VTR,struct_sparse>
                            ::eval(A.get_type(),A.rows(),A.cols());
            ret = matcl::Matrix(tmp,false);
            return;
        };

        if (D >= 0)
        {
            D   = std::min(D,ld);
            ret = matcl::Matrix(A.make_view(1,A.rows(),A.cols(),fd,D),false);
            return;
        };

        //TODO: make view when generalized views are allowed

        VTR Z       = md::default_value<VTR>(A.get_type()); 

        ret_type res(A.get_type(),A.rows(), A.cols(), fd, 0);

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();

        for (Integer i = fd; i <= 0; ++i)
        {
            VTR* ptr_res    = res.rep_ptr() + res.first_elem_diag(i);
            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(i);
            Integer s		= res.diag_length(i);

            if (i > D)
            {
                for (Integer k = 1; k <= s; ++k)
                {
                    mrd::reset_helper(ptr_res[0],Z);
                    ptr_res += res_ld;
                };
            }
            else
            {
                for (Integer k = 1; k <= s; ++k)
                {
                    mrd::reset_helper(ptr_res[0],ptr_A[0]);
                    ptr_res += res_ld;
                    ptr_A   += A_ld;
                };
            };
        };

        res.set_struct(md::predefined_struct::get_tril(A.get_struct(), D, is_real_matrix(A)));
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class ret_type,class M>
struct eval_triu_helper<ret_type,M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& A, Integer D, bool)
    {	
        if (A.rows() == 0 || A.cols() == 0)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (D <= 0 && A.get_struct().is_triu())
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        if (D <= -1 && A.get_struct().is_hessu())
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        Integer fd  = A.first_diag();
        Integer ld  = A.last_diag();

        if (D <= fd)
        {
            ret = matcl::Matrix(A,false);
            return;
        };

        using VTR   = typename ret_type::value_type;
        using VT    = typename M::value_type;

        if (D > ld)
        {
            auto tmp = matcl::details::zero_matrix<VTR,struct_sparse>
                            ::eval(A.get_type(),A.rows(),A.cols());
            ret = matcl::Matrix(tmp,false);
            return;
        };        

        if (D <= 0)
        {
            D   = std::max(D,fd);
            ret = matcl::Matrix(A.make_view(1,A.rows(),A.cols(),D,ld),false);
            return;
        };

        //TODO: make view when generalized views are allowed

        VTR Z               = md::default_value<VTR>(A.get_type());

        ret_type res(A.get_type(), A.rows(), A.cols(), 0, ld);

        Integer res_ld      = res.ld();
        Integer A_ld        = A.ld();

        for (Integer i = 0; i <= ld; ++i)
        {
            VTR* ptr_res    = res.rep_ptr() + res.first_elem_diag(i);
            const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(i);
            Integer s		= res.diag_length(i);

            if (i < D)
            {
                for (Integer k = 1; k <= s; ++k)
                {
                    mrd::reset_helper(ptr_res[0], Z);
                    ptr_res += res_ld;
                };
            }
            else
            {
                for (Integer k = 1; k <= s; ++k)
                {
                    mrd::reset_helper(ptr_res[0],ptr_A[0]);
                    ptr_res += res_ld;
                    ptr_A   += A_ld;
                };
            };
        };

        res.set_struct(md::predefined_struct::get_triu(A.get_struct(), D, is_real_matrix(A)));
        ret = matcl::Matrix(res,true);
        return;
    };
};

template<class M>
struct eval_select_band_helper<M,struct_banded>
{
    static void eval(matcl::Matrix& ret, const M& mat, Integer fd0, Integer ld0)
    {                
        Integer fd, ld;
        estim_bandwidth::eval(mat, fd0, ld0, fd, ld);

        if (ld < fd)
        {
            using val_type      = typename M::value_type;
            using Mat_sparse    = raw::Matrix<val_type, struct_sparse>;
            Mat_sparse res(mat.get_type(), mat.rows(),mat.cols());
            ret = matcl::Matrix(res, false);
            return;
        }

        if (fd <= 0 && ld >= 0)
        {
            ret = matcl::Matrix(mat.make_view(1,mat.rows(),mat.cols(),fd,ld),false);
            return;
        };

        //TODO: make view
        Integer r           = mat.rows();
        Integer c           = mat.cols();

        using val_type      = typename M::value_type;

        M res(mat.get_type(),r,c,fd,ld);

        const val_type* ptr = mat.rep_ptr();

        val_type* ptr_ret   = res.rep_ptr();
        Integer mat_ld      = mat.ld();
        Integer res_ld      = res.ld();

        if (fd > 0)
        {
            val_type Z      = md::default_value<val_type>(res.get_type());

            for (Integer d = 0; d < fd; ++d)
            {
                Integer s               = res.diag_length(d);
                val_type* ptr_rd        = ptr_ret + res.first_elem_diag(d);

                for (Integer i = 0; i < s; ++i)
                {
                    mrd::reset_helper(ptr_rd[0], Z);
                    ptr_rd              += res_ld;
                };
            };
        };
      
        if (ld < 0)
        {
            val_type Z      = md::default_value<val_type>(res.get_type());

            for (Integer d = ld + 1; d <= 0; ++d)
            {
                Integer s               = res.diag_length(d);
                val_type* ptr_rd        = ptr_ret + res.first_elem_diag(d);

                for (Integer i = 0; i < s; ++i)
                {
                    mrd::reset_helper(ptr_rd[0], Z);
                    ptr_rd              += res_ld;
                };
            };
        };

        for (Integer d = fd; d <= ld; ++d)
        {
            Integer s               = res.diag_length(d);
            val_type* ptr_rd        = ptr_ret + res.first_elem_diag(d);
            const val_type* ptr_d   = ptr + mat.first_elem_diag(d);

            for (Integer i = 0; i < s; ++i)
            {
                mrd::reset_helper(ptr_rd[0], ptr_d[0]);
                ptr_rd              += res_ld;
                ptr_d               += mat_ld;
            };
        };

        struct_flag ms  = mat.get_struct();
        struct_flag sf  = ms;
        bool is_mat_sq  = mat.rows() == mat.cols();
        sf              = md::predefined_struct::get_tril(sf, ld, is_real_matrix(mat));
        sf              = md::predefined_struct::get_triu(sf, fd, is_real_matrix(mat));

        if (ms.is_symmetric(is_mat_sq, is_real_matrix(mat)) && fd <= 0 && -fd == ld)
            sf.add(predefined_struct_type::sym);

        if (ms.is_hermitian(is_mat_sq, is_real_matrix(mat)) && fd <= 0 && -fd == ld)
            sf.add(predefined_struct_type::her);

        res.set_struct(sf);
        ret = matcl::Matrix(res,false);
    }
};

template<class M>
void
details::manip_tr_helper<M>::eval_tril(matcl::Matrix& ret, const M& A, Integer d, bool rvalue)
{
    using struct_type = typename M::struct_type;
    return eval_tril_helper<ret_type_tril,M,struct_type>::eval(ret, A,d, rvalue);
};

template<class M>
void
details::manip_tr_helper<M>::eval_triu(matcl::Matrix& ret, const M& A, Integer d, bool rvalue)
{
    using struct_type = typename M::struct_type;
    return eval_triu_helper<ret_type_triu,M,struct_type>::eval(ret, A,d, rvalue);
};

template<class M>
void
details::manip_tr_helper<M>::eval_select_band(matcl::Matrix& ret, const M& A, Integer fd, Integer ld)
{
    using struct_type = typename M::struct_type;
    return eval_select_band_helper<M,struct_type>::eval(ret, A, fd, ld);
};

template<class V>
Integer get_ld(const Matrix<V,struct_sparse>& A, Integer min, bool use_flag)
{
    if (A.nnz() == 0)
        return 0;

    if (use_flag == true)
    {
        struct_flag::diag_type dt = A.get_struct().get_ldiags();
        if (dt == struct_flag::zero)
            return 0;
    };

    Integer r = A.rows(), c = A.cols();
    using val_type = V;

    const details::sparse_ccs<val_type>& Ad = A.rep();
    const Integer * Ad_c = Ad.ptr_c();
    const Integer * Ad_r = Ad.ptr_r();
    const val_type* Ad_x = Ad.ptr_x();

    Integer ld = 0;
    for (Integer j = 0; j < c; ++j)
    {
        Integer k = Ad_c[j+1]-1;
        while(k >= Ad_c[j])
        {
            if (mrd::is_zero(Ad_x[k]))
            {
                --k;
                continue;
            };

            Integer pos = Ad_r[k] - j;
            ld = std::max(ld,pos);

            if (r - 1 - j <= ld)
                return ld;

            break;
        };

        if (min >= 0 && ld > min)
            return ld;
    };

    return ld;
};

template<class V>
Integer get_ud(const Matrix<V,struct_sparse>& A, Integer min, bool use_flag)
{
    Integer c = A.cols();

    if (A.nnz() == 0)
        return 0;

    if (use_flag)
    {
        struct_flag::diag_type dt = A.get_struct().get_udiags();
        if (dt == struct_flag::zero)
            return 0;
    };

    using val_type = V;

    const details::sparse_ccs<val_type>& Ad = A.rep();
    const Integer * Ad_c = Ad.ptr_c();
    const Integer * Ad_r = Ad.ptr_r();
    const val_type* Ad_x = Ad.ptr_x();

    Integer ud = 0;
    for (Integer j = c-1; j >= 0; --j)
    {
        Integer k = Ad_c[j];
        while(k < Ad_c[j+1])
        {
            if (mrd::is_zero(Ad_x[k]))
            {
                ++k;
                continue;
            };
            Integer pos = j - Ad_r[k];
            ud = std::max(ud,pos);

            if (j <= ud)
                return ud;

            break;
        };

        if (min >= 0 && ud > min)
            return ud;
    };
    return ud;
};

template<class V>
Integer get_ld(const Matrix<V,struct_dense>& A, Integer min, bool use_flag)
{
    Integer r = A.rows();
    Integer c = A.cols();

    using val_type = V;

    bool isq = false;

    if (use_flag)
    {
        struct_flag::diag_type dt = A.get_struct().get_ldiags();
        
        if (dt == struct_flag::zero)
            return 0;

        if (dt != struct_flag::general)
            isq = true;
    };

    const val_type* ptr = A.ptr();
    Integer A_ld        = A.ld();

    if (isq)
    {
        for (Integer j = 0; j < c; ++j)
        {
            Integer i = j+1;
            if (i >= r)
                return 0;

            if (mrd::is_zero(ptr[i]) == false)
                return 1;

            ptr += A_ld;
        };

        return 0;
    };

    Integer ld = 0;
    for (Integer j = 0; j < c; ++j)
    {
        Integer k = r-1;
        while(k >= 0)
        {
            if (mrd::is_zero(ptr[k]))
            {
                --k;
                continue;
            };

            Integer pos = k - j;
            ld = std::max(ld,pos);

            if (r - 1 - j <= ld)
                return ld;

            break;
        };

        ptr += A_ld;

        if (min >= 0 && ld > min)
            return ld;
    };

    return ld;
};

template<class V>
Integer get_ud(const Matrix<V,struct_dense>& A, Integer min, bool use_flag)
{
    Integer r       = A.rows();
    Integer c       = A.cols();
    Integer A_ld    = A.ld();
    bool isq        = false;

    using val_type  = V;

    if (use_flag)
    {
        struct_flag::diag_type dt = A.get_struct().get_udiags();
        if (dt == struct_flag::zero)
            return 0;

        if (dt != struct_flag::general)
            isq = true;
    };

    if (isq)
    {
        const val_type* ptr = A.ptr() + A_ld;

        for (Integer j = 1; j < c; ++j)
        {
            Integer i = j-1;
            if (i >= r)
                return 0;

            if (mrd::is_zero(ptr[i]) == false)
                return 1;

            ptr += A_ld;
        };
        return 0;
    };

    const val_type* ptr = A.ptr() + imult(c-1,A_ld);

    Integer ud = 0;
    for (Integer j = c-1; j >= 0; --j)
    {
        Integer k = 0;
        while(k < r)
        {
            if (mrd::is_zero(ptr[k]))
            {
                ++k;
                continue;
            };

            Integer pos = j - k;
            ud = std::max(ud,pos);

            if (j <= ud)
                return ud;

            break;
        };

        ptr -= A_ld;
        if (min >= 0 && ud > min)
            return ud;
    };

    return ud;
};

template<class V>
Integer get_ld(const Matrix<V,struct_banded>& A, Integer min, bool use_flag)
{    
    (void)min;

    Integer r       = A.rows();
    Integer c       = A.cols();
    Integer ldn     = A.number_subdiagonals();
    Integer A_ld    = A.ld();

    if (ldn == 0 || r == 0 || c == 0)
        return 0;

    if (use_flag)
    {
        struct_flag::diag_type dt = A.get_struct().get_ldiags();
        if (dt == struct_flag::zero)
            return 0;

        if (dt != struct_flag::general)
            ldn = 1;
    };

    while(ldn > 0)
    {
        Integer	s       = A.diag_length(-ldn);
        const V* ptr    = A.rep_ptr() + A.first_elem_diag(-ldn);;

        for (Integer i = 0; i < s; ++i)
        {
            if (mrd::is_zero(ptr[0]) == false)
                return ldn;

            ptr         += A_ld;
        };

        --ldn;
    };

    return ldn;
};

template<class V>
Integer get_ud(const Matrix<V,struct_banded>& A, Integer min, bool use_flag)
{
    (void)min;

    Integer udn     = A.number_superdiagonals();
    Integer r       = A.rows();
    Integer c       = A.cols();
    Integer A_ld    = A.ld();

    if (udn == 0 || r == 0 || c == 0)
        return 0;

    if (use_flag)
    {
        struct_flag::diag_type dt = A.get_struct().get_udiags();
        if (dt == struct_flag::zero)
            return 0;

        if (dt != struct_flag::general)
            udn = 1;
    };

    while(udn > 0)
    {
        Integer	s       = A.diag_length(udn);
        const V* ptr    = A.rep_ptr() + A.first_elem_diag(udn);

        for (Integer i = 0; i < s; ++i)
        {
            if (mrd::is_zero(ptr[0]) == false)
                return udn;

            ptr += A_ld;
        };

        --udn;
    };

    return udn;
};

template<class V>
Integer raw::nnz_total(const Matrix<V,struct_dense>& A)
{
    Integer r       = A.rows();
    Integer c       = A.cols();
    Integer A_ld    = A.ld();
    const V* ptr    = A.ptr();

    Integer nz      = 0;

    for (Integer j = 0; j < c; ++j)
    {
        for (Integer i = 0; i < r; ++i)
        {
            if (mrd::is_zero<V>(ptr[i]) == false)
                ++nz;
        };

        ptr         += A_ld;
    };

    return nz;
};

template<class V>
Integer raw::nnz_total(const Matrix<V,struct_sparse>& A)
{
    if (A.nnz() == 0)
        return 0;

    Integer c = A.cols();

    const details::sparse_ccs<V>& Ad = A.rep();
    const Integer * Ad_c    = Ad.ptr_c();
    const V* Ad_x           = Ad.ptr_x();

    Integer k0              = Ad_c[0];
    Integer k1              = Ad_c[c];
    Integer nz              = 0;

    for (Integer i = k0; i < k1; ++i)
    {
        if (mrd::is_zero<V>(Ad_x[i]) == false)
            ++nz;
    };

    return nz;
};

template<class V>
Integer raw::nnz_total(const Matrix<V,struct_banded>& A)
{
    Integer r       = A.rows();
    Integer c       = A.cols();
    Integer A_ld    = A.ld();
    const V* ptr    = A.rep_ptr();

    if (r == 0 || c == 0)
        return 0;

    Integer nz      = 0;

    for (Integer j = 0; j < c; ++j)
    {
        Integer fr      = A.first_row(j);
        Integer lr      = A.last_row(j);
        Integer pos     = A.first_elem_pos(j);

        for(Integer i = fr; i <= lr; ++i, ++pos)
        {
            if (mrd::is_zero<V>(ptr[pos]) == false)
                ++nz;
        };

        ptr             += A_ld;
    };

    return nz;
};

template<class V> struct is_real_type   { static const bool value = false; };
template<> struct is_real_type<Integer> { static const bool value = true; };
template<> struct is_real_type<Real>    { static const bool value = true; };
template<> struct is_real_type<Float>   { static const bool value = true; };

template<class T>
struct neq_helper_tol{};

template<>
struct neq_helper_tol<Integer>
{
    static const bool eval(Integer a, Integer b, Real)
    {
        return a != b;
    };
};

template<>
struct neq_helper_tol<Float>
{
    static const bool eval(Float a, Float b, Real tol0)
    {
        if (tol0 > 0.0)
        {
            Float tol = matcl::eps(a) * Float(tol0);
            return  mrd::abs_helper<Float>::eval(a - b) > tol;
        }
        else
        {
            return a != b;
        }
    };
};

template<>
struct neq_helper_tol<Real>
{
    static const bool eval(Real a, Real b, Real tol0)
    {
        if (tol0 > 0)
        {
            Real tol = matcl::eps(a) * tol0;
            return mrd::abs_helper<Real>::eval(a - b) > tol;
        }
        else
        {
            return a != b;
        }
    };
};

template<>
struct neq_helper_tol<Complex>
{
    static const bool eval(const Complex& a, const Complex& b, Real tol0)
    {
        if (tol0 > 0)
        {
            Real tol    = matcl::eps(a) * tol0;
            Complex dif = md::minus_c(a,b);
            return mrd::abs_helper<Complex>::eval(dif) > tol;
        }
        else
        {
            return a != b;
        }
    };
};

template<>
struct neq_helper_tol<Float_complex>
{
    static const bool eval(const Float_complex& a, const Float_complex& b, Real tol0)
    {
        if (tol0 > 0.0)
        {
            Real tol    = matcl::eps(a) * tol0;
            Float_complex dif = md::minus_c(a,b);
            return mrd::abs_helper<Float_complex>::eval(dif) > tol;
        }
        else
        {
            return a != b;
        }
    };
};

template<>
struct neq_helper_tol<Object>
{
    static const bool eval(const Object& a, const Object& b, Real tol0)
    {
        if (tol0 > 0)
        {
            Object tol  = dynamic::operator*(matcl::eps(a), Object(tol0));
            Object dif  = dynamic::operator-(a,b); 
            Object adif = mrd::abs_helper<Object>::eval(dif);
            return (bool)dynamic::operator>(adif, tol);
        }
        else
        {
            return (bool)dynamic::operator!=(a, b);
        }
    };
};

template<class V>
bool is_sym(const Matrix<V,struct_dense>& A, Real tol, bool use_flag)
{
    if (A.rows() != A.cols())
        return false;

    if (use_flag)
    {
        if (A.get_struct().is_symmetric(true, is_real_matrix(A)) == true)
            return true;
    };    

    Integer s       = A.rows();
    const V* ptr_0  = A.ptr();
    const V* ptr    = ptr_0;
    Integer A_ld    = A.ld();

    for (Integer i = 0; i < s; ++i)
    {
        const V* ptr_s = ptr_0 + i + i*A_ld;
        for (Integer k = i; k < s; ++k)
        {
            const V& val1 = ptr[k];
            const V& val2 = ptr_s[0];

            if (neq_helper_tol<V>::eval(val1,val2, tol) 
                && (bool)mrd::isnan_helper<V>::eval(val1) == false)
            {
                return false;
            };

            ptr_s += A_ld;
        };

        ptr += A_ld;
    };
    return true;
};

template<class V>
bool is_her(const Matrix<V,struct_dense>& A, Real tol, bool use_flag)
{
    if (A.rows() != A.cols())
        return false;

    if (use_flag)
    {
        if (A.get_struct().is_hermitian(true, is_real_matrix(A)))
            return true;
    };

    Integer s       = A.rows();
    const V* ptr_0  = A.ptr();
    const V* ptr    = ptr_0;
    Integer A_ld    = A.ld();

    for (Integer i = 0; i < s; ++i)
    {
        const V* ptr_s = ptr_0 + i + i*A_ld;
        for (Integer k = i; k < s; ++k)
        {
            const V& val1   = ptr[k];
            const V& val2   = ptr_s[0];
            V val2_conj     = mrd::conj_helper<V>::eval(val2);

            if (neq_helper_tol<V>::eval(val1 , val2_conj, tol)
                && (bool)mrd::isnan_helper<V>::eval(val1) == false)
            {
                return false;
            };

            ptr_s += A_ld;
        };

        ptr += A_ld;
    };
    return true;
};

template<class V>
bool is_sym(const Matrix<V,struct_banded>& A, Real tol, bool use_flag)
{
    if (A.rows() != A.cols())
        return false;

    if (use_flag)
    {
        if (A.get_struct().is_symmetric(true, is_real_matrix(A)) == true)
            return true;
    };

    Integer ld      = raw::get_ld(A,-1,use_flag);
    Integer ud      = raw::get_ud(A,-1,use_flag);
    Integer A_ld    = A.ld();

    if (ld != ud)
        return false;

    for (Integer i = 1; i <= ld; ++i)
    {
        Integer s       = A.diag_length(i);
        const V* ptr_l  = A.rep_ptr() + A.first_elem_diag(-i);
        const V* ptr_u  = A.rep_ptr() + A.first_elem_diag(i);

        for (Integer j = 0; j < s; ++j)
        {
            const V& val1 = ptr_l[0];
            const V& val2 = ptr_u[0];

            if (neq_helper_tol<V>::eval(val1 , val2, tol) 
                && (bool)mrd::isnan_helper<V>::eval(val1) == false)
            {
                return false;
            };

            ptr_l += A_ld;
            ptr_u += A_ld;
        };
    };

    return true;
};

template<class V>
bool is_her(const Matrix<V,struct_banded>& A, Real tol, bool use_flag)
{
    if (A.rows() != A.cols())
        return false;

    if (use_flag)
    {
        if (A.get_struct().is_hermitian(true, is_real_matrix(A)))
            return true;
    };

    Integer ld      = raw::get_ld(A,-1,use_flag);
    Integer ud      = raw::get_ud(A,-1,use_flag);

    Integer A_ld    = A.ld();

    if (ld != ud)
        return false;

    for (Integer i = 1; i <= ld; ++i)
    {
        Integer s       = A.diag_length(i);
        const V* ptr_l  = A.rep_ptr() + A.first_elem_diag(-i);
        const V* ptr_u  = A.rep_ptr() + A.first_elem_diag(i);

        for (Integer j = 0; j < s; ++j)
        {
            const V& val1   = ptr_l[0];
            const V& val2   = ptr_u[0];
            V val2_conj     = mrd::conj_helper<V>::eval(val2);

            if (neq_helper_tol<V>::eval(val1 , val2_conj, tol) 
                && (bool)mrd::isnan_helper<V>::eval(val1) == false)
            {
                return false;
            };

            ptr_l += A_ld;
            ptr_u += A_ld;
        };
    };

    if (md::is_complex<V>::value == false)
        return true;

    //check diagonal
    {
        Integer s       = A.diag_length(0);
        const V* ptr_l  = A.rep_ptr() + A.first_elem_diag(0);

        for (Integer j = 0; j < s; ++j)
        {
            const V& val1   = ptr_l[0];
            
            if (mrd::is_zero(mrd::imag_helper<V>::eval(val1)) == false
                    && (bool)mrd::isnan_helper<V>::eval(val1) == false)
            {
                return false;
            };

            ptr_l += A_ld;
        };
    };

    return true;
};

template<class V>
bool is_sym(const Matrix<V,struct_sparse>& A, Real tol, bool use_flag)
{
    if (A.rows() != A.cols())
        return false;

    if (use_flag)
    {
        if (A.get_struct().is_diag() == true)
            return true;
        if (A.get_struct().is_symmetric(true, is_real_matrix(A)) == true)
            return true;

        if (A.get_struct().is_hermitian(true, is_real_matrix(A)) == true)
        {
            bool is_real = is_real_type<V>::value;
            if (is_real == true)
            {
                return true;
            };
        }
    };

    if (A.nnz() == 0)
        return true;

    const Integer* ptr_c    = A.rep().ptr_c();
    const Integer* ptr_r    = A.rep().ptr_r();
    const V* ptr_x          = A.rep().ptr_x();
    Integer c               = A.cols();

    for (Integer j = 0; j < c; ++j)
    {
        for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
        {
            Integer r       = ptr_r[k];
            const V& val1   = ptr_x[k];

            V val2(A.rep().get_elem(j,r));

            if (neq_helper_tol<V>::eval(val1,val2, tol) 
                    && (bool)mrd::isnan_helper<V>::eval(val1) == false)
            {
                return false;
            };
        };
    };

    return true;
};

template<class V>
bool is_her(const Matrix<V,struct_sparse>& A, Real tol, bool use_flag)
{
    if (A.rows() != A.cols())
        return false;

    if (use_flag)
    {
        if (A.get_struct().is_hermitian(true, is_real_matrix(A)))
            return true;
    };

    if (A.nnz() == 0)
        return true;

    const Integer* ptr_c    = A.rep().ptr_c();
    const Integer* ptr_r    = A.rep().ptr_r();
    const V* ptr_x          = A.rep().ptr_x();
    Integer c               = A.cols();

    for (Integer j = 0; j < c; ++j)
    {
        for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
        {
            Integer r       = ptr_r[k];
            const V& val1   = ptr_x[k];

            if (r == j)
            {
                V val_conj  = mrd::conj_helper<V>::eval(val1);
                if (neq_helper_tol<V>::eval(val1, val_conj, tol) 
                    && (bool)mrd::isnan_helper<V>::eval(val1) == false)
                {
                    return false;
                };
                continue;
            };
            
            V val2      = A.rep().get_elem(j,r);
            V val2_conj = mrd::conj_helper<V>::eval(val2);

            if (neq_helper_tol<V>::eval(val1 , val2_conj, tol))
                return false;
        };
    };

    return true;
};

template<class Mat>
bool all_finite_helper<Mat>::eval(const Mat& mat)
{
    using V = typename Mat::value_type;
    using S = typename Mat::struct_type;
    return mrd::all_finite_helper_impl<V,S>::eval(mat);
};

template<class Val>
void details::manip_trans_helper_mat<Val>::eval_trans(Val* ret_ptr, Integer ld, const MP& m)
{
    return eval_trans_helper<MP,MP,struct_dense>::eval_mat(ret_ptr, ld, m);
}

template<class Val>
void details::manip_trans_helper_mat<Val>::eval_ctrans(Val* ret_ptr, Integer ld, const MP& m)
{
    using struct_type = struct_dense;
    return eval_ctrans_helper<MP,MP,true,struct_type>::eval_mat(ret_ptr,ld,m);
};

template<class Val, class Val_ret>
void details::manip_trans_reshaper_helper<Val,Val_ret>::eval_trans(ret_type& ret, const in_type& m, 
                            Integer max_row, Integer max_col, Integer ret_rows, Integer ret_cols)
{
    return eval_trans_helper<ret_type,in_type,struct_sparse>
        ::eval_reshape(ret, m, max_row, max_col,ret_rows, ret_cols);
};

template<class Val, class Val_ret>
void details::manip_trans_reshaper_helper<Val,Val_ret>::eval_ctrans(ret_type& ret, const in_type& m, 
                        Integer max_row, Integer max_col, Integer ret_rows, Integer ret_cols)
{
    static const bool is_compl  = md::is_complex<Val>::value;
    static const bool is_obj    = std::is_same<Val,Object>::value;

    return eval_ctrans_helper<ret_type,in_type,is_compl || is_obj,struct_sparse>
        ::eval_reshape(ret, m, max_row, max_col,ret_rows, ret_cols);
};

template integer_dense full(const matcl::raw::integer_sparse&);
template real_dense full(const matcl::raw::real_sparse&);
template float_dense full(const matcl::raw::float_sparse&);
template complex_dense full(const matcl::raw::complex_sparse&);
template float_complex_dense full(const matcl::raw::float_complex_sparse&);
template object_dense full(const matcl::raw::object_sparse&);

template Matrix<Integer,struct_sparse> sparse(const Matrix<Integer,struct_dense>&);
template Matrix<Real,struct_sparse> sparse(const Matrix<Real,struct_dense>&);
template Matrix<Float,struct_sparse> sparse(const Matrix<Float,struct_dense>&);
template Matrix<Complex,struct_sparse> sparse(const Matrix<Complex,struct_dense>&);
template Matrix<Float_complex,struct_sparse> sparse(const Matrix<Float_complex,struct_dense>&);
template Matrix<Object,struct_sparse> sparse(const Matrix<Object,struct_dense>&);

template Matrix<Integer,struct_sparse> sparse(const Matrix<Integer,struct_banded>&);
template Matrix<Real,struct_sparse> sparse(const Matrix<Real,struct_banded>&);
template Matrix<Float,struct_sparse> sparse(const Matrix<Float,struct_banded>&);
template Matrix<Complex,struct_sparse> sparse(const Matrix<Complex,struct_banded>&);
template Matrix<Float_complex,struct_sparse> sparse(const Matrix<Float_complex,struct_banded>&);
template Matrix<Object,struct_sparse> sparse(const Matrix<Object,struct_banded>&);

template Matrix<Integer,struct_banded> band(const Matrix<Integer,struct_dense>&);
template Matrix<Real,struct_banded> band(const Matrix<Real,struct_dense>&);
template Matrix<Float,struct_banded> band(const Matrix<Float,struct_dense>&);
template Matrix<Complex,struct_banded> band(const Matrix<Complex,struct_dense>&);
template Matrix<Float_complex,struct_banded> band(const Matrix<Float_complex,struct_dense>&);
template Matrix<Object,struct_banded> band(const Matrix<Object,struct_dense>&);

template Matrix<Integer,struct_banded> band(const Matrix<Integer,struct_sparse>&);
template Matrix<Real,struct_banded> band(const Matrix<Real,struct_sparse>&);
template Matrix<Float,struct_banded> band(const Matrix<Float,struct_sparse>&);
template Matrix<Complex,struct_banded> band(const Matrix<Complex,struct_sparse>&);
template Matrix<Float_complex,struct_banded> band(const Matrix<Float_complex,struct_sparse>&);
template Matrix<Object,struct_banded> band(const Matrix<Object,struct_sparse>&);

template struct details::manip_trans_converter_helper<struct_dense, Integer, Integer>;
template struct details::manip_trans_converter_helper<struct_dense, Integer, Real>;
template struct details::manip_trans_converter_helper<struct_dense, Integer, Complex>;
template struct details::manip_trans_converter_helper<struct_dense, Float, Float>;
template struct details::manip_trans_converter_helper<struct_dense, Float, Real>;
template struct details::manip_trans_converter_helper<struct_dense, Float, Complex>;
template struct details::manip_trans_converter_helper<struct_dense, Float, Float_complex>;
template struct details::manip_trans_converter_helper<struct_dense, Real, Real>;
template struct details::manip_trans_converter_helper<struct_dense, Real, Complex>;
template struct details::manip_trans_converter_helper<struct_dense, Complex, Complex>;
template struct details::manip_trans_converter_helper<struct_dense, Float_complex, Float_complex>;
template struct details::manip_trans_converter_helper<struct_dense, Float_complex, Complex>;
template struct details::manip_trans_converter_helper<struct_dense, Object, Object>;

template struct details::manip_trans_converter_helper<struct_sparse, Integer, Integer>;
template struct details::manip_trans_converter_helper<struct_sparse, Integer, Real>;
template struct details::manip_trans_converter_helper<struct_sparse, Integer, Complex>;
template struct details::manip_trans_converter_helper<struct_sparse, Float, Float>;
template struct details::manip_trans_converter_helper<struct_sparse, Float, Real>;
template struct details::manip_trans_converter_helper<struct_sparse, Float, Complex>;
template struct details::manip_trans_converter_helper<struct_sparse, Float, Float_complex>;
template struct details::manip_trans_converter_helper<struct_sparse, Real, Real>;
template struct details::manip_trans_converter_helper<struct_sparse, Real, Complex>;
template struct details::manip_trans_converter_helper<struct_sparse, Complex, Complex>;
template struct details::manip_trans_converter_helper<struct_sparse, Float_complex, Float_complex>;
template struct details::manip_trans_converter_helper<struct_sparse, Float_complex, Complex>;
template struct details::manip_trans_converter_helper<struct_sparse, Object, Object>;

template struct details::manip_trans_converter_helper<struct_banded, Integer, Integer>;
template struct details::manip_trans_converter_helper<struct_banded, Integer, Real>;
template struct details::manip_trans_converter_helper<struct_banded, Integer, Complex>;
template struct details::manip_trans_converter_helper<struct_banded, Float, Float>;
template struct details::manip_trans_converter_helper<struct_banded, Float, Real>;
template struct details::manip_trans_converter_helper<struct_banded, Float, Complex>;
template struct details::manip_trans_converter_helper<struct_banded, Float, Float_complex>;
template struct details::manip_trans_converter_helper<struct_banded, Real, Real>;
template struct details::manip_trans_converter_helper<struct_banded, Real, Complex>;
template struct details::manip_trans_converter_helper<struct_banded, Complex, Complex>;
template struct details::manip_trans_converter_helper<struct_banded, Float_complex, Float_complex>;
template struct details::manip_trans_converter_helper<struct_banded, Float_complex, Complex>;
template struct details::manip_trans_converter_helper<struct_banded, Object, Object>;

template struct details::manip_trans_reshaper_helper<Integer, Integer>;
template struct details::manip_trans_reshaper_helper<Integer, Real>;
template struct details::manip_trans_reshaper_helper<Integer, Complex>;
template struct details::manip_trans_reshaper_helper<Float, Float>;
template struct details::manip_trans_reshaper_helper<Float, Real>;
template struct details::manip_trans_reshaper_helper<Float, Complex>;
template struct details::manip_trans_reshaper_helper<Float, Float_complex>;
template struct details::manip_trans_reshaper_helper<Real, Real>;
template struct details::manip_trans_reshaper_helper<Real, Complex>;
template struct details::manip_trans_reshaper_helper<Complex, Complex>;
template struct details::manip_trans_reshaper_helper<Float_complex, Float_complex>;
template struct details::manip_trans_reshaper_helper<Float_complex, Complex>;
template struct details::manip_trans_reshaper_helper<Object, Object>;

};};

MACRO_INSTANTIATE_G_1(matcl::raw::details::manip_reshape_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::manip_reshape_helper)

MACRO_INSTANTIATE_G_1(matcl::raw::details::manip_trans_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::manip_trans_helper)

MACRO_INSTANTIATE_G_1(matcl::raw::details::manip_tr_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::details::manip_tr_helper)

MACRO_INSTANTIATE_G_1(matcl::raw::all_finite_helper)
MACRO_INSTANTIATE_S_1(matcl::raw::all_finite_helper)

MACRO_INSTANTIATE_SCAL_1(matcl::raw::details::manip_trans_helper_mat)

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::integer_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::float_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::real_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::float_complex_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::complex_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::object_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::integer_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::float_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::real_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::complex_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::float_complex_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::object_sparse&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ld(const matcl::raw::integer_dense&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ld(const matcl::raw::float_dense&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ld(const matcl::raw::real_dense&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ld(const matcl::raw::float_complex_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::complex_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::object_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::integer_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::float_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::real_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::complex_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::float_complex_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ud(const matcl::raw::object_dense&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::integer_band&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::float_band&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::real_band&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::complex_band&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::float_complex_band&, Integer, bool);

template MATCL_MATREP_EXPORT
matcl::Integer matcl::raw::get_ld(const matcl::raw::object_band&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ud(const matcl::raw::integer_band&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ud(const matcl::raw::float_band&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ud(const matcl::raw::real_band&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ud(const matcl::raw::complex_band&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ud(const matcl::raw::float_complex_band&, Integer, bool);

template MATCL_MATREP_EXPORT 
matcl::Integer matcl::raw::get_ud(const matcl::raw::object_band&, Integer, bool);


template bool matcl::raw::is_sym(const matcl::raw::integer_sparse&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::float_sparse&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::real_sparse&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::float_complex_sparse&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::complex_sparse&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::object_sparse&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::integer_sparse&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::float_sparse&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::real_sparse&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::float_complex_sparse&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::complex_sparse&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::object_sparse&, Real, bool);

template bool matcl::raw::is_sym(const matcl::raw::integer_dense&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::float_dense&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::real_dense&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::float_complex_dense&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::complex_dense&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::object_dense&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::integer_dense&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::float_dense&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::real_dense&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::float_complex_dense&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::complex_dense&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::object_dense&, Real, bool);

template bool matcl::raw::is_sym(const matcl::raw::integer_band&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::float_band&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::real_band&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::float_complex_band&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::complex_band&, Real, bool);
template bool matcl::raw::is_sym(const matcl::raw::object_band&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::integer_band&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::float_band&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::real_band&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::complex_band&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::float_complex_band&, Real, bool);
template bool matcl::raw::is_her(const matcl::raw::object_band&, Real, bool);

template matcl::Integer matcl::raw::nnz_total(const matcl::raw::integer_dense&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::float_dense&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::real_dense&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::complex_dense&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::float_complex_dense&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::object_dense&);

template matcl::Integer matcl::raw::nnz_total(const matcl::raw::integer_sparse&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::float_sparse&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::real_sparse&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::complex_sparse&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::float_complex_sparse&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::object_sparse&);

template matcl::Integer matcl::raw::nnz_total(const matcl::raw::integer_band&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::float_band&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::real_band&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::complex_band&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::float_complex_band&);
template matcl::Integer matcl::raw::nnz_total(const matcl::raw::object_band&);

#pragma warning( pop )
