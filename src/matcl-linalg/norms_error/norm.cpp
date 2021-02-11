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

#include "matcl-linalg/norms_error/norm.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/decompositions/svd.h"
#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-core/utils/workspace.h"

namespace matcl
{

namespace details
{

//static inline Real abs2(Real val)       { return val*val; };
//static inline Real abs2(Complex val)    { return abs2(real(val)) + abs2(imag(val)); };

template<class T>
Real norm_diag_impl(const dense_matrix<T>& rep, Integer N, Real p)
{
    using constants::inf;

    const T* ptr = rep.ptr();

    if (p == -1)
    {
        Real value = 0;
        for (Integer i = 0; i < N; ++i)
        {
            if (matcl::is_nan(ptr[i]) == true)
                return matcl::constants::nan();

            value += abs2(ptr[i]);
        };
        return sqrt(value);
    };

    Real value = 0;
    for (Integer i = 0; i < N; ++i)
    {
        if (matcl::is_nan(ptr[i]) == true)
            return matcl::constants::nan();

        value = max(value,abs(ptr[i]));
    };
    return value;
};

Real norm_diag(const Matrix& A, Real p)
{
    switch(A.get_value_code())
    {
        case value_code::v_integer:
            return norm_diag_impl(dense_matrix<Integer>(A),A.length(),p);
        case value_code::v_float:
            return norm_diag_impl(dense_matrix<Float>(A), A.length(),p);
        case value_code::v_real:
            return norm_diag_impl(dense_matrix<Real>(A), A.length(),p);
        case value_code::v_float_complex:
            return norm_diag_impl(dense_matrix<Float_complex>(A), A.length(),p);
        case value_code::v_complex:
            return norm_diag_impl(dense_matrix<Complex>(A), A.length(),p);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("norm");
        default:
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
    };
};

template<class V, class S> struct norm_str {};

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V> 
struct norm_str<V,struct_dense> 
{
    using Mat       = raw::Matrix<V,struct_dense>;
    using VR        = typename md::real_type<V>::type;

    static Real eval(const Mat& A, Real p)
    {
        if (A.size() == 0)
            return 0;

        using VR = typename matcl::details::real_type<V>::type;
        using constants::inf;

        VR * typed_nullptr = nullptr;

        if (p == 1)     return lapack::lange("1",A.rows(),A.cols(),lap(A.ptr()),A.ld(),typed_nullptr);
        if (p == -1)    return lapack::lange("F",A.rows(),A.cols(),lap(A.ptr()),A.ld(),typed_nullptr);
        if (p == -2)    return lapack::lange("M",A.rows(),A.cols(),lap(A.ptr()),A.ld(),typed_nullptr);

        if (p == inf())
        {            
            using workspace         = matcl::pod_workspace<VR>;
            workspace work          = workspace(A.rows());
            VR* ptr_work            = work.ptr();

            return lapack::lange("I",A.rows(),A.cols(),lap(A.ptr()),A.ld(),lap(ptr_work));
        };

        //otherwise second norm

        if (A.rows() == 1 || A.cols() == 1)
        {
            //for vectors second norm is equal to frobenius norm
            return lapack::lange("F",A.rows(),A.cols(),lap(A.ptr()),A.ld(),typed_nullptr);
        };
        
        Matrix S = svd1(Matrix(A,false),svd_algorithm::dc);
        S = S(1);
        return S.get_scalar<Real>();
    };

    static Real eval_vec_all(const Mat& A, basic_vector_norm p)
    {
        switch(p)
        {
            case basic_vector_norm::norm_2:     return eval(A, -1.0);
            case basic_vector_norm::norm_inf:   return eval(A, -2.0);
        };

        //norm_1
        return eval_vec_1_all(A);
    }

    static void eval_vec(Matrix& ret, const Mat& A, basic_vector_norm p, Integer d)
    {
        switch(p)
        {
            case basic_vector_norm::norm_1:
                return eval_vec_1(ret, A, d);
            case basic_vector_norm::norm_2:
                return eval_vec_2(ret, A, d);
            case basic_vector_norm::norm_inf:
                return eval_vec_inf(ret, A, d);
            default:
                //invalid case
                return eval_vec_1(ret, A, d);
        }
    };

    static Real eval_vec_1_all(const Mat& A)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();
            
        if (N == 0 || M == 0)
            return 0.0;

        const V* ptr_A  = A.ptr();
        Integer A_ld    = A.ld();

        Real sum        = 0.0;

        for (Integer i = 0; i < N; ++i)
        {        
            for (Integer j = 0; j < M; ++j)
                sum += abs(ptr_A[j]);

            ptr_A       += A_ld;
        };

        return sum;
    }

    static void eval_vec_1(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            const V* ptr_A  = A.ptr();
            Integer A_ld    = A.ld();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR sum  = VR(0.0);
                
                for (Integer j = 0; j < M; ++j)
                    sum += abs(ptr_A[j]);

                ptr_res[i]  = sum;
                ptr_A       += A_ld;
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();
            const V* ptr_A  = A.ptr();
            Integer A_ld    = A.ld();

            for (Integer j = 0; j < M; ++j)
                ptr_res[j]  = abs(ptr_A[j]);

            ptr_A           += A_ld;

            for (Integer i = 1; i < N; ++i)
            {                
                for (Integer j = 0; j < M; ++j)
                {
                    ptr_res[j]  += abs(ptr_A[j]);
                };

                ptr_A       += A_ld;
            };

            ret = Matrix(res, false);
            return;
        };
    }

    static void eval_vec_2(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            const V* ptr_A  = A.ptr();
            Integer A_ld    = A.ld();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR sum      = VR(0.0);
                
                for (Integer j = 0; j < M; ++j)
                    sum     += abs2(ptr_A[j]);

                ptr_res[i]  = sqrt(sum);
                ptr_A       += A_ld;
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();
            const V* ptr_A  = A.ptr();
            Integer A_ld    = A.ld();

            for (Integer j = 0; j < M; ++j)
                ptr_res[j]  = abs2(ptr_A[j]);

            ptr_A           += A_ld;

            for (Integer i = 1; i < N; ++i)
            {                
                for (Integer j = 0; j < M; ++j)
                {
                    ptr_res[j]  += abs2(ptr_A[j]);
                };

                ptr_A       += A_ld;
            };

            for (Integer j = 0; j < M; ++j)
                ptr_res[j]  = sqrt(ptr_res[j]);

            ret = Matrix(res, false);
            return;
        };
    }

    static void eval_vec_inf(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            const V* ptr_A  = A.ptr();
            Integer A_ld    = A.ld();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR max      = VR(0.0);
                
                for (Integer j = 0; j < M; ++j)
                    max     = std::max(max, abs(ptr_A[j]));

                ptr_res[i]  = max;
                ptr_A       += A_ld;
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();
            const V* ptr_A  = A.ptr();
            Integer A_ld    = A.ld();

            for (Integer j = 0; j < M; ++j)
                ptr_res[j]  = abs(ptr_A[j]);

            ptr_A           += A_ld;

            for (Integer i = 1; i < N; ++i)
            {                
                for (Integer j = 0; j < M; ++j)
                {
                    ptr_res[j]  = std::max(ptr_res[j], abs(ptr_A[j]));
                };

                ptr_A       += A_ld;
            };

            ret = Matrix(res, false);
            return;
        };
    }
};

//--------------------------------------------------------------------------------
//                                  SPARSE
//--------------------------------------------------------------------------------
template<class V> 
struct norm_str<V,struct_sparse> 
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    using VR    = typename md::real_type<V>::type;

    static Real eval(const Mat& A, Real p)
    {
        if (A.nnz() == 0)
            return 0;

        using constants::inf;

        if (p == 1)     return eval_1(A);
        if (p == -1)    return eval_F(A);
        if (p == -2)    return eval_M(A);
        if (p == inf()) return eval_I(A);

        //otherwise second norm

        if (A.rows() == 1 || A.cols() == 1)
        {
            //for vectors second norm is equal to frobenius norm
            return eval_F(A);
        };
        
        Matrix S = svd1(Matrix(A,false),svd_algorithm::dc);
        S = S(1);
        return S.get_scalar<Real>();
    };

    static Real eval_vec_all(const Mat& A, basic_vector_norm p)
    {
        switch(p)
        {
            case basic_vector_norm::norm_2:     return eval(A, -1.0);
            case basic_vector_norm::norm_inf:   return eval(A, -2.0);
        };

        //norm_1
        return eval_vec_1_all(A);
    }

    static void eval_vec(Matrix& ret, const Mat& A, basic_vector_norm p, Integer d)
    {
        switch(p)
        {
            case basic_vector_norm::norm_1:
                return eval_vec_1(ret, A, d);
            case basic_vector_norm::norm_2:
                return eval_vec_2(ret, A, d);
            case basic_vector_norm::norm_inf:
                return eval_vec_inf(ret, A, d);
            default:
                //invalid case
                return eval_vec_1(ret, A, d);
        }
    };

    static Real eval_vec_1_all(const Mat& A)
    {
        Integer N = A.cols();
        const raw::details::sparse_ccs<V>& d = A.rep();

        const Integer* d_c  = d.ptr_c();
        const V* d_x        = d.ptr_x();

        Real acc = 0.0;

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                acc += abs(d_x[k]);
        }

        return acc;
    };

    static Real eval_1(const Mat& A)
    {
        Integer N = A.cols();
        const raw::details::sparse_ccs<V>& d = A.rep();

        const Integer* d_c  = d.ptr_c();
        const V* d_x        = d.ptr_x();

        Real maxs = 0;
        for (Integer j = 0; j < N; ++j)
        {
            Real s = 0.;
            for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
            {
                s += abs(d_x[k]);
            };
            if (s > maxs)
            {
                maxs = s;
            };
        }
        return maxs;
    };
    static Real eval_M(const Mat& A)
    {
        Integer N = A.cols();
        const raw::details::sparse_ccs<V>& d = A.rep();

        const Integer* d_c  = d.ptr_c();
        const V* d_x        = d.ptr_x();

        Real maxs = 0;
        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
            {
                maxs = max(maxs,abs(d_x[k]));
            };
        }
        return maxs;
    };
    static Real eval_I(const Mat& A)
    {
        Integer M = A.rows();
        using VR = typename real_type<V>::type;

        using workspace         = matcl::pod_workspace<VR>;
        workspace work          = workspace(M, VR(0.0));
        VR* sv                  = work.ptr();

        const raw::details::sparse_ccs<V>& d = A.rep();

        const Integer* d_r  = d.ptr_r() + d.offset();
        const V* d_x        = d.ptr_x() + d.offset();

        Integer N = d.nnz();
        for (Integer k = 0; k < N; ++k)
        {
            Integer i   = d_r[k];
            sv[i]       = sv[i] + abs(d_x[k]);            
        };
        Real maxs = 0.;
        for(Integer i = 0; i < M; ++i)
        {
            if (sv[i] > maxs)
                maxs = sv[i];
        };

        return maxs;
    };
    static Real eval_F(const Mat& A)
    {
        Integer N = A.cols();
        const raw::details::sparse_ccs<V>& d = A.rep();

        const Integer* d_c  = d.ptr_c();
        const V* d_x        = d.ptr_x();

        Real maxs = 0;
        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = d_c[j]; k < d_c[j+1]; ++k)
                maxs += abs2(d_x[k]);
        }

        return sqrt(maxs);
    };

    static void eval_vec_1(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        const raw::details::sparse_ccs<V>& rep = A.rep();

        const Integer* d_c  = rep.ptr_c();
        const Integer* d_r  = rep.ptr_r();
        const V* d_x        = rep.ptr_x();

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR acc      = VR(0.0);
                
                for (Integer k = d_c[i]; k < d_c[i+1]; ++k)
                    acc     += abs(d_x[k]);

                ptr_res[i]  = acc;
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), VR(0.0), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();

            for (Integer i = 0; i < N; ++i)
            {                
                for (Integer k = d_c[i]; k < d_c[i+1]; ++k)
                {
                    Integer r   = d_r[k];
                    ptr_res[r]  += abs(d_x[k]);
                };
            };

            ret = Matrix(res, false);
            return;
        };
    }

    static void eval_vec_2(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        const raw::details::sparse_ccs<V>& rep = A.rep();

        const Integer* d_c  = rep.ptr_c();
        const Integer* d_r  = rep.ptr_r();
        const V* d_x        = rep.ptr_x();

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR acc      = VR(0.0);
                
                for (Integer k = d_c[i]; k < d_c[i+1]; ++k)
                    acc     += abs2(d_x[k]);

                ptr_res[i]  = sqrt(acc);
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), VR(0.0), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();

            for (Integer i = 0; i < N; ++i)
            {                
                for (Integer k = d_c[i]; k < d_c[i+1]; ++k)
                {
                    Integer r   = d_r[k];
                    ptr_res[r]  += abs2(d_x[k]);
                };
            };

            for (Integer i = 0; i < M; ++i)
                ptr_res[i]  = sqrt(ptr_res[i]);

            ret = Matrix(res, false);
            return;
        };
    }

    static void eval_vec_inf(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        const raw::details::sparse_ccs<V>& rep = A.rep();

        const Integer* d_c  = rep.ptr_c();
        const Integer* d_r  = rep.ptr_r();
        const V* d_x        = rep.ptr_x();

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR acc      = VR(0.0);
                
                for (Integer k = d_c[i]; k < d_c[i+1]; ++k)
                    acc     = std::max(acc, abs(d_x[k]));

                ptr_res[i]  = acc;
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), VR(0.0), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();

            for (Integer i = 0; i < N; ++i)
            {                
                for (Integer k = d_c[i]; k < d_c[i+1]; ++k)
                {
                    Integer r   = d_r[k];
                    ptr_res[r]  = std::max(ptr_res[r], abs(d_x[k]));
                };
            };

            ret = Matrix(res, false);
            return;
        };
    }
};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class V> 
struct norm_str<V,struct_banded> 
{
    using Mat   = raw::Matrix<V,struct_banded>;
    using VR    = typename md::real_type<V>::type;

    static Real eval(const Mat& A, Real p)
    {
        if (A.size() == 0)
            return 0;

        using constants::inf;

        if (p == 1)     return eval_1(A);
        if (p == -1)    return eval_F(A);
        if (p == -2)    return eval_M(A);
        if (p == inf()) return eval_I(A);

        //otherwise second norm        
        if (A.rows() == 1 || A.cols() == 1)
        {
            //for vectors second norm is equal to frobenius norm
            return eval_F(A);
        };

        Matrix S = svd1(Matrix(A,false),svd_algorithm::dc);
        S = S(1);
        return S.get_scalar<Real>();
    };

    static Real eval_vec_all(const Mat& A, basic_vector_norm p)
    {
        switch(p)
        {
            case basic_vector_norm::norm_2:     return eval(A, -1.0);
            case basic_vector_norm::norm_inf:   return eval(A, -2.0);
        };

        return eval_vec_1_all(A);
    }

    static void eval_vec(Matrix& ret, const Mat& A, basic_vector_norm p, Integer d)
    {
        switch(p)
        {
            case basic_vector_norm::norm_1:
                return eval_vec_1(ret, A, d);
            case basic_vector_norm::norm_2:
                return eval_vec_2(ret, A, d);
            case basic_vector_norm::norm_inf:
                return eval_vec_inf(ret, A, d);
            default:
                //invalid case
                return eval_vec_1(ret, A, d);
        }
    };
    static Real eval_M(const Mat& A)
    {
        Integer N = A.cols();

        const V* ptr = A.rep_ptr();
        Real value  = 0;
        for (Integer J = 0; J < N; ++J, ptr += A.ld())
        {
            Integer fr = A.first_row(J);
            Integer lr = A.last_row(J);
            Integer ii = A.first_elem_pos(J);

            for (Integer i = fr; i <= lr; ++i, ++ii)
                value = max(value, abs(ptr[ii]));
        };

        return value;
    };
    static Real eval_1(const Mat& A)
    {
        Integer N = A.cols();

        const V* ptr = A.rep_ptr();
        Real value = 0;

        for (Integer J = 0; J < N; ++J, ptr += A.ld())
        {
            Integer fr = A.first_row(J);
            Integer lr = A.last_row(J);
            Integer ii = A.first_elem_pos(J);

            Real sum = 0.;
            for (Integer i = fr; i <= lr; ++i, ++ii)
                sum = sum + abs(ptr[ii]);

            value = max(value,sum);
        };
        
        return value;
    };
    static Real eval_I(const Mat& A)
    {
        Integer M = A.rows(), N = A.cols();

        using VR = typename matcl::details::real_type<V>::type;

        using workspace         = matcl::pod_workspace<VR>;
        workspace work          = workspace(M, VR(0.0));
        VR* sv                  = work.ptr();

        const V* ptr = A.rep_ptr();
        for (Integer J = 0; J < N; ++J, ptr += A.ld())
        {
            Integer fr = A.first_row(J);
            Integer lr = A.last_row(J);
            Integer ii = A.first_elem_pos(J);

            for (Integer i = fr; i <= lr; ++i, ++ii)
                sv[i] = sv[i] + abs(ptr[ii]);
        };

        Real maxs = 0.;
        for(Integer i = 0; i < M; ++i)
        {
            if (sv[i] > maxs)
                maxs = sv[i];
        };

        return maxs;
    };

    static Real eval_F(const Mat& A)
    {
        Integer N = A.cols();

        const V* ptr = A.rep_ptr();
        Real value = 0;
        for (Integer J = 0; J < N; ++J, ptr += A.ld())
        {
            Integer fr = A.first_row(J);
            Integer lr = A.last_row(J);
            Integer ii = A.first_elem_pos(J);
    
            for (Integer i = fr; i <= lr; ++i, ++ii)
                value += abs2(ptr[ii]);
        };

        return sqrt(value);
    };

    static Real eval_vec_1_all(const Mat& A)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        const V* ptr_A  = A.rep_ptr();
        Integer A_ld    = A.ld();

        if (N == 0 || M == 0)
            return 0.0;

        Real acc        = 0.0;
            
        for (Integer i = 0; i < N; ++i)
        {
            Integer fr = A.first_row(i);
            Integer lr = A.last_row(i);
            Integer ii = A.first_elem_pos(i);

            for (Integer k = fr; k <= lr; ++k, ++ii)
                acc     += abs(ptr_A[ii]);

            ptr_A       += A_ld;
        };

        return acc;
    }

    static void eval_vec_1(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        const V* ptr_A  = A.rep_ptr();
        Integer A_ld    = A.ld();

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR acc      = VR(0.0);

                Integer fr = A.first_row(i);
                Integer lr = A.last_row(i);
                Integer ii = A.first_elem_pos(i);

                for (Integer k = fr; k <= lr; ++k, ++ii)
                    acc     += abs(ptr_A[ii]);

                ptr_A       += A_ld;
                ptr_res[i]  = acc;
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), VR(0.0), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();

            for (Integer i = 0; i < N; ++i)
            {                
                Integer fr = A.first_row(i);
                Integer lr = A.last_row(i);
                Integer ii = A.first_elem_pos(i);

                for (Integer k = fr; k <= lr; ++k, ++ii)
                    ptr_res[k]  += abs(ptr_A[ii]);

                ptr_A       += A_ld;
            };

            ret = Matrix(res, false);
            return;
        };
    }

    static void eval_vec_2(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        const V* ptr_A  = A.rep_ptr();
        Integer A_ld    = A.ld();

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR acc      = VR(0.0);

                Integer fr = A.first_row(i);
                Integer lr = A.last_row(i);
                Integer ii = A.first_elem_pos(i);

                for (Integer k = fr; k <= lr; ++k, ++ii)
                    acc     += abs2(ptr_A[ii]);

                ptr_A       += A_ld;
                ptr_res[i]  = sqrt(acc);
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), VR(0.0), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();

            for (Integer i = 0; i < N; ++i)
            {
                Integer fr = A.first_row(i);
                Integer lr = A.last_row(i);
                Integer ii = A.first_elem_pos(i);

                for (Integer k = fr; k <= lr; ++k, ++ii)
                    ptr_res[k]  += abs2(ptr_A[ii]);

                ptr_A       += A_ld;
            };

            for (Integer i = 0; i < M; ++i)
                ptr_res[i]  = sqrt(ptr_res[i]);

            ret = Matrix(res, false);
            return;
        };
    }

    static void eval_vec_inf(Matrix& ret, const Mat& A, Integer d)
    {
        Integer M   = A.rows();
        Integer N   = A.cols();

        using Mat_R = raw::Matrix<VR,struct_dense>;

        const V* ptr_A  = A.rep_ptr();
        Integer A_ld    = A.ld();

        if (d == 1)
        {
            Mat_R res(A.get_type(), 1, N);
            
            if (N == 0)
                return;

            VR* ptr_res     = res.ptr();
            
            for (Integer i = 0; i < N; ++i)
            {
                VR acc      = VR(0.0);

                Integer fr = A.first_row(i);
                Integer lr = A.last_row(i);
                Integer ii = A.first_elem_pos(i);

                for (Integer k = fr; k <= lr; ++k, ++ii)
                    acc     = std::max(acc, abs(ptr_A[ii]));

                ptr_A       += A_ld;
                ptr_res[i]  = acc;
            };

            ret = Matrix(res, false);
            return;
        }
        else
        {
            Mat_R res(A.get_type(), VR(0.0), M, 1);

            if (M == 0)
                return;

            VR* ptr_res     = res.ptr();

            for (Integer i = 0; i < N; ++i)
            {                
                Integer fr = A.first_row(i);
                Integer lr = A.last_row(i);
                Integer ii = A.first_elem_pos(i);

                for (Integer k = fr; k <= lr; ++k, ++ii)
                    ptr_res[k]  = std::max(ptr_res[k],abs(ptr_A[ii]));

                ptr_A       += A_ld;
            };

            ret = Matrix(res, false);
            return;
        };
    }
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct norm_impl
{
    using M = raw::Matrix<V,S>;

    static Real eval(const M& A, Real p)
    {
        return norm_str<V,S>::eval(A,p);
    };
};

template<class S>
struct norm_impl<Integer,S>
{
    using M = raw::Matrix<Integer,S>;
    static Real eval(const M& A, Real p)
    {
        using MC    = raw::Matrix<Real,S>;
        MC AC       = raw::converter<MC,M>::eval(A);

        return norm_impl<Real,S>::eval(AC,p);
    };
};

template<class V, class S>
struct norm_vec_impl
{
    using M = raw::Matrix<V,S>;

    static void eval(Matrix& ret, const M& A, basic_vector_norm p, Integer d)
    {
        return norm_str<V,S>::eval_vec(ret, A,p,d);
    };
    static Real eval_all(const M& A, basic_vector_norm p)
    {
        return norm_str<V,S>::eval_vec_all(A,p);
    };
};

template<class S>
struct norm_vec_impl<Integer,S>
{
    using M = raw::Matrix<Integer,S>;
    static void eval(Matrix& ret, const M& A, basic_vector_norm p, Integer d)
    {
        using MC    = raw::Matrix<Real,S>;
        MC AC       = raw::converter<MC,M>::eval(A);

        return norm_vec_impl<Real,S>::eval(ret, AC,p,d);
    };
    static Real eval_all(const M& A, basic_vector_norm p)
    {
        using MC    = raw::Matrix<Real,S>;
        MC AC       = raw::converter<MC,M>::eval(A);

        return norm_vec_impl<Real,S>::eval_all(AC,p);
    };
};

struct unary_visitor_norm : public extract_type_switch<Real,unary_visitor_norm,true>
{
    template<class T>
    static Real eval(const Matrix&, const T& mat, Real p)
    {        
        using V = typename T::value_type;
        using S = typename T::struct_type;

        test_nan_remember_inf test_object;
        bool has_no_nan = test_range<V,S,test_nan_remember_inf>::eval(mat, test_object);
        bool has_no_inf = test_object.get_has_no_inf();

        if (has_no_nan == false)
            return matcl::constants::nan();

        if (has_no_inf == false)
            return matcl::constants::inf();
        
        return details::norm_impl<V,S>::eval(mat,p);
    };

    template<class T>
    static Real eval_scalar(const Matrix&, const T& mat, Real)
    {
        return abs(mat);
    };

    template<class S>
    static Real eval(const Matrix&, const raw::Matrix<Object,S>&, Real)
    {        
        throw error::object_value_type_not_allowed("norm_vec");
    };

    static Real eval_scalar(const Matrix&, const Object&, Real )
    {
        throw error::object_value_type_not_allowed("norm");
    };
};

struct unary_visitor_norm_vec : public extract_type_switch<void,unary_visitor_norm_vec,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, basic_vector_norm p, Integer d)
    {        
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::norm_vec_impl<V,S>::eval(ret, mat, p, d);
    };

    template<class T>
    static void eval_scalar(const Matrix&, const T& mat, Matrix& ret, basic_vector_norm, Integer)
    {
        ret = abs(mat);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, basic_vector_norm, Integer)
    {        
        throw error::object_value_type_not_allowed("norm_vec");
    };

    static void eval_scalar(const Matrix&, const Object&, Matrix&, basic_vector_norm, Integer)
    {
        throw error::object_value_type_not_allowed("norm_vec");
    };
};
struct unary_visitor_norm_vec_all : public extract_type_switch<Real,unary_visitor_norm_vec_all,true>
{
    template<class T>
    static Real eval(const Matrix&, const T& mat, basic_vector_norm p)
    {        
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::norm_vec_impl<V,S>::eval_all(mat, p);
    };

    template<class T>
    static Real eval_scalar(const Matrix&, const T& mat, basic_vector_norm)
    {
        return abs(mat);
    };

    template<class S>
    static Real eval(const Matrix&, const raw::Matrix<Object,S>&, basic_vector_norm)
    {        
        throw error::object_value_type_not_allowed("norm_vec_all");
    };

    static Real eval_scalar(const Matrix&, const Object&, basic_vector_norm)
    {
        throw error::object_value_type_not_allowed("norm_vec_all");
    };
};

template<class Val, bool Is_real = md::is_float_real_scalar<Val>::value>
struct normest_it
{
    using VR    = typename md::real_type<Val>::type;
    using Mat   = raw::Matrix<Val, struct_dense>;

    static Real eval(const linear_operator& f, Real p)
    {
        Integer M   = f.rows();
        Integer N   = f.cols();
        Integer MN  = std::max(M,N);

        if (M == 0 || N == 0)
            return 0.0;
        
        bool norm_1             = (p == 1.0);

        if (norm_1 == false)
            std::swap(M, N);

        using VTR_pod           = matcl::pod_type<Val>;
        using workspace         = matcl::pod_workspace<VTR_pod>;
        using iworkspace        = matcl::pod_workspace<Integer>;

        workspace work_V        = workspace(MN);
        iworkspace work_ISGN    = iworkspace(Is_real? MN : 0);

        Mat X(ti::ti_type<Val>(), MN, 1);

        Val* ptr_V          = reinterpret_cast<Val*>(work_V.ptr());
        Val* ptr_X          = X.ptr();
        Integer* ptr_ISGN   = work_ISGN.ptr();

        VR EST;
        Integer KASE        = 0;
        Integer ISAVE[3];

        lapack::lacn2(MN, lap(ptr_V), lap(ptr_X), ptr_ISGN, &EST, &KASE, ISAVE );

        trans_type t1       = norm_1? trans_type::no_trans  : trans_type::conj_trans;
        trans_type t2       = norm_1? trans_type::conj_trans : trans_type::no_trans;

        while ( KASE != 0 )
        {
            if ( KASE == 1 )
            {
                Matrix mat_x    = Matrix(X.make_view(1, N),false);
                Matrix tmp      = f.mmul_right(std::move(mat_x), t1);

                assign_res(ptr_X, M, N, tmp, false);
            }
            else
            {
                Matrix mat_x    = Matrix(X.make_view(1, M),false);
                Matrix tmp      = f.mmul_right(std::move(mat_x), t2);

                assign_res(ptr_X, M, N, tmp, true);
            };

            lapack::lacn2(MN, lap(ptr_V), lap(ptr_X), ptr_ISGN, &EST, &KASE, ISAVE );
        };

        return EST;
    };

    static void assign_res(Val* ptr_X, Integer M, Integer N, const Matrix& tmp, bool trans)
    {
        Integer exp_size    = (trans == false)? M : N;
        Integer MN          = std::max(M,N);
        
        const Val* ptr_res  = tmp.get_array<Val>();
        Integer res_M       = tmp.rows();
        Integer res_N       = tmp.cols();

        check_res_size(res_M, res_N, exp_size);

        if (ptr_X != ptr_res)
        {
            //evaluation was not done inplace; copy is needed
            for (Integer i = 0; i < exp_size; ++i)
                ptr_X[i]    = ptr_res[i];
        };

        for (Integer i = exp_size; i < MN; ++i)
            ptr_X[i]        = Val(0.0);
    };

    static void check_res_size(Integer res_M, Integer res_N, Integer exp_size)
    {
        if (res_M != exp_size || res_N != 1)
            throw error::invalid_size2(res_M, res_N, exp_size, 1);
    };
};

static Real normest_impl(const linear_operator& f, Real p)
{
    value_code vc = f.get_value_code();

    switch (vc)
    {
        case value_code::v_integer:
        {
            linear_operator fc = f.convert(value_code::v_real);
            return details::normest_it<Real>::eval(fc, p);
        }
        case value_code::v_real:
            return details::normest_it<Real>::eval(f, p);
        case value_code::v_float:
            return details::normest_it<Float>::eval(f, p);
        case value_code::v_complex:
            return details::normest_it<Complex>::eval(f, p);
        case value_code::v_float_complex:
            return details::normest_it<Float_complex>::eval(f, p);
        case value_code::v_object:
            throw error::object_value_type_not_allowed("normest");
        default:
        {
            matcl_assert(0,"invalid value code");
            throw error::error_general("invalid value code");
        }
    };
}

};

Real matcl::norm(const Matrix& A, Real p)
{
    if (p == 1 || p == 2 || p == -1 || p == -2 || p == constants::inf())
    {
    }
    else
    {
        throw error::error_norm_type();
    };

    if (A.structural_nnz() == 0)
        return 0.;

    if (A.get_struct().is_diag())
        return details::norm_diag(full(get_diag(A)),p);

    if (is_unitary(A.get_struct()))
    {
        if (p==2)
            return 1.;
    }

    return details::unary_visitor_norm::make<const Matrix&>(A,p);
};

Matrix matcl::norm_vec(const Matrix& A, basic_vector_norm p, Integer d)
{
    error::check_dim(d);

    if (A.structural_nnz() == 0)
    {
        if (d == 1)
            return zeros(1, A.cols());
        else
            return zeros(A.rows(), 1);
    };    

    if (A.get_struct().is_diag())
    {
        Matrix E = abs(get_diag(A));
        
        if (d == 1)
        {
            Integer N   = A.cols();
            E.resize(N, 1);
            return reshape(E, 1, N);
        }
        else
        {
            Integer M   = A.rows();
            E.resize(M, 1);
            return reshape(E, M, 1);
        };
    };

    Matrix ret;
    details::unary_visitor_norm_vec::make<const Matrix&>(A,ret,p,d);
    return ret;
};

Real matcl::norm_vec_all(const Matrix& A, basic_vector_norm p)
{
    if (A.structural_nnz() == 0)
        return 0.0;

    return details::unary_visitor_norm_vec_all::make<const Matrix&>(A,p);
};

Real matcl::norm(const Matrix& A, basic_vector_norm p)
{
    switch (p)
    {
        case basic_vector_norm::norm_1:     return norm(A,1.0);
        case basic_vector_norm::norm_2:     return norm(A,2.0);
        case basic_vector_norm::norm_inf:   return norm(A,constants::inf());

        default:
        {
            matcl_assert(0,"invalid case");
            return norm(A, 1.0);
        }
    };
}

Real matcl::normest_1(const linear_operator& f)
{
    return details::normest_impl(f, 1);
}

Real matcl::normest_inf(const linear_operator& f)
{
    return details::normest_impl(f, constants::inf());
}

Real matcl::normest_2(const linear_operator& f)
{
    //svds cannot solve very small problems
    static const Integer min_size = 5;

    if (f.rows() < min_size || f.cols() < min_size)
    {
        Matrix I    = speye(f.cols(), f.cols(), f.get_value_code());
        Matrix A    = f.mmul_right(I, trans_type::no_trans);
        return norm(A, 2.0);
    };

    options opts;

    //usuall only few iterations are required
    opts.set(opt::speigs::maxit(10000));
    opts.set(opt::speigs::tol(1e-1));

    Matrix E;
    bool conv;
    tie(E, conv) = svds1(f, 1, cluster_type::LM, opts);

    if (conv == true)
        return E.get_scalar<Real>();

    throw error::iterative_solver_not_converged();
};

Real matcl::abs_normest_r_inf(const linear_operator& f, const Matrix& X)
{
    trans_type t    = trans_type::no_trans;
    Integer N       = X.cols();

    error::check_mul(f.rows(), f.cols(), X.rows(), X.cols(), t, t);    

    if (N > 1)
        throw error::vector_required(X.rows(), X.cols()); 

    if (f.rows() == 0 || f.cols() == 0 || N == 0)
        return 0.0;

    linear_operator D   = bdiag(abs(X));
    linear_operator AD  = linop_mmul(f, D);
    Real n              = normest_inf(AD);

    return n;
};
Real matcl::abs_normest_l_inf(const linear_operator& f, const Matrix& X)
{
    trans_type t    = trans_type::no_trans;
    Integer N       = X.rows();

    error::check_mul(X.rows(), X.cols(), f.rows(), f.cols(), t, t);    

    if (N > 1)
        throw error::vector_required(X.rows(), X.cols()); 

    if (f.rows() == 0 || f.cols() == 0 || N == 0)
        return 0.0;

    linear_operator D   = bdiag(abs(X));
    linear_operator AD  = linop_mmul(D, f);
    Real n              = normest_1(AD);

    return n;
};

Matrix matcl::abs_normest_r_vec_inf(const linear_operator& f, const Matrix& X)
{
    trans_type t    = trans_type::no_trans;
    Integer N       = X.cols();

    error::check_mul(f.rows(), f.cols(), X.rows(), X.cols(), t, t);    

    if (f.rows() == 0 || f.cols() == 0 || N == 0)
        return zeros(1,N);

    using Mat_D = raw::Matrix<Real, struct_dense>;

    Mat_D norm(ti::ti_empty(), 1, N);
    Real* ptr_n     = norm.ptr();

    for (Integer i = 0; i < N; ++i)
    {       
        linear_operator D   = bdiag(abs(X(colon(), i+1)));
        linear_operator AD  = linop_mmul(f, D);
        Real n              = normest_inf(AD);
        ptr_n[i]            = n;
    };

    return Matrix(norm,false);
};
Matrix matcl::abs_normest_l_vec_inf(const linear_operator& f, const Matrix& X)
{
    trans_type t    = trans_type::no_trans;
    Integer N       = X.rows();

    error::check_mul(X.rows(), X.cols(), f.rows(), f.cols(), t, t);    

    if (f.rows() == 0 || f.cols() == 0 || N == 0)
        return zeros(N,1);

    using Mat_D = raw::Matrix<Real, struct_dense>;

    Mat_D norm(ti::ti_empty(), N, 1);
    Real* ptr_n     = norm.ptr();

    for (Integer i = 0; i < N; ++i)
    {       
        linear_operator D   = bdiag(abs(X(i+1, colon())));
        linear_operator AD  = linop_mmul(D, f);
        Real n              = normest_1(AD);
        ptr_n[i]            = n;
    };

    return Matrix(norm,false);
};

};