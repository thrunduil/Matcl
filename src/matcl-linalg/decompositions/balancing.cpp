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

#include "matcl-linalg/decompositions/balancing.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-core/utils/workspace.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"

namespace matcl { namespace details
{

//---------------------------------------------------------------
//                  HELPERS
//---------------------------------------------------------------

template<class V, bool Is_compl = details::is_complex<V>::value>
struct ldexp_helper
{};

template<>
struct ldexp_helper<Float,false>
{
    static Float eval(const Float& val, Integer e)
    {
        return ldexp(val,e);
    };
};
template<>
struct ldexp_helper<Real,false>
{
    static Real eval(const Real& val, Integer e)
    {
        return ldexp(val,e);
    };
};

template<class V>
struct ldexp_helper<V,true>
{
    static V eval(const V& val, Integer e)
    {
        using VR = typename details::real_type<V>::type;

        VR re   = ldexp(real(val),e);
        VR im   = ldexp(imag(val),e);
        return V(re,im);
    };
};

static void round_to_pow2(Matrix& D, bool pow2)
{
    if (pow2 == false)
        return;

    D   = log2(D);
    D   = round(D);
    D   = matcl::pow2(D);
    return;
};
//---------------------------------------------------------------
//                  matrix_balance_func
//---------------------------------------------------------------
template<class Mat, class ST = typename Mat::struct_type>
struct matrix_balance_func
{};

template<class Mat>
struct matrix_balance_func<Mat,struct_dense>
{
    using VT    = typename Mat::value_type;
    using VR    = typename md::real_type<VT>::type;

    static VR abs_diag(const Mat& A, Integer j)
    {
        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        return abs(ptr_A[j + j * A_ld]);
    };
    static VR maxabs_superdiags(const Mat& A, Integer j, const VR* ptr_D)
    {
        Integer A_ld    = A.ld();
        const VT* ptr_A = A.ptr() + j * A_ld;        

        VR d = VR(0.0);
        for (Integer i = 0; i < j; ++i)
        {
            VR val  = abs(ptr_A[i]) * ptr_D[i];
            d       = std::max(d, val);
        };

        return d;
    };

    static void maxabs_rows(const Mat& A, VR* ptr)
    {
        Integer N       = A.cols();
        Integer M       = A.rows();

        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < M; ++i)
                ptr[i]  = std::max(ptr[i], abs(ptr_A[i]));

            ptr_A       += A_ld;;
        };

        return;
    };
    static void maxabs_cols_scaled(const Mat& A, VR* ptr, const VR* ptr_scale)
    {
        Integer N       = A.cols();
        Integer M       = A.rows();

        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < M; ++i)
                ptr[j]  = std::max(ptr[j], abs(ptr_A[i]) * ptr_scale[i]);

            ptr_A       += A_ld;;
        };

        return;
    };

    static void maxabs_rows_scaled2(const Mat& A, VR* ptr, const VR* ptr_scal_rows, 
                                    const VR* ptr_scal_cols)
    {
        Integer N       = A.cols();
        Integer M       = A.rows();

        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < M; ++i)
            {
                VR tmp  = abs(ptr_A[i]) * ptr_scal_rows[i] * ptr_scal_cols[j];
                ptr[i]  = std::max(ptr[i], tmp);
            }

            ptr_A       += A_ld;;
        };

        return;
    };

    static void maxabs_cols_scaled2(const Mat& A, VR* ptr, const VR* ptr_scal_rows, 
                                    const VR* ptr_scal_cols)
    {
        Integer N       = A.cols();
        Integer M       = A.rows();

        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < M; ++i)
            {
                VR tmp  = abs(ptr_A[i]) * ptr_scal_rows[i] * ptr_scal_cols[j];
                ptr[j]  = std::max(ptr[j], tmp);
            }

            ptr_A       += A_ld;;
        };

        return;
    };

    static void sum_rows(const Mat& A, VR* ptr)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
                ptr[i]  += abs2(ptr_A[i]);

            ptr_A       += A_ld;;
        };
        return;
    };

    static void sum_cols(const Mat& A, VR* ptr)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
                ptr[j]  += abs2(ptr_A[i]);

            ptr_A       += A_ld;;
        };
        return;
    };

    static void sum_rows_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
            {
                if (i != j)
                    ptr[i]  += abs(ptr_A[i]);
                else
                    ptr[i]  += abs(ptr_A[i]) * scale_diag;
            }

            ptr_A       += A_ld;;
        };
        return;
    };
    static void max_rows_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
            {
                if (i != j)
                    ptr[i]  = std::max(ptr[i], abs(ptr_A[i]));
                else
                    ptr[i]  = std::max(ptr[i], abs(ptr_A[i]) * scale_diag);
            }

            ptr_A       += A_ld;;
        };
        return;
    };

    static void sum_cols_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
            {
                if (i != j)
                    ptr[j]  += abs(ptr_A[i]);
                else
                    ptr[j]  += abs(ptr_A[i]) * scale_diag;
            }

            ptr_A       += A_ld;;
        };
        return;
    };
    static void max_cols_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
            {
                if (i != j)
                    ptr[j]  = std::max(ptr[j], abs(ptr_A[i]));
                else
                    ptr[j]  = std::max(ptr[j], abs(ptr_A[i]) * scale_diag);
            }

            ptr_A       += A_ld;;
        };
        return;
    };

    static void scale_rows(Mat& A, const Integer* ptr, Integer mult)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
                ptr_A[i]    = ldexp_helper<VT>::eval(ptr_A[i],ptr[i]*mult);

            ptr_A       += A_ld;;
        };
        return;
    };
    static void scale_cols(Mat& A, const Integer* ptr, Integer mult)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer scale   = ptr[j]*mult;

            for (Integer i = 0; i < N; ++i)
                ptr_A[i]    = ldexp_helper<VT>::eval(ptr_A[i],scale);

            ptr_A       += A_ld;;
        };
        return;
    };

    static void scale_rows(Mat& A, const VR* ptr)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer i = 0; i < N; ++i)
                ptr_A[i]    = ptr_A[i] * ptr[i];

            ptr_A       += A_ld;;
        };
        return;
    };

    static void scale_cols_inv(Mat& A, const VR* ptr)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            VR scale    = VR(1.0)/ptr[j];

            for (Integer i = 0; i < N; ++i)
                ptr_A[i]    = ptr_A[i] * scale;

            ptr_A       += A_ld;;
        };
        return;
    };

};

template<class Mat>
struct matrix_balance_func<Mat,struct_sparse>
{
    using VT    = typename Mat::value_type; 
    using VR    = typename md::real_type<VT>::type;

    static VR abs_diag(const Mat& A, Integer j)
    {
        const auto& sprep   = A.rep();

        Integer k;
        bool has            = sprep.has_element(j,j,k);

        if (has == false)
            return VR(0.0);

        const VT* ptr_x     = sprep.ptr_x();
        return abs(ptr_x[k]);
    };

    static VR maxabs_superdiags(const Mat& A, Integer j, const VR* ptr_D)
    {
        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();

        VR d            = VR(0.0);

        for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
        {
            Integer r       = ptr_r[k];
            const VT& v     = ptr_x[k];

            if (r >= j)
                return d;

            VR val  = abs(v) * ptr_D[r];
            d       = std::max(d, val);
        };

        return d;
    };
    static void maxabs_rows(const Mat& A, VR* ptr)
    {
        Integer N       = A.cols();

        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                ptr[r]          = std::max(ptr[r], abs(ptr_x[k]));
            };
        };
    };

    static void maxabs_cols_scaled(const Mat& A, VR* ptr, const VR* ptr_scale)
    {
        Integer N       = A.cols();

        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                ptr[j]          = std::max(ptr[j], abs(ptr_x[k]) * ptr_scale[r]);
            };
        };
    };

    static void maxabs_rows_scaled2(const Mat& A, VR* ptr, const VR* ptr_scal_rows, 
                                    const VR* ptr_scal_cols)
    {
        Integer N       = A.cols();

        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                VR tmp          = abs(ptr_x[k]) * ptr_scal_rows[r] * ptr_scal_cols[j];
                ptr[r]          = std::max(ptr[r], tmp);
            };
        };
    };

    static void maxabs_cols_scaled2(const Mat& A, VR* ptr, const VR* ptr_scal_rows, 
                                    const VR* ptr_scal_cols)
    {
        Integer N       = A.cols();

        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                VR tmp          = abs(ptr_x[k]) * ptr_scal_rows[r] * ptr_scal_cols[j];
                ptr[j]          = std::max(ptr[j], tmp);
            };
        };
    };

    static void sum_rows(const Mat& A, VR* ptr)
    {
        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                ptr[r]          += abs2(ptr_x[k]);
            };
        };
    };

    static void sum_cols(const Mat& A, VR* ptr)
    {
        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const VT* ptr_x         = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
                ptr[j]          += abs2(ptr_x[k]);
        };
    };

    static void sum_rows_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];

                if (r != j)
                    ptr[r]      += abs(ptr_x[k]);
                else
                    ptr[r]      += abs(ptr_x[k]) * scale_diag;
            };
        };
    };
    static void max_rows_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];

                if (r != j)
                    ptr[r]      = std::max(ptr[r], abs(ptr_x[k]));
                else
                    ptr[r]      = std::max(ptr[r], abs(ptr_x[k]) * scale_diag);
            };
        };
    };

    static void sum_cols_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                if (r != j)
                    ptr[j]      += abs(ptr_x[k]);
                else
                    ptr[j]      += abs(ptr_x[k]) * scale_diag;
            }
        };
    };
    static void max_cols_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        const auto& sprep       = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        const VT* ptr_x         = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                if (r != j)
                    ptr[j]      = std::max(ptr[j], abs(ptr_x[k]));
                else
                    ptr[j]      = std::max(ptr[j], abs(ptr_x[k]) * scale_diag);
            }
        };
    };

    static void scale_rows(Mat& A, const Integer* ptr, Integer mult)
    {
        auto sprep              = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        VT* ptr_x               = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                ptr_x[k]        = ldexp_helper<VT>::eval(ptr_x[k],ptr[r]*mult);
            };
        };
    };
    static void scale_cols(Mat& A, const Integer* ptr, Integer mult)
    {
        auto sprep              = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        VT* ptr_x               = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            Integer scale       = ptr[j] * mult;

            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
                ptr_x[k]        = ldexp_helper<VT>::eval(ptr_x[k], scale);
        };
    };

    static void scale_rows(Mat& A, const VR* ptr)
    {
        auto sprep              = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        const Integer* ptr_r    = sprep.ptr_r();
        VT* ptr_x               = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
            {
                Integer r       = ptr_r[k];
                ptr_x[k]        = ptr_x[k] * ptr[r];
            };
        };
    };
    static void scale_cols_inv(Mat& A, const VR* ptr)
    {
        auto sprep              = A.rep();
        const Integer* ptr_c    = sprep.ptr_c();
        VT* ptr_x               = sprep.ptr_x();
        Integer N               = A.rows();

        for (Integer j = 0; j < N; ++j)
        {
            VR scale            = VR(1.0)/ptr[j];

            for (Integer k = ptr_c[j]; k < ptr_c[j+1]; ++k)
                ptr_x[k]        = ptr_x[k] * scale;
        };
    };

};

template<class Mat>
struct matrix_balance_func<Mat,struct_banded>
{
    using VT    = typename Mat::value_type; 
    using VR    = typename md::real_type<VT>::type;

    static VR abs_diag(const Mat& A, Integer j)
    {
        if (A.has_diag(0) == false)
            return VR(0.0);

        Integer A_ld    = A.ld();
        const VT* ptr_A = A.rep_ptr() + A.first_elem_diag(0);

        return abs(ptr_A[j * A_ld]);
    };

    static VR maxabs_superdiags(const Mat& A, Integer j, const VR* ptr_D)
    {
        Integer A_ld    = A.ld();
        const VT* ptr_A = A.rep_ptr() + j * A_ld;        

        Integer fr      = A.first_row(j);
        Integer lr      = A.last_row(j);
        Integer pos     = A.first_elem_pos(j);

        lr              = std::min(lr, j - 1);
        VR d            = VR(0.0);

        for (Integer i = fr; i <= lr; ++i, ++pos)
        {
            VR val  = abs(ptr_A[pos]) * ptr_D[i];
            d       = std::max(d, val);
        };

        return d;
    };

    static void maxabs_rows(const Mat& A, VR* ptr)
    {
        Integer N       = A.cols();

        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr[i]  = std::max(ptr[i], abs(ptr_A[pos]));

            ptr_A       += A_ld;
        };
    };

    static void maxabs_cols_scaled(const Mat& A, VR* ptr, const VR* ptr_scale)
    {
        Integer N       = A.cols();

        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr[j]  = std::max(ptr[j], abs(ptr_A[pos]) * ptr_scale[i]);

            ptr_A       += A_ld;
        };
    };

    static void maxabs_rows_scaled2(const Mat& A, VR* ptr, const VR* ptr_scal_rows, 
                                    const VR* ptr_scal_cols)
    {
        Integer N       = A.cols();

        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
            {
                VR tmp  = abs(ptr_A[pos]) * ptr_scal_rows[fr] * ptr_scal_cols[j];
                ptr[i]  = std::max(ptr[i], tmp);
            }

            ptr_A       += A_ld;
        };
    };

    static void maxabs_cols_scaled2(const Mat& A, VR* ptr, const VR* ptr_scal_rows, 
                                    const VR* ptr_scal_cols)
    {
        Integer N       = A.cols();

        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
            {
                VR tmp  = abs(ptr_A[pos]) * ptr_scal_rows[fr] * ptr_scal_cols[j];
                ptr[j]  = std::max(ptr[j], tmp);
            }

            ptr_A       += A_ld;
        };
    };

    static void sum_rows(const Mat& A, VR* ptr)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr[i]  += abs2(ptr_A[pos]);

            ptr_A       += A_ld;
        };
    };
    static void sum_cols(const Mat& A, VR* ptr)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr[j]  += abs2(ptr_A[pos]);

            ptr_A       += A_ld;
        };
    };

    static void sum_rows_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
            {
                if (i != j)
                    ptr[i]  += abs(ptr_A[pos]);
                else
                    ptr[i]  += abs(ptr_A[pos]) * scale_diag;
            }

            ptr_A       += A_ld;
        };
    };
    static void max_rows_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
            {
                if (i != j)
                    ptr[i]  = std::max(ptr[i], abs(ptr_A[pos]));
                else
                    ptr[i]  = std::max(ptr[i], abs(ptr_A[pos]) * scale_diag);
            }

            ptr_A       += A_ld;
        };
    };

    static void sum_cols_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
            {
                if (i != j)
                    ptr[j]  += abs(ptr_A[pos]);
                else
                    ptr[j]  += abs(ptr_A[pos]) * scale_diag;
            }

            ptr_A       += A_ld;
        };
    };
    static void max_cols_off(const Mat& A, VR* ptr, VR scale_diag)
    {
        Integer N       = A.rows();
        const VT* ptr_A = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
            {
                if (i != j)
                    ptr[j]  = std::max(ptr[j], abs(ptr_A[pos]));
                else
                    ptr[j]  = std::max(ptr[j], abs(ptr_A[pos]) * scale_diag);
            }

            ptr_A       += A_ld;
        };
    };

    static void scale_rows(Mat& A, const Integer* ptr, Integer mult)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr_A[pos]  = ldexp_helper<VT>::eval(ptr_A[pos],ptr[i]*mult);

            ptr_A       += A_ld;
        };
    };
    static void scale_cols(Mat& A, const Integer* ptr, Integer mult)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            Integer scale   = ptr[j] *mult;

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr_A[pos]  = ldexp_helper<VT>::eval(ptr_A[pos],scale);

            ptr_A       += A_ld;
        };
    };

    static void scale_rows(Mat& A, const VR* ptr)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr_A[pos]  = ptr_A[pos] * ptr[i];

            ptr_A       += A_ld;
        };
    };
    static void scale_cols_inv(Mat& A, const VR* ptr)
    {
        Integer N       = A.rows();
        VT* ptr_A       = A.rep_ptr();
        Integer A_ld    = A.ld();

        for (Integer j = 0; j < N; ++j)
        {
            Integer fr  = A.first_row(j);
            Integer lr  = A.last_row(j);
            Integer pos = A.first_elem_pos(j);

            VR scale    = VR(1.0)/ptr[j];

            for (Integer i = fr; i <= lr; ++i, ++pos)
                ptr_A[pos]  = ptr_A[pos] * scale;

            ptr_A       += A_ld;
        };
    };

};

//---------------------------------------------------------------
//                  balance_posdef
//---------------------------------------------------------------
template<class V, class S, bool Is_compl = md::is_complex<V>::value>
struct balance_posdef_str
{
    using Mat   = raw::Matrix<V,S>;
    static void eval(Matrix& ret, const Mat& mat, Real tol_sing, Real& d_max, Real& d_cond)
    {
        using Mat_D = raw::Matrix<V,struct_dense>;
        Mat_D Ac    = raw::converter<Mat_D, Mat>::eval(mat);
        return balance_posdef_str<V,struct_dense>::eval(ret, Ac, tol_sing, d_max, d_cond);
    };
};

template<class V>
struct balance_posdef_str<V,struct_dense, true>
{    
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,struct_dense>;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(Matrix& ret, const Mat& mat, Real tol_sing, Real& d_max, Real& d_cond)
    {
        Mat_R mat_r = raw::converter<Mat_R, Mat>::eval(mat);
        return balance_posdef_str<VR, struct_dense>::eval(ret, mat_r, tol_sing, d_max, d_cond);
    };
};

template<class V>
struct balance_posdef_str<V,struct_dense, false>
{
    using Mat       = raw::Matrix<V,struct_dense>;
    static void eval(matcl::Matrix& ret, const Mat& mat, Real tol_sing, Real& d_max, Real& d_cond)
    {        
        Integer M       = mat.rows();
        Integer N       = mat.cols();

        if (N == 0 || M == 0)
        {
            d_max   = 0;
            d_cond  = 1;
            ret     = Matrix(mat,false);
            return;
        };

        if (N != 1)
            throw error::vector_required(M, N);

        Mat D           = mat.make_unique();
        V* ptr_D        = D.ptr();

        d_max           = ptr_D[0];
        Real d_min      = d_max;
        tol_sing        = std::max(tol_sing, 0.0);
        V tol_t         = (V)tol_sing;

        for (Integer j = 0; j < N; ++j)
        {
            V d             = ptr_D[j];

            if (d <= tol_t)
                d           = 0;

            d_min           = std::min(d_min, Real(d));
            d_max           = std::max(d_max, Real(d));

            if (d == 0)
                ptr_D[j]    = 1.0;
            else
                ptr_D[j]    = V(1) / sqrt(ptr_D[j]);
        };

        d_cond  = div_1(d_min, d_max);
        ret     = matcl::Matrix(D,false);
        return;
    };
};

//---------------------------------------------------------------
//                  balance_gen
//---------------------------------------------------------------

template<class V, class S>
struct balance_gen_str
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,S>;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(const Mat& A, Matrix& ret_R, Matrix& ret_C, Real tol_sing)
    {
        Integer M       = A.rows();
        Integer N       = A.cols();

        Mat_R R(A.get_type(), M, 1);
        Mat_R C(A.get_type(), N, 1);

        VR* ptr_R       = R.ptr();
        VR* ptr_C       = C.ptr();

        if (M == 0 || N == 0)
        {
            for (Integer i = 0; i < M; ++i)
                ptr_R[i]    = VR(1.0);

            for (Integer i = 0; i < N; ++i)
                ptr_C[i]    = VR(1.0);

            ret_R   = Matrix(R,false);
            ret_C   = Matrix(C,false);
            return;
        };

        tol_sing        = std::max(tol_sing, 0.0);
        VR tol_t        = (VR)tol_sing;

        for (Integer i = 0; i < M; ++i)
            ptr_R[i]    = VR(0.0);

        for (Integer i = 0; i < N; ++i)
            ptr_C[i]    = VR(0.0);

        //Find the maximum element in each row.
        matrix_balance_func<Mat>::maxabs_rows(A, ptr_R);

        const VR SMLNUM = lapack::lamch<VR>("S");
        const VR BIGNUM = VR(1.0) / SMLNUM;

        // Invert the scale factors.
        for (Integer i = 0; i < M; ++i)
        {
            if (ptr_R[i] <= tol_t)
                ptr_R[i]    = VR(1.0);
            else
                ptr_R[i]    = VR(1.0) / std::min( std::max(ptr_R[i], SMLNUM), BIGNUM );
        };

        // Find the maximum element in each column, assuming the row scaling computed above.
        matrix_balance_func<Mat>::maxabs_cols_scaled(A, ptr_C, ptr_R);

        // Invert the scale factors.
        for (Integer i = 0; i < M; ++i)
        {
            if (ptr_C[i] <= tol_t)
                ptr_C[i]    = VR(1.0);
            else
                ptr_C[i]    = VR(1.0) / std::min( std::max(ptr_C[i], SMLNUM), BIGNUM );
        };

        ret_R   = Matrix(R,false);
        ret_C   = Matrix(C,false);
        return;
    };
};

//---------------------------------------------------------------
//                  balance_inf
//---------------------------------------------------------------

template<class V, class S>
struct balance_inf_str
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,S>;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(const Mat& A, Matrix& ret_R, Matrix& ret_C, Real tol_sing, const options& opts)
    {
        Integer M       = A.rows();
        Integer N       = A.cols();

        const Integer max_it  = opts.get_option<Integer>(opt::linsolve::balancing_iter());

        Mat_R R(A.get_type(), M, 1);
        Mat_R C(A.get_type(), N, 1);

        VR* ptr_R       = R.ptr();
        VR* ptr_C       = C.ptr();

        for (Integer i = 0; i < M; ++i)
            ptr_R[i]    = VR(1.0);

        for (Integer i = 0; i < N; ++i)
            ptr_C[i]    = VR(1.0);

        if (M == 0 || N == 0 || max_it <= 0)
        {
            ret_R   = Matrix(R,false);
            ret_C   = Matrix(C,false);
            return;
        };

        Mat_R work(A.get_type(), N + M, 1);
        VR* ptr_WR  = work.ptr();
        VR* ptr_WC  = work.ptr() + M;
        
        tol_sing        = std::max(tol_sing, 0.0);
        VR tol_s        = sqrt((VR)tol_sing);

        VR tol_min      = constants::eps<VR>();
        VR tol          = (VR)opts.get_option<Real>(opt::linsolve::balancing_tol());

        if (tol < VR(0.0))
            tol         = VR(1e-3);
        else
            tol         = std::min(tol, tol_min);

        Integer it      = 1;
        VR eps_r, eps_c;

        for (;; ++it)
        {
            for (Integer i = 0; i < M; ++i)
                ptr_WR[i]   = VR(0.0);

            for (Integer i = 0; i < N; ++i)
                ptr_WC[i]   = VR(0.0);

            //Find the maximum element in each row and columns.
            matrix_balance_func<Mat>::maxabs_rows_scaled2(A, ptr_WR, ptr_R, ptr_C);
            matrix_balance_func<Mat>::maxabs_cols_scaled2(A, ptr_WC, ptr_R, ptr_C);

            //update scaling
            for (Integer i = 0; i < M; ++i)
            {
                VR val      = sqrt(ptr_WR[i]);

                if (val <= tol_s)
                    val     = VR(1.0);
                else
                    val     = VR(1.0) / val;

                ptr_R[i]    = ptr_R[i] * val;
            };

            for (Integer i = 0; i < N; ++i)
            {
                VR val      = sqrt(ptr_WC[i]);

                if (val <= tol_s)
                    val     = VR(1.0);
                else
                    val     = VR(1.0) / val;

                ptr_C[i]    = ptr_C[i] * val;
            };

            // check convergence
            eps_r           = VR(0.0);
            for (Integer i = 0; i < M; ++i)
                eps_r       = std::max(eps_r, abs(ptr_WR[i] - VR(1.0)));

            eps_c           = VR(0.0);
            for (Integer i = 0; i < N; ++i)
                eps_c       = std::max(eps_c, abs(ptr_WC[i] - VR(1.0)));

            if (eps_r <= tol && eps_c <= tol)
                break;
            if (it >= max_it)
                break;
        };

        ret_R   = Matrix(R,false);
        ret_C   = Matrix(C,false);
        return;
    };
};

//---------------------------------------------------------------
//                  balance_sym
//---------------------------------------------------------------

template<class V, class S>
struct balance_sym_str
{
    using VR    = typename md::real_type<V>::type;
    using Mat   = raw::Matrix<V,S>;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(const Mat& A, Matrix& ret_D, Real tol_sing)
    {
        Integer N       = A.cols();

        Mat_R D(A.get_type(), N, 1);

        VR* ptr_D       = D.ptr();

        if (N == 0)
        {
            ret_D   = Matrix(D,false);
            return;
        };
        
        tol_sing    = std::max(tol_sing, 0.0);
        VR tol_t    = sqrt((VR)tol_sing);

        for (Integer j = 0; j < N; ++j)
        {
            VR Tii          = matrix_balance_func<Mat>::abs_diag(A, j);   
            VR Tij          = matrix_balance_func<Mat>::maxabs_superdiags(A, j, ptr_D);
            VR d            = std::max(sqrt(Tii), Tij);

            if (d <= tol_t)
                d           = VR(1.0);
            else
                d           = VR(1.0) / d;

            ptr_D[j]        = d;
        };

        ret_D   = Matrix(D,false);
        return;
    };
};

//---------------------------------------------------------------
//                  balance_eig
//---------------------------------------------------------------

template<class VT, class Mat1>
struct balance2_eig_mat_impl
{
    static void eval(const Mat1& A0, mat_tup_2& ret)
    {
        using VR        = typename md::real_type<VT>::type;
        using Mat_DR    = raw::Matrix<VR,struct_dense>;

        Integer N       = A0.rows();

        Mat_DR D        = Mat_DR(A0.get_type(), VR(1.0), N, 1);

        using workspace     = matcl::pod_workspace<VR>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        workspace WORK      = workspace(2*N);
        iworkspace IWORK    = iworkspace(N);
        VR* ptr_WORK_r      = WORK.ptr();
        VR* ptr_WORK_c      = ptr_WORK_r + N;
        Integer* ptr_IWORK  = IWORK.ptr();

        VR* SCALE           = D.ptr();

        Mat1 A              = A0.make_unique();
        A.get_struct().reset_value();

        Integer max_iter    = 10;

        // uses scaling technique proposed by Parlett-Reinsch

        for (Integer iter = 1; iter <= max_iter; ++iter)
        {
            Integer emax    = 0;
            Integer emin    = 0;

            for (Integer i = 0; i < N; ++i)
            {
                ptr_WORK_r[i] = VR(0.0);
                ptr_WORK_c[i] = VR(0.0);
            }

            matrix_balance_func<Mat1>::sum_rows(A, ptr_WORK_r);
            matrix_balance_func<Mat1>::sum_cols(A, ptr_WORK_c);

            for (Integer i = 0; i < N; ++i)
            {
                VR dr       = ptr_WORK_r[i];
                VR dc       = ptr_WORK_c[i];

                if (dr == VR(0.0))
                    continue;

                VR d        = dc/dr;
                Integer e   = lround(std::log2(d)/4);
                ptr_IWORK[i]= e;

                emax        = std::max(emax,e);
                emin        = std::min(emin,e);

                SCALE[i]    = ldexp_helper<VR>::eval(SCALE[i], -e);
            };        

            matrix_balance_func<Mat1>::scale_rows(A, ptr_IWORK, 1);
            matrix_balance_func<Mat1>::scale_cols(A, ptr_IWORK, -1);

            // Stop if norms are all between 1/2 and 2
            if (emax <= 1 && emin >= -1)
                break;
        };

        ret = mat_tup_2(Matrix(A,true), Matrix(D,false));
        return;
    };

    static void eval_off(const Mat1& A0, mat_tup_2& ret)
    {
        using VR        = typename md::real_type<VT>::type;
        using Mat_DR    = raw::Matrix<VR,struct_dense>;

        Integer N       = A0.rows();

        Mat_DR D        = Mat_DR(A0.get_type(), VR(1.0), N, 1);

        using workspace     = matcl::pod_workspace<VR>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        workspace WORK      = workspace(2*N);

        VR* ptr_WORK_r      = WORK.ptr();
        VR* ptr_WORK_c      = ptr_WORK_r + N;

        VR* SCALE           = D.ptr();

        //avoid explosion of scaling; 
        VR eps              = VR(1e-2) / VR(N);
        eps                 = eps * eps;
        VR tol              = VR(1e-1);
        VR inf              = constants::inf<VR>();

        Mat1 A              = A0.make_unique();
        A.get_struct().reset_value();

        Integer max_iter    = 10;

        for (Integer iter = 1; iter <= max_iter; ++iter)
        {
            VR emax         = VR(0.0);
            VR emin         = inf;

            for (Integer i = 0; i < N; ++i)
            {
                ptr_WORK_r[i] = VR(0.0);
                ptr_WORK_c[i] = VR(0.0);
            }

            matrix_balance_func<Mat1>::max_rows_off(A, ptr_WORK_r, VR(0.0));
            matrix_balance_func<Mat1>::max_cols_off(A, ptr_WORK_c, VR(0.0));

            for (Integer i = 0; i < N; ++i)
            {
                VR dr       = ptr_WORK_r[i];
                VR dc       = ptr_WORK_c[i];

                if (dr == VR(0.0) || dc == VR(0.0))
                {
                    VR d    =  matrix_balance_func<Mat1>::abs_diag(A, i);

                    if (iter == 1)
                    {
                        if (d == VR(0.0))
                            d   = std::max(dr, dc);

                        if (d == VR(0.0))
                            d   = VR(1.0);

                        if (dr == VR(0.0))
                            dr  = eps * d;
                        if (dc == VR(0.0))
                            dc  = eps * d;
                    }
                    else
                    {
                        dr  = VR(1.0);
                        dc  = VR(1.0);
                    };
                };

                VR d            = dc/dr;
                VR e            = sqrt(d);
                ptr_WORK_r[i]   = e;

                emax        = std::max(emax,e);
                emin        = std::min(emin,e);

                SCALE[i]    = SCALE[i] / e;
            };        

            matrix_balance_func<Mat1>::scale_rows(A, ptr_WORK_r);
            matrix_balance_func<Mat1>::scale_cols_inv(A, ptr_WORK_r);

            // Stop if norms are all between 1/2 and 2
            if (abs(emax - VR(1.0)) <= tol && abs(emin - VR(1.0)) <= tol)
                break;
        };

        ret = mat_tup_2(Matrix(A,true), Matrix(D,false));
        return;
    };
};

template<class Mat1>
struct balance2_eig_mat_impl<Object,Mat1>
{
    static void eval(const Mat1&, mat_tup_2&)
    {
        throw error::object_value_type_not_allowed("balance_eig");
    };
    static void eval_off(const Mat1&, mat_tup_2&)
    {
        throw error::object_value_type_not_allowed("balance_eig_off");
    };
};

template<class VT, class Mat1, class Mat2>
struct balance_eig_mat_impl
{
    static void eval(const Mat1& A0, const Mat2& B0, mat_tup_4& ret)
    {
        using VR        = typename md::real_type<VT>::type;
        using Mat_DR    = raw::Matrix<VR,struct_dense>;

        Integer N       = A0.rows();

        Mat_DR Dl       = Mat_DR(A0.get_type(), VR(1.0), N, 1);
        Mat_DR Dr       = Mat_DR(A0.get_type(), VR(1.0), N, 1);

        using workspace     = matcl::pod_workspace<VR>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        workspace WORK      = workspace(N);
        iworkspace IWORK    = iworkspace(N);
        VR* ptr_WORK        = WORK.ptr();
        Integer* ptr_IWORK  = IWORK.ptr();

        VR* LSCALE          = Dl.ptr();
        VR* RSCALE          = Dr.ptr();

        Mat1 A              = A0.make_unique();
        Mat2 B              = B0.make_unique();

        A.get_struct().reset_value();
        B.get_struct().reset_value();

        Integer max_iter    = 10;

        // uses scaling technique proposed by Lemonnier and Van Dooren
        //
        //  Lemonnier, Damien, and Paul Van Dooren. "Balancing regular matrix pencils." 
        //  SIAM journal on matrix analysis and applications 28.1 (2006): 253-263.

        for (Integer iter = 1; iter <= max_iter; ++iter)
        {
            Integer emax    = 0;
            Integer emin    = 0;

            for (Integer i = 0; i < N; ++i)
                ptr_WORK[i] = VR(0.0);

            // scale the rows of abs2(A) + abs2(B) to have approximate row sum 1
            matrix_balance_func<Mat1>::sum_rows(A, ptr_WORK);
            matrix_balance_func<Mat2>::sum_rows(B, ptr_WORK);

            for (Integer i = 0; i < N; ++i)
            {
                VR d        = ptr_WORK[i];

                if (d == VR(0.0))
                    continue;

                Integer e   = -lround(std::log2(d)/2);
                ptr_IWORK[i]= e;

                emax        = std::max(emax,e);
                emin        = std::min(emin,e);

                LSCALE[i]   = ldexp_helper<VR>::eval(LSCALE[i], e);
            };        

            matrix_balance_func<Mat1>::scale_rows(A, ptr_IWORK,1);
            matrix_balance_func<Mat2>::scale_rows(B, ptr_IWORK,1);

            // scale the columns of abs2(A) + abs2(B) to have approximate column sum 1

            for (Integer i = 0; i < N; ++i)
                ptr_WORK[i] = VR(0.0);

            matrix_balance_func<Mat1>::sum_cols(A, ptr_WORK);
            matrix_balance_func<Mat2>::sum_cols(B, ptr_WORK);

            for (Integer i = 0; i < N; ++i)
            {
                VR d        = ptr_WORK[i];

                if (d == VR(0.0))
                    continue;

                Integer e   = -lround(std::log2(d)/2);
                ptr_IWORK[i]= e;

                emax        = std::max(emax,e);
                emin        = std::min(emin,e);

                RSCALE[i]   = ldexp_helper<VR>::eval(RSCALE[i], e);
            };

            matrix_balance_func<Mat1>::scale_cols(A, ptr_IWORK,1);
            matrix_balance_func<Mat2>::scale_cols(B, ptr_IWORK,1);

            // Stop if norms are all between 1/2 and 2
            if (emax <= 1 && emin >= -1)
                break;
        };

        ret = mat_tup_4(Matrix(A,true), Matrix(B,true), Matrix(Dl,false), Matrix(Dr,false));
        return;
    };
};

template<class Mat1, class Mat2>
struct balance_eig_mat_impl<Object, Mat1, Mat2>
{
    static void eval(const Mat1&, const Mat2&, mat_tup_4&)
    {
        throw error::object_value_type_not_allowed("balance_eig");
    };
};

//---------------------------------------------------------------
//                  CONVERSIONS
//---------------------------------------------------------------

template<class V, class S>
struct balance_posdef_val
{
    using Mat   = raw::Matrix<V,S>;
    static void eval(Matrix& ret, const Mat& mat, Real tol_sing, Real& d_max, Real& d_cond)
    {
        return balance_posdef_str<V,S>::eval(ret, mat, tol_sing, d_max, d_cond);
    };
};

template<class S>
struct balance_posdef_val<Integer,S>
{
    using Mat   = raw::Matrix<Integer,S>;
    static void eval(Matrix& ret, const Mat& mat, Real tol_sing, Real& d_max, Real& d_cond)
    {
        using Mat_D = raw::Matrix<Real, struct_dense>;
        Mat_D Ac    = raw::converter<Mat_D, Mat>::eval(mat);

        return balance_posdef_val<Real,struct_dense>::eval(ret, Ac, tol_sing, d_max, d_cond);
    };
};

template<class V, class S>
struct balance_gen_val
{
    using Mat   = raw::Matrix<V,S>;
    static void eval(const Mat& mat, Matrix& R, Matrix& C, Real tol_sing)
    {
        return balance_gen_str<V,S>::eval(mat, R, C, tol_sing);
    };
};

template<class S>
struct balance_gen_val<Integer,S>
{
    using Mat   = raw::Matrix<Integer,S>;
    static void eval(const Mat& mat, Matrix& R, Matrix& C, Real tol_sing)
    {
        using Mat_D = raw::Matrix<Real, S>;
        Mat_D Ac    = raw::converter<Mat_D, Mat>::eval(mat);

        return balance_gen_val<Real,S>::eval(Ac, R, C, tol_sing);
    };
};

template<class V, class S>
struct balance_inf_val
{
    using Mat   = raw::Matrix<V,S>;
    static void eval(const Mat& mat, Matrix& R, Matrix& C, Real tol_sing, const options& opts)
    {
        return balance_inf_str<V,S>::eval(mat, R, C, tol_sing, opts);
    };
};

template<class S>
struct balance_inf_val<Integer,S>
{
    using Mat   = raw::Matrix<Integer,S>;
    static void eval(const Mat& mat, Matrix& R, Matrix& C, Real tol_sing, const options& opts)
    {
        using Mat_D = raw::Matrix<Real, S>;
        Mat_D Ac    = raw::converter<Mat_D, Mat>::eval(mat);

        return balance_inf_val<Real,S>::eval(Ac, R, C, tol_sing, opts);
    };
};

template<class V, class S>
struct balance_sym_val
{
    using Mat   = raw::Matrix<V,S>;
    static void eval(const Mat& mat, Matrix& D, Real tol_sing)
    {
        return balance_sym_str<V,S>::eval(mat, D, tol_sing);
    };
};

template<class S>
struct balance_sym_val<Integer,S>
{
    using Mat   = raw::Matrix<Integer,S>;
    static void eval(const Mat& mat, Matrix& D, Real tol_sing)
    {
        using Mat_D = raw::Matrix<Real, S>;
        Mat_D Ac    = raw::converter<Mat_D, Mat>::eval(mat);

        return balance_sym_val<Real,S>::eval(Ac, D, tol_sing);
    };
};

//---------------------------------------------------------------
//                  visitors
//---------------------------------------------------------------
struct balance_posdef_vis : public extract_type_switch<void,balance_posdef_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, Real tol_sing, Real& d_max, Real& d_cond)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::balance_posdef_val<V,S>::eval(ret, mat, tol_sing, d_max, d_cond);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix& h, const T& A, Matrix& ret, Real tol_sing, 
                            Real& d_max, Real& d_cond)
    {
        using Mat   = raw::Matrix<T,struct_dense>;
        Mat Ac      = raw::converter<Mat,T>::eval(A);
        return eval<Mat>(h, Ac, ret, tol_sing, d_max, d_cond);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, Real, Real&, Real&)
    {
        throw error::object_value_type_not_allowed("balance_posdef");
    };

    static void eval_scalar(const matcl::Matrix&, const Object&, Matrix&, Real, Real&, Real&)
    {
        throw error::object_value_type_not_allowed("balance_posdef");
    };
};

struct balance_gen_vis : public extract_type_switch<void,balance_gen_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& R, Matrix& C, Real tol)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::balance_gen_val<V,S>::eval(mat, R, C, tol);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix& h, const T& A, Matrix& R, Matrix& C, Real tol)
    {
        using Mat   = raw::Matrix<T,struct_dense>;
        Mat Ac      = raw::converter<Mat,T>::eval(A);
        return eval<Mat>(h, Ac, R, C, tol);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, Matrix&, Real)
    {
        throw error::object_value_type_not_allowed("balance_gen");
    };

    static void eval_scalar(const matcl::Matrix&, const Object&, Matrix&, Matrix&, Real)
    {
        throw error::object_value_type_not_allowed("balance_gen");
    };
};

struct balance_inf_vis : public extract_type_switch<void,balance_inf_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& R, Matrix& C, Real tol,
                     const options& opts)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::balance_inf_val<V,S>::eval(mat, R, C, tol, opts);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix& h, const T& A, Matrix& R, Matrix& C,
                            Real tol, const options& opts)
    {
        using Mat   = raw::Matrix<T,struct_dense>;
        Mat Ac      = raw::converter<Mat,T>::eval(A);
        return eval<Mat>(h, Ac, R, C, tol, opts);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, Matrix&, 
                     Real, const options&)
    {
        throw error::object_value_type_not_allowed("balance_inf");
    };

    static void eval_scalar(const matcl::Matrix&, const Object&, Matrix&, Matrix&, 
                            Real, const options&)
    {
        throw error::object_value_type_not_allowed("balance_inf");
    };
};

struct balance_sym_vis : public extract_type_switch<void,balance_sym_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& D, Real min_diag)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return details::balance_sym_val<V,S>::eval(mat, D, min_diag);
    };

    template<class T>
    static void eval_scalar(const matcl::Matrix& h, const T& A, Matrix& D, Real min_diag)
    {
        using Mat   = raw::Matrix<T,struct_dense>;
        Mat Ac      = raw::converter<Mat,T>::eval(A);
        return eval<Mat>(h, Ac, D, min_diag);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, Real)
    {
        throw error::object_value_type_not_allowed("balance_sym");
    };

    static void eval_scalar(const matcl::Matrix&, const Object&, Matrix&, Real)
    {
        throw error::object_value_type_not_allowed("balance_sym");
    };
};

struct balance_eig_vis : public extract_type2_switch<void,balance_eig_vis,mr::val_type_corrector_int>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, mat_tup_4& ret)
    {
        using V1    = typename T1::value_type;
        using V2    = typename T2::value_type;
        using S1    = typename T1::struct_type;
        using S2    = typename T2::struct_type;

        using VT    = typename md::unify_types2<V1,V2,Float>::type;
        using Mat1  = raw::Matrix<VT, S1>;
        using Mat2  = raw::Matrix<VT, S2>;

        const Mat1& Ac  = raw::converter<Mat1,T1>::eval(A);
        const Mat2& Bc  = raw::converter<Mat2,T2>::eval(B);

        return balance_eig_mat_impl<VT, Mat1, Mat2>::eval(Ac, Bc, ret);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, mat_tup_4& ret)
    {
        using DM    = raw::Matrix<T2,struct_dense>;
        DM BD       = DM(ti::get_ti(B), B, 1,1);
        return balance_eig_vis::eval_mat_mat<T1,DM>(A, BD, ret);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, mat_tup_4& ret)
    {
        using DM    = raw::Matrix<T1,struct_dense>;
        DM AD       = DM(ti::get_ti(A), A, 1,1);
        return balance_eig_vis::eval_mat_mat<DM,T2>(AD, B, ret);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, mat_tup_4& ret)
    {
        using DM1   = raw::Matrix<T1,struct_dense>;
        using DM2   = raw::Matrix<T2,struct_dense>;
        DM1 AD      = DM1(ti::get_ti(A), A, 1,1);
        DM2 BD      = DM2(ti::get_ti(B), B, 1,1);

        return balance_eig_vis::eval_mat_mat<DM1,DM2>(AD, BD, ret);
    };
};

struct balance2_eig_vis : public extract_type_switch<void,balance2_eig_vis,true>
{
    template<class T1>
    static void eval(const Matrix&, const T1& A, mat_tup_2& ret)
    {
        using V1    = typename T1::value_type;
        using S1    = typename T1::struct_type;

        using VT    = typename md::unify_types<V1,Float>::type;
        using Mat1  = raw::Matrix<VT, S1>;

        const Mat1& Ac  = raw::converter<Mat1,T1>::eval(A);

        return balance2_eig_mat_impl<VT, Mat1>::eval(Ac, ret);
    };

    template<class T1>
    static void eval_scalar(const Matrix& h, const T1& A, mat_tup_2& ret)
    {
        using DM1   = raw::Matrix<T1,struct_dense>;
        DM1 AD      = DM1(ti::get_ti(A), A, 1,1);

        return balance2_eig_vis::eval<DM1>(h, AD, ret);
    };
};

struct balance2_eig_off_vis : public extract_type_switch<void,balance2_eig_off_vis,true>
{
    template<class T1>
    static void eval(const Matrix&, const T1& A, mat_tup_2& ret)
    {
        using V1    = typename T1::value_type;
        using S1    = typename T1::struct_type;

        using VT    = typename md::unify_types<V1,Float>::type;
        using Mat1  = raw::Matrix<VT, S1>;

        const Mat1& Ac  = raw::converter<Mat1,T1>::eval(A);

        return balance2_eig_mat_impl<VT, Mat1>::eval_off(Ac, ret);
    };

    template<class T1>
    static void eval_scalar(const Matrix& h, const T1& A, mat_tup_2& ret)
    {
        using DM1   = raw::Matrix<T1,struct_dense>;
        DM1 AD      = DM1(ti::get_ti(A), A, 1,1);

        return balance2_eig_vis::eval<DM1>(h, AD, ret);
    };
};

//---------------------------------------------------------------
//                  TOP LEVEL IMPL
//---------------------------------------------------------------

static void balance_eig_impl(const Matrix& A, const Matrix& B, mat_tup_4& ret)
{
    if (!A.is_square() || !B.is_square() || A.rows() != B.rows())
        throw error::error_size_geig(A.rows(),A.cols(), B.rows(), B.cols());

    bool is_obj         = (A.get_value_code() == value_code::v_object)
                        ||(B.get_value_code() == value_code::v_object);

    if (is_obj == true)
        throw error::object_value_type_not_allowed("balance_eig");

    Integer N           = A.rows();

    bool isv            = A.all_finite() && B.all_finite();

    value_code vt0      = matrix_traits::unify_value_types(A.get_value_code(), B.get_value_code());
    value_code vt       = matrix_traits::unify_value_types(vt0, value_code::v_float);

    if (isv == false)
    {
        Matrix AA       = details::make_nan_matrix(N, N, vt);
        Matrix BB       = AA;
        Matrix Dl       = details::make_nan_matrix(N, 1, vt);
        Matrix Dr       = details::make_nan_matrix(N, 1, vt);

        ret             = mat_tup_4(AA,BB,Dl, Dr);
        return;
    }

    if (A.structural_nnz() + B.structural_nnz() == 0 
        || ( A.get_struct(). is_id() || is_unitary(A.get_struct()))
                && ( B.get_struct(). is_id() || is_unitary(B.get_struct()))
       )
    {
        Matrix AA       = convert(A, matrix_traits::get_matrix_type(vt,A.get_struct_code()));
        Matrix BB       = convert(B, matrix_traits::get_matrix_type(vt,B.get_struct_code()));
        Matrix Dl       = ones(N,1,vt);
        Matrix Dr       = Dl;
    
        ret             = mat_tup_4(AA,BB,Dl, Dr);
        return;
    }

    return details::balance_eig_vis::make(A,B,ret);
};

static void balance_eig_impl(const Matrix& A, mat_tup_2& ret)
{
    if (!A.is_square())
        throw error::error_size_eig(A.rows(),A.cols());

    bool is_obj         = (A.get_value_code() == value_code::v_object);

    if (is_obj == true)
        throw error::object_value_type_not_allowed("balance_eig");

    Integer N           = A.rows();
    
    bool isv            = A.all_finite();

    value_code vt0      = A.get_value_code();
    value_code vt       = matrix_traits::unify_value_types(vt0, value_code::v_float);

    if (isv == false)
    {
        Matrix AA       = details::make_nan_matrix(N, N, vt);
        Matrix D        = details::make_nan_matrix(N, 1, vt);

        ret             = mat_tup_2(AA,D);
        return;
    }

    if ( A.structural_nnz() == 0 ||  A.get_struct(). is_id() || is_unitary(A.get_struct()))
    {
        Matrix AA       = convert(A, matrix_traits::get_matrix_type(vt,A.get_struct_code()));
        Matrix D        = ones(N,1,vt);
    
        ret             = mat_tup_2(AA,D);
        return;
    }

    return details::balance2_eig_vis::make<const Matrix&>(A,ret);
};

static void balance_eig_off_impl(const Matrix& A, mat_tup_2& ret)
{
    if (!A.is_square())
        throw error::error_size_eig(A.rows(),A.cols());

    bool is_obj         = (A.get_value_code() == value_code::v_object);

    if (is_obj == true)
        throw error::object_value_type_not_allowed("balance_eig");

    Integer N           = A.rows();
    
    bool isv            = A.all_finite();

    value_code vt0      = A.get_value_code();
    value_code vt       = matrix_traits::unify_value_types(vt0, value_code::v_float);

    if (isv == false)
    {
        Matrix AA       = details::make_nan_matrix(N, N, vt);
        Matrix D        = details::make_nan_matrix(N, 1, vt);

        ret             = mat_tup_2(AA,D);
        return;
    }

    if ( A.structural_nnz() == 0 ||  A.get_struct(). is_id() || is_unitary(A.get_struct()))
    {
        Matrix AA       = convert(A, matrix_traits::get_matrix_type(vt,A.get_struct_code()));
        Matrix D        = ones(N,1,vt);
    
        ret             = mat_tup_2(AA,D);
        return;
    }

    return details::balance2_eig_off_vis::make<const Matrix&>(A,ret);
};

};};

namespace matcl
{

//--------------------------------------------------------------------
//              LINEAR EQUATIONS
//--------------------------------------------------------------------
mat_tup_2 matcl::balance_posdef2(const Matrix& A, bool pow2, Real tol_sing)
{
    Matrix D    = balance_posdef(A, pow2, tol_sing);
    Matrix AA   = scale_rowscols(A, D, D);
    return mat_tup_2(AA,D);
};

mat_tup_2 matcl::balance_posdef2(Matrix&& A, bool pow2, Real tol_sing)
{
    Matrix D    = balance_posdef(A, pow2, tol_sing);
    Matrix AA   = scale_rowscols(std::move(A), D, D);
    return mat_tup_2(AA,D);
};

Matrix matcl::balance_posdef(const Matrix& A, bool pow2, Real tol_sing)
{
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    Matrix D = abs(get_diag(A));

    Matrix ret;
    Real d_max, d_cond;
    details::balance_posdef_vis::make<const Matrix&>(D, ret, tol_sing, d_max, d_cond);
    round_to_pow2(ret, pow2);
    return ret;
};

mat_tup_2 matcl::balance_sym2(const Matrix& A, bool pow2, Real tol_sing)
{
    Matrix D    = balance_sym(A, pow2, tol_sing);
    Matrix AA   = scale_rowscols(A, D, D);
    return mat_tup_2(AA,D);
};

mat_tup_2 matcl::balance_sym2(Matrix&& A, bool pow2, Real tol_sing)
{
    Matrix D    = balance_sym(A, pow2, tol_sing);
    Matrix AA   = scale_rowscols(std::move(A), D, D);
    return mat_tup_2(AA,D);
};

Matrix matcl::balance_sym(const Matrix& A, bool pow2, Real tol_sing)
{
    if (A.rows() != A.cols())
        throw error::square_matrix_required(A.rows(), A.cols());

    Matrix D;
    details::balance_sym_vis::make<const Matrix&>(A, D, tol_sing);
    round_to_pow2(D, pow2);
    return D;
};

mat_tup_2 matcl::balance_gen(const Matrix& A, bool pow2, Real tol_sing)
{
    Matrix R, C;
    details::balance_gen_vis::make<const Matrix&>(A, R, C, tol_sing);
    round_to_pow2(R, pow2);
    round_to_pow2(C, pow2);
    return mat_tup_2(R,C);
};

mat_tup_3 matcl::balance_gen2(const Matrix& A, bool pow2, Real tol_sing)
{
    Matrix R, C;
    tie(R,C)    = balance_gen(A, pow2, tol_sing);
    Matrix B    = scale_rowscols(A, R, C);
    return mat_tup_3(B,R,C);
};

mat_tup_3 matcl::balance_gen2(Matrix&& A, bool pow2, Real tol_sing)
{
    Matrix R, C;
    tie(R,C)    = balance_gen(A, pow2, tol_sing);
    Matrix B    = scale_rowscols(std::move(A), R, C);
    return mat_tup_3(B,R,C);
};

mat_tup_2 matcl::balance_inf(const Matrix& A, bool pow2, Real tol_sing, const options& opts)
{
    Matrix R, C;
    details::balance_inf_vis::make<const Matrix&>(A, R, C, tol_sing, opts);
    round_to_pow2(R, pow2);
    round_to_pow2(C, pow2);
    return mat_tup_2(R,C);
};

mat_tup_3 matcl::balance_inf2(const Matrix& A, bool pow2, Real tol_sing, const options& opts)
{
    Matrix R, C;
    tie(R,C)    = balance_inf(A,pow2,tol_sing, opts);
    Matrix B    = scale_rowscols(A, R, C);
    return mat_tup_3(B,R,C);
};

mat_tup_3 matcl::balance_inf2(Matrix&& A, bool pow2, Real tol_sing, const options& opts)
{
    Matrix R, C;
    tie(R,C)    = balance_inf(A,pow2,tol_sing,opts);
    Matrix B    = scale_rowscols(std::move(A), R, C);
    return mat_tup_3(B,R,C);
};

//--------------------------------------------------------------------
//                      EIGENVALUES
//--------------------------------------------------------------------
mat_tup_4 matcl::balance_eig(const Matrix& A0, const Matrix& B0)
{
    Matrix A(A0);
    Matrix B(B0);

    mat_tup_4 ret;
    details::balance_eig_impl(A,B,ret);
    return ret;
};
mat_tup_4 matcl::balance_eig(Matrix&& A0, const Matrix& B0)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    mat_tup_4 ret;
    details::balance_eig_impl(A,B,ret);
    return ret;
};
mat_tup_4 matcl::balance_eig(const Matrix& A0, Matrix&& B0)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    mat_tup_4 ret;
    details::balance_eig_impl(A,B,ret);
    return ret;
};
mat_tup_4 matcl::balance_eig(Matrix&& A0, const Matrix&& B0)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    mat_tup_4 ret;
    details::balance_eig_impl(A,B,ret);
    return ret;
};

mat_tup_2 matcl::balance_eig(const Matrix& A0)
{
    Matrix A(A0);

    mat_tup_2 ret;
    details::balance_eig_impl(A,ret);
    return ret;
};

mat_tup_2 matcl::balance_eig(Matrix&& A0)
{
    Matrix A(std::move(A0));

    mat_tup_2 ret;
    details::balance_eig_impl(A,ret);
    return ret;
};

mat_tup_2 matcl::balance_offdiag(const Matrix& A0)
{
    Matrix A(A0);

    mat_tup_2 ret;
    details::balance_eig_off_impl(A,ret);
    return ret;
};

mat_tup_2 matcl::balance_offdiag(Matrix&& A0)
{
    Matrix A(std::move(A0));

    mat_tup_2 ret;
    details::balance_eig_off_impl(A,ret);
    return ret;
};

};
