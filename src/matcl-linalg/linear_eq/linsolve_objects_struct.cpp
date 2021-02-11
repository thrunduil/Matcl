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

#include "linsolve_objects_struct.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-linalg/decompositions/qr.h"
#include "linsolve_utils.inl"
#include "matcl-linalg/utils/linalg_utils.h"

#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-linalg/norms_error/norm.h"

namespace matcl { namespace details
{

template<class Func>
void eval_dense(const Matrix& A, Func& f)
{
    value_code vc = A.get_value_code();

    switch(vc)
    {
        case value_code::v_integer:
        {
            using type      = Integer;
            using Mat       = raw::Matrix<type,struct_dense>;
            const Mat& rep  = A.impl<Mat>();
            f.eval(rep);
            return;
        }
        case value_code::v_real:
        {
            using type      = Real;
            using Mat       = raw::Matrix<type,struct_dense>;
            const Mat& rep  = A.impl<Mat>();
            f.eval(rep);
            return;
        }
        case value_code::v_float:
        {
            using type      = Float;
            using Mat       = raw::Matrix<type,struct_dense>;
            const Mat& rep  = A.impl<Mat>();
            f.eval(rep);
            return;
        }
        case value_code::v_complex:
        {
            using type      = Complex;
            using Mat       = raw::Matrix<type,struct_dense>;
            const Mat& rep  = A.impl<Mat>();
            f.eval(rep);
            return;
        }
        case value_code::v_float_complex:
        {
            using type      = Float_complex;
            using Mat       = raw::Matrix<type,struct_dense>;
            const Mat& rep  = A.impl<Mat>();
            f.eval(rep);
            return;
        }
        case value_code::v_object:
        {
            using type      = Object;
            using Mat       = raw::Matrix<type,struct_dense>;
            const Mat& rep  = A.impl<Mat>();
            f.eval(rep);
            return;
        }
        default:
        {
            matcl_assert(0,"invalid case");
            throw error::error_general("invalid case");
        }
    };
};

struct logdet_diag_functor
{
    Real m_det;

    logdet_diag_functor()
        :m_det(0.0)
    {};

    template<class Val>
    void eval(const raw::Matrix<Val,struct_dense>& A)
    {
        const Val* ptr  = A.ptr();
        Integer N       = A.cols();
        Integer M       = A.rows();
        Integer A_ld    = A.ld();

        Real ld         = 0.0;

        for (Integer i = 0; i < N; ++i)
        {
            for (Integer j = 0; j < M; ++j);
                ld      += log(abs(ptr[i]));

            ptr         += A_ld;
        };

        m_det           = ld;
    };

    void eval(const raw::Matrix<Object,struct_dense>&)
    {
        throw error::object_value_type_not_allowed("linsolve_obj_diag::log_det");
    };
};

static Real log_det_diag(const Matrix& d)
{
    logdet_diag_functor f;
    eval_dense(d, f);
    return f.m_det;
};

//--------------------------------------------------------------------------------
//                  inverse diagonal 2x2
//--------------------------------------------------------------------------------

template<class V>
struct inv_22
{
    using VR = typename md::real_type<V>::type;
    static void eval(V& d_11, V& d_12, V& d_21, V& d_22, VR& sig_max, VR& sig_min)
    {
        return matcl::inv_22(d_11, d_12, d_21, d_22, sig_max, sig_min);
    };
};

template<class V>
struct inv_diag22_cont
{
    using VR = typename md::real_type<V>::type;

    static void eval(Integer M, const V* d_0, const V* d_p1, const V* d_m1, 
                     V* rd_0, V* rd_p1, V* rd_m1, Integer A_ld, Integer res_ld, Real& logdet,
                     Real& sig_max_r, Real& sig_min_r, Real diag_tol, bool& modif)
    {
        VR  sig_max = VR(0.0);
        VR  sig_min = constants::inf<VR>();
        logdet      = 0.0;
        modif       = false;

        for (Integer i = 0; i < M; ++i)
        {
            bool block_2    = false;

            V d_00          = d_0[0];

            if (i < M - 1)
            {
                if (d_m1[0] != V(0.0) || d_p1[0] != V(0.0))
                    block_2 = true;
            };

            if (abs(d_00) < diag_tol)
            {
                d_00    = V(VR(diag_tol)) * (abs(d_00) == VR(0.0) ? VR(1.0) : sign(d_00));
                modif   = true;
            }

            if (block_2 == false)
            {
                if (d_00 == V(0.0))
                    throw error::error_singular();

                rd_0[0] = V(1.0) / d_00;
                logdet  += matcl::log(matcl::abs(d_00));

                sig_max = std::max(sig_max, matcl::abs(d_00));
                sig_min = std::min(sig_min, matcl::abs(d_00));

                if (i < M - 1)
                {
                    rd_m1[0]    = V(0.0);
                    rd_p1[0]    = V(0.0);
                };

                d_0             += A_ld;
                d_p1            += A_ld;
                d_m1            += A_ld;

                rd_0            += res_ld;
                rd_p1           += res_ld;
                rd_m1           += res_ld;
            }
            else
            {
                V d_01          = d_p1[0];
                V d_10          = d_m1[0];
                V d_11          = d_0[1*A_ld];

                VR sig_max_l, sig_min_l;

                //invert 2x2 block
                inv_22<V>::eval(d_00, d_01, d_10, d_11, sig_max_l, sig_min_l);

                if (sig_min_l < diag_tol)
                {
                    sig_min_l   = VR(diag_tol);
                    modif       = true;
                }
                if (sig_max_l < diag_tol)
                {
                    sig_max_l   = VR(diag_tol);
                    modif       = true;
                };

                if (sig_min_l == VR(0.0))
                    throw error::error_singular();

                sig_max         = std::max(sig_max, sig_max_l);
                sig_min         = std::min(sig_min, sig_min_l);

                logdet          += matcl::log(sig_max_l);
                logdet          += matcl::log(sig_min_l);

                rd_0[0]         = d_00;
                rd_m1[0]        = d_10;
                rd_p1[0]        = d_11;
                rd_0[1*res_ld]  = d_11;

                d_0             += A_ld;
                d_p1            += A_ld;
                d_m1            += A_ld;

                rd_0            += res_ld;
                rd_p1           += res_ld;
                rd_m1           += res_ld;

                i               += 1;

                if (i < M - 1)
                {
                    rd_m1[0]    = V(0.0);
                    rd_p1[0]    = V(0.0);

                    if (d_m1[0] != V(0.0) || d_p1[0] != V(0.0))
                        throw error::invalid_diagonal_22();
                };

                d_0             += A_ld;
                d_p1            += A_ld;
                d_m1            += A_ld;

                rd_0            += res_ld;
                rd_p1           += res_ld;
                rd_m1           += res_ld;
            };
        };

        if (sig_min == VR(0.0))
            throw error::error_singular();

        sig_max_r   = sig_max;
        sig_min_r   = sig_min;
    };
};

template<class V, class S>
struct linsolve_obj_diag_22_inv_str
{};

template<class V>
struct linsolve_obj_diag_22_inv_str<V,struct_sparse>
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    using Mat_B = raw::Matrix<V,struct_banded>;

    static void eval(Matrix& ret, const Mat& A, Real& logdet, Real& sig_max, Real& sig_min,
                     Real diag_tol, bool& modif)
    {
        Mat_B Ac    = raw::converter<Mat_B,Mat>::eval(A);
        return linsolve_obj_diag_22_inv_str<V, struct_banded>
                    ::eval(ret, Ac, logdet, sig_max, sig_min, diag_tol, modif);
    };
};

template<class V>
struct linsolve_obj_diag_22_inv_str<V,struct_dense>
{
    using Mat       = raw::Matrix<V,struct_dense>;
    using Mat_B     = raw::Matrix<V,struct_banded>;

    static void eval(Matrix& ret, const Mat& A, Real& logdet, Real& sig_max, Real& sig_min, 
                     Real diag_tol, bool& modif)
    {
        Integer M   = A.rows();

        Mat_B res(A.get_type(), M, M, -1, 1);

        const V* d_0    = A.ptr();
        const V* d_p1   = A.ptr() + 1;
        const V* d_m1   = A.ptr() + 1 * A.ld();

        V* rd_0         = res.rep_ptr() + res.first_elem_diag(0);
        V* rd_p1        = res.rep_ptr() + res.first_elem_diag(1);
        V* rd_m1        = res.rep_ptr() + res.first_elem_diag(-1);

        Integer A_ld    = A.ld() + 1;
        Integer res_ld  = res.ld();

        inv_diag22_cont<V>::eval(M, d_0, d_p1, d_m1, rd_0, rd_p1, rd_m1, A_ld, res_ld, 
                                 logdet, sig_max, sig_min, diag_tol, modif);

        ret = Matrix(res,false);
    };
};

template<class V>
struct linsolve_obj_diag_22_inv_str<V,struct_banded>
{
    using Mat   = raw::Matrix<V,struct_banded>;

    static void eval(Matrix& ret, const Mat& A, Real& logdet, Real& sig_max, Real& sig_min, 
                     Real diag_tol, bool& modif)
    {
        Integer M   = A.rows();

        Mat res(A.get_type(), M, M, -1, 1);

        const V* d_0    = A.rep_ptr() + A.first_elem_diag(0);
        const V* d_p1   = A.rep_ptr() + A.first_elem_diag(1);
        const V* d_m1   = A.rep_ptr() + A.first_elem_diag(-1);

        V* rd_0         = res.rep_ptr() + res.first_elem_diag(0);
        V* rd_p1        = res.rep_ptr() + res.first_elem_diag(1);
        V* rd_m1        = res.rep_ptr() + res.first_elem_diag(-1);

        Integer A_ld    = A.ld();
        Integer res_ld  = res.ld();

        inv_diag22_cont<V>::eval(M, d_0, d_p1, d_m1, rd_0, rd_p1, rd_m1, A_ld, res_ld, logdet, 
                                 sig_max, sig_min, diag_tol, modif);

        ret = Matrix(res,false);
    };
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_base
//--------------------------------------------------------------------------------

linsolve_obj_base::linsolve_obj_base(Integer M, Integer N, value_code vc, const ti::ti_object& ti,
                                     bool modif)
    :m_M(M), m_N(N), m_vc(vc), m_ti(ti), m_modified(modif)
{};

linsolve_obj_base::~linsolve_obj_base()
{};

Integer linsolve_obj_base::rows() const
{
    return m_M;
};

Integer linsolve_obj_base::cols() const
{
    return m_N;
};

value_code linsolve_obj_base::get_value_code() const
{
    return m_vc;
};

bool linsolve_obj_base::is_modified() const
{
    return m_modified;
};

bool linsolve_obj_base::is_direct() const
{
    return true;
};

ti::ti_object linsolve_obj_base::get_type() const
{
    return m_ti;
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_empty
//--------------------------------------------------------------------------------

linsolve_obj_empty::linsolve_obj_empty(value_code vc, const ti::ti_object& ti)
    :linsolve_obj_base(0, 0, vc, ti, false)
{};

linsolve_obj_empty::~linsolve_obj_empty()
{};

linsolve_obj_empty::data_ptr linsolve_obj_empty::convert(value_code new_val_code) const
{
    return linsolve_obj_empty::data_ptr(new linsolve_obj_empty(new_val_code, m_ti));
};

matcl::Matrix linsolve_obj_empty::inv() const
{
    return matcl::zeros(0, 0, m_vc);
};

matcl::Matrix linsolve_obj_empty::base_matrix() const
{
    return matcl::zeros(0, 0, m_vc);
};

matcl::Matrix linsolve_obj_empty::solve(const Matrix& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve(0, 0, X.rows(), X.cols());
    return matcl::zeros(0, X.cols(), m_vc);
}

matcl::Matrix linsolve_obj_empty::solve(Matrix&& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve(0, 0, X.rows(), X.cols());
    return matcl::zeros(0, X.cols(), m_vc);
};

matcl::Matrix linsolve_obj_empty::solve_rev(const Matrix& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve_rev(0, 0, X.rows(), X.cols());
    return matcl::zeros(X.rows(), 0, m_vc);
}

matcl::Matrix linsolve_obj_empty::solve_rev(Matrix&& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve_rev(0, 0, X.rows(), X.cols());
    return matcl::zeros(X.rows(), 0, m_vc);
};

Matrix linsolve_obj_empty::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(0, 0, X.rows(), X.cols(), t, trans_type::no_trans);
    return details::convert_value(X, m_vc);
}

Matrix linsolve_obj_empty::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(0, 0, X.rows(), X.cols(), t, trans_type::no_trans);
    return details::convert_value(X, m_vc);
}

Matrix linsolve_obj_empty::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), 0, 0, trans_type::no_trans, t);
    return details::convert_value(X, m_vc);
}

Matrix linsolve_obj_empty::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), 0, 0, trans_type::no_trans, t);
    return details::convert_value(X, m_vc);
}

//--------------------------------------------------------------------------------
//                  linsolve_obj_nan
//--------------------------------------------------------------------------------
linsolve_obj_nan::linsolve_obj_nan(Integer N, value_code vc, const ti::ti_object& ti)
    :linsolve_obj_base(N, N, vc, ti, false)
{};

linsolve_obj_nan::~linsolve_obj_nan()
{};

linsolve_obj_nan::data_ptr linsolve_obj_nan::convert(value_code new_val_code) const
{
    return linsolve_obj_nan::data_ptr(new linsolve_obj_nan(m_N, new_val_code, m_ti));
};

matcl::Matrix linsolve_obj_nan::inv() const
{
    return md::make_nan_matrix(m_N, m_N, m_vc);
};

matcl::Matrix linsolve_obj_nan::base_matrix() const
{
    return md::make_nan_matrix(m_N, m_N, m_vc);
};

matcl::Matrix linsolve_obj_nan::solve(const Matrix& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve(m_N, m_N, X.rows(), X.cols());
    return md::make_nan_matrix(m_N, X.cols(), m_vc);
}

matcl::Matrix linsolve_obj_nan::solve(Matrix&& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve(m_N, m_N, X.rows(), X.cols());
    return md::make_nan_matrix(m_N, X.cols(), m_vc);
};

matcl::Matrix linsolve_obj_nan::solve_rev(const Matrix& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve_rev(m_N, m_N, X.rows(), X.cols());
    return md::make_nan_matrix(X.rows(), m_N, m_vc);
}

matcl::Matrix linsolve_obj_nan::solve_rev(Matrix&& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve_rev(m_N, m_N, X.rows(), X.cols());
    return md::make_nan_matrix(X.rows(), m_N, m_vc);
};

Matrix linsolve_obj_nan::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(m_N, m_N, X.rows(), X.cols(), t, trans_type::no_trans);
    return md::make_nan_matrix(m_N, X.cols(), m_vc);
}

Matrix linsolve_obj_nan::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(m_N, m_N, X.rows(), X.cols(), t, trans_type::no_trans);
    return md::make_nan_matrix(m_N, X.cols(), m_vc);
}

Matrix linsolve_obj_nan::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_N, m_N, trans_type::no_trans, t);
    return md::make_nan_matrix(X.rows(), m_N, m_vc);
}

Matrix linsolve_obj_nan::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_N, m_N, trans_type::no_trans, t);
    return md::make_nan_matrix(X.rows(), m_N, m_vc);
}

//--------------------------------------------------------------------------------
//                  linsolve_obj_id
//--------------------------------------------------------------------------------

linsolve_obj_id::linsolve_obj_id(Integer N, value_code vc, const ti::ti_object& ti)
    :linsolve_obj_base(N, N, vc, ti, false)
{};

linsolve_obj_id::~linsolve_obj_id()
{};

linsolve_obj_id::data_ptr linsolve_obj_id::convert(value_code new_val_code) const
{
    return linsolve_obj_id::data_ptr(new linsolve_obj_id(m_N, new_val_code, m_ti));
};

matcl::Matrix linsolve_obj_id::inv() const
{
    return matcl::speye(m_N, m_N, m_vc);
};

matcl::Matrix linsolve_obj_id::solve(const Matrix& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve(m_N, m_N, X.rows(), X.cols());
    return details::convert_value(X, m_vc);
}

matcl::Matrix linsolve_obj_id::solve(Matrix&& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve(m_N, m_N, X.rows(), X.cols());
    return details::convert_value(X, m_vc);
};

matcl::Matrix linsolve_obj_id::solve_rev(const Matrix& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve_rev(m_N, m_N, X.rows(), X.cols());
    return details::convert_value(X, m_vc);
}

matcl::Matrix linsolve_obj_id::solve_rev(Matrix&& X, trans_type tA) const
{
    (void)tA;
    error::check_lsolve_rev(m_N, m_N, X.rows(), X.cols());
    return details::convert_value(X, m_vc);
};

Matrix linsolve_obj_id::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return details::convert_value(X, m_vc);
}

Matrix linsolve_obj_id::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return details::convert_value(X, m_vc);
}

Matrix linsolve_obj_id::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return details::convert_value(X, m_vc);
}

Matrix linsolve_obj_id::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);
    return details::convert_value(X, m_vc);
}

matcl::Matrix linsolve_obj_id::base_matrix() const
{
    return beye(m_N,m_N,0,0,m_vc);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_scalar
//--------------------------------------------------------------------------------
linsolve_obj_scalar::linsolve_obj_scalar(const Matrix& A, const options& opts)
    :linsolve_obj_base(1,1,A.get_value_code(), A.get_type(), false), m_A(A)
{
    if (A.is_scalar() == false)
        throw error::scalar_required(A.rows(), A.cols());

    Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());

    if (A == 0 && tol >= 0)
        throw error::error_singular();

    if (tol < 0)
    {
        tol             = -tol;

        if (abs(A) < tol)
        {
            Matrix tol_m    = tol;
            tol_m           = details::convert_value(tol_m, A.get_value_code());
            m_A             = tol_m * (A == 0 ? 1.0f : sign(A));
            m_modified      = true;
        };
    };
};

linsolve_obj_scalar::linsolve_obj_scalar(const Matrix& A, bool modif, from_inv)
    :linsolve_obj_base(1,1,A.get_value_code(), A.get_type(), modif), m_A(A)
{};

linsolve_obj_scalar::~linsolve_obj_scalar()
{};

bool linsolve_obj_scalar::is_hermitian() const
{
    bool is_real = matrix_traits::is_float_real(get_value_code());
    return is_real? true : false;
};

bool linsolve_obj_scalar::is_posdef() const
{
    bool is_real = matrix_traits::is_float_real(get_value_code());
    return is_real? bool(m_A > 0) : false;
};

Real linsolve_obj_scalar::log_det() const
{
    return matcl::log(matcl::abs(m_A)).get_scalar<Real>();
};

Real linsolve_obj_scalar::normest_1() const
{
    return matcl::abs(m_A).get_scalar<Real>();
};

Real linsolve_obj_scalar::normest_inf() const
{
    return matcl::abs(m_A).get_scalar<Real>();
};
Real linsolve_obj_scalar::normest_2() const
{
    return matcl::abs(m_A).get_scalar<Real>();
};

Real linsolve_obj_scalar::mat_normest_1() const
{
    return 1.0 / normest_1();
};
Real linsolve_obj_scalar::mat_normest_2() const
{
    return 1.0 / normest_2();
};

Real linsolve_obj_scalar::mat_normest_inf() const
{
    return 1.0 / normest_inf();
};

linsolve_obj_scalar::data_ptr 
linsolve_obj_scalar::convert(value_code new_val_code) const
{
    Matrix Ai = details::convert_value(m_A,new_val_code);
    return linsolve_obj_scalar::data_ptr(new linsolve_obj_scalar(Ai, m_modified, from_inv()));
};

matcl::Matrix linsolve_obj_scalar::inv() const
{
    return 1.0f / m_A;
};

matcl::Matrix linsolve_obj_scalar::base_matrix() const
{
    return m_A;
};

matcl::Matrix linsolve_obj_scalar::solve(const Matrix& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return inv() * X;
        case trans_type::conj_trans:
            return conj(inv()) * X;
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};

matcl::Matrix linsolve_obj_scalar::solve(Matrix&& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return inv() * std::move(X);
        case trans_type::conj_trans:
            return conj(inv()) * std::move(X);
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}

matcl::Matrix linsolve_obj_scalar::solve_rev(const Matrix& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return X * inv();
        case trans_type::conj_trans:
            return X * conj(inv());
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};
matcl::Matrix linsolve_obj_scalar::solve_rev(Matrix&& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return std::move(X) * inv();
        case trans_type::conj_trans:
            return std::move(X) * conj(inv());
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};

Matrix linsolve_obj_scalar::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(1, 1, X.rows(), X.cols(), t, trans_type::no_trans);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return m_A * X;
        case trans_type::conj_trans:
            return conj(m_A) * X;
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};

Matrix linsolve_obj_scalar::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(1, 1, X.rows(), X.cols(), t, trans_type::no_trans);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return m_A * std::move(X);
        case trans_type::conj_trans:
            return conj(m_A) * std::move(X);
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}
Matrix linsolve_obj_scalar::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), 1, 1, trans_type::no_trans, t);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return X * m_A;
        case trans_type::conj_trans:
            return X * conj(m_A);
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};
Matrix linsolve_obj_scalar::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), 1, 1, trans_type::no_trans, t);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return std::move(X) * m_A;
        case trans_type::conj_trans:
            return std::move(X) * conj(m_A);
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_diag
//--------------------------------------------------------------------------------
linsolve_obj_diag::linsolve_obj_diag(const Matrix& A, const options& opts)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), false)
{
    initialize(A, opts);
};

linsolve_obj_diag::linsolve_obj_diag(const Matrix& Di, from_vec_inv)
    :linsolve_obj_base(Di.length(), Di.length(), Di.get_value_code(), Di.get_type(), false)
{
    initialize_vec_inv(Di);
};

linsolve_obj_diag::linsolve_obj_diag(const Matrix& D, const Matrix& Di, bool pd, bool modif, from_inv)
    :linsolve_obj_base(D.length(), D.length(), D.get_value_code(), D.get_type(), modif)
    , m_D(D), m_Di(Di), m_is_posdef(pd)
{};

void linsolve_obj_diag::initialize(const Matrix& A, const options& opts)
{
    m_D             = get_diag(A);
    Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());
    
    Matrix D2;
    bool sing, modif;
    tie(D2, sing, modif)    = make_nonsingular(A, tol);

    if (sing == true)
        throw error::error_singular();
    if (modif == true)
        m_modified  = true;

    m_Di            = div(1.f, D2);
    
    initialize_posdef(m_D);
};

void linsolve_obj_diag::initialize_vec_inv(const Matrix& Di)
{
    m_Di                = Di;

    bool is_nsing       = all_vec(m_Di);

    if (is_nsing == false)
        throw error::error_singular();

    m_D                 = div(1.f, m_Di);

    initialize_posdef(m_D);
};

void linsolve_obj_diag::initialize_posdef(const Matrix& A)
{
    bool is_real = matrix_traits::is_float_real(A.get_value_code());

    if (is_real == false)
    {
        m_is_posdef = false;
        return;
    };

    m_is_posdef = all_vec(A > 0.0f);
};

linsolve_obj_diag::~linsolve_obj_diag()
{};

bool linsolve_obj_diag::is_hermitian() const
{
    bool is_real = matrix_traits::is_float_real(get_value_code());
    return is_real? true : false;
};

bool linsolve_obj_diag::is_posdef() const
{
    return m_is_posdef;
};

Real linsolve_obj_diag::log_det() const
{
    return log_det_diag(m_D);
};
Real linsolve_obj_diag::normest_1() const
{
    return norm_vec_all(m_Di, basic_vector_norm::norm_inf);
};
Real linsolve_obj_diag::normest_2() const
{
    return norm_vec_all(m_Di, basic_vector_norm::norm_inf);
}

Real linsolve_obj_diag::normest_inf() const
{
    return norm_vec_all(m_Di, basic_vector_norm::norm_inf);
};

Real linsolve_obj_diag::mat_normest_1() const
{
    return norm_vec_all(m_D, basic_vector_norm::norm_inf);
};
Real linsolve_obj_diag::mat_normest_2() const
{
    return norm_vec_all(m_D, basic_vector_norm::norm_inf);
};
Real linsolve_obj_diag::mat_normest_inf() const
{
    return norm_vec_all(m_D, basic_vector_norm::norm_inf);
};

linsolve_obj_diag::data_ptr linsolve_obj_diag::convert(value_code new_val_code) const
{
    Matrix Di   = details::convert_value(m_Di,new_val_code);
    Matrix D    = details::convert_value(m_D,new_val_code);
    return linsolve_obj_diag::data_ptr(new linsolve_obj_diag(D, Di, m_is_posdef, m_modified, from_inv()));
};

matcl::Matrix linsolve_obj_diag::inv() const
{
    Matrix ret = bdiag(m_Di);

    if (is_hermitian() == true)
    {
        ret.add_struct(predefined_struct_type::her);

        if (is_posdef() == true)
            ret.add_struct(posdef_flag());
    };

    return ret;
};

matcl::Matrix linsolve_obj_diag::base_matrix() const
{
    Matrix ret = bdiag(m_D);

    if (is_hermitian() == true)
    {
        ret.add_struct(predefined_struct_type::her);

        if (is_posdef() == true)
            ret.add_struct(posdef_flag());
    };

    return ret;
};

matcl::Matrix linsolve_obj_diag::solve(const Matrix& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_rows(X,m_Di);    
        case trans_type::conj_trans:
            return scale_rows(X,conj(m_Di));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}
matcl::Matrix linsolve_obj_diag::solve(Matrix&& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_rows(std::move(X),m_Di);        
        case trans_type::conj_trans:
            return scale_rows(std::move(X), conj(m_Di));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}
matcl::Matrix linsolve_obj_diag::solve_rev(const Matrix& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_cols(X,m_Di);        
        case trans_type::conj_trans:
            return scale_cols(X,conj(m_Di));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}
matcl::Matrix linsolve_obj_diag::solve_rev(Matrix&& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_cols(std::move(X),m_Di);
        case trans_type::conj_trans:
            return scale_cols(std::move(X),conj(m_Di));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}

Matrix linsolve_obj_diag::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_rows(X, m_D);
        case trans_type::conj_trans:
            return scale_rows(X, conj(m_D));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}

Matrix linsolve_obj_diag::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_rows(std::move(X), m_D);
        case trans_type::conj_trans:
            return scale_rows(std::move(X), conj(m_D));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}

Matrix linsolve_obj_diag::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_cols(X, m_D);
        case trans_type::conj_trans:
            return scale_cols(X, conj(m_D));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}

Matrix linsolve_obj_diag::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);

    switch(t)
    {
        case trans_type::no_trans:
        case trans_type::trans:
            return scale_cols(std::move(X), m_D);
        case trans_type::conj_trans:
            return scale_cols(std::move(X), conj(m_D));
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}

//--------------------------------------------------------------------------------
//                  linsolve_obj_unitary
//--------------------------------------------------------------------------------

linsolve_obj_unitary::linsolve_obj_unitary(const Matrix& A)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), false)
    , m_A(A)
{}
linsolve_obj_unitary::linsolve_obj_unitary(const Matrix& A, from_inv)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), false)
    , m_A(A)
{};

linsolve_obj_unitary::~linsolve_obj_unitary()
{};

Real linsolve_obj_unitary::normest_1() const
{
    return mat_normest_inf();
};
Real linsolve_obj_unitary::normest_inf() const
{
    return mat_normest_1();
};
Real linsolve_obj_unitary::mat_normest_1() const
{
    return norm(m_A, 1.0);
};
Real linsolve_obj_unitary::mat_normest_inf() const
{
    return norm(m_A, constants::inf());
};

matcl::Matrix linsolve_obj_unitary::inv() const
{
    return ctrans(m_A);
};

matcl::Matrix linsolve_obj_unitary::base_matrix() const
{
    return m_A;
};

linsolve_obj_unitary::data_ptr linsolve_obj_unitary::convert(value_code new_val_code) const
{
    Matrix Ai = details::convert_value(m_A,new_val_code);
    return linsolve_obj_diag::data_ptr(new linsolve_obj_unitary(Ai, from_inv()));
};

Matrix linsolve_obj_unitary::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(m_A.rows(), m_A.cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return mmul(m_A, X, t);
};
Matrix linsolve_obj_unitary::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(m_A.rows(), m_A.cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return mmul(m_A, std::move(X), t);
};

Matrix linsolve_obj_unitary::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_A.rows(), m_A.cols(), trans_type::no_trans, t);
    return mmul(X, m_A, trans_type::no_trans, t);
}
Matrix linsolve_obj_unitary::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_A.rows(), m_A.cols(), trans_type::no_trans, t);
    return mmul(std::move(X), m_A, trans_type::no_trans, t);
}

matcl::Matrix linsolve_obj_unitary::solve(const Matrix& X, trans_type t) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    trans_type_ext te   = details::trans_manip::link_trans(
                            details::trans_manip::convert_trans(t), 
                            trans_type_ext::conj_trans);

    return mmul(m_A, X, te, trans_type_ext::no_trans);
};

matcl::Matrix linsolve_obj_unitary::solve(Matrix&& X, trans_type t) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());
    trans_type_ext te   = details::trans_manip::link_trans(
                            details::trans_manip::convert_trans(t), 
                            trans_type_ext::conj_trans);

    return mmul(m_A, std::move(X), te, trans_type_ext::no_trans);
};

matcl::Matrix linsolve_obj_unitary::solve_rev(const Matrix& X, trans_type t) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    trans_type_ext te   = details::trans_manip::link_trans(
                            details::trans_manip::convert_trans(t), 
                            trans_type_ext::conj_trans);

    return mmul(X, m_A, trans_type_ext::no_trans, te);
};
matcl::Matrix linsolve_obj_unitary::solve_rev(Matrix&& X, trans_type t) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    trans_type_ext te   = details::trans_manip::link_trans(
                            details::trans_manip::convert_trans(t), 
                            trans_type_ext::conj_trans);

    return mmul(std::move(X), m_A, trans_type_ext::no_trans, te);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_unitary_mat
//--------------------------------------------------------------------------------

linsolve_obj_unitary_mat::linsolve_obj_unitary_mat(const unitary_matrix& A)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), false)
    ,m_A(A)
{}

linsolve_obj_unitary_mat::~linsolve_obj_unitary_mat()
{};

matcl::Matrix linsolve_obj_unitary_mat::inv() const
{
    return ctrans(m_A).to_matrix();
};

matcl::Matrix linsolve_obj_unitary_mat::base_matrix() const
{
    return m_A.to_matrix();
};

linsolve_obj_unitary_mat::data_ptr linsolve_obj_unitary_mat::convert(value_code new_val_code) const
{
    unitary_matrix Ai = m_A.convert(new_val_code);
    return linsolve_obj_diag::data_ptr(new linsolve_obj_unitary_mat(Ai));
};

matcl::Matrix linsolve_obj_unitary_mat::solve(const Matrix& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
            return ctrans(m_A) * X;
        case trans_type::trans:
            return ctrans(trans(m_A)) * X;
        case trans_type::conj_trans:
            return m_A * X;
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};

matcl::Matrix linsolve_obj_unitary_mat::solve(Matrix&& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
            return ctrans(m_A) * std::move(X);
        case trans_type::trans:
            return ctrans(trans(m_A)) * std::move(X);
        case trans_type::conj_trans:
            return m_A * std::move(X);
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
}

matcl::Matrix linsolve_obj_unitary_mat::solve_rev(const Matrix& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
            return X * ctrans(m_A);
        case trans_type::trans:
            return X * ctrans(trans(m_A));
        case trans_type::conj_trans:
            return X * m_A;
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};
matcl::Matrix linsolve_obj_unitary_mat::solve_rev(Matrix&& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    switch(tA)
    {
        case trans_type::no_trans:
            return std::move(X) * ctrans(m_A);
        case trans_type::trans:
            return X * ctrans(trans(m_A));
        case trans_type::conj_trans:
            return std::move(X) * m_A;
        default:
        {
            matcl_assert(0,"unknown case");
            throw error::error_general("invalid case");
        }
    };
};

Matrix linsolve_obj_unitary_mat::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(m_A.rows(), m_A.cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return trans(m_A, t) * X;
};
Matrix linsolve_obj_unitary_mat::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(m_A.rows(), m_A.cols(), X.rows(), X.cols(), t, trans_type::no_trans);
    return trans(m_A, t) * std::move(X);
};

Matrix linsolve_obj_unitary_mat::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_A.rows(), m_A.cols(), trans_type::no_trans, t);
    return X * trans(m_A, t);
}
Matrix linsolve_obj_unitary_mat::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_A.rows(), m_A.cols(), trans_type::no_trans, t);
    return std::move(X) * trans(m_A, t);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_triang
//--------------------------------------------------------------------------------

linsolve_obj_triang::linsolve_obj_triang(const Matrix& A, const options& opts)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), false)
    , m_LU(A), m_p(permvec::identity(A.rows())), m_q(permvec::identity(A.cols()))
{
    if (matcl::get_ld(A,0) == 0)
    {
        m_tril  = false;
        m_LU.add_struct(predefined_struct_type::triu);
    }
    else
    {
        m_tril  = true;
        m_LU.add_struct(predefined_struct_type::tril);
    };

    correct_singular(m_LU, opts);
    test_singular();
};
linsolve_obj_triang::linsolve_obj_triang(const Matrix& A, const permvec& p, 
                                        const permvec& q, const options& opts)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), false)
    , m_LU(A), m_p(p), m_q(q)
{
    if (matcl::get_ld(A,0) == 0)
    {
        m_LU.add_struct(predefined_struct_type::triu);
        m_tril  = false;
    }
    else
    {
        m_LU.add_struct(predefined_struct_type::tril);
        m_tril  = true;
    };

    correct_singular(m_LU, opts);
    test_singular();
};
linsolve_obj_triang::linsolve_obj_triang(const Matrix& A, const permvec& p, 
                                        const permvec& q, bool modif, from_inv)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), modif)
    , m_LU(A), m_p(p), m_q(q)
{
    if (matcl::get_ld(A,0) == 0)
    {
        m_LU.add_struct(predefined_struct_type::triu);
        m_tril  = false;
    }
    else
    {
        m_LU.add_struct(predefined_struct_type::tril);
        m_tril  = true;
    };

    test_singular();
};

void linsolve_obj_triang::correct_singular(Matrix& A, const options& opts)
{
    Real tol        = opts.get_option<Real>(opt::linsolve::tol_sing());
    
    bool sing, modif;
    tie(A, sing, modif) = make_nonsingular(A, tol);

    if (sing == true)
        throw error::error_singular();
    if (modif == true)
        m_modified  = true;
};

linsolve_obj_triang::~linsolve_obj_triang()
{};

void linsolve_obj_triang::test_singular() const
{
    matcl::Matrix Ad    = get_diag(m_LU);
    bool is_nsing       = all_vec(Ad);

    if (is_nsing == false)
        throw error::error_singular();
};

linsolve_obj_triang::data_ptr linsolve_obj_triang::convert(value_code new_val_code) const
{
    Matrix Uc = details::convert_value(m_LU,new_val_code);
    return linsolve_obj_triang::data_ptr(new linsolve_obj_triang(Uc,m_p,m_q, m_modified, from_inv()));
};

template<class Fac>
struct triang_inv_vis : public extract_type_switch<void, triang_inv_vis<Fac>,true>
{
    template<class Mat>
    static void eval(const Matrix&, const Mat&, Matrix& ret, const Fac& f)
    {
        using V     = typename Mat::value_type;
        using S     = typename Mat::struct_type;
        using VR    = typename md::unify_types<V,Float>::type;
        using Mat_R = raw::Matrix<VR,S>;

        return f.inv_impl<Mat_R>(ret);
    };

    template<class V>
    static void eval_scalar(const Matrix&, const V& LU, Matrix& ret, const Fac& f)
    {
        using VR    = typename md::unify_types<V,Float>::type;
        return f.inv_scal_impl<VR>(ret, LU);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, const Fac&)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::inv");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, const Fac&)
    {
        throw error::object_value_type_not_allowed("linsolve_obj::inv");
    };

};

template<class V, class S>
struct eval_inv_triang
{
    template<class Fac>
    static void eval(Matrix& ret, const Matrix& LU, const permvec&, const permvec&, 
                     bool, const Fac& f)
    {
        Integer N   = LU.rows();
        Matrix X    = speye(N, N, LU.get_value_code());

        ret = f.solve(std::move(X), trans_type::no_trans);
        return;
    };
};

template<class V>
struct eval_inv_triang<V,struct_dense>
{
    using Mat   = raw::Matrix<V,struct_dense>;

    template<class Fac>
    static void eval(Matrix& ret, const Matrix& LU, const permvec& p, const permvec& q, 
                     bool tril, const Fac&)
    {
        const char* UPLO    = tril? "L" : "U";
        const char* DIAG    = "N";
        Integer N           = LU.rows();

        Matrix LU_s(LU);
        Mat A               = LU_s.impl_unique<Mat>();

        Integer info;
        lapack::trtri(UPLO, DIAG, N, lap(A.ptr()), A.ld(), info);

        if (info > 0)
            throw error::error_singular();

        Matrix Am           = Matrix(std::move(A),false);

        if (p.is_id() == false && q.is_id() == false)
            Am              = std::move(Am)(p,q);

        ret = Am;
    };
};

template<class Mat>
void linsolve_obj_triang::inv_impl(Matrix& ret) const
{
    using V     = typename Mat::value_type;
    using S     = typename Mat::struct_type;

    return eval_inv_triang<V,S>::eval(ret, m_LU, m_p, m_q, m_tril,*this);
};

template<class V>
void linsolve_obj_triang::inv_scal_impl(Matrix& ret, const V& scal) const
{
    if (scal == V(0.0))
        throw error::error_singular();

    ret = V(1.0) / scal;
};

matcl::Matrix linsolve_obj_triang::inv() const
{
    using Fac = linsolve_obj_triang;

    Matrix ret;
    triang_inv_vis<Fac>::make<const Matrix&>(m_LU, ret, *this);
    return ret;
};

matcl::Matrix linsolve_obj_triang::base_matrix() const
{
    if (m_p.is_id())
    {
        if (m_q.is_id())
            return m_LU;
        else
            return m_LU(colon(), m_q.invperm());
    }
    else
    {
        if (m_q.is_id())
            return m_LU(m_p.invperm(), colon());
        else
            return m_LU(m_p.invperm(), m_q.invperm());
    }
};

matcl::Matrix linsolve_obj_triang::solve(const Matrix& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());
    return linsolve(m_LU, m_p, m_q, X, tA);
}
matcl::Matrix linsolve_obj_triang::solve(Matrix&& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());
    return linsolve(m_LU, m_p, m_q, std::move(X), tA);
}
matcl::Matrix linsolve_obj_triang::solve_rev(const Matrix& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());
    return linsolve_rev2(m_LU, m_p, m_q, X, tA);
}
matcl::Matrix linsolve_obj_triang::solve_rev(Matrix&& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());
    return linsolve_rev2(m_LU, m_p, m_q, std::move(X), tA);
}

Real linsolve_obj_triang::log_det() const
{
    matcl::Matrix Ad    = get_diag(m_LU);
    return log_det_diag(Ad);
};

Matrix linsolve_obj_triang::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(m_LU.rows(), m_LU.cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    permvec p = m_p;
    permvec q = m_q;

    if (t != trans_type::no_trans)
        std::swap(p, q);

    //LU(p^-1,q^-1) * X = (LU * X(q,:))(p^-1,:)
    //op(LU(p^-1,q^-1)) * X = (op(LU) * X(p,:))(q^-1,:)
    Matrix LUX;
    if (q.is_id() == true)
        LUX = mmul(m_LU, X, t);
    else
        LUX = mmul(m_LU, X(q,colon()), t);

    if (p.is_id() == true)
        return LUX;
    else
        return std::move(LUX)(p.invperm(), colon());
};

Matrix linsolve_obj_triang::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(m_LU.rows(), m_LU.cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    permvec p = m_p;
    permvec q = m_q;

    if (t != trans_type::no_trans)
        std::swap(p, q);

    //LU(p^-1,q^-1) * X = (LU * X(q,:))(p^-1,:)
    //op(LU(p^-1,q^-1)) * X = (op(LU) * X(p,:))(q^-1,:)
    Matrix LUX;
    if (q.is_id() == true)
        LUX = mmul(m_LU, std::move(X), t);
    else
        LUX = mmul(m_LU, std::move(X)(q,colon()), t);

    if (p.is_id() == true)
        return LUX;
    else
        return std::move(LUX)(p.invperm(), colon());
};

Matrix linsolve_obj_triang::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_LU.rows(), m_LU.cols(), trans_type::no_trans, t);

    permvec p = m_p;
    permvec q = m_q;

    if (t != trans_type::no_trans)
        std::swap(p, q);

    //X * LU(p^-1,q^-1) = (X(:,p) * LU)(:,q^-1)
    //X * op(LU(p^-1,q^-1)) = (X(:,q) * op(LU))(:,p^-1)

    Matrix LUX;
    if (p.is_id() == true)
        LUX = mmul(X, m_LU, trans_type::no_trans, t);
    else
        LUX = mmul(m_LU, X(colon(), p), trans_type::no_trans, t);

    if (q.is_id() == true)
        return LUX;
    else
        return std::move(LUX)(colon(), q.invperm());
};

Matrix linsolve_obj_triang::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_LU.rows(), m_LU.cols(), trans_type::no_trans, t);

    permvec p = m_p;
    permvec q = m_q;

    if (t != trans_type::no_trans)
        std::swap(p, q);

    //X * LU(p^-1,q^-1) = (X(:,p) * LU)(:,q^-1)
    //X * op(LU(p^-1,q^-1)) = (X(:,q) * op(LU))(:,p^-1)

    Matrix LUX;
    if (p.is_id() == true)
        LUX = mmul(std::move(X), m_LU, trans_type::no_trans, t);
    else
        LUX = mmul(m_LU, std::move(X)(colon(), p), trans_type::no_trans, t);

    if (q.is_id() == true)
        return LUX;
    else
        return std::move(LUX)(colon(), q.invperm());
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_diag_22
//--------------------------------------------------------------------------------
linsolve_obj_diag_22::linsolve_obj_diag_22(const Matrix& A, const options& opts)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), false)
    , m_A(A)
{
    initialize(A, m_logdet, m_sig_max, m_sig_min, opts);

    bool is_real    = matrix_traits::is_float_real(get_value_code());
    m_hermitian     = A.get_struct().is_hermitian(A.is_square(), is_real);
    m_posdef        = (m_hermitian == false) ? false : matcl::is_posdef(A.get_struct());

    if (is_hermitian() == true)
    {
        m_Ai.add_struct(predefined_struct_type::her);

        if (is_posdef() == true)
            m_Ai.add_struct(posdef_flag());
    };
};
linsolve_obj_diag_22::linsolve_obj_diag_22(const Matrix& A, const Matrix& Ai, bool herm, bool pd,
                                           Real logdet, Real sig_max, Real sig_min, bool modif, from_inv)
    :linsolve_obj_base(A.rows(), A.cols(), A.get_value_code(), A.get_type(), modif)
    , m_A(A), m_Ai(Ai), m_hermitian(herm), m_logdet(logdet), m_posdef(pd)
    , m_sig_max(sig_max), m_sig_min(sig_min)
{};

template<class V, class S>
struct linsolve_obj_diag_22_inv_val
{
    using Mat   = raw::Matrix<V,S>;

    static void eval(Matrix& ret, const Mat& A, Real& logdet, Real& sig_max, Real& sig_min, 
                     Real diag_tol, bool& modif)
    {
        return linsolve_obj_diag_22_inv_str<V,S>::eval(ret, A, logdet, sig_max, sig_min, diag_tol, modif);
    }
};

template<class S>
struct linsolve_obj_diag_22_inv_val<Integer,S>
{
    using Mat       = raw::Matrix<Integer,S>;
    using Mat_B     = raw::Matrix<Real,struct_banded>;

    static void eval(Matrix& ret, const Mat& A, Real& logdet, Real& sig_max, Real& sig_min,
                     Real diag_tol, bool& modif)
    {
        Mat_B Ac    = raw::converter<Mat_B,Mat>::eval(A);
        return linsolve_obj_diag_22_inv_val<Real, struct_banded>
                    ::eval(ret, Ac, logdet, sig_max, sig_min, diag_tol, modif);
    };
};

struct linsolve_obj_diag_22_inv : public extract_type_switch<void, linsolve_obj_diag_22_inv,true>
{
    template<class T>
    static void eval(const Matrix&, const T& A, Matrix& ret, Real& logdet, Real& sig_max, Real& sig_min,
                     Real diag_tol, bool& modif)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return linsolve_obj_diag_22_inv_val<V,S>::eval(ret, A, logdet, sig_max, sig_min, diag_tol, modif);
    }

    template<class T>
    static void eval_scalar(const Matrix&, const T& A0, Matrix& ret, Real& logdet, Real& sig_max, Real& sig_min,
                            Real tol, bool& modif)
    {
        if (A0 == T(0.0) && tol == 0.0)
            throw error::error_singular();

        using VR    = typename md::real_type<T>::type;
        modif       = false;

        T A;
        if (abs(A0) < tol)
        { 
            A       = (VR)tol * (A0 == T(0.0) ? VR(1.0) : sign(A0));
            modif   = true;
        }
        else
        {
            A       = A0;
        }
        ret         = T(1.0) / A;
        logdet      = matcl::log(matcl::abs(A));
        sig_max     = abs(A);
        sig_min     = abs(A);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, Real&, Real&, Real&, Real, bool&)
    {
        throw error::object_value_type_not_allowed("linsolve_diag_22");
    };

    static void eval_scalar(const Matrix&, const Object&, Matrix&, Real&, Real&, Real&, Real, bool&)
    {
        throw error::object_value_type_not_allowed("linsolve_diag_22");
    };
};

void linsolve_obj_diag_22::initialize(const Matrix& A, Real& logdeg, Real& sig_max, Real& sig_min,
                                      const options& opts)
{
    Real tol            = opts.get_option<Real>(opt::linsolve::tol_sing());
    
    // calculate small pivot tolerance
    if (tol > 0.0)
    {
        Real max_val    = norm_vec_all(A, basic_vector_norm::norm_inf);

        if (max_val == 0.0)
            throw error::error_singular();

        tol             = pow(constants::eps(A.get_value_code()), tol) * max_val;
    }
    else if (tol < 0)
    {
        tol             = -tol;
    };

    Matrix Ai;
    bool modif;
    linsolve_obj_diag_22_inv::make<const Matrix&>(A, Ai, logdeg, sig_max, sig_min, tol, modif);

    if (modif == true)
        m_modified  = true;

    m_Ai    = Ai;
};

linsolve_obj_diag_22::~linsolve_obj_diag_22()
{};

bool linsolve_obj_diag_22::is_hermitian() const
{
    return m_hermitian;
};
bool linsolve_obj_diag_22::is_posdef() const
{
    return m_posdef;
};

linsolve_obj_diag_22::data_ptr linsolve_obj_diag_22::convert(value_code new_val_code) const
{
    Matrix Ai = details::convert_value(m_Ai,new_val_code);
    Matrix A  = details::convert_value(m_A,new_val_code);
    return linsolve_obj_diag_22::data_ptr(new linsolve_obj_diag_22(A, Ai, m_hermitian, m_posdef, m_logdet, 
                                                m_sig_max, m_sig_min, m_modified, from_inv()));
};

Real linsolve_obj_diag_22::log_det() const
{
    return m_logdet;
};

matcl::Matrix linsolve_obj_diag_22::inv() const
{
    return m_Ai;
};

matcl::Matrix linsolve_obj_diag_22::base_matrix() const
{
    return m_A;
};

Real linsolve_obj_diag_22::normest_2() const
{
    return 1.0/m_sig_min;
}
Real linsolve_obj_diag_22::mat_normest_2() const
{
    return m_sig_max;
};

matcl::Matrix linsolve_obj_diag_22::solve(const Matrix& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    return mmul(m_Ai, X, tA);
}
matcl::Matrix linsolve_obj_diag_22::solve(Matrix&& X, trans_type tA) const
{
    error::check_lsolve(this->rows(), this->cols(), X.rows(), X.cols());

    return mmul(m_Ai, std::move(X), tA);
}
matcl::Matrix linsolve_obj_diag_22::solve_rev(const Matrix& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());

    return mmul(X, m_Ai, trans_type::no_trans, tA);
}
matcl::Matrix linsolve_obj_diag_22::solve_rev(Matrix&& X, trans_type tA) const
{
    error::check_lsolve_rev(this->rows(), this->cols(), X.rows(), X.cols());
    return mmul(std::move(X), m_Ai, trans_type::no_trans, tA);
}

Matrix linsolve_obj_diag_22::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(m_A.rows(), m_A.cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    return mmul(m_A, X, t);
};
Matrix linsolve_obj_diag_22::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(m_A.rows(), m_A.cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    return mmul(m_A, std::move(X), t);
};

Matrix linsolve_obj_diag_22::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_A.rows(), m_A.cols(), trans_type::no_trans, t);
    return mmul(X, m_A, trans_type::no_trans, t);
};
Matrix linsolve_obj_diag_22::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), m_A.rows(), m_A.cols(), trans_type::no_trans, t);
    return mmul(std::move(X), m_A, trans_type::no_trans, t);
};

//--------------------------------------------------------------------------------
//                  linsolve_obj_seq_2
//--------------------------------------------------------------------------------
linsolve_obj_seq_2::linsolve_obj_seq_2()
    :linsolve_obj_base(0, 0, value_code::v_real, ti::ti_object_type<Real>(), false)
    , m_has_base(false)
{};

linsolve_obj_seq_2::linsolve_obj_seq_2(const linsolve_obj& A1, const linsolve_obj& A2)
    :linsolve_obj_base(A1.rows(), A2.cols(), A1.get_value_code(), A1.get_type(), false)
    , m_A1(A1), m_A2(A2), m_has_base(false)
{
    m_vc        = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, value_code::v_float);  
    m_modified  = A1.is_modified() || A2.is_modified();

    m_A1    = m_A1.convert(m_vc);
    m_A2    = m_A2.convert(m_vc);

    m_ti    = m_A1.get_type();
};
linsolve_obj_seq_2::linsolve_obj_seq_2(const Matrix& base, const linsolve_obj& A1, const linsolve_obj& A2)
    :linsolve_obj_base(A1.rows(), A2.cols(), base.get_value_code(), base.get_type(), false)
    , m_A1(A1), m_A2(A2), m_has_base(true), m_base_mat(base)
{
    m_modified  = A1.is_modified() || A2.is_modified();
    m_vc        = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, base.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, value_code::v_float);  

    m_base_mat  = details::convert_value(m_base_mat, m_vc);
    m_A1        = m_A1.convert(m_vc);
    m_A2        = m_A2.convert(m_vc);

    m_ti    = m_A1.get_type();
};

bool linsolve_obj_seq_2::is_direct() const
{
    return m_A1.is_direct() && m_A2.is_direct();
};

void linsolve_obj_seq_2::reset(const linsolve_obj& A1, const linsolve_obj& A2)
{
    m_vc        = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, value_code::v_float);  

    m_A1        = A1.convert(m_vc);
    m_A2        = A2.convert(m_vc);

    m_ti        = m_A1.get_type();
    m_M         = m_A1.rows();
    m_N         = m_A2.cols();
    m_has_base  = false;
    m_modified  = A1.is_modified() || A2.is_modified();
};
void linsolve_obj_seq_2::reset(const Matrix& base, const linsolve_obj& A1, const linsolve_obj& A2)
{
    m_vc        = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, base.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, value_code::v_float);  

    m_base_mat  = details::convert_value(base, m_vc);
    m_A1        = A1.convert(m_vc);
    m_A2        = A2.convert(m_vc);

    m_ti        = m_A1.get_type();
    m_M         = m_A1.rows();
    m_N         = m_A2.cols();
    m_has_base  = true;
    m_modified  = A1.is_modified() || A2.is_modified();
};

linsolve_obj_seq_2::~linsolve_obj_seq_2()
{};

Real linsolve_obj_seq_2::log_det() const
{
    return m_A1.log_det() + m_A2.log_det();
}

bool linsolve_obj_seq_2::is_hermitian() const
{
    if (m_has_base == false)
        return false;

    bool is_real    = matrix_traits::is_float_real(m_base_mat.get_value_code());
    return m_base_mat.get_struct().is_hermitian(m_base_mat.is_square(), is_real);
};

bool linsolve_obj_seq_2::is_posdef() const
{
    if (m_has_base == false)
        return false;

    return matcl::is_posdef(m_base_mat.get_struct());
};

linsolve_obj_seq_2::data_ptr 
linsolve_obj_seq_2::convert(value_code vc) const
{
    if (m_has_base == false)
    {
        return data_ptr(new linsolve_obj_seq_2(m_A1.convert(vc), m_A2.convert(vc)));
    }
    else
    {
        Matrix Ac = details::convert_value(m_base_mat, vc);
        return data_ptr(new linsolve_obj_seq_2(Ac, m_A1.convert(vc), m_A2.convert(vc)));
    }
};

matcl::Matrix linsolve_obj_seq_2::solve(const Matrix& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //(A1*A2)Y = X => Y = inv(A2)*inv(A1)*X
        Matrix Y    = m_A1.solve(X, trans_type::no_trans);
        Y           = m_A2.solve(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // op(A1*A2)Y = X => oop(A2)*op(A1)Y = X = > 
        // Y = op(op(A1))*inv(op(A2))*X
        Matrix Y    = m_A2.solve(X, tA);
        Y           = m_A1.solve(std::move(Y), tA);
        return Y;
    };
};
matcl::Matrix linsolve_obj_seq_2::solve(Matrix&& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //(A1*A2)Y = X => Y = inv(A2)*inv(A1)*X
        Matrix Y    = m_A1.solve(std::move(X), trans_type::no_trans);
        Y           = m_A2.solve(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // op(A1*A2)Y = X => oop(A2)*op(A1)Y = X = > 
        // Y = op(op(A1))*inv(op(A2))*X
        Matrix Y    = m_A2.solve(std::move(X), tA);
        Y           = m_A1.solve(std::move(Y), tA);
        return Y;
    };
};

/// solve Y * op(A) = X, where A is the matrix represented by this object
matcl::Matrix linsolve_obj_seq_2::solve_rev(const Matrix& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //Y*(A1*A2) = X => Y = X*inv(A2)*inv(A1)
        Matrix Y    = m_A2.solve_rev(X, trans_type::no_trans);
        Y           = m_A1.solve_rev(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // Y*op(A1*A2) = X => Y*op(A2)*op(A1) = X = > 
        // Y = X*op(op(A1))*inv(op(A2))
        Matrix Y    = m_A1.solve_rev(X, tA);
        Y           = m_A2.solve_rev(std::move(Y), tA);
        return Y;
    };
};
matcl::Matrix linsolve_obj_seq_2::solve_rev(Matrix&& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //Y*(A1*A2) = X => Y = X*inv(A2)*inv(A1)
        Matrix Y    = m_A2.solve_rev(std::move(X), trans_type::no_trans);
        Y           = m_A1.solve_rev(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // Y*op(A1*A2) = X => Y*op(A2)*op(A1) = X = > 
        // Y = X*op(op(A1))*inv(op(A2))
        Matrix Y    = m_A1.solve_rev(std::move(X), tA);
        Y           = m_A2.solve_rev(std::move(Y), tA);
        return Y;
    };
};

matcl::Matrix linsolve_obj_seq_2::inv() const
{
    Matrix X    = m_A1.inv();
    X           = m_A2.solve(std::move(X), trans_type::no_trans);

    if (is_hermitian() == true)
    {
        X.add_struct(predefined_struct_type::her);

        if (is_posdef() == true)
            X.add_struct(posdef_flag());
    };

    return X;
};

matcl::Matrix linsolve_obj_seq_2::base_matrix() const
{
    if (m_has_base)
        return m_base_mat;
    else
        return m_A1.base_matrix() * m_A2.base_matrix();
};

Matrix linsolve_obj_seq_2::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    if (m_has_base)
        return mmul(m_base_mat, X, t);

    if (t == trans_type::no_trans)
    {
        //(A1*A2)Y
        Matrix Y    = m_A2.mmul_right(X, trans_type::no_trans);
        Y           = m_A1.mmul_right(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //op(A1*A2) * Y = op(A2)*op(A1) * Y
        Matrix Y    = m_A1.mmul_right(X, t);
        Y           = m_A2.mmul_right(std::move(Y), t);
        return Y;
    };
}

Matrix linsolve_obj_seq_2::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    if (m_has_base)
        return mmul(m_base_mat, std::move(X), t);

    if (t == trans_type::no_trans)
    {
        //(A1*A2)Y
        Matrix Y    = m_A2.mmul_right(std::move(X), trans_type::no_trans);
        Y           = m_A1.mmul_right(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //op(A1*A2) * Y = op(A2)*op(A1) * Y
        Matrix Y    = m_A1.mmul_right(std::move(X), t);
        Y           = m_A2.mmul_right(std::move(Y), t);
        return Y;
    };
}

Matrix linsolve_obj_seq_2::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);

    if (m_has_base)
        return mmul(X, m_base_mat, trans_type::no_trans, t);

    if (t == trans_type::no_trans)
    {
        //Y*(A1*A2)
        Matrix Y    = m_A1.mmul_left(X, trans_type::no_trans);
        Y           = m_A2.mmul_left(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //Y * op(A1*A2) = Y * op(A2)*op(A1)
        Matrix Y    = m_A2.mmul_left(X, t);
        Y           = m_A1.mmul_left(std::move(Y), t);
        return Y;
    };
}

Matrix linsolve_obj_seq_2::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);

    if (m_has_base)
        return mmul(std::move(X), m_base_mat, trans_type::no_trans, t);

    if (t == trans_type::no_trans)
    {
        //Y*(A1*A2)
        Matrix Y    = m_A1.mmul_left(std::move(X), trans_type::no_trans);
        Y           = m_A2.mmul_left(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //Y * op(A1*A2) = Y * op(A2)*op(A1)
        Matrix Y    = m_A2.mmul_left(std::move(X), t);
        Y           = m_A1.mmul_left(std::move(Y), t);
        return Y;
    };
}

Real linsolve_obj_seq_2::mat_normest_1() const
{
    if (m_has_base)
        return norm(m_base_mat, 1.0);
    else
        return linsolve_obj_base::mat_normest_1();
}
Real linsolve_obj_seq_2::mat_normest_inf() const
{
    if (m_has_base)
        return norm(m_base_mat, constants::inf());
    else
        return linsolve_obj_base::mat_normest_inf();
}

//--------------------------------------------------------------------------------
//                  linsolve_obj_seq_3
//--------------------------------------------------------------------------------
linsolve_obj_seq_3::linsolve_obj_seq_3()
    :linsolve_obj_base(0, 0, value_code::v_real, ti::ti_object_type<Real>(), false)
    ,m_has_base(false)
{};

linsolve_obj_seq_3::linsolve_obj_seq_3(const linsolve_obj& A1, const linsolve_obj& A2, 
                                             const linsolve_obj& A3)
    :linsolve_obj_base(A1.rows(), A3.cols(), A1.get_value_code(), A1.get_type(), false)
    , m_A1(A1), m_A2(A2), m_A3(A3),m_has_base(false)
{
    m_modified  = A1.is_modified() || A2.is_modified() || A3.is_modified();

    m_vc    = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc    = matrix_traits::unify_value_types(m_vc, A3.get_value_code());  
    m_vc    = matrix_traits::unify_value_types(m_vc, value_code::v_float);  

    m_A1    = m_A1.convert(m_vc);
    m_A2    = m_A2.convert(m_vc);
    m_A3    = m_A3.convert(m_vc);

    m_ti    = m_A1.get_type();
};
linsolve_obj_seq_3::linsolve_obj_seq_3(const Matrix& base, const linsolve_obj& A1, const linsolve_obj& A2, 
                                             const linsolve_obj& A3)
    :linsolve_obj_base(A1.rows(), A3.cols(), A1.get_value_code(), A1.get_type(), false)
    , m_A1(A1), m_A2(A2), m_A3(A3),m_has_base(true), m_base_mat(base)
{
    m_modified  = A1.is_modified() || A2.is_modified() || A3.is_modified();

    m_vc    = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc    = matrix_traits::unify_value_types(m_vc, A3.get_value_code());  
    m_vc    = matrix_traits::unify_value_types(m_vc, m_base_mat.get_value_code());  
    m_vc    = matrix_traits::unify_value_types(m_vc, value_code::v_float);  

    m_A1        = m_A1.convert(m_vc);
    m_A2        = m_A2.convert(m_vc);
    m_A3        = m_A3.convert(m_vc);
    m_base_mat  = details::convert_value(m_base_mat, m_vc);

    m_ti    = m_A1.get_type();
};

void linsolve_obj_seq_3::reset(const linsolve_obj& A1, const linsolve_obj& A2, const linsolve_obj& A3)
{
    m_vc        = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, A3.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, value_code::v_float);  

    m_A1        = A1.convert(m_vc);
    m_A2        = A2.convert(m_vc);
    m_A3        = A3.convert(m_vc);

    m_ti        = m_A1.get_type();
    m_M         = m_A1.rows();
    m_N         = m_A3.cols();
    m_has_base  = false;
    m_modified  = A1.is_modified() || A2.is_modified() || A3.is_modified();
};

void linsolve_obj_seq_3::reset(const Matrix& base, const linsolve_obj& A1, const linsolve_obj& A2, 
                               const linsolve_obj& A3)
{
    m_vc        = matrix_traits::unify_value_types(A1.get_value_code(), A2.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, A3.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, base.get_value_code());  
    m_vc        = matrix_traits::unify_value_types(m_vc, value_code::v_float);  

    m_A1        = A1.convert(m_vc);
    m_A2        = A2.convert(m_vc);
    m_A3        = A3.convert(m_vc);
    m_base_mat  = details::convert_value(base, m_vc);

    m_ti        = m_A1.get_type();
    m_M         = m_A1.rows();
    m_N         = m_A3.cols();
    m_has_base  = true;
    m_modified  = A1.is_modified() || A2.is_modified() || A3.is_modified();
};

linsolve_obj_seq_3::~linsolve_obj_seq_3()
{};

bool linsolve_obj_seq_3::is_direct() const
{
    return m_A1.is_direct() && m_A2.is_direct() && m_A3.is_direct();
};

Real linsolve_obj_seq_3::log_det() const
{
    return m_A1.log_det() + m_A2.log_det() + m_A3.log_det();
};

Real linsolve_obj_seq_3::mat_normest_1() const
{
    if (m_has_base)
        return norm(m_base_mat, 1.0);
    else
        return linsolve_obj_base::mat_normest_1();
}
Real linsolve_obj_seq_3::mat_normest_inf() const
{
    if (m_has_base)
        return norm(m_base_mat, constants::inf());
    else
        return linsolve_obj_base::mat_normest_inf();
}

bool linsolve_obj_seq_3::is_hermitian() const
{
    if (m_has_base == false)
        return false;

    bool is_real    = matrix_traits::is_float_real(m_base_mat.get_value_code());
    return m_base_mat.get_struct().is_hermitian(m_base_mat.is_square(), is_real);
};

bool linsolve_obj_seq_3::is_posdef() const
{
    if (m_has_base == false)
        return false;

    return matcl::is_posdef(m_base_mat.get_struct());
};

linsolve_obj_seq_3::data_ptr 
linsolve_obj_seq_3::convert(value_code vc) const
{
    if (m_has_base == false)
        return data_ptr(new linsolve_obj_seq_3(m_A1.convert(vc), m_A2.convert(vc), m_A3.convert(vc)));

    Matrix Ac = details::convert_value(m_base_mat, vc);
    return data_ptr(new linsolve_obj_seq_3(Ac, m_A1.convert(vc), m_A2.convert(vc), m_A3.convert(vc)));
};

matcl::Matrix linsolve_obj_seq_3::solve(const Matrix& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //(A1*A2*A3)Y = X => Y = inv(A3)*inv(A2)*inv(A1)*X
        Matrix Y    = m_A1.solve(X, trans_type::no_trans);
        Y           = m_A2.solve(std::move(Y), trans_type::no_trans);
        Y           = m_A3.solve(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // op(A1*A2*A3)Y = X => op(A3)*op(A2)*op(A1)Y = X = > 
        // Y = op(op(A1))*inv(op(A2))*inv(op(A3))*X
        Matrix Y    = m_A3.solve(X, tA);
        Y           = m_A2.solve(std::move(Y), tA);
        Y           = m_A1.solve(std::move(Y), tA);
        return Y;
    };
};
matcl::Matrix linsolve_obj_seq_3::solve(Matrix&& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //(A1*A2*A3)Y = X => Y = inv(A3)*inv(A2)*inv(A1)*X
        Matrix Y    = m_A1.solve(std::move(X), trans_type::no_trans);
        Y           = m_A2.solve(std::move(Y), trans_type::no_trans);
        Y           = m_A3.solve(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // op(A1*A2*A3)Y = X => op(A3)*op(A2)*op(A1)Y = X = > 
        // Y = op(op(A1))*inv(op(A2))*inv(op(A3))*X
        Matrix Y    = m_A3.solve(std::move(X), tA);
        Y           = m_A2.solve(std::move(Y), tA);
        Y           = m_A1.solve(std::move(Y), tA);
        return Y;
    };
};

/// solve Y * op(A) = X, where A is the matrix represented by this object
matcl::Matrix linsolve_obj_seq_3::solve_rev(const Matrix& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //Y*(A1*A2*A3) = X => Y = X*inv(A3)*inv(A2)*inv(A1)
        Matrix Y    = m_A3.solve_rev(X, trans_type::no_trans);
        Y           = m_A2.solve_rev(std::move(Y), trans_type::no_trans);
        Y           = m_A1.solve_rev(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // Y*op(A1*A2*A3) = X => Y*op(A3)*op(A2)*op(A1) = X = > 
        // Y = X*op(op(A1))*inv(op(A2))*inv(op(A3))
        Matrix Y    = m_A1.solve_rev(X, tA);
        Y           = m_A2.solve_rev(std::move(Y), tA);
        Y           = m_A3.solve_rev(std::move(Y), tA);
        return Y;
    };
};
matcl::Matrix linsolve_obj_seq_3::solve_rev(Matrix&& X, trans_type tA) const
{
    if (tA == trans_type::no_trans)
    {
        //Y*(A1*A2*A3) = X => Y = X*inv(A3)*inv(A2)*inv(A1)
        Matrix Y    = m_A3.solve_rev(std::move(X), trans_type::no_trans);
        Y           = m_A2.solve_rev(std::move(Y), trans_type::no_trans);
        Y           = m_A1.solve_rev(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        // Y*op(A1*A2*A3) = X => Y*op(A3)*op(A2)*op(A1) = X = > 
        // Y = X*op(op(A1))*inv(op(A2))*inv(op(A3))
        Matrix Y    = m_A1.solve_rev(std::move(X), tA);
        Y           = m_A2.solve_rev(std::move(Y), tA);
        Y           = m_A3.solve_rev(std::move(Y), tA);
        return Y;
    };
};

matcl::Matrix linsolve_obj_seq_3::inv() const
{
    Matrix X    = m_A1.inv();
    X           = m_A2.solve(std::move(X), trans_type::no_trans);
    X           = m_A3.solve(std::move(X), trans_type::no_trans);

    if (is_hermitian() == true)
    {
        X.add_struct(predefined_struct_type::her);

        if (is_posdef() == true)
            X.add_struct(posdef_flag());
    };

    return X;
};

matcl::Matrix linsolve_obj_seq_3::base_matrix() const
{
    if (m_has_base)
        return m_base_mat;
    else
        return chain_mult(m_A1.base_matrix(), m_A2.base_matrix(), m_A3.base_matrix());
};

Matrix linsolve_obj_seq_3::mmul_right(const Matrix& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    if (m_has_base)
        return mmul(m_base_mat, X, t);

    if (t == trans_type::no_trans)
    {
        //(A1*A2*A3)Y
        Matrix Y    = m_A3.mmul_right(X, trans_type::no_trans);
        Y           = m_A2.mmul_right(std::move(Y), trans_type::no_trans);
        Y           = m_A1.mmul_right(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //op(A1*A2*A3) * Y = op(A3)*op(A2)*op(A1) * Y
        Matrix Y    = m_A1.mmul_right(X, t);
        Y           = m_A2.mmul_right(std::move(Y), t);
        Y           = m_A3.mmul_right(std::move(Y), t);
        return Y;
    };
}

Matrix linsolve_obj_seq_3::mmul_right(Matrix&& X, trans_type t) const
{
    error::check_mul(this->rows(), this->cols(), X.rows(), X.cols(), t, trans_type::no_trans);

    if (m_has_base)
        return mmul(m_base_mat, std::move(X), t);

    if (t == trans_type::no_trans)
    {
        //(A1*A2*A3)Y
        Matrix Y    = m_A3.mmul_right(std::move(X), trans_type::no_trans);
        Y           = m_A2.mmul_right(std::move(Y), trans_type::no_trans);
        Y           = m_A1.mmul_right(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //op(A1*A2*A3) * Y = op(A3)*op(A2)*op(A1) * Y
        Matrix Y    = m_A1.mmul_right(X, t);
        Y           = m_A2.mmul_right(std::move(Y), t);
        Y           = m_A3.mmul_right(std::move(Y), t);
        return Y;
    };
}

Matrix linsolve_obj_seq_3::mmul_left(const Matrix& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);

    if (m_has_base)
        return mmul(X, m_base_mat, trans_type::no_trans, t);

    if (t == trans_type::no_trans)
    {
        //Y*(A1*A2*A3)
        Matrix Y    = m_A1.mmul_left(X, trans_type::no_trans);
        Y           = m_A2.mmul_left(std::move(Y), trans_type::no_trans);
        Y           = m_A3.mmul_left(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //Y * op(A1*A2*A3) = Y * op(A3)*op(A2)*op(A1)
        Matrix Y    = m_A3.mmul_left(X, t);
        Y           = m_A2.mmul_left(std::move(Y), t);
        Y           = m_A1.mmul_left(std::move(Y), t);
        return Y;
    };
}

Matrix linsolve_obj_seq_3::mmul_left(Matrix&& X, trans_type t) const
{
    error::check_mul(X.rows(), X.cols(), this->rows(), this->cols(), trans_type::no_trans, t);

    if (m_has_base)
        return mmul(std::move(X), m_base_mat, trans_type::no_trans, t);

    if (t == trans_type::no_trans)
    {
        //Y*(A1*A2*A3)
        Matrix Y    = m_A1.mmul_left(std::move(X), trans_type::no_trans);
        Y           = m_A2.mmul_left(std::move(Y), trans_type::no_trans);
        Y           = m_A3.mmul_left(std::move(Y), trans_type::no_trans);
        return Y;
    }
    else
    {
        //Y * op(A1*A2*A3) = Y * op(A3)*op(A2)*op(A1)
        Matrix Y    = m_A3.mmul_left(std::move(X), t);
        Y           = m_A2.mmul_left(std::move(Y), t);
        Y           = m_A1.mmul_left(std::move(Y), t);
        return Y;
    };
}

//--------------------------------------------------------------------------------
//                  linsolve_obj_uhess
//--------------------------------------------------------------------------------

linsolve_obj_uhess::linsolve_obj_uhess(const Matrix& A, const options& opts)
    :linsolve_obj_seq_2()
{
    unitary_matrix Q;
    Matrix U;

    if (A.get_struct_code() == struct_code::struct_dense)
        tie(Q, U)   = qr_givens(A);
    else
        tie(Q, U)   = qr2(A);

    linsolve_obj_seq_2::reset(A, linsolve_unitary(Q), linsolve_triang(U, opts));
};

linsolve_obj_uhess::~linsolve_obj_uhess()
{};

//--------------------------------------------------------------------------------
//                  linsolve_obj_lhess_dense
//--------------------------------------------------------------------------------
linsolve_obj_lhess_dense::linsolve_obj_lhess_dense(const Matrix& A, const options& opts)
    :linsolve_obj_seq_2()
{
    Matrix L;
    unitary_matrix Q;    

    tie(Q, L)       = lq_givens(A);

    linsolve_obj_seq_2::reset(A, linsolve_triang(L, opts), linsolve_unitary(Q));
};

linsolve_obj_lhess_dense::~linsolve_obj_lhess_dense()
{};

};};