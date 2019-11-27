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

#include "matcl-linalg/decompositions/qr.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-linalg/utils/optim_params.h"
#include "matcl-internals/func/inplace.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-blas-lapack/blas/blas.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-linalg/decompositions/qr_band.h"
#include "matcl-linalg/decompositions/quern/matcl_quern_solver.h"
#include "matcl-linalg/decompositions/householder_q.h"
#include "matcl-linalg/decompositions/givens.h"

namespace matcl
{
    enum class qr_type
    {
        qr, ql, lq, rq
    };
};

namespace matcl { namespace details
{

template<class V, class S>
struct qr_givens_str
{};

template<class V>
struct qr_givens_str<V,struct_dense>
{
    using Mat = raw::Matrix<V,struct_dense>;
    static void eval(qr2_return& ret, const Mat& A0, qr_type qt, Integer ldiags, Integer udiags)
    {
        using VR        = typename md::real_type<V>::type;
        using Mat_R     = raw::Matrix<VR,struct_dense>;
        using Mat_I     = raw::Matrix<Integer,struct_dense>;

        Mat A           = A0.make_unique();
        A.get_struct().reset();

        Integer M       = A.rows();
        Integer N       = A.cols();

        V* ptr_A        = A.ptr();
        Integer A_ld    = A.ld();

        Integer SLEN;

        if (qt == qr_type::qr || qt == qr_type::rq)
            SLEN        = ldiags * M;
        else
            SLEN        = udiags * N;

        Mat_R C         = Mat_R(ti::ti_type<VR>(), SLEN, 1);
        Mat S           = Mat(ti::ti_type<V>(), SLEN, 1);
        Mat_I IND       = Mat_I(ti::ti_type<Integer>(), SLEN, 2);

        VR* ptr_C       = C.ptr();
        V* ptr_S        = S.ptr();
        Integer* ptr_I1 = IND.ptr();
        Integer* ptr_I2 = IND.ptr() + IND.ld();

        Integer SLEN_ALL    = 0;
        Integer NU          = 0;
        bool from_left      = false;

        if (qt == qr_type::qr)
        {
            //anihilate subdiagonal
            Integer INFO;
            lapack::huundl(M, N, ldiags, lap(ptr_A), A_ld, lap(ptr_C), lap(ptr_S), ptr_I1, 
                            ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid parameter passed to huundl");

            SLEN_ALL    += SLEN;
            ptr_C       += SLEN;
            ptr_S       += SLEN;
            ptr_I1      += SLEN;
            ptr_I2      += SLEN;

            A.get_struct().add(predefined_struct_type::triu);

            NU              = M;
            from_left       = true;
        }
        else if (qt == qr_type::rq)
        {
            //it should already be checked, that M == N

            //anihilate subdiagonal
            Integer INFO;
            lapack::huundr(N, ldiags, 0, lap(ptr_A), A_ld, lap(ptr_C), lap(ptr_S), ptr_I1, 
                            ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid parameter passed to huundr");

            SLEN_ALL    += SLEN;
            ptr_C       += SLEN;
            ptr_S       += SLEN;
            ptr_I1      += SLEN;
            ptr_I2      += SLEN;

            A.get_struct().add(predefined_struct_type::triu);

            NU              = N;
            from_left       = false;
        }
        if (qt == qr_type::lq)
        {
            //anihilate subdiagonal
            Integer INFO;
            lapack::hlundr(M, N, udiags, lap(ptr_A), A_ld, lap(ptr_C), lap(ptr_S), ptr_I1, 
                            ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid parameter passed to hlundr");

            SLEN_ALL    += SLEN;
            ptr_C       += SLEN;
            ptr_S       += SLEN;
            ptr_I1      += SLEN;
            ptr_I2      += SLEN;

            A.get_struct().add(predefined_struct_type::tril);

            NU              = N;
            from_left       = false;
        }
        else if (qt == qr_type::ql)
        {
            //it should already be checked, that M == N

            //anihilate subdiagonal
            Integer INFO;
            lapack::hlundl(N, udiags, lap(ptr_A), A_ld, lap(ptr_C), lap(ptr_S), ptr_I1, 
                            ptr_I2, SLEN, INFO);

            if (INFO != 0)
                throw error::error_general("invalid parameter passed to hlundl");

            SLEN_ALL    += SLEN;
            ptr_C       += SLEN;
            ptr_S       += SLEN;
            ptr_I1      += SLEN;
            ptr_I2      += SLEN;

            A.get_struct().add(predefined_struct_type::tril);

            NU              = M;
            from_left       = true;
        }

        C.assign_to_fresh(C.resize(SLEN_ALL,1));
        S.assign_to_fresh(S.resize(SLEN_ALL,1));
        IND.assign_to_fresh(IND.resize(SLEN_ALL,2));

        unitary_matrix Q = givens_to_unitary(NU, Matrix(C,false), Matrix(S,false), 
                                             Matrix(IND,false), from_left);
        Q                = ctrans(Q);

        ret = qr2_return(Q, Matrix(A,true));
        return;
    };
};

template<class V>
struct qr_givens_str<V,struct_sparse>
{
    using Mat   = raw::Matrix<V,struct_sparse>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval(qr2_return& ret, const Mat& A, qr_type qt, Integer ldiags, Integer udiags)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return qr_givens_str<V,struct_dense>::eval(ret, Ac,qt,ldiags,udiags);
    };
};

template<class V>
struct qr_givens_str<V,struct_banded>
{
    using Mat   = raw::Matrix<V,struct_banded>;
    using Mat_D = raw::Matrix<V,struct_dense>;

    static void eval(qr2_return& ret, const Mat& A, qr_type qt, Integer ldiags, Integer udiags)
    {
        Mat_D Ac    = raw::converter<Mat_D,Mat>::eval(A);
        return qr_givens_str<V,struct_dense>::eval(ret, Ac,qt,ldiags,udiags);
    };
};

template<class V, class S>
struct qr_givens_val
{
    using Mat = raw::Matrix<V,S>;
    static void eval(qr2_return& ret, const Mat& A, qr_type qt, Integer ldiags, Integer udiags)
    {
        return qr_givens_str<V,S>::eval(ret, A, qt, ldiags, udiags);
    };
};

template<class S>
struct qr_givens_val<Integer,S>
{
    using Mat   = raw::Matrix<Integer,S>;
    using Mat_R = raw::Matrix<Real,S>;
    static void eval(qr2_return& ret, const Mat& A, qr_type qt, Integer ldiags, Integer udiags)
    {
        Mat_R Ac    = raw::converter<Mat_R,Mat>::eval(A);
        return qr_givens_val<Real, S>::eval(ret, Ac, qt, ldiags, udiags);
    };
};

struct qr_givens_vis : public extract_type_switch<void, qr_givens_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, qr2_return& ret, qr_type qt, Integer ldiags, Integer udiags)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;

        return qr_givens_val<V, S>::eval(ret, mat, qt, ldiags, udiags);
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T&, qr2_return& ret, 
                            qr_type, Integer, Integer)
    {
        using VR = typename md::real_type_int_real<T>::type;
        ret = qr2_return(unitary_matrix(VR(1.0),false),Matrix(handle));
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, qr2_return&, 
                     qr_type, Integer, Integer)
    {
        throw error::object_value_type_not_allowed("qr_givens");
    };

    static void eval_scalar(const Matrix&, const Object&, qr2_return&, qr_type, Integer, Integer)
    {
        throw error::object_value_type_not_allowed("qr_givens");
    };
};

static void qr_givens_impl(qr2_return& ret, const Matrix& A, qr_type qt)
{
    matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
    matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

    Integer M   = A.rows();
    Integer N   = A.cols();

    if (qt == qr_type::rq || qt == qr_type::ql)
    {
        if (M != N)
            throw error::invalid_size(M,N);
    };

    if (A.structural_nnz() == 0)
    {
        if (qt == qr_type::rq || qt == qr_type::qr)
            A.add_struct(predefined_struct_type::triu);
        else
            A.add_struct(predefined_struct_type::tril);

        Matrix Q = speye(A.rows(),A.rows(), vt);
        ret = qr2_return(unitary_matrix(Q,false), A);
        return;
    };

    if (is_unitary(A.get_struct()) == true)
    {
        ret = qr2_return(unitary_matrix(A,false), speye(A.rows(), A.rows(), vt));
        return;
    }

    Integer ldiags  = matcl::get_ld(A, -1);
    Integer udiags  = matcl::get_ud(A, -1);

    if (ldiags == 0)
    {
        A.add_struct(predefined_struct_type::triu);

        if (qt == qr_type::qr || qt == qr_type::rq)
        {            
            Matrix Q = speye(A.rows(),A.rows(), vt);
            ret = qr2_return(unitary_matrix(Q,false), A);
            return;
        };
    }
    if (udiags == 0)
    {
        A.add_struct(predefined_struct_type::tril);

        if (qt == qr_type::ql || qt == qr_type::lq)
        {            
            Matrix Q = speye(A.rows(),A.rows(), vt);
            ret = qr2_return(unitary_matrix(Q,false), A);
            return;
        };
    };

    return details::qr_givens_vis::make<const Matrix&>(A, ret, qt, ldiags, udiags);
};

}};

namespace matcl
{

qr2_return matcl::qr_givens(const Matrix& A0)
{
    Matrix A(A0);

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::qr);
    return ret;
};
qr2_return matcl::qr_givens(Matrix&& A0)
{
    Matrix A(std::move(A0));

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::qr);
    return ret;
};

qr2_return matcl::ql_givens(const Matrix& A0)
{
    Matrix A(A0);

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::ql);
    return ret;
};
qr2_return matcl::ql_givens(Matrix&& A0)
{
    Matrix A(std::move(A0));

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::ql);
    return ret;
};

qr2_return matcl::rq_givens(const Matrix& A0)
{
    Matrix A(A0);

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::rq);
    return ret;
};
qr2_return matcl::rq_givens(Matrix&& A0)
{
    Matrix A(std::move(A0));

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::rq);
    return ret;
};

qr2_return matcl::lq_givens(const Matrix& A0)
{
    Matrix A(A0);

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::lq);
    return ret;
};
qr2_return matcl::lq_givens(Matrix&& A0)
{
    Matrix A(std::move(A0));

    qr2_return ret;
    details::qr_givens_impl(ret, A, qr_type::lq);
    return ret;
};

};
