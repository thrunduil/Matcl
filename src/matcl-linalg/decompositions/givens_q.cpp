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

#include "matcl-linalg/decompositions/givens_q.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/decompositions/qr_utils.h"

namespace matcl { namespace details
{

//------------------------------------------------------------------------
//                          UTILS
//------------------------------------------------------------------------
const char* get_trans_char(trans_type t)
{
    switch (t)
    {
        case trans_type::no_trans:      return "N";
        case trans_type::trans:         return "T";
        case trans_type::conj_trans:    return "C";
        default:                        return "N"; //should not happened
    };
}

template<class Val1, class Val2, class Struct>
struct givens_mult_struct{};

//------------------------------------------------------------------------
//                          DENSE
//------------------------------------------------------------------------
template<class Val1, class Val2>
struct givens_mult_struct<Val1, Val2, struct_dense>
{
    using Mat           = raw::Matrix<Val2,struct_dense>;
    using umatrix       = givens_q<Val1>;

    static const bool is_real_umatrix   = md::is_float_real_scalar<Val1>::value;

    static void eval_right(Matrix& ret, const Mat& X, const umatrix& data, trans_type t_unitary)
    {
        using VR        = typename md::real_type<Val1>::type;
        using VTR       = typename md::unify_types<Val1,Val2>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;

        Integer ret_N   = X.cols();
        Integer ret_M   = data.m_mat_size;

        ret_type Y      = make_copy<VTR, Val2>::eval(X, ret_M, ret_N);
        VTR* ptr_Y      = Y.ptr();
        Integer Y_ld    = Y.ld();

        const char* TYPE    = data.m_from_left ? "L" : "R";
        const char* SIDE    = "L";
        const char* TRANS   = get_trans_char(t_unitary);
        Integer K           = data.seq_size();

        const VR* ptr_C     = data.m_C.get().ptr();
        const Val1* ptr_S   = data.m_S.get().ptr();
        const Integer*ptr_I = data.m_I.get().ptr();
        Integer I_ld        = data.m_I.get().ld();

        const Integer* ptr_I1   = ptr_I;
        const Integer* ptr_I2   = ptr_I + I_ld;
        Integer LDIAGS      = matcl::get_ld(Matrix(X,false), -1);
        Integer UDIAGS      = matcl::get_ud(Matrix(X,false), -1);

        Integer INFO;

        using VL            = typename md::lapack_value_type<Val1>::type;
        lapack::rotseq<VL>(TYPE, SIDE, TRANS, K, lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, X.rows(), X.cols(),
                       LDIAGS, UDIAGS, lap(ptr_Y), Y_ld, INFO);
       
        if (INFO != 0)
            throw error::error_general("invalid parameter passed to rotseq");

        struct_flag su;
        su.set_user(unitary_flag());

        bool is_square  = Y.rows() == Y.cols();
        bool is_sq_X    = X.rows() == X.cols();
        Y.set_struct(predefined_struct_ext::mult_struct(su,X.get_struct(),t_unitary,trans_type::no_trans,
                            is_real_umatrix, is_real_matrix(X), true, is_sq_X, is_square));

        if (LDIAGS == 0)
            Y.add_struct(predefined_struct_type::triu);
        else if (LDIAGS == 1)
            Y.add_struct(predefined_struct_type::hessu);

        if (UDIAGS == 0)
            Y.add_struct(predefined_struct_type::tril);
        else if (UDIAGS == 1)
            Y.add_struct(predefined_struct_type::hessl);

        ret = Matrix(Y,true);
        return;
    };

    static void eval_left(Matrix& ret, const Mat& X, const umatrix& data, trans_type t_unitary)
    {
        using VR        = typename md::real_type<Val1>::type;
        using VTR       = typename md::unify_types<Val1,Val2>::type;
        using ret_type  = raw::Matrix<VTR,struct_dense>;

        Integer ret_N   = data.m_mat_size;
        Integer ret_M   = X.rows();

        ret_type Y      = make_copy<VTR, Val2>::eval(X, ret_M, ret_N);

        VTR* ptr_Y      = Y.ptr();
        Integer Y_ld    = Y.ld();

        const char* TYPE    = data.m_from_left ? "L" : "R";
        const char* SIDE    = "R";
        const char* TRANS   = get_trans_char(t_unitary);
        Integer K           = data.seq_size();

        const VR* ptr_C     = data.m_C.get().ptr();
        const Val1* ptr_S   = data.m_S.get().ptr();
        const Integer*ptr_I = data.m_I.get().ptr();
        Integer I_ld        = data.m_I.get().ld();

        const Integer* ptr_I1   = ptr_I;
        const Integer* ptr_I2   = ptr_I + I_ld;
        Integer LDIAGS      = matcl::get_ld(Matrix(X,false), -1);
        Integer UDIAGS      = matcl::get_ud(Matrix(X,false), -1);

        Integer INFO;

        using VL            = typename md::lapack_value_type<Val1>::type;
        lapack::rotseq<VL>(TYPE, SIDE, TRANS, K, lap(ptr_C), lap(ptr_S), ptr_I1, ptr_I2, X.rows(), X.cols(),
                       LDIAGS, UDIAGS, lap(ptr_Y), Y_ld, INFO);
       
        if (INFO != 0)
            throw error::error_general("invalid parameter passed to rotseq");

        struct_flag su;
        su.set_user(unitary_flag());

        bool is_square  = Y.rows() == Y.cols();
        bool is_sq_X    = X.rows() == X.cols();
        Y.set_struct(predefined_struct_ext::mult_struct(X.get_struct(),su,trans_type::no_trans,t_unitary,
                            is_real_matrix(X), is_real_umatrix, is_sq_X, true, is_square));

        if (LDIAGS == 0)
            Y.add_struct(predefined_struct_type::triu);
        else if (LDIAGS == 1)
            Y.add_struct(predefined_struct_type::hessu);

        if (UDIAGS == 0)
            Y.add_struct(predefined_struct_type::tril);
        else if (UDIAGS == 1)
            Y.add_struct(predefined_struct_type::hessl);

        ret = Matrix(Y,true);
        return;
    };
};

//------------------------------------------------------------------------
//                          SPARSE
//------------------------------------------------------------------------
template<class Val1, class Val2>
struct givens_mult_struct<Val1, Val2, struct_sparse>
{
    using Mat           = raw::Matrix<Val2,struct_sparse>;
    using umatrix       = givens_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return givens_mult_struct<Val1, Val2, struct_dense>::eval_right(ret,md,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return givens_mult_struct<Val1, Val2, struct_dense>::eval_left(ret,md,data,t_unitary);
    };
};

//------------------------------------------------------------------------
//                          BAND
//------------------------------------------------------------------------
template<class Val1, class Val2>
struct givens_mult_struct<Val1, Val2, struct_banded>
{
    using Mat           = raw::Matrix<Val2,struct_banded>;
    using umatrix       = givens_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return givens_mult_struct<Val1, Val2, struct_dense>::eval_right(ret,md,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using DM        = raw::Matrix<Val2, struct_dense>;
        DM md           = raw::converter<DM,Mat>::eval(mat);
        return givens_mult_struct<Val1, Val2, struct_dense>::eval_left(ret,md,data,t_unitary);
    };
};

//------------------------------------------------------------------------
//                          givens_helper
//------------------------------------------------------------------------
template<class Val1, class Val2, class Struct>
struct givens_mult_impl
{
    using Mat           = raw::Matrix<Val2,Struct>;
    using umatrix       = givens_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        return givens_mult_struct<Val1,Val2,Struct>::eval_right(ret,mat,data,t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        return givens_mult_struct<Val1,Val2,Struct>::eval_left(ret,mat,data,t_unitary);
    };
};

template<class Val1, class Struct>
struct givens_mult_impl<Val1, Integer, Struct>
{
    using Mat           = raw::Matrix<Integer,Struct>;
    using umatrix       = givens_q<Val1>;

    static void eval_right(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using Mat_D     = raw::Matrix<Real,struct_dense>;
        Mat_D md        = raw::converter<Mat_D,Mat>::eval(mat);
        return givens_mult_impl<Val1,Real,struct_dense>::eval_right(ret, md, data, t_unitary);
    };
    static void eval_left(Matrix& ret, const Mat& mat, const umatrix& data, trans_type t_unitary)
    {
        using Mat_D     = raw::Matrix<Real,struct_dense>;
        Mat_D md        = raw::converter<Mat_D,Mat>::eval(mat);
        return givens_mult_impl<Val1,Real,struct_dense>::eval_left(ret, md, data, t_unitary);
    };
};

template<class Val>
struct givens_mult_right : public extract_type_switch<void, givens_mult_right<Val>,true>
{
    using umatrix   = givens_q<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using VM    = typename T::value_type;
        using ST    = typename T::struct_type;

        return givens_mult_impl<Val,VM,ST>::eval_right(ret, mat, data, t_unitary);
    };
    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, t_unitary, data);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_right");
    };
};

template<class Val>
struct givens_mult_left : public extract_type_switch<void, givens_mult_left<Val>,true>
{
    using umatrix   = givens_q<Val>;

    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using VM    = typename T::value_type;
        using ST    = typename T::struct_type;

        return givens_mult_impl<Val,VM, ST>::eval_left(ret, mat, data, t_unitary);
    };
    template<class T>
    static void eval_scalar(const Matrix& h, const T& mat, Matrix& ret, trans_type t_unitary, 
                            const umatrix& data)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM full_mat = DM(ti::get_ti<T>(mat),mat,1,1);
        return eval(h,full_mat, ret, t_unitary, data);
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, trans_type, 
                            const umatrix&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::mult_left");
    };
};

template<class Val>
struct givens_q_convert_vis : public extract_type_switch<void, givens_q_convert_vis<Val>,true>
{
    using umatrix   = givens_q<Val>;
    using data_ptr  = typename unitary_matrix_data::data_ptr;

    template<class T>
    static void eval(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        using VM    = typename T::value_type;

        ret = data.convert_impl<VM>();
    };
    template<class T>
    static void eval_scalar(const Matrix&, const T&, const umatrix& data, data_ptr& ret)
    {
        ret = data.convert_impl<T>();
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
    static void eval_scalar(const Matrix&, const Object&, const umatrix&, data_ptr&)
    {
        throw error::object_value_type_not_allowed("unitary_matrix::convert");
    };
};

//------------------------------------------------------------------------
//                          givens_helper
//------------------------------------------------------------------------
template<class Val>
struct givens_helper
{
    using Mat       = raw::Matrix<Val,struct_dense>;
    using umatrix   = givens_q<Val>;

    static void to_matrix(Matrix& ret, const umatrix& data);
    static void eval_mult_right(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                                const umatrix& data);
    static void eval_mult_left(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                                const umatrix& data);
};

template<class Val>
void givens_helper<Val>::to_matrix(Matrix& ret, const umatrix& data)
{
    Matrix I = speye(data.rows(), data.cols(), data.get_value_code());
    return eval_mult_right(ret, I, trans_type::no_trans, data);
};

template<class Val>
void givens_helper<Val>::eval_mult_right(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                            const umatrix& data)
{
    return givens_mult_right<Val>::make<const Matrix&>(mat, ret, t_unitary, data);
};

template<class Val>
void givens_helper<Val>::eval_mult_left(Matrix& ret, const matcl::Matrix& mat, trans_type t_unitary, 
                            const umatrix& data)
{
    return givens_mult_left<Val>::make<const Matrix&>(mat, ret, t_unitary, data);
};

//------------------------------------------------------------------------
//                          givens_q
//------------------------------------------------------------------------
template<class Val>
givens_q<Val>::givens_q()
    :m_mat_size(0), m_C(ti::ti_type<VR>()), m_S(ti::ti_type<Val>()), m_I(ti::ti_type<Integer>())
    ,m_from_left(false)
{};

template<class Val>
givens_q<Val>::givens_q(Integer N, const Mat_R& C, const Mat& S, const Mat_I& I, bool from_left)
    :m_mat_size(N), m_C(C), m_S(S), m_I(I), m_from_left(from_left)
{};

template<class Val>
givens_q<Val>::~givens_q()
{};

template<class Val>
serialization_helper<unitary_matrix_data>*
givens_q<Val>::get_serialization_helper() const
{
    std::ostringstream msg;
    msg << "givens_q" << "<" << get_type_name<Val>::eval() << ">";

    return serialization_helper<unitary_matrix_data>
            ::get<givens_q<Val>>(msg.str());
};

template<class Val>
void givens_q<Val>::save(oarchive& os) const
{
    //general parameters
    os << this->m_mat_size;
    os << this->m_from_left;

    os.get() << m_C.get();
    os.get() << m_S.get();
    os.get() << m_I.get();
};

template<class Val>
unitary_matrix_data* givens_q<Val>::load(iarchive& ar)
{
    using ptr_type  = std::unique_ptr<givens_q<Val>>;
    ptr_type ret    = ptr_type(new givens_q<Val>());

    Integer fl;
    //general parameters
    ar >> ret->m_mat_size;
    ar >> fl;

    ret->m_from_left    = fl ? true : false;

    Mat_R C_loc(ret->m_C.get().get_type());
    Mat   S_loc(ret->m_S.get().get_type());
    Mat_I I_loc(ret->m_I.get().get_type());

    ar.get() >> C_loc;
    ar.get() >> S_loc;
    ar.get() >> I_loc;

    ret->m_C.rebind(C_loc);
    ret->m_S.rebind(S_loc);
    ret->m_I.rebind(I_loc);

    return ret.release();
};

template<class Val>
void givens_q<Val>::save(std::ostream& os) const
{
    os << " ";

    //general parameters
    os << this->m_mat_size << " ";
    os << this->m_from_left << " ";

    os << Matrix(m_C.get(), false);
    os << Matrix(m_S.get(), false);
    os << Matrix(m_I.get(), false);
};

template<class Val>
unitary_matrix_data* givens_q<Val>::load(std::istream& is)
{
    using ptr_type = std::unique_ptr<givens_q<Val>>;
    ptr_type ret = ptr_type(new givens_q<Val>());

    //general parameters
    is >> ret->m_mat_size;
    is >> ret->m_from_left;

    Matrix C, S, I;
    is >> C;
    is >> S;
    is >> I;

    Matrix Cc   = matcl::convert(C, Mat_R::matrix_code);
    Matrix Sc   = matcl::convert(S, Mat::matrix_code);
    Matrix Ic   = matcl::convert(C, Mat_I::matrix_code);

    ret->m_C.rebind(Cc.get_impl<Mat_R>());
    ret->m_S.rebind(Sc.get_impl<Mat>());
    ret->m_I.rebind(Ic.get_impl<Mat_I>());

    return ret.release();
};

template<class Val>
void givens_q<Val>::mult_right(matcl::Matrix& ret, const matcl::Matrix& X, 
                                     trans_type t_unit) const
{
    if (X.is_scalar())
    {
        this->to_matrix(ret);
        ret = mmul(ret,X,t_unit);
        return;
    };

    error::check_mul(m_mat_size, m_mat_size, X.rows(), X.cols(), t_unit, trans_type::no_trans);

    Integer M   = this->m_mat_size;
    Integer N   = X.cols();

    if (X.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(X.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        ret = spzeros(M, N, 0, vt);
        return;
    };

    return givens_helper<Val>::eval_mult_right(ret, X, t_unit, *this);
};

template<class Val>
void givens_q<Val>::mult_left(matcl::Matrix& ret, const matcl::Matrix& X, 
                                     trans_type t_unit) const
{
    if (X.is_scalar())
    {
        this->to_matrix(ret);
        ret = mmul(X,ret,trans_type::no_trans,t_unit);
        return;
    };

    error::check_mul(X.rows(), X.cols(), m_mat_size, m_mat_size, trans_type::no_trans, t_unit);

    Integer M   = X.rows();
    Integer N   = this->m_mat_size;

    if (X.structural_nnz() == 0)
    {
        matcl::value_code vt0 = matrix_traits::real_value_type(X.get_value_code());
        matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

        ret = spzeros(M, N, 0, vt);
        return;
    };

    return givens_helper<Val>::eval_mult_left(ret, X, t_unit, *this);
};


template<class Val>
void givens_q<Val>::to_matrix(matcl::Matrix& ret) const
{
    return givens_helper<Val>::to_matrix(ret, *this);
};

template<class Val>
typename givens_q<Val>::data_ptr 
givens_q<Val>::convert(value_code new_val_code) const
{
    value_code vc   = matrix_traits::unify_value_types(new_val_code, value_code::v_float);
    Matrix v        = zeros(0,0,vc);

    data_ptr ret;
    givens_q_convert_vis<Val>::make<const Matrix&>(v,*this, ret);
    return ret;
};

template<class Val>
template<class T>
typename givens_q<Val>::data_ptr
givens_q<Val>::convert_impl() const
{
    using TR        = typename details::unify_types<T, Float>::type;
    using TRR       = typename details::real_type<TR>::type;

    using Mat_C     = raw::Matrix<TR,struct_dense>;
    using Mat_RC    = raw::Matrix<TRR,struct_dense>;

    const Mat_C& S  = raw::converter<Mat_C, Mat>::eval(m_S.get());
    const Mat_RC& C = raw::converter<Mat_RC, Mat_R>::eval(m_C.get());

    return data_ptr(new givens_q<TR>(m_mat_size, C, S, m_I.get(), m_from_left));
};

template class givens_q<Real>;
template class givens_q<Float>;
template class givens_q<Complex>;
template class givens_q<Float_complex>;

}};