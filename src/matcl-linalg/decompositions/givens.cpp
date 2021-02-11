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

#include "matcl-linalg/decompositions/givens.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/decompositions/givens_q.h"

namespace matcl { namespace details
{

template<class Val>
struct construct_givens_str
{
    using Mat = raw::Matrix<Val,struct_dense>;

    static void eval(unitary_matrix& ret, const Mat& mat_S, Integer N, const Matrix& C, 
                     const Matrix& Ind, bool from_left)
    {
        using VR        = typename md::real_type<Val>::type;
        using Mat_R     = raw::Matrix<VR,struct_dense>;
        using Mat_I     = raw::Matrix<Integer,struct_dense>;

        const Mat_R&  mat_C   = C.impl<Mat_R>();
        const Mat_I&  mat_I   = Ind.impl<Mat_I>();

        using unitary_matrix_data_ptr = unitary_matrix::unitary_matrix_data_ptr;

        bool isv        = mat_S.all_finite() && mat_C.all_finite();

        if (isv == false)
        {
            ret = unitary_matrix::from_nan(N,N,matrix_traits::value_code<Val>::value);
            return;
        };

        unitary_matrix_data_ptr data(new givens_q<Val>(N,mat_C, mat_S, mat_I, from_left) );
        ret = unitary_matrix(data);
        return;
    };
};

template<class Val, class S>
struct construct_givens_val
{
    using Mat   = raw::Matrix<Val,S>;
    using Mat_D = raw::Matrix<Val,struct_dense>;

    static void eval(unitary_matrix& ret, const Mat& matS, Integer N, const Matrix& C, 
                     const Matrix& Ind, bool from_left)
    {
        Mat_D Sc    = raw::converter<Mat_D, Mat>::eval(matS);
        return construct_givens_str<Val>::eval(ret, Sc, N, C, Ind, from_left);
    };
};

template<class S>
struct construct_givens_val<Integer, S>
{
    using Mat   = raw::Matrix<Integer,S>;
    using Mat_R = raw::Matrix<Real,struct_dense>;

    static void eval(unitary_matrix& ret, const Mat& matS, Integer N, const Matrix& C, 
                     const Matrix& Ind, bool from_left)
    {
        Mat_R Sc    = raw::converter<Mat_R, Mat>::eval(matS);
        return construct_givens_val<Real, struct_dense>::eval(ret, Sc, N, C, Ind, from_left);
    };
};

struct givens_to_unitary_vis : public extract_type_switch<void, givens_to_unitary_vis,true>
{
    template<class T>
    static void eval(const Matrix&, const T& matS, unitary_matrix& ret, Integer N, const Matrix& C, 
                     const Matrix& Ind, bool from_left)
    {
        using V = typename T::value_type;
        using S = typename T::struct_type;
        return construct_givens_val<V,S>::eval(ret, matS, N, C, Ind, from_left);
    };

    template<class T>
    static void eval_scalar(const Matrix& h, const T& S, unitary_matrix& ret, Integer N, const Matrix& C, 
                     const Matrix& Ind, bool from_left)
    {
        using DM    = raw::Matrix<T,struct_dense>;
        DM Sc       = DM(ti::get_ti(S), S, 1, 1);
        return eval<DM>(h, Sc, ret, N, C, Ind, from_left); 
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, unitary_matrix&, Integer, const Matrix&, 
                     const Matrix&, bool)
    {
        throw error::object_value_type_not_allowed("givens_to_unitary");
    };

    static void eval_scalar(const Matrix&, const Object&, unitary_matrix&, Integer, const Matrix&, 
                     const Matrix&, bool)
    {
        throw error::object_value_type_not_allowed("givens_to_unitary");
    };
};

static void givens_to_unitary_impl(unitary_matrix& ret, Integer N, const Matrix& C, const Matrix& S, 
                                  const Matrix& Ind, bool from_left)
{
    Integer K   = C.length();

    // check size
    if (C.is_vector() == false)
        throw error::error_invalid_givens_seq_size(C.rows(), C.cols(), S.rows(),
                        S.cols(), Ind.rows(), Ind.cols(), 1);

    if (S.is_vector() == false)
        throw error::error_invalid_givens_seq_size(C.rows(), C.cols(), S.rows(),
                        S.cols(), Ind.rows(), Ind.cols(), 2);

    if (Ind.cols() != 2)
        throw error::error_invalid_givens_seq_size(C.rows(), C.cols(), S.rows(),
                        S.cols(), Ind.rows(), Ind.cols(), 3);

    if (S.length() != K)
        throw error::error_invalid_givens_seq_size(C.rows(), C.cols(), S.rows(),
                        S.cols(), Ind.rows(), Ind.cols(), 2);

    if (Ind.rows() != K)
        throw error::error_invalid_givens_seq_size(C.rows(), C.cols(), S.rows(),
                        S.cols(), Ind.rows(), Ind.cols(), 3);

    // check value types;
    value_code vc_C = C.get_value_code();
    value_code vc_S = S.get_value_code();
    value_code vc_I = Ind.get_value_code();

    if (vc_C == value_code::v_object || vc_S == value_code::v_object || vc_I == value_code::v_object)
        throw error::object_value_type_not_allowed("givens_to_unitary");

    if (vc_I != value_code::v_integer)
        throw error::error_invalid_givens_seq_value_code(vc_I, 3);

    if (matrix_traits::is_float_complex(vc_C) == true)
        throw error::error_invalid_givens_seq_value_code(vc_I, 1);

    // check indices
    using Mat_I = raw::Matrix<Integer,struct_dense>;

    const Mat_I& mat_I      = Ind.impl<Mat_I>();
    const Integer* ptr_I    = mat_I.ptr();
    Integer I_ld            = mat_I.ld();

    for (Integer j = 0; j < 2; ++j)
    {
        for (Integer i = 0; i < K; ++i)
        {
            Integer ind     = ptr_I[i];

            if (ind < 1 || ind > N)
                throw error::error_invalid_givens_seq_index(i + 1, j + 1, ind, N);
        };

        ptr_I               += I_ld;
    };

    // convert value codes
    value_code vc_S2        = matrix_traits::unify_value_types(vc_S, vc_C);
    value_code vc_Sr        = matrix_traits::unify_value_types(vc_S2, value_code::v_float);
    value_code vc_Cr        = matrix_traits::real_value_type(vc_Sr);

    Matrix C_c              = convert(C, matrix_traits::get_matrix_type(vc_Cr,struct_code::struct_dense));
    Matrix S_c              = convert(S, matrix_traits::get_matrix_type(vc_Sr,struct_code::struct_dense));

    // extract type based on S matrix
    return givens_to_unitary_vis::make<const Matrix&>(S_c, ret, N, C_c, Ind, from_left);

};

}};

namespace matcl
{

template<class Val, class Enable>
MATCL_LINALG_EXPORT
void matcl::construct_givens2(const Val& a, const Val& b, typename md::real_type<Val>::type& c,
        Val& s, Val& r)
{
    using matcl::details::lap;

    lapack::lartg(lap(a), lap(b), lap(&c), lap(&s), lap(&r));
};

unitary_matrix matcl::givens_to_unitary(Integer N, const Matrix& C, const Matrix& S, 
                                        const Matrix& Ind, bool from_left)
{
    unitary_matrix ret;
    details::givens_to_unitary_impl(ret, N, C, S, Ind, from_left);
    return ret;
};

template void construct_givens<Real,void>(const Real&, const Real&, Real&, Real&);
template void construct_givens<Float,void>(const Float&, const Float&, Float&, Float&);
template void construct_givens<Complex,void>(const Complex&, const Complex&, Real&, Complex&);
template void construct_givens<Float_complex,void>(const Float_complex&, const Float_complex&,
                                                   Float&, Float_complex&);

template void MATCL_LINALG_EXPORT
matcl::construct_givens2<Real,void>(const Real&, const Real&, Real&, Real&, Real&);

template void MATCL_LINALG_EXPORT 
construct_givens2<Float,void>(const Float&, const Float&, Float&, Float&, Float&);

template void MATCL_LINALG_EXPORT
construct_givens2<Complex,void>(const Complex&, const Complex&, Real&, Complex&, Complex&);

template void MATCL_LINALG_EXPORT
construct_givens2<Float_complex,void>(const Float_complex&, const Float_complex&,
                                                   Float&, Float_complex&, Float_complex&);

};
