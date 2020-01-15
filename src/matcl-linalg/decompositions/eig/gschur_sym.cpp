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

#include "matcl-linalg/decompositions/schur.h"
#include "matcl-linalg/decompositions/gschur_sym.h"
#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-blas-lapack/lapack/lapack.h"
#include "matcl-blas-ext/lapack_ext/lapack_ext.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"
#include "matcl-linalg/decompositions/eig/schur_utils.h"
#include "matcl-linalg/decompositions/qr.h"
#include "matcl-linalg/utils/linalg_utils.h"
#include "matcl-internals/base/utils.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/details/struct_flag_predefined.h"
#include "matcl-matrep/lib_functions/matrix_gen.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant

namespace matcl { namespace details
{

template<class str> struct str_to_code_gs{};
template<> struct str_to_code_gs<struct_banded>	{ static const int value = 1; };
template<> struct str_to_code_gs<struct_sparse>	{ static const int value = 2; };
template<> struct str_to_code_gs<struct_dense>  { static const int value = 2; };

template<int str> struct code_to_str_gs{};
template<> struct code_to_str_gs<1>             {using type = struct_banded; };
template<> struct code_to_str_gs<2>             {using type = struct_dense; };

template<class S1, class S2>
struct unify_struct
{
    static const int SC1    = str_to_code_gs<S1>::value;
    static const int SC2    = str_to_code_gs<S2>::value;
    static const int SCR    = (SC1 > SC2)? SC1 : SC2;
    using type              = typename code_to_str_gs<SCR>::type;
};

template<class Val, class Struct>
struct gschur_sym_str
{};

template<class Val>
struct gschur_sym_str<Val, struct_dense>
{
    using VC    = typename md::complex_type<Val>::type;
    using VR    = typename md::real_type<Val>::type;

    using Mat   = raw::Matrix<Val,struct_dense>;
    using Mat_C = raw::Matrix<VC,struct_dense>;
    using Mat_R = raw::Matrix<VR,struct_dense>;

    static void eval(const Mat& A0, const Mat& B0, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        Mat A   = A0.make_unique();
        Mat B   = B0.make_unique();

        A.set_struct(struct_flag());
        B.set_struct(struct_flag());

        int itype   = 1;
        switch(type)
        {
            case gschur_sym_type::A_B:
                itype    = 1;
                break;
            case gschur_sym_type::AB:
                itype    = 2;
                break;
            case gschur_sym_type::BA:
                itype    = 3;
                break;
        };

        const char* jobz    = with_V? "V" : "N";
        const char* uplo    = "U";
        Integer info        = 0;

        Integer N           = A.rows();

        Mat_R w(A.get_type(), N, 1);

        Val work_query;
        VR  rwork_query;
        Integer iwork_query;

        lapack::hegvd(itype, jobz, uplo, N, lap(A.ptr()), A.ld(), lap(B.ptr()), B.ld(),
                          lap(w.ptr()), lap(&work_query), -1, lap(&rwork_query), -1,
                          lap(&iwork_query), -1, &info);
        
        Integer lwork       = (Integer) real(work_query);
        Integer lrwork      = (Integer) rwork_query;
        Integer liwork      = iwork_query;

        using VTR_pod       = matcl::pod_type<Val>;
        using workspace     = md::workspace2<Val>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        using rworkspace    = matcl::pod_workspace<VR>;

        workspace WORK      = workspace(lwork);
        iworkspace IWORK    = iworkspace(liwork);
        rworkspace RWORK    = rworkspace(liwork);

        Val* ptr_WORK       = WORK.ptr();
        Integer* ptr_IWORK  = IWORK.ptr();
        VR* ptr_RWORK       = RWORK.ptr();

        lapack::hegvd(itype, jobz, uplo, N, lap(A.ptr()), A.ld(), lap(B.ptr()), B.ld(),
                          lap(w.ptr()), lap(ptr_WORK), lwork ,lap(ptr_RWORK), lrwork,
                          lap(ptr_IWORK), liwork, &info);

        if (info)
        {
            if(info > N)
                throw error::error_nonposdef(false);
            else if (info > 0)
                throw error::error_gschur();
            else
                throw error::error_general("invalid argument passed to hegvd");
        };
        
        gd.m_eig    = Matrix(w,true);
        gd.m_D      = bdiag(gd.m_eig);

        if (with_V)
        {
            gd.m_V  = Matrix(A,true);
        };

        if (gd.test_factors() == false)
            return;
    };    
};

template<class Val>
struct gschur_sym_str<Val, struct_banded>
{
    using VC    = typename md::complex_type<Val>::type;
    using VR    = typename md::real_type<Val>::type;

    using Mat   = raw::Matrix<Val,struct_banded>;
    using Mat_C = raw::Matrix<VC,struct_banded>;
    using Mat_R = raw::Matrix<VR,struct_dense>;
    using Mat_D = raw::Matrix<Val,struct_dense>;

    static void eval(const Mat& A0, const Mat& B0, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        if (type != gschur_sym_type::A_B)
        {
            //only this version is available for band matrices            

            Mat_D A = raw::converter<Mat_D,Mat>::eval(A0);
            Mat_D B = raw::converter<Mat_D,Mat>::eval(B0);
            return gschur_sym_str<Val,struct_dense>::eval(A,B,gd,type,with_V);
        };

        Mat A   = A0.make_unique();
        Mat B   = B0.make_unique();

        A.set_struct(struct_flag());
        B.set_struct(struct_flag());

        const char* jobz    = with_V? "V" : "N";
        const char* uplo    = "U";
        Integer info        = 0;

        if (A.has_diag(0) == false)
            throw error::band_matrix_with_main_diag_required(A.first_diag(), A.last_diag());
        if (B.has_diag(0) == false)
            throw error::band_matrix_with_main_diag_required(B.first_diag(), B.last_diag());

        Integer N           = A.rows();
        Integer KA          = A.number_superdiagonals();
        Integer KB          = B.number_superdiagonals();
        Integer NZ          = with_V? N : 1;

        //hbgvd requires, that KB <= KA
        if (KA < KB)
        {
            Mat Ar          = A.resize(A.rows(), A.cols(), -KB, KB); 
            A.assign_to_fresh(Ar);
            KA              = A.number_superdiagonals();
        };

        Mat_R w(A.get_type(), N, 1);
        Mat_D Z(A.get_type(), NZ, NZ);

        Val work_query;
        VR  rwork_query;
        Integer iwork_query;

        Val* ptr_A          = A.rep_ptr();
        Val* ptr_B          = B.rep_ptr();
        Val* ptr_Z          = Z.ptr();
        Integer ld_Z        = Z.ld();

        lapack::hbgvd(jobz, uplo, N, KA, KB, lap(ptr_A), A.ld(), lap(ptr_B), B.ld(),
                          lap(w.ptr()), lap(ptr_Z), ld_Z, lap(&work_query), -1, lap(&rwork_query), -1,
                          lap(&iwork_query), -1, info);
        
        Integer lwork       = (Integer) real(work_query);
        Integer lrwork      = (Integer) rwork_query;
        Integer liwork      = iwork_query;

        using VTR_pod       = matcl::pod_type<Val>;
        using workspace     = md::workspace2<Val>;
        using iworkspace    = matcl::pod_workspace<Integer>;
        using rworkspace    = matcl::pod_workspace<VR>;

        workspace WORK      = workspace(lwork);
        iworkspace IWORK    = iworkspace(liwork);
        rworkspace RWORK    = rworkspace(liwork);

        Val* ptr_WORK       = WORK.ptr();
        Integer* ptr_IWORK  = IWORK.ptr();
        VR* ptr_RWORK       = RWORK.ptr();

        lapack::hbgvd(jobz, uplo, N, KA, KB, lap(ptr_A), A.ld(), lap(ptr_B), B.ld(),
                          lap(w.ptr()), lap(ptr_Z), ld_Z, lap(ptr_WORK), lwork ,lap(ptr_RWORK), lrwork,
                          lap(ptr_IWORK), liwork, info);

        if (info)
        {
            if(info > N)
                throw error::error_nonposdef(false);
            else if (info > 0)
                throw error::error_gschur();
            else
                throw error::error_general("invalid argument passed to hegvd");
        };
        
        gd.m_eig    = Matrix(w,true);
        gd.m_D      = bdiag(gd.m_eig);

        if (with_V)
        {
            gd.m_V  = Matrix(Z,true);
        };

        if (gd.test_factors() == false)
            return;
    };
};

template<class Val>
struct gschur_sym_str<Val, struct_sparse>
{
    using Mat   = raw::Matrix<Val,struct_sparse>;
    using Mat_D = raw::Matrix<Val,struct_dense>;

    static void eval(const Mat& A, const Mat& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        Mat_D Ac = raw::converter<Mat_D,Mat>::eval(A);
        Mat_D Bc = raw::converter<Mat_D,Mat>::eval(B);
        return gschur_sym_str<Val,struct_dense>::eval(Ac,Bc,gd,type,with_V);
    };
};

template<class V1, class S1, class S2>
struct gschur_sym_impl2
{
    using Mat1  = raw::Matrix<V1,S1>;
    using Mat2  = raw::Matrix<V1,S2>;

    using VR        = V1;
    using SR        = typename md::unify_struct<S1,S2>::type;
    using Mat_R     = raw::Matrix<VR,SR>;

    static void eval(const Mat1& A, const Mat2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        const Mat_R& Ac = raw::converter<Mat_R, Mat1>::eval(A);
        const Mat_R& Bc = raw::converter<Mat_R, Mat2>::eval(B);

        return gschur_sym_impl2<VR,SR,SR>::eval(Ac,Bc,gd,type,with_V);
    };
};

template<class V1, class S1>
struct gschur_sym_impl2<V1,S1,S1>
{
    using Mat1  = raw::Matrix<V1,S1>;

    static void eval(const Mat1& A, const Mat1& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        return gschur_sym_str<V1,S1>::eval(A,B,gd,type,with_V);
    };
};

template<class V1, class V2, class S1, class S2>
struct gschur_sym_impl
{
    static_assert(std::is_same<V1,V2>::value == false, "invalid template specialization"); 

    using Mat1  = raw::Matrix<V1,S1>;
    using Mat2  = raw::Matrix<V2,S2>;

    using VR        = typename md::unify_types<V1,V2>::type;
    using SR        = typename md::unify_struct<S1,S2>::type;
    using Mat_R     = raw::Matrix<VR,SR>;

    static void eval(const Mat1& A, const Mat2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        Mat_R Ac    = raw::converter<Mat_R, Mat1>::eval(A);
        Mat_R Bc    = raw::converter<Mat_R, Mat2>::eval(B);

        return gschur_sym_impl<VR,VR,SR,SR>::eval(Ac,Bc,gd,type,with_V);
    };
};

template<class V1, class S1, class S2>
struct gschur_sym_impl<V1,V1,S1,S2>
{
    using Mat1  = raw::Matrix<V1,S1>;
    using Mat2  = raw::Matrix<V1,S2>;

    static void eval(const Mat1& A, const Mat2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        return gschur_sym_impl2<V1,S1,S2>::eval(A,B,gd,type,with_V);
    };
};

template<class S1, class S2>
struct gschur_sym_impl<Integer,Integer,S1,S2>
{
    using Mat1  = raw::Matrix<Integer,S1>;
    using Mat2  = raw::Matrix<Integer,S2>;

    using VR        = Real;
    using SR        = typename md::unify_struct<S1,S2>::type;
    using Mat_R     = raw::Matrix<VR,SR>;

    static void eval(const Mat1& A, const Mat2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        Mat_R Ac    = raw::converter<Mat_R, Mat1>::eval(A);
        Mat_R Bc    = raw::converter<Mat_R, Mat2>::eval(B);

        return gschur_sym_impl<Real,Real,SR,SR>::eval(Ac,Bc,gd,type,with_V);
    };
};
template<class S1, class S2>
struct gschur_sym_impl<Object,Object,S1,S2>
{
    using Mat1  = raw::Matrix<Object,S1>;
    using Mat2  = raw::Matrix<Object,S2>;

    static void eval(const Mat1&, const Mat2&, gschur_sym_decomposition&, gschur_sym_type, bool)
    {
        throw error::object_value_type_not_allowed("gschur_sym_decomposition");
    };
};

struct gschur_sym_vis : public extract_type2_switch<void, gschur_sym_vis, mr::val_type_corrector_diag>
{
    template<class T1, class T2>
    static void eval_mat_mat(const T1& A, const T2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        using V1    = typename T1::value_type;
        using S1    = typename T1::struct_type;

        using V2    = typename T2::value_type;
        using S2    = typename T2::struct_type;

        return details::gschur_sym_impl<V1,V2,S1,S2>::eval(A,B,gd,type,with_V);
    };

    template<class T1, class T2>
    static void eval_scal_scal(const T1& A, const T2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        using Dense_2 = raw::Matrix<T2,struct_dense>;

        Dense_1 Ac(ti::get_ti(A),A,1,1);
        Dense_2 Bc(ti::get_ti(B),B,1,1);

        return eval_mat_mat<Dense_1,Dense_2>(Ac,Bc, gd, type, with_V);
    };

    template<class T1, class T2>
    static void eval_mat_scal(const T1& A, const T2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        using Dense_2 = raw::Matrix<T2,struct_dense>;
        Dense_2 Bc(ti::get_ti(B),B,1,1);
        return eval_mat_mat<T1,Dense_2>(A,Bc, gd, type, with_V);
    };

    template<class T1, class T2>
    static void eval_scal_mat(const T1& A, const T2& B, gschur_sym_decomposition& gd, gschur_sym_type type, bool with_V)
    {
        using Dense_1 = raw::Matrix<T1,struct_dense>;
        Dense_1 Ac(ti::get_ti(A),A,1,1);
        return eval_mat_mat<Dense_1,T2>(Ac,B, gd, type, with_V);
    };
};

};};

namespace matcl
{

gschur_sym_decomposition::gschur_sym_decomposition()
{
    clear();
};

gschur_sym_decomposition::~gschur_sym_decomposition()
{};

void gschur_sym_decomposition::clear()
{
    m_V         = zeros(0,0);
    m_D         = m_V;
    m_eig       = zeros(0,1);
    m_has_Q     = false;
    m_is_nan    = false;
};

gschur_sym_decomposition::gschur_sym_decomposition(const Matrix &A0, const Matrix &B0, gschur_sym_type type, 
                                                   bool with_V)
{
    Matrix A(A0);
    Matrix B(B0);

    clear();
    compute(A, B, type, with_V);
};
gschur_sym_decomposition::gschur_sym_decomposition(Matrix &&A0, const Matrix &B0, gschur_sym_type type, bool with_V)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    clear();
    compute(A, B, type, with_V);
};
gschur_sym_decomposition::gschur_sym_decomposition(const Matrix &A0, Matrix &&B0, gschur_sym_type type, bool with_V)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    clear();
    compute(A, B, type, with_V);
};
gschur_sym_decomposition::gschur_sym_decomposition(Matrix &&A0, Matrix && B0, gschur_sym_type type, bool with_V)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    clear();
    compute(A, B, type, with_V);
};

gschur_sym_decomposition& 
gschur_sym_decomposition::operator()(const Matrix &A0, const Matrix &B0, gschur_sym_type type, bool with_V)
{
    Matrix A(A0);
    Matrix B(B0);

    clear();
    compute(A, B, type, with_V);
    return *this;
};
gschur_sym_decomposition& 
gschur_sym_decomposition::operator()(Matrix &&A0, const Matrix &B0, gschur_sym_type type, bool with_V)
{
    Matrix A(std::move(A0));
    Matrix B(B0);

    clear();
    compute(A, B, type, with_V);
    return *this;
};
gschur_sym_decomposition& 
gschur_sym_decomposition::operator()(const Matrix &A0, Matrix &&B0, gschur_sym_type type, bool with_V)
{
    Matrix A(A0);
    Matrix B(std::move(B0));

    clear();
    compute(A, B, type, with_V);
    return *this;
};
gschur_sym_decomposition& 
gschur_sym_decomposition::operator()(Matrix && A0, Matrix && B0, gschur_sym_type type, bool with_V)
{
    Matrix A(std::move(A0));
    Matrix B(std::move(B0));

    clear();
    compute(A, B, type, with_V);
    return *this;
};

void gschur_sym_decomposition::compute(const Matrix &A, const Matrix &B, gschur_sym_type type, bool with_V)
{
    if (!A.is_square() || !B.is_square() || A.rows() != B.rows())
        throw error::error_size_geig(A.rows(),A.cols(), B.rows(), B.cols());

    value_code vt_A     = A.get_value_code();
    value_code vt_B     = B.get_value_code();
    value_code vt0      = matrix_traits::unify_value_types(vt_A, vt_B);
    value_code vt       = matrix_traits::unify_value_types(vt0, value_code::v_float);

    bool is_complex_A   = matrix_traits::is_float_complex(vt_A);
    bool is_complex_B   = matrix_traits::is_float_complex(vt_B);

    bool is_obj         = (vt == value_code::v_object);

    m_has_Q             = with_V;

    if (is_obj == true)
        throw error::object_value_type_not_allowed("gschur_sym");

    Integer N           = A.rows();

    // FAST EXITS
    if (N == 0)
    {
        m_V             = speye(0,0, vt);
        m_D             = speye(0,0, vt);
        m_eig           = m_D;
        return;
    }
    if(B.structural_nnz() == 0)
    {
        throw error::error_nonposdef(false);
    }
    
    bool isv            = A.all_finite() && B.all_finite();
    
    if (isv == false)
    {
        m_V         = details::make_nan_matrix(N, N, vt);
        m_D         = m_V;
        m_eig       = details::make_nan_matrix(N, 1, vt);
        m_is_nan    = true;
        return;
    }
    else
    {
        m_is_nan    = false;
    };

    // both A and B trivial:
    if (N == 1)
    {
        if (B > 0)
        {
            Matrix D, V;

            Matrix one      = ones(1,1,vt);

            Matrix Ac       = matcl::convert_value(A, vt);
            Matrix Bc       = matcl::convert_value(B, vt);

            Matrix sq_B     = sqrt(Bc);

            switch (type)
            {
                case gschur_sym_type::A_B:
                {
                    D   = div(Ac, Bc);
                    V   = div(one, sq_B);
                    break;
                }
                case gschur_sym_type::AB:
                {
                    D   = Ac * Bc;
                    V   = div(one, sq_B);
                    break;
                }
                case gschur_sym_type::BA:
                {                    
                    D   = Ac * Bc;
                    V   = sq_B;
                    break;
                };
            }

            m_D     = D;
            m_V     = V;
            m_eig   = D;
            return;
        }
        else
        {
            throw error::error_nonposdef(false);
        }
    }

    // B trivial/easy
    if (B.get_struct().is_id())
    {
        schur_decomposition schur_obj;
        if (is_complex_B == true)
        {
            mat_code mc = matrix_traits::get_matrix_type(vt, A.get_struct_code());
            Matrix A2   = convert(A,mc);
            A2.add_struct(predefined_struct_type::her);
            schur_obj(std::move(A2));
        }
        else if (A.get_struct().is_hermitian(true, is_complex_A == false))
        {
            schur_obj(A);
        }
        else if (A.is_unique())
        {
            A.add_struct(predefined_struct_type::her);
            schur_obj(A);
        }
        else
        {
            Matrix A2 = A;
            A2.make_unique();
            A2.add_struct(predefined_struct_type::her);
            schur_obj(std::move(A));
        };

        m_V     = schur_obj.U();
        m_D     = schur_obj.TA();
        m_eig   = schur_obj.eig();
        return;
    }

    if(is_diag(B) && is_diag(A))
    {
        Matrix diag_A   = get_diag(A);
        Matrix diag_B   = get_diag(B);

        mat_code mc     = matrix_traits::get_matrix_type(vt, struct_code::struct_dense);
        diag_A          = convert(diag_A, mc);
        diag_B          = convert(diag_B, mc);

        if (any_vec(diag_B <= 0))
            throw error::error_nonposdef(false);

        Matrix D;
        Matrix V = ones(1,1,vt);

        switch (type)
        {
            case gschur_sym_type::A_B:
            {
                D       = div(diag_A, diag_B);

                if (with_V)
                    V   = div(V,sqrt(diag_B));

                break;
            }
            case gschur_sym_type::AB:
            {
                D       = mul(diag_A, diag_B);

                if (with_V)
                    V   = div(V,sqrt(diag_B));
                break;
            }
            case gschur_sym_type::BA:
            {
                D       = mul(diag_B, diag_A);

                if (with_V)
                    V   = sqrt(diag_B);
                break;
            }
        }

        m_eig   = D;
        m_D     = bdiag(D);

        if (with_V)
            m_V = bdiag(V);

        return;
    }

    return details::gschur_sym_vis::make(A,B,*this,type,with_V);
};

bool gschur_sym_decomposition::test_factors()
{    
    bool isv_V      = (m_has_Q == true) && m_V.all_finite() || m_has_Q == false;

    bool isv        = m_eig.all_finite() && isv_V;

    Integer N       = m_D.rows();

    if (isv == false)
    {
        value_code vt   = m_D.get_value_code();

        m_D             = details::make_nan_matrix(N, N, vt);
        m_eig           = details::make_nan_matrix(N, 1, vt);

        if (m_has_Q)
            m_V         = m_D;

        m_is_nan        = true;

        return false;
    }
    else
    {
        m_is_nan        = false;
    };
    
    return true;
};

Matrix gschur_sym_decomposition::V() const
{
    if (m_has_Q)
        return m_V;
    else
        throw error::error_schur_U_not_computed();
};
Matrix gschur_sym_decomposition::D() const
{
    return m_D;
};
Matrix gschur_sym_decomposition::eig() const
{
    return m_eig;
};

};

#pragma warning( pop )