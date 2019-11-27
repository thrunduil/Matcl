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
#include "matcl-linalg/linear_eq/linsolve_objects_decomp.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/decompositions/balancing.h"

namespace matcl
{

namespace details
{

template<class V, class S>
struct lapack_xgeqrf_maker{};

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V> 
struct lapack_xgeqrf_maker<V,struct_dense>
{
    using Mat = raw::Matrix<V,struct_dense>;

    static bool can_do_banded_alg(const Mat& Ac, Integer ld, Integer ud)
    {
        Integer dominant_dim    = std::max(Ac.rows(), Ac.cols());
        Integer bandwidth       = ld + ud;
        
        return bandwidth * linalg_optim_params::qrband_band_crossover() < dominant_dim;
    }

    static void eval(Mat& Ac, Mat& tau, Integer& ret_ld)
    {
        Integer ld, ud;
        {
            Matrix tmp(Ac,false);
            ld                  = matcl::get_ld(tmp);
            ud                  = matcl::get_ud(tmp);
        };
        ret_ld = ld;

        if (can_do_banded_alg(Ac,ld,ud) == true)
            eval_banded_in_dense_rep(Ac, tau, ld, ud);
        else
            eval_dense_in_dense_rep(Ac, tau);
    };

    static void eval_dense_in_dense_rep(Mat& Ac, Mat& tau)
    {
        Integer M   = Ac.rows();
        Integer N   = Ac.cols();

        Integer info1 = 0;
        Mat work1(ti::ti_real(), 1, 1);
        
        lapack::geqrf(M, N, lap(Ac.ptr()), Ac.ld(), lap(tau.ptr()),
                        lap(work1.ptr()), -1, lap(&info1));

        Integer lwork1 = (Integer) real(work1.ptr()[0]);
        work1.reset_unique(lwork1, 1);
        
        lapack::geqrf(M, N, lap(Ac.ptr()), Ac.ld(), lap(tau.ptr()), 
                        lap(work1.ptr()), lwork1, lap(&info1));
        
        if (info1 != 0)
        {
            std::ostringstream msg;
            msg << "geqrf returned info = " << info1;
            throw error::error_general(msg.str());
        };
    }

    static void eval_banded_in_dense_rep(Mat& Ac, Mat& tau, Integer ld, Integer ud)
    {
        Integer info1 = 0;

        Integer M   = Ac.rows();
        Integer N   = Ac.cols();
        Integer ML  = ld;
        Integer MU  = ud;

        lapack::i_type acld = Ac.ld();

        if (ML + MU > linalg_optim_params::qrband_block_crossover())
        {
            //block version
            Integer block_size = linalg_optim_params::qrband_block_size(); 
            
            // workspace based on documentation of gebqrf
            Integer work_size = block_size * block_size + 
                                std::min(N, ML + MU) * block_size + 
                                std::min(M, ML + block_size) * block_size + 1; 

            Mat work1(ti::ti_real(), work_size, 1);
            
            lapack::gebqrf(block_size, M, N, ML, MU, lap(Ac.ptr()), acld, lap(tau.ptr()),
                           lap(work1.ptr()), info1);
        }
        else
        {
            //non block version
            Mat work1(ti::ti_real(), min(N,MU+ML), 1);

            lapack::gebqr2(M, N, ML, MU, lap(Ac.ptr()), acld, lap(tau.ptr()), 
                           lap(work1.ptr()), info1);
        }

        if (info1 != 0)
        {
            std::ostringstream msg;
            msg << "gebqrf/gebqr2 returned info = " << info1;

            throw error::error_general(msg.str());
        }
    };
};

template<class V> 
struct geqp3type_selector
{};

template<> 
struct geqp3type_selector<Real>
{
    using Mat       = raw::Matrix<Real, struct_dense>;
    using Mat_int   = raw::Matrix<Integer, struct_dense>;

    static void eval(Mat& Ac, Mat_int& jpvt, Mat& tau, Mat& work, Integer lwork, 
                        Integer* info2)
    {
        Integer M = Ac.rows();
        Integer N = Ac.cols();

        lapack::dgeqp3(M, N, lap(Ac.ptr()), Ac.ld(), lap(jpvt.ptr()), lap(tau.ptr()), 
                        lap(work.ptr()), lwork, lap(info2));
    };
};
template<> 
struct geqp3type_selector<Float>
{
    using Mat       = raw::Matrix<Float, struct_dense>;
    using Mat_int   = raw::Matrix<Integer, struct_dense>;

    static void eval(Mat& Ac, Mat_int& jpvt, Mat& tau, Mat& work, Integer lwork, 
                        Integer* info2)
    {
        Integer M = Ac.rows();
        Integer N = Ac.cols();

        lapack::sgeqp3(M, N, lap(Ac.ptr()), Ac.ld(), lap(jpvt.ptr()), lap(tau.ptr()), 
                        lap(work.ptr()), lwork, lap(info2));
    };
};
template<>
struct geqp3type_selector<Complex>
{
    using Mat       = raw::Matrix<Complex, struct_dense>;
    using Mat_i     = raw::Matrix<Integer, struct_dense>;
    using rwork_Mat = raw::Matrix<Real, struct_dense>;

    static void eval(Mat& Ac, Mat_i& jpvt, Mat& tau, Mat& work, Integer lwork, 
                    Integer* info2)
    {
        Integer M = Ac.rows();
        Integer N = Ac.cols();

        rwork_Mat rwork(ti::ti_real(), 2*N, 1);

        lapack::zgeqp3(M, N, lap(Ac.ptr()), Ac.ld(), lap(jpvt.ptr()), lap(tau.ptr()), 
                        lap(work.ptr()), lwork, lap(rwork.ptr()), lap(info2));
    };
};
template<>
struct geqp3type_selector<Float_complex>
{
    using Mat       = raw::Matrix<Float_complex, struct_dense>;
    using Mat_i     = raw::Matrix<Integer, struct_dense>;
    using rwork_Mat = raw::Matrix<Float, struct_dense>;

    static void eval(Mat& Ac, Mat_i& jpvt, Mat& tau, Mat& work, Integer lwork, 
                    Integer* info2)
    {
        Integer M = Ac.rows();
        Integer N = Ac.cols();

        rwork_Mat rwork(ti::ti_real(), 2*N, 1);

        lapack::cgeqp3(M, N, lap(Ac.ptr()), Ac.ld(), lap(jpvt.ptr()), lap(tau.ptr()), 
                        lap(work.ptr()), lwork, lap(rwork.ptr()), lap(info2));
    };
};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class Val> 
struct lapack_xgeqrf_maker<Val, struct_banded>
{
    using Mat       = raw::Matrix<Val, struct_banded>;
    using dense_Mat = raw::Matrix<Val, struct_dense>;

    static void eval(Mat& Ac_band, dense_Mat& tau)
    {
        matcl_assert(Ac_band.has_diag(0) == true, "this case should already be processed");

        Integer M       = Ac_band.rows();
        Integer N       = Ac_band.cols();
        Integer ML      = Ac_band.number_subdiagonals();
        Integer MU      = Ac_band.number_superdiagonals();
        Integer info1   = 0;
        
        if (ML + MU > linalg_optim_params::qrband_block_crossover())
        {
            Integer block_size = linalg_optim_params::qrband_block_size();

            // workspace based on documentation of gbbqrf
            Integer work_size = block_size * block_size + 
                                min(N, ML + MU) * block_size + 
                                min(M, ML + block_size) * block_size + 1; 

            dense_Mat work1(ti::ti_real(), work_size, 1);
            
            // expand Ac_band representation upwars (more superdiagonals) and
            // rightwise (more columns) to accomodate matrix R
            Mat Ac          = Ac_band.resize(M, max(N, ML + MU + block_size), -ML, 
                                            ML + MU + block_size - 1);
            lapack::i_type acld = Ac.ld();

            lapack::gbbqrf(block_size, M, N, ML, MU, lap(Ac.rep_ptr()), acld, 
                           lap(tau.ptr()), lap(work1.ptr()), info1);

            //Ac_band is allocated on stack, so is fresh
            Ac_band.assign_to_fresh(Ac.resize(M, N, -ML, ML + MU + block_size + 1));
        }
        else 
        {
            // expand Ac_band representation upwars (more superdiagonals) and
            // rightwise (more columns) to accomodate matrix R
            Mat Ac = Ac_band.resize(M, max(N, ML + MU + 1), -ML, ML + MU);            

            lapack::i_type acld = Ac.ld();
            
            dense_Mat work1(ti::ti_real(), min(N,MU+ML), 1);
            
            lapack::gbbqr2(M, N, ML, MU, lap(Ac.rep_ptr()), acld, lap(tau.ptr()), 
                           lap(work1.ptr()), info1);

            //Ac_band is allocated on stack, so is fresh
            Ac_band.assign_to_fresh(Ac.resize(M, N, -ML, MU + ML));
        }

        if (info1 != 0)
        {
            std::ostringstream msg;
            msg << "dgbbqrf/dgbbqr2 returned info = " << info1;

            throw error::error_general(msg.str());
        };
    }
};

template<class V, class S> 
struct qr_str{};

//--------------------------------------------------------------------------------
//                                  DENSE
//--------------------------------------------------------------------------------
template<class V> 
struct qr_str<V,struct_dense>
{
    using Mat   = raw::Matrix<V, struct_dense>;
    using Mat_i = raw::Matrix<Integer, struct_dense>;

    static void eval(Matrix& ret, const Mat& A)
    {
        Matrix tau;
        return eval2(ret,tau,A);
    };

    static void eval2(Matrix& ret, Matrix& Tau, const Mat& A)
    {
        Integer M = A.rows();
        Integer N = A.cols();
        Integer K = min(M, N);

        Mat tau(A.get_type(), K, 1);

        Mat Ac = A.make_unique();

        Integer ld;
        lapack_xgeqrf_maker<V, struct_dense>::eval(Ac, tau, ld);

        Ac.set_struct(struct_flag());
        ret = Matrix (Ac, true);
        Tau = Matrix (tau, false);
    };
    static void eval(const Mat& A, unitary_matrix& ret_Q, Matrix& ret_R, const bool economy)
    {
        Integer M       = A.rows();
        Integer N       = A.cols();
        Integer M_eco   = economy ? min(M,N) : M;
        Integer K       = min(M_eco, N);

        Mat tau(A.get_type(),K,1);

        Mat Ac(A.get_type());

        if (A.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(A);
        else
            Ac.assign_to_fresh(A.copy());
        
        Ac.set_struct(struct_flag());

        Integer ld = 0;
        lapack_xgeqrf_maker<V, struct_dense>::eval(Ac, tau, ld);

        Mat Qc          = Ac.copy();
        Qc.assign_to_fresh(Qc.resize(M, M_eco));

        using householder_impl = std::shared_ptr<householder_q<V>>;

        bool isv        = Qc.all_finite() && tau.all_finite();

        if (isv == false)
        {
            ret_Q   = unitary_matrix::from_nan(Qc.rows(), M_eco, matrix_traits::value_code<V>::value);
        }
        else
        {
            householder_impl Q_impl(new details::householder_q<V>(M_eco, Qc, tau, ld, 0));

            ret_Q   = unitary_matrix(Q_impl);
        };

        Ac.assign_to_fresh(Ac.resize(M_eco, N));
        Ac.set_struct(struct_flag());

        matcl::raw::inplace::make_triu<Mat>::eval(ret_R,Ac,0);
    };

    static void eval(const Mat& A, unitary_matrix& ret_Q, Matrix& ret_R,permvec& ret_E, 
                        bool economy)
    {
        Integer M       = A.rows();
        Integer N       = A.cols();
        Integer M_eco   = economy ? min(M, N) : M;
        Integer K       = min(M_eco, N);

        Mat tau(A.get_type(), K, 1);

        Mat_i jpvt(ti::ti_int(), 0, N, 1);

        Mat Ac(A.get_type());

        if (A.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(A);
        else
            Ac.assign_to_fresh(A.copy());

        {
            Integer info1 = 0;
            Mat work1(A.get_type(), 1, 1);

            geqp3type_selector<V>::eval(Ac, jpvt, tau, work1, -1, &info1);
            
            Integer lwork1 = (Integer) real(work1.ptr()[0]);
            work1.reset_unique(lwork1, 1);
            
            geqp3type_selector<V>::eval(Ac, jpvt, tau, work1, lwork1, &info1);
            
            if (info1 != 0)
            {
                std::ostringstream msg;
                msg << "geqp3 returned info = " << info1;
                throw error::error_general(msg.str());
            };
        }   

        Mat Qc  = Ac.copy();
        Qc.assign_to_fresh(Qc.resize(M, M_eco));

        bool isv    = Qc.all_finite() && tau.all_finite();

        if (isv == false)
        {
            ret_Q   = unitary_matrix::from_nan(Qc.rows(), M_eco, matrix_traits::value_code<V>::value);
        }
        else
        {
            using householder_impl = std::shared_ptr<householder_q<V>>;
            householder_impl Q_impl(new details::householder_q<V>(M_eco, Qc, tau, M, 0));

            ret_Q   = unitary_matrix(Q_impl);
        };

        Ac.assign_to_fresh(Ac.resize(M_eco, N));
        Ac.set_struct(struct_flag());

        matcl::raw::inplace::make_triu<Mat>::eval(ret_R,Ac,0);

        permvec p = permvec::from_matrix(Matrix(jpvt,false));
        ret_E = p;
    };
};

//--------------------------------------------------------------------------------
//                                  BAND
//--------------------------------------------------------------------------------
template<class V> 
struct qr_str<V,struct_banded>
{
    using Mat       = raw::Matrix<V, struct_banded>;
    using Mat_dense = raw::Matrix<V, struct_dense>;
    using Mat_i     = raw::Matrix<Integer, struct_dense>;

    static void eval(Matrix& ret, const Mat& A)
    {
        Matrix tau;
        return eval2(ret,tau,A);
    };
    static void eval2(Matrix& ret, Matrix& Tau, const Mat& A)
    {
        Integer M = A.rows();
        Integer N = A.cols();
        Integer K = min(M, N);

        Mat_dense tau(A.get_type(), K, 1);

        Mat Ac = A.make_unique();

        lapack_xgeqrf_maker<V, struct_banded>::eval(Ac, tau);

        Ac.set_struct(struct_flag());
        ret = Matrix (Ac, true);
        Tau = Matrix(tau,false);
    };
    static void eval(const Mat& A, unitary_matrix& ret_Q, Matrix& ret_R, const bool economy)
    {
        Integer M       = A.rows();
        Integer N       = A.cols();
        Integer M_eco   = economy ? min(M,N) : M;
        Integer K       = min(M_eco, N);

        Mat_dense tau(A.get_type(),K,1);

        Mat Ac(A.get_type());

        if (A.get_refstr()->is_unique() == true)
            Ac.assign_to_fresh(A);
        else
            Ac.assign_to_fresh(A.copy());
        
        Ac.set_struct(struct_flag());

        lapack_xgeqrf_maker<V, struct_banded>::eval(Ac, tau);
        
        Mat Qc  = Ac.copy();
        Qc.assign_to_fresh(Qc.resize(M, M_eco));

        bool isv    = Qc.all_finite();

        if (isv == false)
        {
            ret_Q   = unitary_matrix::from_nan(Qc.rows(), M_eco, matrix_traits::value_code<V>::value);
        }
        else
        {
            using householder_impl = std::shared_ptr<householder_band_q<V>>;
            householder_impl Q_impl(new details::householder_band_q<V>(M_eco, Qc, tau));

            ret_Q   = unitary_matrix(Q_impl);
        };

        Ac.assign_to_fresh(Ac.resize(M_eco, N));
        Ac.set_struct(struct_flag());

        matcl::raw::inplace::make_triu<Mat>::eval(ret_R,Ac,0);
    };
    static void eval(const Mat& A, unitary_matrix& ret_Q, Matrix& ret_R, permvec& ret_E, 
                        bool economy)
    {
        Mat_dense Ac = raw::converter<Mat_dense,Mat>::eval(A);
        return qr_str<V,struct_dense>::eval(Ac,ret_Q,ret_R,ret_E,economy);
    };

};

//--------------------------------------------------------------------------------
//                                  SPARSE
//--------------------------------------------------------------------------------
template<class V> 
struct qr_str<V,struct_sparse>
{
    using S = struct_sparse;

    using Mat       = raw::Matrix<V, S>;
    using dense_Mat = raw::Matrix<V, struct_dense>;
    using Mat_i     = raw::Matrix<Integer, struct_dense>;
    
    static void eval(Matrix& ret, const Mat& A)
    {
        quern_solver qs(matcl::Matrix(A,false),false,false,false);
        ret = qs.get_r();
        return;
    };
    static void eval2(Matrix& ret, Matrix& Tau, const Mat& A)
    {
        dense_Mat Ac = raw::converter<dense_Mat,Mat>::eval(A);
        return qr_str<V,struct_dense>::eval2(ret,Tau,Ac);
    };

    static void eval(const Mat& A, unitary_matrix& ret_Q, Matrix& ret_R, bool economy)
    {
        quern_solver qs(matcl::Matrix(A,false),true,economy,false);
        ret_R = qs.get_r();
        ret_Q = qs.get_unitary_matrix();
        return;
    };
    static void eval(const Mat& A, unitary_matrix& ret_Q, Matrix& ret_R, permvec& ret_E, 
                        bool economy)
    {
        dense_Mat Ac = raw::converter<dense_Mat,Mat>::eval(A);
        return qr_str<V,struct_dense>::eval(Ac,ret_Q,ret_R,ret_E,economy);
    };
};

//--------------------------------------------------------------------------------
//                                  VALUE TYPE
//--------------------------------------------------------------------------------
template<class V, class S>
struct qr_impl
{
    using M = raw::Matrix<V, S>;

    static void eval(Matrix& ret, const M& A)
    {
        return qr_str<V,S>::eval(ret, A);
    };
    static void eval2(Matrix& ret, Matrix& tau, const M& A)
    {
        return qr_str<V,S>::eval2(ret,tau, A);
    };
    static void eval(const M& A, unitary_matrix& ret_Q, Matrix& ret_R, bool economy)
    {
        qr_str<V,S>::eval(A, ret_Q, ret_R, economy);
    };
    static void eval(const M& A, unitary_matrix& ret_Q, Matrix& ret_R, permvec& ret_E, 
                    bool economy)
    {
        qr_str<V,S>::eval(A, ret_Q, ret_R, ret_E, economy);
    };
};

template<class S>
struct qr_impl<Integer, S>
{
    using M = raw::Matrix<Integer, S>;

    static void eval(Matrix& ret, const M& A)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return qr_impl<Real,S>::eval(ret, AC);
    };
    static void eval2(Matrix& ret, Matrix& tau, const M& A)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return qr_impl<Real,S>::eval2(ret, tau, AC);
    };

    static void eval(const M& A, unitary_matrix& ret_Q, Matrix& ret_R, const bool economy)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        qr_impl<Real, S>::eval(AC, ret_Q, ret_R, economy);
    };

    static void eval(const M& A, unitary_matrix& ret_Q, Matrix& ret_R, permvec& ret_E, 
                        const bool economy)
    {
        using MC = raw::Matrix<Real,S>;
        MC AC = raw::converter<MC,M>::eval(A);
        return qr_impl<Real,S>::eval(AC, ret_Q, ret_R, ret_E, economy);
    };
};

template<class S>
struct qr_impl<Object, S>
{
    using M = raw::Matrix<Object, S>;
    static void eval(Matrix&, const M&)
    {
        throw error::object_value_type_not_allowed("qr");
    };

    static void eval2(Matrix&, Matrix&, const M&)
    {
        throw error::object_value_type_not_allowed("qr");
    };
    static void eval(const M& , Matrix& , Matrix&, bool)
    {
        throw error::object_value_type_not_allowed("qr");
    };
    static void eval(const M& , Matrix& , Matrix& , permvec&, bool)
    {
        throw error::object_value_type_not_allowed("qr");
    };
};

struct unary_visitor_qr1 : public extract_type_switch<void,unary_visitor_qr1,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename md::real_type_int_real<V>::type;

        bool isv    = mat.all_finite();

        if (isv == false)
        {
            ret = details::make_nan_matrix<VR>(mat.rows(), mat.cols());
            return;
        }

        return qr_impl<V, S>::eval(ret, mat);
    };
    template<class T>
    static void eval(const Matrix&, const T& mat, Matrix& ret, Matrix& tau)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename md::real_type_int_real<V>::type;

        bool isv    = mat.all_finite();

        if (isv == false)
        {
            ret = details::make_nan_matrix<VR>(mat.rows(), mat.cols());
            return;
        }

        return qr_impl<V, S>::eval2(ret, tau, mat);
    };
    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, Matrix& ret, Matrix& tau)
    {
        static_assert(std::is_same<T,Object>::value == false, "object not allowed");
        using constants::nan;
    
        using VR = typename md::real_type_int_real<T>::type;

        bool isv = is_finite(v);

        if (isv == false)
        {
            ret = nan<VR>();
            tau = nan<VR>();
            return;
        }
        
        ret = handle;
        tau = VR(0.0);
        return;
    };
    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, Matrix& ret)
    {
        static_assert(std::is_same<T,Object>::value == false, "object not allowed");
        using constants::nan;
    
        using VR = typename md::real_type_int_real<T>::type;

        bool isv = is_finite(v);

        if (isv == false)
        {
            ret = nan<VR>();
            return;
        }
        
        ret = handle;
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&)
    {
        throw error::object_value_type_not_allowed("qr");
    };
    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, Matrix&, Matrix&)
    {
        throw error::object_value_type_not_allowed("qr");
    };

    static void eval_scalar(const Matrix&, const Object&, Matrix&)
    {
        throw error::object_value_type_not_allowed("qr");
    };
    static void eval_scalar(const Matrix&, const Object&, Matrix&, Matrix&)
    {
        throw error::object_value_type_not_allowed("qr");
    };
};

struct unary_visitor_qr2 : public extract_type_switch<void, unary_visitor_qr2,true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, qr2_return& ret, const bool economy)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename md::real_type_int_real<V>::type;

        bool isv    = mat.all_finite();        

        if (isv == false)
        {            
            Matrix ret_R;
            Integer N_Q, M_Q;

            if(economy == false)
            {
                M_Q     = mat.rows();
                N_Q     = mat.rows();

                ret_R = details::make_nan_matrix<VR>(mat.rows(), mat.cols());
            }
            else
            {
                Integer K   = min(mat.rows(), mat.cols());
                M_Q         = mat.rows();
                N_Q         = K;

                ret_R = details::make_nan_matrix<VR>(K, mat.cols());
            }

            ret = qr2_return(unitary_matrix::from_nan(M_Q, N_Q, ret_R.get_value_code()), ret_R);
            return;
        };

        unitary_matrix ret_Q;
        Matrix ret_R;

        qr_impl<V, S>::eval(mat, ret_Q, ret_R, economy);
            
        ret_R.add_struct(predefined_struct_type::triu);

        ret = qr2_return(ret_Q, ret_R);
        return;
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, qr2_return& ret, bool  )
    {
        static_assert(std::is_same<T,Object>::value == false, "object not allowed");
        using VR = typename md::real_type_int_real<T>::type;

        using constants::nan;
        
        bool isv = is_finite(mat);
        if (isv == false)
        {
            unitary_matrix Q = unitary_matrix::from_nan(1,1,matrix_traits::value_code<VR>::value);
            ret = qr2_return(Q, nan<VR>());
            return;
        }
        
        ret = qr2_return(unitary_matrix(VR(1.0),false),Matrix(handle)); // (q,r)
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, qr2_return&, bool)
    {
        throw error::object_value_type_not_allowed("qr2");
    };

    static void eval_scalar(const Matrix&, const Object&, qr2_return&, bool)
    {
        throw error::object_value_type_not_allowed("qr2");
    };
};

struct unary_visitor_qr3 : public extract_type_switch<void, unary_visitor_qr3, true>
{
    template<class T>
    static void eval(const Matrix&, const T& mat, qr3_return& ret, const bool economy)
    {
        using V     = typename T::value_type;
        using S     = typename T::struct_type;
        using VR    = typename md::real_type_int_real<V>::type;

        bool isv    = mat.all_finite();
        
        if (isv == false)
        {
            value_code vc   = matrix_traits::value_code<VR>::value;

            permvec p = permvec::identity(mat.cols());

            if (economy == false)
            {
                unitary_matrix Q = unitary_matrix::from_nan(mat.rows(), mat.rows(), vc);

                ret = qr3_return(Q,details::make_nan_matrix<VR>(mat.rows(), mat.cols()),p);
                return;
            }
            else
            {
                Integer K = min(mat.rows(), mat.cols());
                unitary_matrix Q = unitary_matrix::from_nan(mat.rows(), K, vc);

                ret = qr3_return(Q, details::make_nan_matrix<VR>(K, mat.cols()), p);
                return;
            }
        }
        
        unitary_matrix q;
        Matrix r;
        permvec p;
        
        details::qr_impl<V, S>::eval(mat, q, r, p, economy);       

        r.add_struct(predefined_struct_type::triu);
        
        ret = qr3_return(q, r, p) ;
    };

    template<class T>
    static void eval_scalar(const Matrix& handle, const T& v, qr3_return& ret, const bool  )
    {
        static_assert(std::is_same<T,Object>::value == false, "object not allowed");
        using VR = typename md::real_type_int_real<T>::type;

        using constants::nan;
        
        bool isv = is_finite(v);
        if (isv == false)
        {
            unitary_matrix Q = unitary_matrix::from_nan(1,1,matrix_traits::value_code<VR>::value);
            ret = qr3_return(Q, nan<VR>(), permvec::identity(1));
            return;
        }
        
        ret = qr3_return(unitary_matrix(VR(1.0),false), handle, permvec::identity(1));
        return;
    };

    template<class S>
    static void eval(const Matrix&, const raw::Matrix<Object,S>&, qr3_return&, bool)
    {
        throw error::object_value_type_not_allowed("qr3");
    };

    static void eval_scalar(const Matrix&, const Object&, qr3_return&, bool)
    {
        throw error::object_value_type_not_allowed("qr3");
    };
};

static void qr_internal_impl(Matrix& ret, const Matrix& A)
{
    matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
    matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

    if (A.structural_nnz() == 0 || is_triu(A))
    {
        A.add_struct(predefined_struct_type::triu);
        
        ret = A;
        return;
    }

    if (is_unitary(A.get_struct()) == true)
    {
        ret = speye(A.rows(), A.rows(), vt);
        return;
    }

    return details::unary_visitor_qr1::make<const Matrix&>(A, ret);
};
static void qr_internal2_impl(Matrix& ret, Matrix& tau, const Matrix& A)
{
    matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
    matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

    if (A.structural_nnz() == 0 || is_triu(A))
    {
        Integer K   = A.cols();

        A.add_struct(predefined_struct_type::triu);
        
        ret = A;
        tau = matcl::zeros(K,1,vt);
        return;
    }

    return details::unary_visitor_qr1::make<const Matrix&>(A, ret, tau);
};

static void qr2_impl(qr2_return& ret, const Matrix& A, bool economy)
{    
    matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
    matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);


    if (A.structural_nnz() == 0 || is_triu(A))
    {
        if (economy == false)
        {
            A.add_struct(predefined_struct_type::triu);

            Matrix Q = speye(A.rows(),A.rows(), vt);
            ret = qr2_return(unitary_matrix(Q,false), A);
            return;
        }
        else
        {
            Matrix temp = A(colon(1, min(A.rows(), A.cols())), colon());
            temp.add_struct(predefined_struct_type::triu);
    
            Matrix Q = speye(A.rows(), min(A.rows(), A.cols()), vt);
            ret = qr2_return(unitary_matrix(Q,false), temp);
            return;
        }
    };

    if (is_unitary(A.get_struct()) == true)
    {
        // A is square - no economy version & R is square
        ret = qr2_return(unitary_matrix(A,false), speye(A.rows(), A.rows(), vt));
        return;
    }

    return details::unary_visitor_qr2::make<const Matrix&>(A, ret, economy);
};

static void qr3_impl(qr3_return& ret, const Matrix& A, const bool economy)
{
    matcl::value_code vt0 = matrix_traits::real_value_type(A.get_value_code());
    matcl::value_code vt  = matrix_traits::unify_value_types(vt0, value_code::v_float);

    if (A.structural_nnz() == 0 || is_triu(A))
    {
        if(!economy)
        {
            A.add_struct(predefined_struct_type::triu);

            permvec p   = permvec::identity(A.cols());
            Matrix I    = speye(A.rows(), A.rows(), vt);
            ret         = qr3_return(unitary_matrix(I,false), A, p);
            return;
        }
        else
        {
            Matrix temp = A(colon(1, min(A.rows(), A.cols())), colon());
            temp.add_struct(predefined_struct_type::triu);

            permvec p   = permvec::identity(A.cols());
            Matrix I    = speye(A.rows(), min(A.rows(), A.cols()), vt);

            ret = qr3_return(unitary_matrix(I,false) , temp, p );
            return;
        }
    }

    if (is_unitary(A.get_struct()) == true)
    {
        // A is square - no economy version & R is square
        permvec p   = permvec::identity(A.cols());
        ret         = qr3_return(unitary_matrix(A,false), speye(A.rows(), A.rows(), vt), p);
        return;
    }

    return details::unary_visitor_qr3::make<const Matrix&>(A, ret, economy);
};

};

Matrix matcl::qr_internal(const Matrix& A0)
{
    //increase refcount
    Matrix A(A0);

    Matrix ret;
    details::qr_internal_impl(ret, A);
    return ret;
};
Matrix matcl::qr_internal(Matrix&& A0)
{
    Matrix A(std::move(A0));

    Matrix ret;
    details::qr_internal_impl(ret, A);
    return ret;
};

qr_internal2_return matcl::qr_internal2(const Matrix& A0)
{
    Matrix A(std::move(A0));

    Matrix ret, tau;
    details::qr_internal2_impl(ret, tau, A);
    return qr_internal2_return(ret,tau);
};
qr_internal2_return matcl::qr_internal2(Matrix&& A0)
{
    Matrix A(std::move(A0));

    Matrix ret, tau;
    details::qr_internal2_impl(ret,tau,A);
    return qr_internal2_return(ret,tau);
};

Matrix matcl::qr(const Matrix& A)
{
    Matrix ret = triu(qr_internal(A));
    return ret ;
};
Matrix matcl::qr(Matrix&& A)
{
    Matrix ret = triu(qr_internal(std::move(A)));
    return ret ;
};

qr2_return matcl::qr2(const Matrix& A0, bool economy)
{
    //increase refcount
    Matrix A(A0);

    qr2_return ret;
    details::qr2_impl(ret,A,economy);
    return ret;
};
qr2_return matcl::qr2(Matrix&& A0, bool economy)
{
    Matrix A(std::move(A0));

    qr2_return ret;
    details::qr2_impl(ret,A,economy);
    return ret;
};

qr3_return matcl::qr3(const Matrix& A0, bool economy)
{
    //increase refcount
    Matrix A(A0);

    qr3_return ret;
    details::qr3_impl(ret,A,economy);
    return ret;
};
qr3_return matcl::qr3(Matrix&& A0, bool economy)
{
    Matrix A(std::move(A0));

    qr3_return ret;
    details::qr3_impl(ret,A,economy);
    return ret;
};

linsolve_obj matcl::linsolve_qr(const Matrix& A, const Matrix& Q, const Matrix& R, const options& opts)
{
    permvec p = permvec::identity(R.cols());
    return linsolve_qr(A, unitary_matrix(Q,true), R, p, opts);
};
linsolve_obj matcl::linsolve_qr(const Matrix& A, const Matrix& Q, const Matrix& R, 
                                const permvec& p, const options& opts)
{
    return linsolve_qr(A, unitary_matrix(Q,true), R, p, opts);
};
linsolve_obj matcl::linsolve_qr(const Matrix& A, const unitary_matrix& Q, const Matrix& R, const options& opts)
{
    permvec p = permvec::identity(R.cols());
    return linsolve_qr(A, Q, R, p, opts);
};
linsolve_obj matcl::linsolve_qr(const Matrix& A, const unitary_matrix& Q, const Matrix& R, 
                                const permvec& p, const options& opts)
{
    if (A.rows() != Q.rows())
        throw error::invalid_qr_factors();
    if (Q.cols() != R.rows())
        throw error::invalid_qr_factors();

    if (Q.rows() < Q.cols())
        throw error::invalid_qr_factors();
    if (A.cols() != R.cols())
        throw error::invalid_qr_factors();

    if (p.length() != R.cols())
        throw error::invalid_qr_factors();

    if (matcl::is_triu(R) == false)
        throw error::invalid_qr_factors();

    if (R.rows() != R.cols())
        throw error::square_matrix_required(R.rows(), R.cols());

    if (Q.rows() != Q.cols())
        throw error::square_matrix_required(Q.rows(), Q.cols());

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_qr");
    if (Q.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_qr");
    if (R.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_qr");

    Integer N   = R.rows();

    if (N)
    {
        value_code vc   = matrix_traits::unify_value_types(Q.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(R.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, R.get_type())));
    };

    bool isv            = Q.all_finite() && R.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(Q.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(R.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, R.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_qr(A,Q,R,p, opts)));
}

linsolve_obj matcl::linsolve_ql(const Matrix& A, const Matrix& Q, const Matrix& L, const options& opts)
{
    return linsolve_ql(A, unitary_matrix(Q,true), L, opts);
};
linsolve_obj matcl::linsolve_ql(const Matrix& A, const unitary_matrix& Q, const Matrix& L, const options& opts)
{
    if (A.rows() != Q.rows())
        throw error::invalid_qr_factors();
    if (A.cols() != L.cols())
        throw error::invalid_qr_factors();
    if (Q.cols() != L.rows())
        throw error::invalid_qr_factors();

    if (Q.rows() < Q.cols())
        throw error::invalid_qr_factors();

    if (matcl::is_tril(L) == false)
        throw error::invalid_qr_factors();

    if (L.rows() != L.cols())
        throw error::square_matrix_required(L.rows(), L.cols());
    if (Q.rows() != Q.cols())
        throw error::square_matrix_required(Q.rows(), Q.cols());

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_ql");
    if (Q.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_ql");
    if (L.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_ql");

    permvec p = permvec::identity(L.cols());

    Integer N   = L.rows();

    if (N == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(Q.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(L.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, L.get_type())));
    };

    bool isv            = Q.all_finite() && L.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(Q.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(L.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, L.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_qr(A,Q,L,p, opts)));
};

linsolve_obj matcl::linsolve_rq(const Matrix& A, const Matrix& R, const Matrix& Q, const options& opts)
{
    return linsolve_rq(A, R, unitary_matrix(Q,true), opts);
};
linsolve_obj matcl::linsolve_rq(const Matrix& A, const Matrix& R, const unitary_matrix& Q, const options& opts)
{
    if (A.rows() != R.rows())
        throw error::invalid_qr_factors();
    if (A.cols() != Q.cols())
        throw error::invalid_qr_factors();
    if (R.cols() != Q.rows())
        throw error::invalid_qr_factors();

    if (Q.rows() < Q.cols())
        throw error::invalid_qr_factors();

    if (matcl::is_triu(R) == false)
        throw error::invalid_qr_factors();

    if (R.rows() != R.cols())
        throw error::square_matrix_required(R.rows(), R.cols());
    if (Q.rows() != Q.cols())
        throw error::square_matrix_required(Q.rows(), Q.cols());

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_rq");
    if (R.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_rq");
    if (Q.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_rq");

    Integer N   = R.rows();

    if (N == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(R.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(Q.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, R.get_type())));
    };

    bool isv            = R.all_finite() && Q.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(R.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(Q.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, R.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_rq(A, R, Q, opts)));
}

linsolve_obj matcl::linsolve_lq(const Matrix& A, const Matrix& L, const Matrix& Q, const options& opts)
{
    return linsolve_lq(A, L, unitary_matrix(Q,true), opts);
}
linsolve_obj matcl::linsolve_lq(const Matrix& A, const Matrix& L, const unitary_matrix& Q, const options& opts)
{
    if (A.rows() != L.rows())
        throw error::invalid_qr_factors();
    if (A.cols() != Q.cols())
        throw error::invalid_qr_factors();
    if (L.cols() != Q.rows())
        throw error::invalid_qr_factors();

    if (Q.rows() < Q.cols())
        throw error::invalid_qr_factors();

    if (matcl::is_tril(L) == false)
        throw error::invalid_qr_factors();

    if (L.rows() != L.cols())
        throw error::square_matrix_required(L.rows(), L.cols());
    if (Q.rows() != Q.cols())
        throw error::square_matrix_required(Q.rows(), Q.cols());

    if (A.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_lq");
    if (L.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_lq");
    if (Q.get_value_code() == value_code::v_object)
        throw error::object_value_type_not_allowed("linsolve_lq");

    Integer N   = L.rows();

    if (N == 0)
    {
        value_code vc   = matrix_traits::unify_value_types(Q.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(L.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_empty(vc, L.get_type())));
    };

    bool isv            = Q.all_finite() && L.all_finite() && A.all_finite();

    if (isv == false)
    {
        value_code vc   = matrix_traits::unify_value_types(Q.get_value_code(), value_code::v_float);
        vc              = matrix_traits::unify_value_types(L.get_value_code(), vc);
        vc              = matrix_traits::unify_value_types(A.get_value_code(), vc);

        using data_ptr = linsolve_obj::linsolve_data_ptr;
        return linsolve_obj(data_ptr(new details::linsolve_obj_nan(N, vc, L.get_type() )));
    };

    using data_ptr  = linsolve_obj::linsolve_data_ptr;
    return linsolve_obj(data_ptr(new details::linsolve_obj_rq(A, L, Q, opts)));
}

struct linsolve_qr_impl
{
    static linsolve_obj eval(const Matrix& A, bool pivot, const options& opts)
    {
        if (A.rows() != A.cols())
            throw error::square_matrix_required(A.rows(), A.cols());

        bool balance    = false;
        Real tol_sig    = 0.0;

        if ((pivot == false) && opts.get_option<bool>(opt::linsolve::do_balancing()))
        {
            balance     = true;
            tol_sig     = 0.0;
        }
        else if ((pivot == true) && opts.get_option<bool>(opt::linsolve::do_balancing_rr()))
        {
            balance     = true;

            tol_sig     = opts.get_option<Real>(opt::linsolve::tol_sing());
            if (tol_sig < 0.0)
            {
                tol_sig = - tol_sig;
            }
            else if (tol_sig > 0.0)
            {
                tol_sig = norm(A, -2.0) * pow(constants::eps(A.get_value_code()), tol_sig);
            };
        };

        if (pivot == true)
        {
            unitary_matrix Q;
            Matrix R;
            permvec p;

            if (balance == true)
            {
                Matrix B, Dl, Dr;
                tie(B,Dl,Dr)    = balance_gen2(A, true, tol_sig);
                tie(Q,R,p)      = qr3(B);

                linsolve_obj lo_B = linsolve_qr(B, Q,R,p, opts);
                return linsolve_balanced(A, Dl,lo_B,Dr);
            }
            else
            {
                tie(Q,R,p) = qr3(A);
                return linsolve_qr(A, Q,R,p, opts);
            };
        };

        const Real qr_givens_threashold   = details::linalg_optim_params::qr_givens_qr_threashold();

        Integer N       = A.rows();
        Integer max_ld  = Integer(N * qr_givens_threashold);

        Integer ld      = matcl::get_ld(A, max_ld);

        if (ld == 0)
            return linsolve_triang(A, opts);

        Integer ud      = matcl::get_ud(A, max_ld);

        if (ud == 0)
            return linsolve_triang(A, opts);

        bool dense      = (A.get_struct_code() == struct_code::struct_dense);
        bool select_qr  = ld >= max_ld && ud >= max_ld;
        bool select_qrg = ld < max_ld;

        Matrix B, Dl, Dr;
                
        if (balance == true)
            tie(B,Dl,Dr)    = balance_gen2(A, true, tol_sig);

        if (dense == false || select_qr)
        {
            unitary_matrix Q;
            Matrix R;

            if (balance == true)
            {
                tie(Q,R) = qr2(B);
                linsolve_obj lo_B = linsolve_qr(A,Q,R, opts);
                return linsolve_balanced(A, Dl, lo_B, Dr);
            }
            else
            {
                tie(Q,R) = qr2(A);
                return linsolve_qr(A,Q,R, opts);
            };
        }
        else if (select_qrg)
        {
            unitary_matrix Q;
            Matrix R;

            if (balance)
            {
                tie(Q,R) = qr_givens(B);
                linsolve_obj lo_B = linsolve_qr(A,Q,R, opts);
                return linsolve_balanced(A, Dl, lo_B, Dr);
            }
            else
            {
                tie(Q,R) = qr_givens(A);
                return linsolve_qr(A,Q,R, opts);
            };
        }
        else
        {
            unitary_matrix Q;
            Matrix L;

            if (balance)
            {
                tie(Q,L) = lq_givens(A);
                linsolve_obj lo_B = linsolve_lq(A,L,Q, opts);
                return linsolve_balanced(A, Dl, lo_B, Dr);
            }
            {
                tie(Q,L) = lq_givens(A);
                return linsolve_lq(A,L,Q, opts);
            }
        }
    }
};

linsolve_obj matcl::linsolve_qr(const Matrix& A, bool pivot, const options& opts)
{
    return linsolve_qr_impl::eval(A,pivot, opts);
}
linsolve_obj matcl::linsolve_qr(Matrix&& A, bool pivot, const options& opts)
{
    return linsolve_qr_impl::eval(std::move(A),pivot,opts);
}

};