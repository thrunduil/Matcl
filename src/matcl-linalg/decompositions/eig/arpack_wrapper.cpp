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

#include "arpack_wrapper.h"
#include "extern/arpack/arpack.h"
#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-linalg/decompositions/chol.h"
#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-matrep/matcl_matrep.h"

namespace matcl { namespace details
{

// arpack is not thread-safe
using mutex_type    = matcl::default_mutex;

mutex_type arpack_mutex;

//---------------------------------------------------------------------------------
//              arpack_select_eig_helper
//---------------------------------------------------------------------------------
// Helper class to translate enum indicating arpack algorithm to string
class arpack_select_eig_helper
{
    private:
        static std::string select_eig_string(const cluster_type which_in);

    public:
        arpack_select_eig_helper(cluster_type select_type = cluster_type::LM)
        : m_select_type(select_type), m_select_string(select_eig_string(select_type))
        {}
 
        cluster_type        get_select_type() const         { return m_select_type;}
        const std::string&  get_select_type_string() const  { return m_select_string;}

    private:
        cluster_type    m_select_type;
        std::string     m_select_string;
};

std::string arpack_select_eig_helper::select_eig_string(const cluster_type which_in)
{
    switch(which_in)
    {
        case cluster_type::LM: return "LM";
        case cluster_type::SM: return "SM";
        case cluster_type::LR: return "LR";
        case cluster_type::SR: return "SR";
        case cluster_type::LA: return "LA";
        case cluster_type::SA: return "SA";
        case cluster_type::BE: return "BE";
        case cluster_type::LI: return "LI";
        case cluster_type::SI: return "SI";

        default :
            matcl_assert(0,"unknown arpack_select");
            throw error::error_general("unknown arpack_select");
    }
};

//---------------------------------------------------------------------------------
//                          arpack_wrapper
//---------------------------------------------------------------------------------
template<class T>
arpack_wrapper<T>::arpack_wrapper(const linear_operator& A_in, arpack_output out, Integer max_iter,
                                  Integer n_arnoldi, Real tol_in)
    : m_lock(arpack_mutex), m_N(A_in.rows()), m_tol(TR(tol_in)), m_calculate_called(false)
    , m_max_iter(max_iter), m_output(out), m_generalized(false), m_A(A_in), m_converged(false)
{    
    m_kernel_type   = get_kernel_type(A_in);
    m_has_x0        = false;
    m_num_arnoldi   = n_arnoldi;
}

template<class T>
void arpack_wrapper<T>::set_cluster(cluster_type ec)
{
    m_cluster = ec;
};

template<class T>
void arpack_wrapper<T>::set_x0(const Matrix& x0)
{
    m_x0            = x0;
    m_has_x0        = true;
};

template<class T>
arpack_wrapper<T>::arpack_wrapper(const linear_operator& A_in, const linear_operator& B_in, 
                                  bool hermitian, arpack_output out, Integer max_iter, Integer n_arnoldi,
                                  Real tol_in)
    : m_lock(arpack_mutex), m_N(A_in.rows()), m_tol(TR(tol_in)), m_calculate_called(false), m_max_iter(max_iter), 
      m_output(out), m_generalized(true), m_A(A_in), m_B(B_in), m_converged(false)
{
    m_kernel_type   = get_kernel_type(hermitian);
    m_has_x0        = false;
    m_num_arnoldi   = n_arnoldi;
}

template<class T>
arpack_wrapper<T>::~arpack_wrapper()
{};

template<class T>
void arpack_wrapper<T>::calculate(Integer k, bool return_nonconverged)
{
    prepare_dimmensions(k);
    prepare_integer_arrays();
    prepare_arrays();

    Mat_D v     = m_v.get_impl_unique<Mat_D>();
    m_ldv       = v.ld();

    Integer xyaupd_info = do_arpack_loop(v);

    if (xyaupd_info != 0)
    {
        if (xyaupd_info < 0)
        {
            throw error::error_general("invalid parameter passed to yaupd");
        }
        else if (xyaupd_info ==1)
        {
            std::string msg = "maximum number of iterations taken; all possible eigenvalues has been found";
            error::get_global_messanger_linalg()->warning_arpack(msg);
        }
        else if (xyaupd_info == 3)
        {
            std::string msg = "no shifts could be applied during a cycle of the Implicitly restarted Arnoldi iteration; "
                              "one possibility is to increase the number of Arnoldi vectors";
            error::get_global_messanger_linalg()->warning_arpack(msg);
        }
        else
        {
            std::string msg = "unknown error";
            error::get_global_messanger_linalg()->warning_arpack(msg);
        }
    };

    process_xyeupd(xyaupd_info, return_nonconverged, v);

    m_calculate_called = true;
}

template<class T>
void arpack_wrapper<T>::prepare_dimmensions(Integer k)
{
    m_num_eig   = k;
    // using Remark 4 dnaupd, znaupd, dsaupd...
    // NOTE: actually use some possibly overshot value of ncv - to be safer

    if (m_num_arnoldi < 0)
        m_num_vec = std::min(2 * k + 3, m_N);
    else
        m_num_vec = std::min(std::max(k + 3, m_num_arnoldi), m_N);
}

template<class T>
void arpack_wrapper<T>::prepare_integer_arrays()
{
    m_iparam            = make_integer_dense_noinit(11, 1, m_iparam_ptr);
    m_iparam_ptr[0]     = 1; // exact shifts with respect to the reduced tridiagonal matrix T
    //m_iparam_ptr[1]   not used
    m_iparam_ptr[2]     = m_max_iter; // maxiter
    m_iparam_ptr[3]     = 1; // nb, only 1 supported by arpack
    m_iparam_ptr[4]     = 0; // to be filled by arpack
    //m_iparam_ptr[5]   not used

    m_iparam_ptr[6]     = m_generalized == false ? 1 : 2;
    //m_iparam_ptr[7]   not used (by us)
    m_iparam_ptr[8]     = 0; // to be filled by arpack
    m_iparam_ptr[9]     = 0; // to be filled by arpack
    m_iparam_ptr[10]    = 0; // to be filled by arpack

    Integer len         = m_kernel_type == kernel_type::symmetric ? 11 : 14;
    m_ipntr             = make_integer_dense_noinit(len, 1, m_ipntr_ptr);
    m_select            = make_integer_dense_noinit(m_num_vec, 1, m_select_ptr);
}

template<class T>
void arpack_wrapper<T>::prepare_arrays()
{
    if (m_has_x0 == false)
    {
        m_x0            = make_dense_noinit(m_N, 1, get_value_code());
    }

    m_x0_ptr            = m_x0.get_array_unique<T>();

    // NOTE: m_num_vec + 1, not m_num_vec as required by dnaupd - HeapCorruption otherwise!
    // NOTE: in cases other than dnaupd -> this is just precaution
    m_v             = make_dense_noinit(m_N, m_num_vec + 1, get_value_code());

    // take care while handling this matrix to keep its refcount == 1
    m_workd         = make_dense_noinit(3 * m_N, 1, get_value_code());
    
    // taking care... watch out for workd_ptr
    m_workd_ptr     = m_workd.get_array_unique<T>();

    switch(m_kernel_type)
    {
        case kernel_type::symmetric:
            m_lworkl = m_num_vec * m_num_vec + 8 * m_num_vec;
            break;
        case kernel_type::nonsymmetric:
            m_lworkl = 3 * m_num_vec * m_num_vec + 6 * m_num_vec;
            break;
        case kernel_type::complex:
            m_lworkl = 3 * m_num_vec * m_num_vec + 5 * m_num_vec;
            break;
        default:
            matcl_assert(0, "unknown calculator_type in arpack_calculator");
            m_lworkl = 3 * m_num_vec * m_num_vec + 8 * m_num_vec; // safest value just in case...            
            break;
    };

    m_workl         = make_dense_noinit(m_lworkl, 1, get_value_code());
    m_workl_ptr     = m_workl.get_array_unique<T>();

    if (m_kernel_type == kernel_type::complex)
    {
        value_code vtr  = matrix_traits::value_code<TR>::value;
        m_rwork     = make_dense_noinit(m_num_vec,1,vtr);
        m_rwork_ptr = m_rwork.get_array_unique<TR>();
    }

    switch(m_kernel_type)
    {
        case kernel_type::nonsymmetric:
            m_workev        = make_dense_noinit(3 * m_num_vec, 1, get_value_code());
            m_workev_ptr    = m_workev.get_array_unique<T>();
            break;
        case kernel_type::complex:
            m_workev        = make_dense_noinit(2 * m_num_vec, 1, get_value_code());
            m_workev_ptr    = m_workev.get_array_unique<T>();
            break;
        case kernel_type::symmetric:
        default:
            break;
    }
}

template<class T>
Integer arpack_wrapper<T>::do_arpack_loop(Mat_D& v)
{
    Integer info        = m_has_x0? 1 : 0;
    Integer ido         = 0;
    Integer N           = m_N;
    bool do_next_loop   = true;

    while (do_next_loop)
    {
        call_xyaupd(ido, info, v);

        if (info < 0)
        {
            std::ostringstream msg;
            msg << "xyaupd returned code: " << info;
            throw error::error_arpack(msg.str());
        }

        if (ido == -1)
        {
            Integer from    = m_ipntr_ptr[0];
            Integer to      = m_ipntr_ptr[1];

            //eval Y = A * X 
            matcl::Matrix x     = make_dense_foreign(N, 1, m_workd_ptr + from - 1, N);
            matcl::Matrix y     = make_dense_foreign(N, 1, m_workd_ptr + to - 1, N);

            m_A.mmul_right(x, trans_type::no_trans, y);
        }
        else if (ido == 1)
        {
            Integer from    = m_ipntr_ptr[0];
            Integer to      = m_ipntr_ptr[1];

            //eval Y = A * X 
            matcl::Matrix x     = make_dense_foreign(N, 1, m_workd_ptr + from - 1, N);
            matcl::Matrix y     = make_dense_foreign(N, 1, m_workd_ptr + to - 1, N);

            m_A.mmul_right(x, trans_type::no_trans, y);
        }
        else if (ido == 2)
        {
            Integer from    = m_ipntr_ptr[0];
            Integer to      = m_ipntr_ptr[1];

            //eval Y = B * X 
            matcl::Matrix x     = make_dense_foreign(N, 1, m_workd_ptr + from - 1, N);
            matcl::Matrix y     = make_dense_foreign(N, 1, m_workd_ptr + to - 1, N);

            m_B.mmul_right(x, trans_type::no_trans, y);
        }
        else
        {
            do_next_loop = false;
        }
    }

    return info;
}

template<class T>
void arpack_wrapper<T>::process_xyeupd(Integer xyaupd_info, bool return_nonconverged, Mat_D& v)
{
    Integer info    = xyaupd_info;
    Matrix mat_eig;

    switch(m_kernel_type)
    {
        case kernel_type::symmetric:
            mat_eig         = make_dense_noinit(m_num_eig, 1, get_value_code());
            break;
        case kernel_type::nonsymmetric:
            mat_eig         = make_dense_noinit(m_num_eig + 1, 2, get_value_code());
            break;
        case kernel_type::complex:
            mat_eig         = make_dense_noinit(m_num_eig + 1, 1, get_value_code());
            break;
        default:
            throw error::error_arpack("Internal arpack error");
            break;
    }

    Mat_D eig       = mat_eig.get_impl_unique<Mat_D>();
    call_xyeupd(info, v, eig);

    bool info_nconv = (info == -14) ;
    bool none_conv  = m_iparam_ptr[5-1] == 0;
    if (info_nconv == true || none_conv == true && return_nonconverged  == false)
    {
        // XYAUPD  did not find any eigenvalues to sufficient accuracy
        m_V_result  = speye(m_N, 0, get_value_code());
        m_U_result  = m_V_result;
        m_D         = eye(0,1, get_value_code());
        m_T_result  = zeros(0,0, get_value_code());
        m_converged = false;
        return;
    }

    if (info < 0 && info != -14) // other errors
    {
        std::ostringstream msg;
        msg << "xyeupd returned code: " << info;
        throw error::error_arpack(msg.str());
    }

    Integer results_to_return = return_nonconverged ? m_num_eig 
                            : std::min(m_iparam_ptr[4], m_num_eig);
    
    bool comp_V = (m_output == arpack_output::ritz || m_output == arpack_output::schur_ritz);
    bool comp_U = (m_output == arpack_output::schur || m_output == arpack_output::schur_ritz);

    //TODO: check B-orthogonality of Schur vectors

    switch(m_kernel_type)
    {
        case kernel_type::symmetric:
        case kernel_type::complex:
        {
            m_D             = mat_eig(colon(1, results_to_return), 1);

            if (comp_V)
                m_V_result  = m_vec_ritz(colon(), colon(1, results_to_return));

            if (comp_U)
                m_U_result  = m_vec_schur(colon(), colon(1, results_to_return));

            break;
        }
        case kernel_type::nonsymmetric:
        {
            m_D         = mat_eig(colon(1, results_to_return), 1) 
                        + constants::i<TR>() * mat_eig(colon(1, results_to_return), 2);

            if (comp_V)
            {
                mat_row mr;
                for (int i = 1; i <= results_to_return; i++)
                {
                    if (imag(m_D(i)) == 0)
                    {
                        mr, m_vec_ritz(colon(), i);
                    }
                    else if (imag(m_D(i)) > 0)
                    {
                        mr, m_vec_ritz(colon(), i) + constants::i<TR>() * m_vec_ritz(colon(), i + 1);

                        if (i == results_to_return) 
                            break;

                        mr, m_vec_ritz(colon(), i) - constants::i<TR>() * m_vec_ritz(colon(), i + 1);
                        i++;
                    }
                    else
                    {   
                        mr, m_vec_ritz(colon(), i) - constants::i<TR>() * m_vec_ritz(colon(), i + 1);

                        if (i == results_to_return) 
                            break;

                        mr, m_vec_ritz(colon(), i) + constants::i<TR>() * m_vec_ritz(colon(), i + 1);
                        i++;
                    }
                }
                
                m_V_result = mr.to_matrix();
            }

            if (comp_U)
                m_U_result  = m_vec_schur(colon(), colon(1, results_to_return));

            break;
        }
        default:
        {
            throw error::error_arpack("internal arpack error");
            break;
        }
    }

    if (comp_U)
        extract_H(results_to_return);

    m_converged     = xyaupd_info == 0 ? true : false;
}

template<class T>
void arpack_wrapper<T>::extract_H(Integer results_to_return)
{
    if (m_kernel_type == kernel_type::symmetric)
    {
        m_T_result      = bdiag(m_D, 0);
        return;
    };

    Integer offset_H    = m_ipntr_ptr[12-1];
    T* work_H           = m_workl_ptr + offset_H - 1;
    Integer ldh         = m_num_vec;

    m_T_result          = make_dense_noinit(results_to_return, results_to_return, get_value_code());
    Mat_D mat_T         = m_T_result.impl_unique<Mat_D>();
    T* ptr_T            = mat_T.ptr();
    Integer T_ld        = mat_T.ld();

    for (Integer j = 0; j < results_to_return; ++j)
    {
        for (Integer i = 0; i < results_to_return; ++i)
        {
            ptr_T[i]    = work_H[i];
        };

        ptr_T           += T_ld;
        work_H          += ldh;
    };

    m_T_result.set_struct(predefined_struct_type::qtriu);
};

template<class T>
Matrix arpack_wrapper<T>::get_V() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");

    switch (m_output)
    {
        case arpack_output::eig:
        case arpack_output::schur:
            throw error::error_arpack("eigenvectors are not computed");
    }

    return m_V_result;
}

template<class T>
Integer arpack_wrapper<T>::size() const
{
    return m_N;
};

template<class T>
Matrix arpack_wrapper<T>::get_U() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");

    switch (m_output)
    {
        case arpack_output::eig:
        case arpack_output::ritz:
            throw error::error_arpack("schur vectors are not computed");
    }

    return m_U_result;
}

template<class T>
Matrix arpack_wrapper<T>::get_T() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");

    switch (m_output)
    {
        case arpack_output::eig:
        case arpack_output::ritz:
            throw error::error_arpack("schur vectors are not computed");
    }

    return m_T_result;
};

template<class T>
Matrix arpack_wrapper<T>::get_D() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");
    return m_D;
}

template<class T>
bool arpack_wrapper<T>::get_converged() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");
    return m_converged;
}

template<class T>
value_code arpack_wrapper<T>::get_value_code() const
{
    return matrix_traits::value_code<T>::value;
}

template<class T>
Integer arpack_wrapper<T>::number_converged_eigenvalues() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");
    return m_iparam_ptr[5-1];
};

template<class T>
Integer arpack_wrapper<T>::number_Arnoldi_iterations() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");
    return m_iparam_ptr[3-1];
}

template<class T>
Integer arpack_wrapper<T>::number_operations_Ax() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");
    return m_iparam_ptr[9-1];
}

template<class T>
Integer arpack_wrapper<T>::number_operations_Bx() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");
    return m_iparam_ptr[10-1];
}

template<class T>
Integer arpack_wrapper<T>::number_reorthogonalizations() const
{
    matcl_assert(m_calculate_called, "Internal Arpack error, arpack calculations not called");
    return m_iparam_ptr[11-1];
}

template<class T>
kernel_type arpack_wrapper<T>::get_kernel_type(const linear_operator& A)
{
    bool is_complex = md::is_complex<T>::value;

    if (is_complex)
        return kernel_type::complex;
    else if (A.is_hermitian())
        return kernel_type::symmetric;
    else
        return kernel_type::nonsymmetric;
}
template<class T>
kernel_type arpack_wrapper<T>::get_kernel_type(bool hermitian)
{
    bool is_complex = md::is_complex<T>::value;

    if (is_complex)
        return kernel_type::complex;
    else if (hermitian)
        return kernel_type::symmetric;
    else
        return kernel_type::nonsymmetric;
}

template<>
void arpack_wrapper<Real>::call_xyaupd(Integer& ido, Integer& info, Mat_D& v)
{
    arpack_select_eig_helper sel_helper(m_cluster);
    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";

    if (m_kernel_type == kernel_type::symmetric)
    {
        arpack::dsaupd_ (   lap(&ido), bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), 
                            lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), 
                            lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(&info) );
    }
    else 
    {        
        arpack::dnaupd_ (   lap(&ido), bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), 
                            lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), 
                            lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(&info) );
    };
}

template<>
void arpack_wrapper<Float>::call_xyaupd(Integer& ido, Integer& info, Mat_D& v)
{
    arpack_select_eig_helper sel_helper(m_cluster);
    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";

    if (m_kernel_type == kernel_type::symmetric)
    {        
        arpack::ssaupd_ (   lap(&ido), bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), 
                            lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), 
                            lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(&info) );
    }
    else 
    {        
        arpack::snaupd_ (   lap(&ido), bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), 
                            lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), 
                            lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(&info) );
    };
}

template<>
void arpack_wrapper<Complex>::call_xyaupd(Integer& ido, Integer& info, Mat_D& v)
{
    arpack_select_eig_helper sel_helper(m_cluster);
    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";

    arpack::znaupd_ (   lap(&ido), bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), 
                        lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), 
                        lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(m_rwork_ptr), lap(&info) );
}
template<>
void arpack_wrapper<Float_complex>::call_xyaupd(Integer& ido, Integer& info, Mat_D& v)
{
    arpack_select_eig_helper sel_helper(m_cluster);
    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";

    arpack::cnaupd_ (   lap(&ido), bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), 
                        lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), 
                        lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(m_rwork_ptr), lap(&info) );
}

template<>
void arpack_wrapper<Real>::call_xyeupd(Integer& info, Mat_D& v, Mat_D& eig)
{
    arpack_select_eig_helper sel_helper(m_cluster);

    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";
    Integer rvec        = (m_output == arpack_output::eig)? 0 : 1;    

    Mat_D Z             = Mat_D(v.get_type());
    Real* ptr_Z;
    Integer ldz;

    if (m_output != arpack_output::eig)
    {
        Z.assign_to_fresh(v.copy());
        ptr_Z           = Z.ptr();
        ldz             = Z.ld();
    }
    else
    {
        ptr_Z           = nullptr;
        ldz             = 1;
    };

    switch(m_kernel_type)
    {
        case kernel_type::symmetric:
        {
            const char* howmny  = "A";

            arpack::dseupd_((lapack::l_type*)lap(&rvec), howmny, (lapack::l_type*)lap(m_select_ptr), lap(eig.ptr()),
                        lap(ptr_Z), lap(&ldz), nullptr, bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), 
                        lap(m_x0_ptr), lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr),
                        lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(&info) );

            break;
        }
        case kernel_type::nonsymmetric:
        {
            const char* howmny  = (m_output == arpack_output::schur)? "P" : "A";

            arpack::dneupd_((lapack::l_type*)lap(&rvec), howmny, (lapack::l_type*)lap(m_select_ptr), lap(eig.ptr()),
                        lap(eig.ptr() + eig.ld()), lap(ptr_Z), lap(&ldz), nullptr, nullptr, lap(m_workev_ptr), bmat, 
                        lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), lap(&m_num_vec), lap(v.ptr()),
                        lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), lap(m_workd_ptr), lap(m_workl_ptr), 
                        lap(&m_lworkl), lap(&info) );

            break;
        }
        default:
        {
            throw error::error_arpack("internal arpack error");
            break;
        }
    }

    m_vec_ritz  = Matrix(Z,false);
    m_vec_schur = Matrix(v,false);    
}

template<>
void arpack_wrapper<Float>::call_xyeupd(Integer& info, Mat_D& v, Mat_D& eig)
{
    arpack_select_eig_helper sel_helper(m_cluster);
    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";

    Integer rvec        = (m_output == arpack_output::eig)? 0 : 1;    

    Mat_D Z             = Mat_D(v.get_type());
    Float* ptr_Z;
    Integer ldz;

    if (m_output != arpack_output::eig)
    {
        Z.assign_to_fresh(v.copy());
        ptr_Z           = Z.ptr();
        ldz             = Z.ld();
    }
    else
    {
        ptr_Z           = nullptr;
        ldz             = 1;
    };

    switch(m_kernel_type)
    {
        case kernel_type::symmetric:
        {
            const char* howmny  = "A";
            arpack::sseupd_((lapack::l_type*)lap(&rvec), howmny, (lapack::l_type*)lap(m_select_ptr), lap(eig.ptr()),
                        lap(ptr_Z), lap(&ldz), nullptr, bmat, lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), 
                        lap(m_x0_ptr), lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr),
                        lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(&info) );
            break;
        }
        case kernel_type::nonsymmetric:
        {
            const char* howmny  = (m_output == arpack_output::schur)? "P" : "A";

            arpack::sneupd_((lapack::l_type*)lap(&rvec), howmny, (lapack::l_type*)lap(m_select_ptr), lap(eig.ptr()),
                        lap(eig.ptr() + eig.ld()), lap(ptr_Z), lap(&ldz), nullptr, nullptr, lap(m_workev_ptr), bmat, 
                        lap(&m_N), which, lap(&m_num_eig), lap(&m_tol), lap(m_x0_ptr), lap(&m_num_vec), lap(v.ptr()),
                        lap(&m_ldv), lap(m_iparam_ptr), lap(m_ipntr_ptr), lap(m_workd_ptr), lap(m_workl_ptr), 
                        lap(&m_lworkl), lap(&info) );
            break;
        }
        default:
        {
            throw error::error_arpack("internal arpack error");
            break;
        }
    }

    m_vec_ritz  = Matrix(Z,false);
    m_vec_schur = Matrix(v,false);
}

template<>
void arpack_wrapper<Complex>::call_xyeupd(Integer& info, Mat_D& v, Mat_D& eig)
{
    arpack_select_eig_helper sel_helper(m_cluster);
    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";
    Integer rvec        = (m_output == arpack_output::eig)? 0 : 1;
    const char* howmny  = (m_output == arpack_output::schur)? "P" : "A";

    Mat_D Z             = Mat_D(v.get_type());
    Complex* ptr_Z;
    Integer ldz;

    if (m_output != arpack_output::eig)
    {
        Z.assign_to_fresh(v.copy());
        ptr_Z           = Z.ptr();
        ldz             = Z.ld();
    }
    else
    {
        ptr_Z           = nullptr;
        ldz             = 1;
    };

    arpack::zneupd_((lapack::l_type*)lap(&rvec), howmny, (lapack::l_type*)lap(m_select_ptr), lap(eig.ptr()), 
                    lap(ptr_Z), lap(&ldz), nullptr, lap(m_workev_ptr), bmat, lap(&m_N), which, lap(&m_num_eig),
                    lap(&m_tol), lap(m_x0_ptr), lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr),
                    lap(m_ipntr_ptr), lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(m_rwork_ptr), lap(&info) );

    m_vec_ritz  = Matrix(Z,false);
    m_vec_schur = Matrix(v,false);
}
template<>
void arpack_wrapper<Float_complex>::call_xyeupd(Integer& info, Mat_D& v, Mat_D& eig)
{
    arpack_select_eig_helper sel_helper(m_cluster);
    const char * which  = sel_helper.get_select_type_string().c_str();
    const char* bmat    = m_generalized? "G" : "I";

    Integer rvec        = (m_output == arpack_output::eig)? 0 : 1;
    const char* howmny  = (m_output == arpack_output::schur)? "P" : "A";

    Mat_D Z             = Mat_D(v.get_type());
    Float_complex* ptr_Z;
    Integer ldz;

    if (m_output != arpack_output::eig)
    {
        Z.assign_to_fresh(v.copy());
        ptr_Z           = Z.ptr();
        ldz             = Z.ld();
    }
    else
    {
        ptr_Z           = nullptr;
        ldz             = 1;
    };
    arpack::cneupd_((lapack::l_type*)lap(&rvec), howmny, (lapack::l_type*)lap(m_select_ptr), lap(eig.ptr()), 
                    lap(ptr_Z), lap(&ldz), nullptr, lap(m_workev_ptr), bmat, lap(&m_N), which, lap(&m_num_eig),
                    lap(&m_tol), lap(m_x0_ptr), lap(&m_num_vec), lap(v.ptr()), lap(&m_ldv), lap(m_iparam_ptr),
                    lap(m_ipntr_ptr), lap(m_workd_ptr), lap(m_workl_ptr), lap(&m_lworkl), lap(m_rwork_ptr), lap(&info) );

    m_vec_ritz  = Matrix(Z,false);
    m_vec_schur = Matrix(v,false);
}

template class arpack_wrapper<Real>;
template class arpack_wrapper<Float>;
template class arpack_wrapper<Complex>;
template class arpack_wrapper<Float_complex>;

}};
