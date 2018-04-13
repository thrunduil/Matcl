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

#include "matcl-linalg/general/linalg_exception.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-core/lib_functions/constants.h"
#include "matcl-linalg/matrix_eq/lyapunov.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/func/test_inf_nan.h"

#include "matcl-linalg/linear_eq/linsolve.h"
#include "matcl-linalg/decompositions/lu.h"
#include "matcl-linalg/decompositions/eig_functions.h"
#include "matcl-linalg/decompositions/schur.h"
#include "matcl-linalg/norms_error/norm.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include <vector>

//TODO
//#include "extern_lib_src/slicot/include/slicot.h"

#if 0

namespace matcl { namespace details
{

enum class lyapunov_kind
{
    DALE, CALE
};

template<class V>
mat_tup_2 solve_lyapunov_val     (V& A, V& C, const lyapunov_kind kind)
{
    using constants::nan;
    typedef typename V::value_type  vt;
    typedef typename V::struct_type st;

    bool isv =    A.all_finite() && C.all_finite();

    if (isv == false)
    {
        return  mat_tup_2(repmat(constants::nan(), A.rows(), A.cols()), constants::nan());
    }
    Integer N = A.cols();
    
    V u(ti::ti_real(),max(1,N),N);
    V wr(ti::ti_real(),N,1);
    V wi(ti::ti_real(),N,1);

    C.set_struct(struct_flag());

    Integer info = 0;
    char dico = kind == lyapunov_kind::DALE ? 'D' : 'C' ; // discrete or continuous time
    char job = 'B'; // both solution&separation
    char fact = 'N'; // A provided non-factorized
    char trana = 'N';  // no transpose
    lapack::i_type lda = A.ld();
    lapack::i_type ldc = C.ld();
    lapack::i_type ldu = u.ld();
    Real sep = -1;
    Real ferr = -1;
    Real scale = -1;
    raw::integer_dense iwork(ti::ti_int(), N,N);
    Integer ldwork = -1;
    V dwork(A.get_type(),1,1);
    slicot::sb03md_(&dico, &job, &fact, &trana,
                    lap(&N), lap(A.ptr()), lap(&lda), lap(u.ptr()), lap(&ldu), lap(C.ptr()), lap(&ldc),
                    lap(&scale), lap(&sep), lap(&ferr), lap(wr.ptr()), lap(wi.ptr()), lap(iwork.ptr()), 
                    lap(dwork.ptr()), lap(&ldwork), lap(&info));
    ldwork = (Integer) dwork.ptr()[0];
    dwork.reset_unique(ldwork, 1);
    slicot::sb03md_(&dico, &job, &fact, &trana,
                    lap(&N), lap(A.ptr()), lap(&lda), lap(u.ptr()), lap(&ldu), lap(C.ptr()), lap(&ldc),
                    lap(&scale), lap(&sep), lap(&ferr), lap(wr.ptr()), lap(wi.ptr()), lap(iwork.ptr()), 
                    lap(dwork.ptr()), lap(&ldwork), lap(&info));
    if (info != 0 )
    {
        throw error::error_lyapunov();
    }

    return mat_tup_2(div(Matrix(C,true), scale), ferr);
}

mat_tup_2 solve_lyapunov(const Matrix& A, const Matrix& C, const lyapunov_kind kind)
{
    if (!A.is_square() || !C.is_square() || A.rows() != C.rows() )
    {
        throw error::error_size_lyapunov();
    }

    if (C.numel() == 0) return mat_tup_2(C,0.0);

    matcl::value_code vA = A.get_value_code();
    matcl::value_code vC = C.get_value_code();
    matcl::value_code v  = (matcl::value_code) max(Integer(vA),Integer(vC));
    switch (v)
    {
        case value_code::v_integer:
        case value_code::v_float:
        {
            //TODO: impl float
        }
        case value_code::v_real:
        {
            Matrix Ac = convert(A,mat_code::real_dense);
            Matrix Cc = convert(C,mat_code::real_dense);
            typedef matcl::raw::Matrix<Real,struct_dense> DM;
            return details::solve_lyapunov_val<DM>(Ac.get_impl_unique<DM>(), Cc.get_impl_unique<DM>(), kind);
        }
        case value_code::v_float_complex:
        {
            //TODO: impl float
        }
        case value_code::v_complex:
        {
            throw error::error_complex_value_type_not_allowed();
        }
        case value_code::v_object:
        {
            throw error::object_value_type_not_allowed("lyapunov");
        }
    };
    throw error::error_general("impossible type case in solve_dale/cale");
}

/**
 * Class to sparse-solve lyapunov equations (DALE & CALE)
 *
 * DALE: AXA' - X = - GG'
 * CALE: AX + XA' = -GG'
 * Note the position of the transpose!
 */   
class sparse_lyapunov_calculator
{
public:
    virtual ~sparse_lyapunov_calculator()
    {}
    virtual Matrix sparse_solve_lyapunov(const Matrix& A, const Matrix& G, const Integer l_zero_cap);
protected:
    virtual Matrix scalar_solution() const = 0;
    virtual void check_spectrum() const = 0;
    virtual Matrix do_LRCF_ADI_iterations() const = 0;
    virtual Real calculate_s_p_t(const Matrix& P, const Complex t) const = 0;
protected:

    Matrix get_structured_eye() const;

    Matrix D; // part of spectrum to use
    Integer l_zero; // number of ADI-shift parameters to use
    Matrix A, G;
    Matrix set_P;

    Integer l_zero_cap; // maximum number of ADI-shifts parameters to consider
private:
    void calculate_spectrum();
    void build_ADI_shift_parameters_set();
};

class dale_lyapunov : public sparse_lyapunov_calculator
{
protected:
    virtual Matrix scalar_solution() const;
    virtual void check_spectrum() const;
    virtual Matrix do_LRCF_ADI_iterations() const;
    virtual Real calculate_s_p_t(const Matrix& P, const Complex t) const;
};
class cale_lyapunov : public sparse_lyapunov_calculator
{
protected:
    virtual Matrix scalar_solution() const;
    virtual void check_spectrum() const;
    virtual Matrix do_LRCF_ADI_iterations() const;
    virtual Real calculate_s_p_t(const Matrix& P, const Complex t) const;
public:
    cale_lyapunov(const Real V_update_tol_in)
        : V_update_tol(V_update_tol_in)
    {}
private:
    const Real V_update_tol;
};

Matrix sparse_lyapunov_calculator::sparse_solve_lyapunov  (const Matrix& A_in, const Matrix& G_in,
                                                           const Integer l_zero_cap_in)
{
    A = A_in; G = G_in; l_zero_cap = l_zero_cap_in;
    if (!A.is_square() || A.rows() != G.rows() )
    {
        throw error::error_size_lyapunov();
    }
    if (A.numel() == 0) return A;

    matcl::value_code vA = A.get_value_code();
    matcl::value_code vG = G.get_value_code();
    matcl::value_code v  = (matcl::value_code) max(Integer(vA),Integer(vG));
    switch (v)
    {
        case value_code::v_integer:
        case value_code::v_float:
        {
            //TODO: impl float
        }
        case value_code::v_real:
            break; // OK
        case value_code::v_float_complex:
        {
            //TODO: impl float
        }
        case value_code::v_complex:
        {
            throw error::error_complex_value_type_not_allowed();
        }
        case value_code::v_object:
        {
            throw error::object_value_type_not_allowed("lyapunov");
        }
        default:
            throw error::error_general("impossible type case in solve_dale/cale");
    }
    Integer n = A.rows();
    if (n == 1) 
    {
        return scalar_solution();
    }
    if (any(any(is_nan(A), 1), 2) || any(any(is_inf(A), 1), 2) ||
        any(any(is_nan(G), 1), 2) || any(any(is_inf(G), 1), 2) )
    {
        return  repmat(constants::nan(), n, n);
    }
    // algo: suboptimal set of ADI-shift parameters
    calculate_spectrum();
    check_spectrum();
    build_ADI_shift_parameters_set();
    // algo: LRCF_ADI iterations
    return do_LRCF_ADI_iterations();
}

Matrix  cale_lyapunov::scalar_solution() const
{
    if ((A >= -1e-08)) throw error::error_lyapunov();
    return sqrt(div(-G*trans(G), 2 * A));
}
Matrix  dale_lyapunov::scalar_solution() const
{
    if ((abs(A) >= 1 - 1e-08)) throw error::error_lyapunov();
    return sqrt(div(-G*trans(G), (A*A - 1)));
}

void    sparse_lyapunov_calculator::calculate_spectrum()
{
    Integer n = A.rows();
    if (n < 120)
    {
        D = schur_decomposition(A).eig();
        l_zero = n;
    }
    else
    {
        Integer k_plus =    std::min(n/2 + 1, 50);
        Integer k_minus =   std::min(n/2, 25);
        Matrix E_plus;
        bool conv;
        
        options opts{opt::speigs::return_nonconvergent(false), opt::speigs::tol(1e-8)};

        tie(E_plus,conv) = eigs(A, k_plus, cluster_type::LM, opts);
        Matrix E_minus;
        tie(E_minus,conv) = eigs(A, k_minus, cluster_type::SM, opts);
        D =  vec(vertcat(E_plus, E_minus));
        l_zero = std::min(k_plus + k_minus, l_zero_cap);
    }
    D = convert(D, mat_code::complex_dense); // in case D turns out Real
}

void    cale_lyapunov::check_spectrum() const
{
    const Matrix nonnegative_real = find(real(D) >= -1e-08);
    if (nonnegative_real.rows() > 0) 
    {
        throw error::error_lyapunov();
    }
}
void    dale_lyapunov::check_spectrum() const
{
    const Matrix outside_unit = find(abs(D) >= 1.0 - 1e-08);
    if (outside_unit.rows() > 0) 
    {
        throw error::error_lyapunov();
    }
}

void    sparse_lyapunov_calculator::build_ADI_shift_parameters_set()
{
    Matrix max_sp = zeros(D.rows(), 1);
    for (int i = 1; i <= D.rows(); i++)
    {
        Matrix P = D(i);
        for (int j = 1; j <= D.rows(); j ++)
        {
            max_sp(i) = max(max_sp(i), calculate_s_p_t(P, D(j).to_matrix().get_scalar<Complex>()));
        }
    }
    Matrix sorted_max, order;
    tie(sorted_max, order) = sort2(max_sp);
    if ( (imag(D(order(1))) == 0.0))
    {
        set_P = D(order(1));
    }
    else
    {
        set_P = vertcat(D(order(1)), ctrans(D(order(1))));
    }
    set_P.reserve(l_zero,1);
    while( set_P.rows() < l_zero )
    {
        Matrix sp = zeros(D.rows(), 1);
        for (int j = 1; j <= D.rows(); j ++)
        {
            sp(j) = calculate_s_p_t(set_P, D(j).to_matrix().get_scalar<Complex>());
        }
        Matrix sorted_sp, order_sp;
        tie(sorted_sp, order_sp) = sort2(sp, 1, false); // find max
        if ( (imag(D(order_sp(1))) == 0.0) )
        {
            set_P = vertcat(set_P, D(order_sp(1)));
        }
        else
        {
            set_P = (mat_col(), set_P, D(order_sp(1)), ctrans(D(order_sp(1))));
        }
    }
}

Matrix  cale_lyapunov::do_LRCF_ADI_iterations() const
{
    Matrix I_n = get_structured_eye();

    Matrix V = sqrt( -2 * real(set_P(1))) * linsolve(A + set_P(1) * I_n, G);
    mat_row Z;
    Z, V;
    for (int i = 2; i <= set_P.rows(); i++)
    {
        V = sqrt ( div(real(set_P(i)), real(set_P(i-1)) ))
            * (V - (set_P(i) + ctrans(set_P(i-1))) * linsolve(A + set_P(i) * I_n, V));
        Z, V;
        if (norm(V, -2) < V_update_tol) return Z;
    }
    return Z;
}
Matrix  dale_lyapunov::do_LRCF_ADI_iterations() const
{
    Matrix I_n = get_structured_eye();
    Matrix start_M_hat = linsolve(I_n - ctrans(set_P(1)) * A, sqrt(1 - pow(abs(set_P(1)), 2)) * A * G);
    mat_row Z;
    Z, G, start_M_hat;
    for (int i = 2; i <= set_P.rows(); i++)
    {
        // if shift parameter is real, things may work faster in Real domain
        Matrix shift_parameter, c_shift_parameter;
        if((imag(set_P(i)) == 0.0))
        {
            c_shift_parameter = shift_parameter = real(set_P(i));
        }
        else
        {
            shift_parameter = set_P(i);
            c_shift_parameter = ctrans(shift_parameter);
        }
        Matrix lu_L, lu_U; permvec lu_p,lu_q;
        std::tie(lu_L, lu_U, lu_p, lu_q) = lu(I_n - c_shift_parameter * A);
        permvec id_permutation = permvec::identity(A.rows());

        // NOTE: M_hat is real while M_wave usually complex!
        Matrix M_hat  = linsolve(lu_U, id_permutation, lu_q, linsolve(lu_L, lu_p, id_permutation, (sqrt(1 - pow(abs(shift_parameter), 2)) * A * G)));
        Matrix M_wave = linsolve(lu_U, id_permutation, lu_q, linsolve(lu_L, lu_p, id_permutation, A)) * (A - shift_parameter * I_n) * Z;
        
        Z = mat_row(), G, M_hat, M_wave;
    }
    return Z;
}


Real cale_lyapunov::calculate_s_p_t(const Matrix& P, const Complex t) const
{
    Real numerator = abs(prod(t - P)).get_scalar<Real>();
    Real denominator = abs(prod(t + P)).get_scalar<Real>();
    return numerator / denominator;
}
Real dale_lyapunov::calculate_s_p_t(const Matrix& P, const Complex t) const
{
    if (t == 0.0) return 0.0;
    Real numerator = abs(prod(t - P)).get_scalar<Real>();
    Real denominator = abs(prod(div(1.f,t) - P)).get_scalar<Real>();
    return numerator / denominator;
}

Matrix sparse_lyapunov_calculator::get_structured_eye() const
{
    Integer n = A.rows();
    switch (A.get_struct_code())
    {
    case struct_code::struct_banded:
        return beye(n,n,0,0);
    case struct_code::struct_sparse:
    case struct_code::struct_dense:
    default:
        return speye(n,n);
    }
}

} // end namespace details


mat_tup_2 solve_dale            (const Matrix& A, const Matrix& C)
{
    return details::solve_lyapunov(A, C, details::lyapunov_kind::DALE);
}

mat_tup_2 solve_cale            (const Matrix& A, const Matrix& C)
{
    return details::solve_lyapunov(A, C, details::lyapunov_kind::CALE);
}

Matrix    sparse_dale           (const Matrix& A, const Matrix& G, const Integer l_zero_cap)
{
    boost::scoped_ptr<details::sparse_lyapunov_calculator> 
        impl(new details::dale_lyapunov());
    return impl->sparse_solve_lyapunov(trans(A), G, l_zero_cap);
}
Matrix    sparse_cale           (const Matrix& A, const Matrix& G, const Integer l_zero_cap, const Real V_update_tol)
{
    boost::scoped_ptr<details::sparse_lyapunov_calculator> 
        impl(new details::cale_lyapunov(V_update_tol));
    return impl->sparse_solve_lyapunov(trans(A), G, l_zero_cap);
}
}

#endif