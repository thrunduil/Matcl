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

#include "matcl-linalg/decompositions/lu/lu_struct.h"
#include "matcl-internals/func/lapack_utils.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/base/pv_constructor.h"
#include "matcl-matrep/details/struct_flag_predefined.h"

namespace matcl { namespace details
{

template<class value_type, class struct_type>
void lu_diag<value_type,struct_type>::eval(lu_return_type& ret, const Mat& A, opt::linsolve::pivot_type piv, VR tol)
{
    using VL        = typename lapack_value_type<value_type>::type;
    using Mat_D     = raw::Matrix<value_type,struct_dense>;

    Integer M       = A.rows();
    Integer N       = A.cols();    

    value_code vc   = matrix_traits::value_code<value_type>::value;
    Mat_D Ub        = raw::converter<Mat_D,decltype(A.get_diag())>
                        ::eval(A.get_diag()).make_explicit().make_unique();

    if (piv == opt::linsolve::pivot_type::partial)
    {
        Integer K   = std::min(M,N);
        Matrix L    = speye(M, K, vc);        
        Matrix U    = bdiags(Matrix(Ub,false), 0, K, N); 

        permvec P   = permvec::identity(A.rows());
        permvec Q   = permvec::identity(A.cols());

        ret = lu_return_type(L,U,P,Q);
        return;
    };

    Matrix P        = matcl::irange(1, M);
    Matrix Q        = matcl::irange(1, N);

    Integer* P_ptr  = P.get_array_unique<Integer>();
    Integer* Q_ptr  = Q.get_array_unique<Integer>();

    value_type* U_ptr = Ub.ptr();

    Integer K       = std::min(M, N);
    Integer JM      = lapack::amax(K, lap(U_ptr), 1);
    VR vmax         = abs(U_ptr[JM-1]);

    VR TOL          = tol;

    if (TOL < VR(0.0))
        TOL         = VR(10.0) * constants::eps<VR>() * vmax;
    
    Integer rank    = K;

    if (vmax == VR(0.0))
    {
        rank        = 0;
        goto lab_exit;
    };

    for (Integer i = 0; i < K; ++i)
    {
        Integer L = lapack::amax(K-i, lap(U_ptr+i), 1) - 1 + i;
        value_type amax     = U_ptr[L];

        if (abs(amax) <= TOL)
        {
            rank    = i;
            goto lab_exit;
        };

        if (i != L)
        {
            {
                value_type tmp = U_ptr[i];
                U_ptr[i]    = U_ptr[L];
                U_ptr[L]    = tmp;
            };
            {
                Integer tmp = P_ptr[i];
                P_ptr[i]    = P_ptr[L];
                P_ptr[L]    = tmp;
            };
            {
                Integer tmp = Q_ptr[i];
                Q_ptr[i]    = Q_ptr[L];
                Q_ptr[L]    = tmp;
            };
        };
    };

  lab_exit:

    for (Integer j = rank; j < K; ++j)
        U_ptr[j]    = value_type(0.0);

    matcl::Matrix L = speye(M, K, vc); 
    matcl::Matrix U = bdiags(Matrix(Ub,false), 0, K, N); 

    struct_flag so  = struct_flag(predefined_struct_type::diag);
    so.add(U.get_struct());
    U.add_struct(so);

    ret = lu_return_type(L,U, details::pv_constructor::make(P), details::pv_constructor::make(Q));
    return;
};

void details::lu_triu(lu_return_type& ret, const matcl::Matrix& A)
{
    Integer M       = A.rows();
    Integer N       = A.cols();
    Integer K       = std::min(M, N);

    matcl::Matrix L = speye(M, K, A.get_value_code()); 
    matcl::Matrix U = A(colon(1,K),colon());

    struct_flag so  = struct_flag(predefined_struct_type::triu);
    U.add_struct(so);

    permvec P = permvec::identity(A.rows());
    permvec Q = permvec::identity(A.cols());

    ret = lu_return_type(L,U,P,Q);
};

void details::lu_tril_nonunit(lu_return_type& ret, const matcl::Matrix& A)
{
    Integer M       = A.rows();
    Integer N       = A.cols();
    Integer K       = std::min(M, N);

    matcl::Matrix L = A(colon(), colon(1,K));
    matcl::Matrix U = speye(K, N, A.get_value_code()); 

    struct_flag so  = struct_flag(predefined_struct_type::tril);
    L.add_struct(so);

    permvec P = permvec::identity(M);
    permvec Q = permvec::identity(N);

    ret = lu_return_type(L,U,P,Q);
};

};};

template struct matcl::details::lu_diag<matcl::Real,matcl::struct_dense>;
template struct matcl::details::lu_diag<matcl::Float,matcl::struct_dense>;
template struct matcl::details::lu_diag<matcl::Complex,matcl::struct_dense>;
template struct matcl::details::lu_diag<matcl::Float_complex,matcl::struct_dense>;
template struct matcl::details::lu_diag<matcl::Real,matcl::struct_sparse>;
template struct matcl::details::lu_diag<matcl::Float,matcl::struct_sparse>;
template struct matcl::details::lu_diag<matcl::Complex,matcl::struct_sparse>;
template struct matcl::details::lu_diag<matcl::Float_complex,matcl::struct_sparse>;
template struct matcl::details::lu_diag<matcl::Real,matcl::struct_banded>;
template struct matcl::details::lu_diag<matcl::Float,matcl::struct_banded>;
template struct matcl::details::lu_diag<matcl::Complex,matcl::struct_banded>;
template struct matcl::details::lu_diag<matcl::Float_complex,matcl::struct_banded>;
