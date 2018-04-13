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

#include "matcl-linalg/decompositions/cholmod.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/decompositions/cholmod/cholmod.h"
#include "matcl-linalg/decompositions/cholmod_impl.inl"
#include "matcl-matrep/details/extract_type_switch.h" 

namespace matcl
{

namespace details
{
struct visitor_cholmod : public extract_type_switch<void,visitor_cholmod,true>
{
    template<class T>
    static void eval(const matcl::Matrix& handle, const T& mat, cholmod_return_type& ret, bool upper,
                        correction_alg::type alg, correction_type::type corr, Real tol)
    {
        return cholmod_impl<T>::eval(ret,handle,mat,upper,alg,corr,tol);
    };
    template<class T>
    static void eval_scalar(const Matrix& handle, const T& mat, cholmod_return_type& ret, bool upper, 
                        correction_alg::type alg, correction_type::type corr, Real tol)
    {
        using FullMatrix = raw::Matrix<T,struct_dense>;
        return eval<FullMatrix>(handle,FullMatrix(ti::get_ti<T>(mat),mat,1,1),ret,upper,alg,corr,tol);
    };
};
};

cholmod_return_type matcl::cholmod(const Matrix& A0, bool upper, 
                                   correction_alg::type alg,correction_type::type corr, Real tol)
{
    //increase refcount
    Matrix A(A0);
    A = full(A);

    cholmod_return_type ret;
    details::visitor_cholmod::make<const Matrix&>(A, ret, upper, alg, corr, tol);
    return ret;
};
cholmod_return_type matcl::cholmod(Matrix&& A0, bool upper, 
                                   correction_alg::type alg,correction_type::type corr, Real tol)
{
    Matrix A(std::move(A0));
    A = full(A);

    cholmod_return_type ret;
    details::visitor_cholmod::make<const Matrix&>(A, ret, upper, alg, corr, tol);
    return ret;
};

};