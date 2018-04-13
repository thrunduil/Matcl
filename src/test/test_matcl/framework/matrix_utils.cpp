/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe� Kowal 2018
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

#include "matrix_utils.h"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-internals/base/utils.h"
#include "matcl-matrep/details/extract_type_switch.h" 
#include "matcl-scalar/details/scalfunc_helpers.h" 
#include "matcl-scalar/details/matfunc_helpers.h" 
#include "matcl-internals/base/optim_params.h"
#include "matcl-internals/func/converter.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-matrep/details/matrix.inl"

namespace matcl
{

namespace md = matcl::details;
namespace mrd = matcl::raw::details;
namespace mdyp = matcl::dynamic::functions;

template<class VT, class ST, class T>
struct norm_1_impl
{};

template<class VT, class T>
struct norm_1_impl<VT,struct_dense,T>
{
    using real_type = typename details::object_or_type<VT,void, Real>::type;

    static Real eval(const T& mat)
    {
        Integer r           = mat.rows();
        Integer c           = mat.cols();
        ti::ti_type<VT> ti  = ti::get_return_ti<ti::ti_type<VT>>(mdyp::real::eval(),ti::get_ti(mat));
        real_type out       = md::default_value<real_type>(ti);    

        const VT* ptr_m     = mat.ptr();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = 0; i < r; ++i)
            {
                real_type val = matcl::raw::details::abs_helper<VT>::eval(ptr_m[i]);
                if ((bool)mrd::gt_helper<real_type,real_type>::eval(val, out) 
                              && (bool)mrd::isnan_helper<real_type>::eval(val) == false
                              && (bool)mrd::isinf_helper<real_type>::eval(val) == false)
                {
                    out = val;
                };
            }
            ptr_m += mat.ld();
        };
        return to_real(out);
    };

    template<class RT>
    static Real to_real(const RT& out)
    {
        return out;
    };
    static Real to_real(const Object& out)
    {
        return cast_real(out);
    }
};

template<class VT, class T>
struct norm_1_impl<VT,struct_sparse,T>
{
    using real_type = typename details::object_or_type<VT,void, Real>::type;

    static Real eval(const T& mat)
    {
        Integer c           = mat.cols();
        ti::ti_type<VT> ti  = ti::get_return_ti<ti::ti_type<VT>>(mdyp::real::eval(),ti::get_ti(mat));
        real_type out       = md::default_value<real_type>(ti);    

        const raw::details::sparse_ccs<VT>& rep = mat.rep();

        const Integer * Ad_c	= rep.ptr_c();
        const VT * Ad_x	        = rep.ptr_x();

        for (Integer j = 0; j < c; ++j)
        {
            for (Integer i = Ad_c[j]; i < Ad_c[j+1]; ++i)
            {
                real_type val = matcl::raw::details::abs_helper<VT>::eval(Ad_x[i]);
                if ((bool)mrd::gt_helper<real_type,real_type>::eval(val, out) 
                              && (bool)mrd::isnan_helper<real_type>::eval(val) == false
                              && (bool)mrd::isinf_helper<real_type>::eval(val) == false)
                {
                    out = val;
                };
            }
        };
        return to_real(out);
    };
    template<class RT>
    static Real to_real(const RT& out)
    {
        return out;
    };
    static Real to_real(const Object& out)
    {
        return cast_real(out);
    }
};
template<class VT, class T>
struct norm_1_impl<VT,struct_banded,T>
{
    static Real eval(const T& mat)
    {
        Real density = Real(mat.nnz())/(mat.rows()+1.)/(mat.cols()+1.);

        if (density < optim_params::max_sparse_density_min)
        {
            using SparseMatrix = raw::Matrix<VT,struct_sparse>;
            SparseMatrix Bc = raw::converter<SparseMatrix,T>::eval(mat);
            return norm_1_impl<VT,struct_sparse,SparseMatrix>::eval(Bc);
        }
        else
        {
            using DenseMatrix = raw::Matrix<VT,struct_dense>;
            DenseMatrix Bc = raw::converter<DenseMatrix,T>::eval(mat);

            return norm_1_impl<VT,struct_dense,DenseMatrix>::eval(Bc);
        };
    };
};

struct unary_visitor_norm_1 : public matcl::details::extract_type_switch<Real,unary_visitor_norm_1,true>
{
    template<class T>
    static Real eval(const Matrix&, const T& mat)
    {        
        using V = typename T::value_type;
        using S = typename T::struct_type;

        return norm_1_impl<V,S,T>::eval(mat);
    };
    template<class T>
    static Real eval_scalar(const Matrix& m, const T& mat)
    {
        (void)mat;

        Matrix f = full(m);
        using DM = matcl::raw::Matrix<T,struct_dense>;
        return eval<DM>(f,f.impl<DM>());
    };
};

Real matcl::norm_1(const Matrix& A)
{
    if (A.structural_nnz() == 0)
        return 0;

    return unary_visitor_norm_1::make<const Matrix&>(A);
};

};