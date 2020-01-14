/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2018
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

#include "matcl-matrep/container/matrix2.inl"
#include "matcl-matrep/container/matrix_container.inl"
#include "matcl-matrep/lib_functions/manip.h"
#include "matcl-scalar/object.h"
#include "matcl-internals/base/instantiate.h"
#include "matcl-matrep/objects/details/type_info_object.h"
#include "matcl-matrep/details/details_manip.h"

#pragma warning( push )
#pragma warning(disable:4127)	// conditional expression is constant
#pragma warning(disable:4702)	// unreachable code

namespace matcl
{

namespace details
{
    template<class T>
    MATCL_MATREP_EXPORT const T* get_rep_functor<T>::eval_array_c(const Matrix& m)
    {
        using DM                = raw::Matrix<T,struct_dense>;
        matcl::mat_code ret_mt  = matrix_traits::mat_type_info_type<DM>::matrix_code;

        if (m.get_matrix_code() != ret_mt)
            details::matrix_data_accesser::assign_to_const_mat(m,matcl::convert(m,ret_mt));

        const DM& mat = m.get_impl<DM>();        
 
        if (mat.is_explicit() == false)
        {
            details::matrix_data_accesser::assign_to_const_mat(m,Matrix(mat.make_explicit(),false));

            const DM& mat2 = m.get_impl<DM>();
            return mat2.ptr();
        }
        else
        {
            return mat.ptr();
        };
    };

    template<class T>
    T* get_rep_functor<T>::eval_array_nc(Matrix& m)
    {
        using DM                = raw::Matrix<T,struct_dense>;
        matcl::mat_code ret_mt  = matrix_traits::mat_type_info_type<DM>::matrix_code;

        if (m.get_matrix_code() != ret_mt)
            m = matcl::convert(m,ret_mt);

        DM& mat = m.get_impl_unique<DM>();
        m.set_struct(matcl::struct_flag());

        if (mat.is_explicit() == false)
        {
            m           = Matrix(mat.make_explicit(),false);
            DM& mat2    = m.get_impl_unique<DM>();

            return mat2.ptr();
        }
        else
        {
            return mat.ptr();
        };
    };

    template<class M, class T>
    struct get_scalar
    {
        static M eval_get(const T& val)
        {
            return matcl::details::convert_scalar<M,T>(val);
        };

        static M& eval_get_unique(T& )
        {
            matcl::mat_code ret = matrix_traits::mat_type_info_type<M>::matrix_code;
            matcl::mat_code in  = matrix_traits::mat_type_info_type<T>::matrix_code;
            throw error::invalid_type_get(ret,in);
        };
    };

    template<class M>
    struct get_scalar<M,M>
    {
        static const M& eval_get(const M& val)
        {
            return val;
        };

        static M eval_get(M&& val)
        {
            return val;
        };

        static M& eval_get(M& val)
        {
            return val;
        };

        static M& eval_get_unique(M& val)
        {
            return val;
        };
    };

    template<class M>
    M get_scalar_functor<M>::eval_get(const matrix_base& mat)
    {
        switch(mat.m_type)
        {
            case mat_code::integer_scalar:
            {
                return get_scalar<M,Integer>::eval_get(mat.m_value.val_int);
            }
            case mat_code::real_scalar:
            {
                return get_scalar<M,Real>::eval_get(mat.m_value.val_real);
            }
            case mat_code::float_scalar:
            {
                return get_scalar<M,Float>::eval_get(mat.m_value.val_float);
            }
            case mat_code::complex_scalar:
            {
                return get_scalar<M,Complex>::eval_get(*reinterpret_cast<const Complex*>
                                                       (&mat.m_value.val_complex));
            }
            case mat_code::float_complex_scalar:
            {
                const Float_complex& tmp = *reinterpret_cast<const Float_complex*>(&mat.m_value.val_fcomplex);
                return get_scalar<M,Float_complex>::eval_get(tmp);
            }
            case mat_code::object_scalar:
            {
                return M(get_scalar<M,Object>::eval_get(mat.get_object()));
            }
            default:
            {
                return eval_get(matrix_data_accesser::get_base(mat.m_value.m_mat.m_mat_ptr->get_scalar()));
            }
        };
    };

    template<class M>
    M& get_scalar_functor<M>::eval_get_unique(matcl::Matrix& mat)
    {
        matcl::mat_code ret_mt = matrix_traits::mat_type_info_type_2<M,struct_scalar>::matrix_code;

        if (static_cast<int>(mat.get_value_code()) != static_cast<int>(ret_mt))
            mat = matcl::convert(mat,ret_mt);

        details::matrix_base& mb = details::matrix_data_accesser::get_base(mat);
        switch(ret_mt)
        {
            case mat_code::integer_scalar:
            {
                return get_scalar<M,Integer>::eval_get_unique(mb.m_value.val_int);
            }
            case mat_code::real_scalar:
            {
                return get_scalar<M,Real>::eval_get_unique(mb.m_value.val_real);
            }
            case mat_code::float_scalar:
            {
                return get_scalar<M,Float>::eval_get_unique(mb.m_value.val_float);
            }
            case mat_code::complex_scalar:
            {
                return get_scalar<M,Complex>::eval_get_unique(*reinterpret_cast<Complex*>
                                                              (&mb.m_value.val_complex));
            }
            case mat_code::float_complex_scalar:
            {
                return get_scalar<M,Float_complex>::eval_get_unique(*reinterpret_cast<Float_complex*>
                                                                    (&mb.m_value.val_fcomplex));
            }
            case mat_code::object_scalar:
            {
                return get_scalar<M,Object>::eval_get_unique(mb.get_object());
            }
            default:
            {
                matcl_assert(0,"invalid matrix type");
                throw;
            }
        };
    };

    template<class T,class TS>
    struct to_scalar_impl
    {
        static TS eval(const T& val)
        {
            return TS(val(1,1));
        };
    };

    template<class T>
    struct to_scalar
    {
        using TS = typename T::value_type;
        static TS eval(const T& val)
        {
            return to_scalar_impl<T,TS>::eval(val);
        };
    };

    template<class T>
    struct constructor_helper_scal
    {};

    template<>
    struct constructor_helper_scal<Integer>
    {
        static void eval(Integer val,matrix_base& mat)
        {
            mat.m_value.val_int = val;
            mat.m_type = mat_code::integer_scalar; 
        };
    };
    
    template<>
    struct constructor_helper_scal<Real>
    {
        static void eval(Real val,matrix_base& mat)
        {
            mat.m_value.val_real    = val;
            mat.m_type              = mat_code::real_scalar; 
        };
    };

    template<>
    struct constructor_helper_scal<Float>
    {
        static void eval(Float val,matrix_base& mat)
        {
            mat.m_value.val_float   = val;
            mat.m_type              = mat_code::float_scalar; 
        };
    };

    template<>
    struct constructor_helper_scal<Float_complex>
    {
        static void eval(const Float_complex& val,matrix_base& mat)
        {
            mat.m_value.val_fcomplex[0] = matcl::real(val);
            mat.m_value.val_fcomplex[1] = matcl::imag(val);
            mat.m_type                  = mat_code::float_complex_scalar; 
        };
    };
    
    template<>
    struct constructor_helper_scal<Object>
    {
        static void eval(const Object& val,matrix_base& mat)
        {
            new(&mat.get_object()) Object(val);
            mat.m_type = mat_code::object_scalar; 
        };
    };

    template<>
    struct constructor_helper_scal<Complex>
    {
        static void eval(const Complex& val,matrix_base& mat)
        {
            mat.m_value.val_complex[0]  = matcl::real(val);
            mat.m_value.val_complex[1]  = matcl::imag(val);
            mat.m_type                  = mat_code::complex_scalar; 
        };
    };

    template<class T,bool is_scal>
    struct constructor_helper_impl
    {        
        static bool eval(const T& val,bool allow_conv, matrix_base& mat)
        {            
            if (allow_conv)                
            {
                if (val.rows() == 1 && val.cols() == 1)
                {
                    constructor_helper_scal<typename T::value_type>::eval(to_scalar<T>::eval(val),mat);
                    return false;
                }
                else
                {
                    Matrix tmp = val.fast_optim();
                    mat.reset(matrix_data_accesser::get_base(tmp));
                    return true;
                };
            }
            else
            {
                mat.m_value.m_mat.m_refcount= val.get_refstr();
                mat.m_value.m_mat.m_mat_ptr = create_container<T>::eval(val);
                mat.m_type			        = type_to_code<typename create_container<T>::matrix_type>::value;
                return true;
            };
        };

        static bool eval(T&& val,bool allow_conv, matrix_base& mat)
        {
            if (allow_conv)                
            {
                if (val.rows() == 1 && val.cols() == 1)
                {
                    constructor_helper_scal<typename T::value_type>::eval(to_scalar<T>::eval(val),mat);
                    return false;
                };
            }

            mat.m_value.m_mat.m_refcount= val.get_refstr();
            mat.m_value.m_mat.m_mat_ptr = create_container<T>::eval(std::move(val));            
            mat.m_type			        = type_to_code<typename create_container<T>::matrix_type>::value;

            return true;
        };
    };

    template<class T>
    struct constructor_helper_impl<raw::sparse_matrix_base<T>, false>
    {        
        using arg_type  = raw::sparse_matrix_base<T>;
        using mat_type  = raw::Matrix<T, struct_sparse>;

        static bool eval(const arg_type& val, bool allow_conv, matrix_base& mat)
        {            
            // we can safely create a (shallow) copy; returned object modification
            // is controlled by Matrix, which should create independent copy when
            // object is modified and is not unique
            mat_type vm   = mat_type(val, mat_type::copy_is_safe());

            return constructor_helper_impl<mat_type, false>::eval(vm, allow_conv, mat);
        };

        static bool eval(arg_type&& val,bool allow_conv, matrix_base& mat)
        {
            mat_type vm   = mat_type(std::move(val));

            return constructor_helper_impl<mat_type, false>::eval(std::move(vm), allow_conv, mat);
        };
    };

    template<class T>
    struct constructor_helper_impl<T,true>
    {
        static bool eval(const T& val,bool ,matrix_base& mat)
        {
            using scalar_type = typename details::promote_scalar<T>::type;
            constructor_helper_scal<scalar_type>::eval(static_cast<scalar_type>(val),mat);
            return false;
        };
    };

    template<class T>
    bool constructor_helper<T>::eval(const T& val,bool allow_conv,matrix_base& mat)
    {
        return constructor_helper_impl<T,is_scalar<T>::value>::eval(val,allow_conv,mat);
    };
    
    template<class T>
    bool constructor_helper<T>::eval(T&& val,bool allow_conv,matrix_base& mat)
    {
        return constructor_helper_impl<T,is_scalar<T>::value>::eval(std::move(val),allow_conv,mat);
    };
};

template<class val_type, class str_type>
void details::matrix_container<val_type,str_type>::disp_impl(const disp_stream_ptr& os) const	
{ 
    raw::disp(os,m_matrix); 
};

bool details::matrix_base::is_same_matrix(const details::matrix_base& A1, const details::matrix_base& A2)
{
    if (A1.m_type < mat_code::integer_dense || A2.m_type < mat_code::integer_dense)
        return false;

    //A1 and A2 are matrices

    // if reference counters are different, then A1 and A2 have own data pointers
    if (A1.m_value.m_mat.m_refcount != A2.m_value.m_mat.m_refcount)
        return false;    

    if (A1.m_value.m_mat.m_mat_ptr == A2.m_value.m_mat.m_mat_ptr)
        return true;

    //A1 and A2 share the same data pointers, but may be different views of the same matrix
    return A1.m_value.m_mat.m_mat_ptr->is_same_matrix(*A2.m_value.m_mat.m_mat_ptr);
};

bool details::matrix_base::check_effective_unique() const
{
    return m_value.m_mat.m_mat_ptr->is_effective_unique();
};

template void details::matrix_container<Integer,struct_dense>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Integer,struct_sparse>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Integer,struct_banded>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Real,struct_dense>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Real,struct_sparse>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Real,struct_banded>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Float,struct_dense>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Float,struct_sparse>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Float,struct_banded>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Complex,struct_dense>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Complex,struct_sparse>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Complex,struct_banded>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Float_complex,struct_dense>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Float_complex,struct_sparse>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Float_complex,struct_banded>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Object,struct_dense>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Object,struct_sparse>::disp_impl(const disp_stream_ptr& os) const;
template void details::matrix_container<Object,struct_banded>::disp_impl(const disp_stream_ptr& os) const;

template const raw::Matrix<Integer,struct_dense>& 	Matrix::get_impl<raw::Matrix<Integer,struct_dense> >() const;
template const raw::Matrix<Real,struct_dense>&		Matrix::get_impl<raw::Matrix<Real,struct_dense> >() const;
template const raw::Matrix<Float,struct_dense>&		Matrix::get_impl<raw::Matrix<Float,struct_dense> >() const;
template const raw::Matrix<Complex,struct_dense>&	Matrix::get_impl<raw::Matrix<Complex,struct_dense> >() const;
template const raw::Matrix<Float_complex,struct_dense>&	Matrix::get_impl<raw::Matrix<Float_complex,struct_dense> >() const;
template const raw::Matrix<Object,struct_dense>&	Matrix::get_impl<raw::Matrix<Object,struct_dense> >() const;

template const raw::Matrix<Integer,struct_sparse>&	Matrix::get_impl<raw::Matrix<Integer,struct_sparse> >() const;
template const raw::Matrix<Real,struct_sparse>&		Matrix::get_impl<raw::Matrix<Real,struct_sparse>>() const;
template const raw::Matrix<Float,struct_sparse>&	Matrix::get_impl<raw::Matrix<Float,struct_sparse>>() const;
template const raw::Matrix<Complex,struct_sparse>&	Matrix::get_impl<raw::Matrix<Complex,struct_sparse> >() const;
template const raw::Matrix<Float_complex,struct_sparse>& Matrix::get_impl<raw::Matrix<Float_complex,struct_sparse> >() const;
template const raw::Matrix<Object,struct_sparse>&	Matrix::get_impl<raw::Matrix<Object,struct_sparse> >() const;

template const raw::Matrix<Integer,struct_banded>&	Matrix::get_impl<raw::Matrix<Integer,struct_banded> >() const;
template const raw::Matrix<Real,struct_banded>&		Matrix::get_impl<raw::Matrix<Real,struct_banded> >() const;
template const raw::Matrix<Float,struct_banded>&	Matrix::get_impl<raw::Matrix<Float,struct_banded> >() const;
template const raw::Matrix<Complex,struct_banded>&	Matrix::get_impl<raw::Matrix<Complex,struct_banded> >() const;
template const raw::Matrix<Float_complex,struct_banded>& Matrix::get_impl<raw::Matrix<Float_complex,struct_banded> >() const;
template const raw::Matrix<Object,struct_banded>&	Matrix::get_impl<raw::Matrix<Object,struct_banded> >() const;

// get_functor

template struct details::get_functor<raw::Matrix<Integer,struct_dense> >;
template struct details::get_functor<raw::Matrix<Real,struct_dense> >;
template struct details::get_functor<raw::Matrix<Float,struct_dense> >;
template struct details::get_functor<raw::Matrix<Complex,struct_dense> >;
template struct details::get_functor<raw::Matrix<Float_complex,struct_dense> >;
template struct details::get_functor<raw::Matrix<Object,struct_dense> >;

template struct details::get_functor<raw::Matrix<Integer,struct_sparse> >;
template struct details::get_functor<raw::Matrix<Real,struct_sparse>>;
template struct details::get_functor<raw::Matrix<Float,struct_sparse>>;
template struct details::get_functor<raw::Matrix<Complex,struct_sparse> >;
template struct details::get_functor<raw::Matrix<Float_complex,struct_sparse> >;
template struct details::get_functor<raw::Matrix<Object,struct_sparse> >;

template struct details::get_functor<raw::Matrix<Integer,struct_banded> >;
template struct details::get_functor<raw::Matrix<Real,struct_banded> >;
template struct details::get_functor<raw::Matrix<Float,struct_banded> >;
template struct details::get_functor<raw::Matrix<Complex,struct_banded> >;
template struct details::get_functor<raw::Matrix<Float_complex,struct_banded> >;
template struct details::get_functor<raw::Matrix<Object,struct_banded> >;

// assign_functor

template struct details::assign_functor<raw::Matrix<Integer,struct_dense>, Integer >;
template struct details::assign_functor<raw::Matrix<Real,struct_dense>, Integer >;
template struct details::assign_functor<raw::Matrix<Float,struct_dense>, Integer >;
template struct details::assign_functor<raw::Matrix<Complex,struct_dense>, Integer >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_dense>, Integer >;
template struct details::assign_functor<raw::Matrix<Object,struct_dense>, Integer >;

template struct details::assign_functor<raw::Matrix<Integer,struct_sparse>, Integer >;
template struct details::assign_functor<raw::Matrix<Real,struct_sparse>, Integer>;
template struct details::assign_functor<raw::Matrix<Float,struct_sparse>, Integer>;
template struct details::assign_functor<raw::Matrix<Complex,struct_sparse>, Integer >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_sparse>, Integer >;
template struct details::assign_functor<raw::Matrix<Object,struct_sparse>, Integer >;

template struct details::assign_functor<raw::Matrix<Integer,struct_banded>, Integer >;
template struct details::assign_functor<raw::Matrix<Real,struct_banded>, Integer >;
template struct details::assign_functor<raw::Matrix<Float,struct_banded>, Integer >;
template struct details::assign_functor<raw::Matrix<Complex,struct_banded>, Integer >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_banded>, Integer >;
template struct details::assign_functor<raw::Matrix<Object,struct_banded>, Integer >;

template struct details::assign_functor<raw::Matrix<Integer,struct_dense>, Real >;
template struct details::assign_functor<raw::Matrix<Real,struct_dense>, Real >;
template struct details::assign_functor<raw::Matrix<Float,struct_dense>, Real >;
template struct details::assign_functor<raw::Matrix<Complex,struct_dense>, Real >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_dense>, Real >;
template struct details::assign_functor<raw::Matrix<Object,struct_dense>, Real >;

template struct details::assign_functor<raw::Matrix<Integer,struct_sparse>, Real >;
template struct details::assign_functor<raw::Matrix<Real,struct_sparse>, Real>;
template struct details::assign_functor<raw::Matrix<Float,struct_sparse>, Real>;
template struct details::assign_functor<raw::Matrix<Complex,struct_sparse>, Real >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_sparse>, Real >;
template struct details::assign_functor<raw::Matrix<Object,struct_sparse>, Real >;

template struct details::assign_functor<raw::Matrix<Integer,struct_banded>, Real >;
template struct details::assign_functor<raw::Matrix<Real,struct_banded>, Real >;
template struct details::assign_functor<raw::Matrix<Float,struct_banded>, Real >;
template struct details::assign_functor<raw::Matrix<Complex,struct_banded>, Real >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_banded>, Real >;
template struct details::assign_functor<raw::Matrix<Object,struct_banded>, Real >;

template struct details::assign_functor<raw::Matrix<Integer,struct_dense>, Float >;
template struct details::assign_functor<raw::Matrix<Real,struct_dense>, Float >;
template struct details::assign_functor<raw::Matrix<Float,struct_dense>, Float >;
template struct details::assign_functor<raw::Matrix<Complex,struct_dense>, Float >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_dense>, Float >;
template struct details::assign_functor<raw::Matrix<Object,struct_dense>, Float >;

template struct details::assign_functor<raw::Matrix<Integer,struct_sparse>, Float >;
template struct details::assign_functor<raw::Matrix<Real,struct_sparse>, Float>;
template struct details::assign_functor<raw::Matrix<Float,struct_sparse>, Float>;
template struct details::assign_functor<raw::Matrix<Complex,struct_sparse>, Float >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_sparse>, Float >;
template struct details::assign_functor<raw::Matrix<Object,struct_sparse>, Float >;

template struct details::assign_functor<raw::Matrix<Integer,struct_banded>, Float >;
template struct details::assign_functor<raw::Matrix<Real,struct_banded>, Float >;
template struct details::assign_functor<raw::Matrix<Float,struct_banded>, Float >;
template struct details::assign_functor<raw::Matrix<Complex,struct_banded>, Float >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_banded>, Float >;
template struct details::assign_functor<raw::Matrix<Object,struct_banded>, Float >;

template struct details::assign_functor<raw::Matrix<Integer,struct_dense>, Complex >;
template struct details::assign_functor<raw::Matrix<Real,struct_dense>, Complex >;
template struct details::assign_functor<raw::Matrix<Float,struct_dense>, Complex >;
template struct details::assign_functor<raw::Matrix<Complex,struct_dense>, Complex >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_dense>, Complex >;
template struct details::assign_functor<raw::Matrix<Object,struct_dense>, Complex >;

template struct details::assign_functor<raw::Matrix<Integer,struct_sparse>, Complex >;
template struct details::assign_functor<raw::Matrix<Real,struct_sparse>, Complex>;
template struct details::assign_functor<raw::Matrix<Float,struct_sparse>, Complex>;
template struct details::assign_functor<raw::Matrix<Complex,struct_sparse>, Complex >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_sparse>, Complex >;
template struct details::assign_functor<raw::Matrix<Object,struct_sparse>, Complex >;

template struct details::assign_functor<raw::Matrix<Integer,struct_banded>, Complex >;
template struct details::assign_functor<raw::Matrix<Real,struct_banded>, Complex >;
template struct details::assign_functor<raw::Matrix<Float,struct_banded>, Complex >;
template struct details::assign_functor<raw::Matrix<Complex,struct_banded>, Complex >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_banded>, Complex >;
template struct details::assign_functor<raw::Matrix<Object,struct_banded>, Complex >;

template struct details::assign_functor<raw::Matrix<Integer,struct_dense>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Real,struct_dense>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Float,struct_dense>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Complex,struct_dense>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_dense>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Object,struct_dense>, Float_complex >;

template struct details::assign_functor<raw::Matrix<Integer,struct_sparse>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Real,struct_sparse>, Float_complex>;
template struct details::assign_functor<raw::Matrix<Float,struct_sparse>, Float_complex>;
template struct details::assign_functor<raw::Matrix<Complex,struct_sparse>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_sparse>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Object,struct_sparse>, Float_complex >;

template struct details::assign_functor<raw::Matrix<Integer,struct_banded>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Real,struct_banded>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Float,struct_banded>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Complex,struct_banded>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Float_complex,struct_banded>, Float_complex >;
template struct details::assign_functor<raw::Matrix<Object,struct_banded>, Float_complex >;

template struct details::assign_functor_scal<Object, Integer >;
template struct details::assign_functor_scal<Object, Real >;
template struct details::assign_functor_scal<Object, Float >;
template struct details::assign_functor_scal<Object, Complex >;
template struct details::assign_functor_scal<Object, Float_complex >;
template struct details::assign_functor_scal<Object, Object >;
template struct details::assign_functor_scal<Complex, Integer >;
template struct details::assign_functor_scal<Complex, Real >;
template struct details::assign_functor_scal<Complex, Float >;
template struct details::assign_functor_scal<Complex, Object >;
template struct details::assign_functor_scal<Complex, Complex >;
template struct details::assign_functor_scal<Complex, Float_complex >;
template struct details::assign_functor_scal<Float_complex, Integer >;
template struct details::assign_functor_scal<Float_complex, Real >;
template struct details::assign_functor_scal<Float_complex, Float >;
template struct details::assign_functor_scal<Float_complex, Object >;
template struct details::assign_functor_scal<Float_complex, Complex >;
template struct details::assign_functor_scal<Float_complex, Float_complex >;
template struct details::assign_functor_scal<Real, Integer >;
template struct details::assign_functor_scal<Real, Real >;
template struct details::assign_functor_scal<Real, Float >;
template struct details::assign_functor_scal<Real, Object >;
template struct details::assign_functor_scal<Real, Complex >;
template struct details::assign_functor_scal<Real, Float_complex >;
template struct details::assign_functor_scal<Float, Integer >;
template struct details::assign_functor_scal<Float, Real >;
template struct details::assign_functor_scal<Float, Float >;
template struct details::assign_functor_scal<Float, Object >;
template struct details::assign_functor_scal<Float, Complex >;
template struct details::assign_functor_scal<Float, Float_complex >;
template struct details::assign_functor_scal<Integer, Integer >;
template struct details::assign_functor_scal<Integer, Real >;
template struct details::assign_functor_scal<Integer, Float >;
template struct details::assign_functor_scal<Integer, Object >;
template struct details::assign_functor_scal<Integer, Complex >;
template struct details::assign_functor_scal<Integer, Float_complex >;

// get_rep_functor

template struct details::get_rep_functor<Integer>;
template struct details::get_rep_functor<Real>;
template struct details::get_rep_functor<Float>;
template struct details::get_rep_functor<Complex>;
template struct details::get_rep_functor<Float_complex>;

}

#define MACRO_INST_CONSTRUCTOR(type,arg,arg2)	\
    template struct matcl::details::constructor_helper<type>;

MACRO_FOREACH_CODE_ALL(MACRO_INST_CONSTRUCTOR,,)

//MACRO_INST_CONSTRUCTOR(bool,,)
MACRO_INST_CONSTRUCTOR(unsigned char,,)
MACRO_INST_CONSTRUCTOR(unsigned short,,)
MACRO_INST_CONSTRUCTOR(unsigned int,,)
MACRO_INST_CONSTRUCTOR(unsigned long,,)
MACRO_INST_CONSTRUCTOR(signed char,,)
MACRO_INST_CONSTRUCTOR(signed short,,)
MACRO_INST_CONSTRUCTOR(signed int,,)
MACRO_INST_CONSTRUCTOR(signed long,,)
MACRO_INST_CONSTRUCTOR(float,,)
MACRO_INST_CONSTRUCTOR(double,,)
MACRO_INST_CONSTRUCTOR(long double,,)
MACRO_INST_CONSTRUCTOR(matcl::Complex,,)
MACRO_INST_CONSTRUCTOR(matcl::Float_complex,,)
#ifdef _MSC_VER
    MACRO_INST_CONSTRUCTOR(matcl::Integer,,)
    MACRO_INST_CONSTRUCTOR(matcl::Real,,)
    MACRO_INST_CONSTRUCTOR(matcl::Float,,)
#endif
MACRO_INST_CONSTRUCTOR(matcl::Object,,)

namespace matcl
{
    template struct details::get_scalar_functor<Integer>;
    template struct details::get_scalar_functor<Real>;
    template struct details::get_scalar_functor<Float>;
    template struct details::get_scalar_functor<Complex>;
    template struct details::get_scalar_functor<Float_complex>;
    template struct details::get_scalar_functor<Object>;
    
    template MATCL_MATREP_EXPORT const Integer*            Matrix::get_array<Integer>() const;
    template MATCL_MATREP_EXPORT const Real*               Matrix::get_array<Real>() const;
    template MATCL_MATREP_EXPORT const Float*              Matrix::get_array<Float>() const;
    template MATCL_MATREP_EXPORT const Complex*            Matrix::get_array<Complex>() const;
    template MATCL_MATREP_EXPORT const Float_complex*      Matrix::get_array<Float_complex>() const;
    template MATCL_MATREP_EXPORT const Object*             Matrix::get_array<Object>() const;

    template MATCL_MATREP_EXPORT Integer*                  Matrix::get_array_unique<Integer>();
    template MATCL_MATREP_EXPORT Real*                     Matrix::get_array_unique<Real>();
    template MATCL_MATREP_EXPORT Float*                    Matrix::get_array_unique<Float>();
    template MATCL_MATREP_EXPORT Complex*                  Matrix::get_array_unique<Complex>();
    template MATCL_MATREP_EXPORT Float_complex*            Matrix::get_array_unique<Float_complex>();
    template MATCL_MATREP_EXPORT Object*                   Matrix::get_array_unique<Object>();

    template MATCL_MATREP_EXPORT md::constructor_helper<raw::sparse_matrix_base<Integer>>;
    template MATCL_MATREP_EXPORT md::constructor_helper<raw::sparse_matrix_base<Float>>;
    template MATCL_MATREP_EXPORT md::constructor_helper<raw::sparse_matrix_base<Real>>;
    template MATCL_MATREP_EXPORT md::constructor_helper<raw::sparse_matrix_base<Complex>>;
    template MATCL_MATREP_EXPORT md::constructor_helper<raw::sparse_matrix_base<Float_complex>>;
    template MATCL_MATREP_EXPORT md::constructor_helper<raw::sparse_matrix_base<Object>>;
};

#pragma warning( pop )
