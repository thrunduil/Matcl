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

#pragma once

#include "matcl-matrep/general/config.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/func/raw/disp.h"
#include "matcl-core/matrix/matrix_traits.h"
#include "matcl-matrep/details/debug.h"
#include "matcl-matrep/container/container_pool.h"
#include "matcl-core/IO/archive.h"
#include "matcl-internals/container/mat_b.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-matrep/matrix/matrix.h"

namespace matcl { namespace details
{

class matrix_container_base : public refcount_str_default
{
    protected:
        using refcount_str      = matcl::details::refcount_str_default;

    private:
        matcl::mat_code         m_type;

    protected:
        matrix_container_base(matcl::mat_code type);        

    private:
        matrix_container_base(const matrix_container_base&);
        matrix_container_base& operator=(const matrix_container_base&);

        friend class matcl::Matrix;

    public:
        virtual void            destroy_data() = 0;
        virtual void            mark_unique(bool unique) = 0;
        virtual bool            is_effective_unique() const = 0;

        matcl::mat_code         get_type() const        { return m_type; };

        Integer                 rows() const; 
        Integer                 cols() const;
        Integer                 structural_nnz() const;
        Integer                 structural_ldiags(bool use_flags) const;
        Integer                 structural_udiags(bool use_flags) const;
        int_tup_2               size() const;
        bool                    is_scalar_true() const;
        ti::ti_object           get_type_info() const;
        matcl::Matrix           get_scalar() const;

        const struct_flag&      get_struct() const;
        struct_flag&            get_struct();
        void                    set_struct(const struct_flag&) const;
        void                    add_struct(const struct_flag&) const;

        bool                    is_same_matrix(const matrix_container_base& other) const;        

        matcl::Matrix           reserve(Integer r, Integer c) const;
        matcl::Matrix           resize(Integer r, Integer c) const;
        matcl::Matrix           reserve_band(Integer r, Integer c, Integer fd, Integer ld) const;
        matcl::Matrix           resize_band(Integer r, Integer c, Integer fd, Integer ld) const;

        void                    disp(const disp_stream_ptr& os) const;

        void                    serialize(oarchive_impl&, const unsigned int ) const    {};
        void                    serialize(iarchive_impl&, const unsigned int )          { m_refcount = 0; };   
};

template<class val_type, class str_type>
class matrix_container : public matrix_container_base
{
    public:
        using Matrix = raw::Matrix<val_type,str_type>;

    private:
        Matrix    m_matrix;    

        matrix_container()
            :matrix_container_base(matrix_traits::mat_type_info_type_2<val_type,str_type>::matrix_code)
            ,m_matrix(ti::get_object_ti<val_type>())
        {};

        friend class matcl::Matrix;
        friend class boost::serialization::access;

    public:        
        matrix_container(const Matrix& mat);
        matrix_container(Matrix&& mat);

        matrix_container(const matrix_container&) = delete;
        matrix_container& operator=(const matrix_container&) = delete;

        virtual void            destroy_data() override     { m_matrix.destroy_data(); };
        virtual void            mark_unique(bool unique) override
                                                            { m_matrix.mark_unique(unique); };
        virtual bool            is_effective_unique() const override
                                                            { return m_matrix.is_unique(); };
        void                    call_destructor()           { m_matrix.~Matrix(); };

        Matrix&                 get()                       { return m_matrix;};
        const Matrix&           get() const                 { return m_matrix;};    

        Integer                 rows_impl() const           { return m_matrix.rows(); }; 
        Integer                 cols_impl() const           { return m_matrix.cols(); };
        Integer                 struct_ldiags_impl(bool use_flag) const;
        Integer                 struct_udiags_impl(bool use_flag) const;
        int_tup_2               size_impl() const           { return int_tup_2(rows(),cols()); };        
        Integer                 struct_nnz_impl() const     { return m_matrix.nnz(); };         
        ti::ti_type<val_type>   get_ti_impl() const         { return m_matrix.get_type(); };
        const struct_flag&      get_struct_impl() const     { return m_matrix.get_struct(); };
        struct_flag&            get_struct_impl()           { return m_matrix.get_struct(); };

        bool                    is_scalar_true_impl() const;
        void                    disp_impl(const disp_stream_ptr& os) const;
        void                    set_struct_impl(const struct_flag& fl) const;
        void                    add_struct_impl(const struct_flag& fl) const;

        bool                    is_same_matrix_impl(const matrix_container_base* other) const; 

        matcl::Matrix           get_elem(Integer i) const;
        matcl::Matrix           get_elem(Integer i,Integer j) const;
        matcl::Matrix           get_scalar_impl() const;

        Matrix                  reserve_impl(Integer r, Integer c) const;
        Matrix                  resize_impl(Integer r, Integer c) const;

        void *                  operator new(size_t size);
        void                    operator delete(void *p, size_t size);

        void serialize(oarchive_impl & ar, const unsigned int ) const
        {
            matrix_container& tmp = const_cast<matrix_container&>(*this);
            ar << boost::serialization::base_object<matrix_container_base>(tmp);
            ar << m_matrix;
        };

        void serialize(iarchive_impl & ar, const unsigned int )
        {
            ar >> boost::serialization::base_object<matrix_container_base>(*this);
            ar >> m_matrix;
        };
};

template<class value_type,class struct_type>
inline matrix_container_base* 
create_matrix_container(const raw::Matrix<value_type,struct_type>& val)
{   
    using container_type    = matrix_container<value_type,struct_type>;
    using matrix_type       = raw::Matrix<value_type,struct_type>;

    void* ptr = container_pool::malloc<(int)type_to_code<matrix_type>::value>();
    ::new (ptr) container_type(val);

    return static_cast<matrix_container_base*>(ptr);
};

template<class value_type,class struct_type>
inline matrix_container_base* 
create_matrix_container(raw::Matrix<value_type,struct_type>&& val)
{   
    using container_type    = matrix_container<value_type,struct_type>;
    using matrix_type       = raw::Matrix<value_type,struct_type>;

    void* ptr = container_pool::malloc<(int)type_to_code<matrix_type>::value>();
    ::new (ptr) container_type(std::move(val));

    return static_cast<matrix_container_base*>(ptr);
};

inline void free_container(matrix_container_base* ptr)
{
    container_pool::free(ptr,(int)ptr->get_type());
};

inline void destroy_container(matrix_container_base* ptr)
{
    ptr->destroy_data();
    free_container(ptr);
};

template<class T,bool is_mat>
struct matrix_container_type
{
    using type = typename details::promote_scalar<T>::type;
};

template<class T>
struct matrix_container_type<T,true>
{
    using type = T;
};

template<class T, bool is_mat>
struct create_container_caller_proxy;

template<class T>
class create_container
{
    public:
        using matrix_type = typename matrix_container_type
                            <T,details::is_matrix<T>::value>::type;

        static details::matrix_container_base* eval(const T& val)
        {
            return eval_impl<details::is_matrix<T>::value>(val);
        };

        static details::matrix_container_base* eval(T&& val)
        {
            return eval_impl<details::is_matrix<T>::value>(std::move(val));
        };

    private:

        template<bool is_mat>
        static details::matrix_container_base* eval_impl(const T& val)
        {
            return create_container_caller_proxy<T, is_mat>::call_helper(val);
        }
    
        template<bool is_mat>
        static details::matrix_container_base* eval_impl(T&& val)
        {
            return create_container_caller_proxy<T, is_mat>::call_helper(std::move(val));
        }

        static details::matrix_container_base* eval_impl_false(const T& val)
        {
            using value_type    = details::promote_scalar<T>::type;
            using Matrix        = matcl::raw::Matrix<value_type,struct_dense>;
            return eval_matrix(Matrix(value_type(val)));
        };

        static details::matrix_container_base* eval_impl_true(const T& val)
        {
            return eval_matrix(val);
        };
        
        static details::matrix_container_base* eval_impl_true(T&& val)
        {
            return eval_matrix(std::move(val));
        };

        template<class value_type, class struct_type>
        static details::matrix_container_base* eval_matrix(const raw::Matrix<value_type,struct_type>& val)
        {
            return details::create_matrix_container<value_type,struct_type>(val);
        };
        
        template<class value_type, class struct_type>
        static details::matrix_container_base* eval_matrix(raw::Matrix<value_type,struct_type>&& val)
        {
            return details::create_matrix_container<value_type,struct_type>(std::move(val));
        };

        template<class K, bool is_mat> 
        friend struct create_container_caller_proxy;
};

template<class T, bool is_mat>
struct create_container_caller_proxy
{
    static details::matrix_container_base * call_helper( const T & val )
    {
        return create_container< T >::eval_impl_false( val );
    }

    static details::matrix_container_base * call_helper(T&& val )
    {
        return create_container< T >::eval_impl_false( val );
    }
};

template<class T>
struct create_container_caller_proxy<T, true>
{
    static details::matrix_container_base * call_helper( const T & val )
    {
        return create_container< T >::eval_impl_true( val );
    }
    
    static details::matrix_container_base * call_helper(T&& val )
    {
        return create_container< T >::eval_impl_true(std::move(val));
    }
};

};};