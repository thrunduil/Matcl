/*
 *  This file is a part of Matrix Computation Library (MATCL)
 *
 *  Copyright (c) Pawe³ Kowal 2017 - 2021
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
#include "matcl-matrep/details/refcount.h"
#include "matcl-core/matrix/complex_type.h"
#include "matcl-scalar/object.h"
#include "matcl-matrep/objects/details/type_info_object.h"

namespace matcl { namespace details
{
    template<class M>   struct get_functor;

    template<class Mat>
    class rvalue_holder
    {
        public:
            using matrix_type       = Mat;

        private:
	        using container_type    =  matrix_container_base;
            friend struct get_functor<Mat>;

        private:
            container_type*         m_container;

        private:
            rvalue_holder(container_type* cont)             : m_container(cont) {};
            rvalue_holder(const rvalue_holder&)             = delete;
            rvalue_holder& operator=(const rvalue_holder&)  = delete;

        public:
            ~rvalue_holder();

            rvalue_holder(rvalue_holder&& other);
            rvalue_holder& operator=(rvalue_holder&&);

            //can be casted to rvalue reference. Reference is invalidated when rvalue_holder
            //is destroyed.
            matrix_type&                get();
    };

    template<class M>
    struct get_functor
    {
        static MATCL_MATREP_EXPORT const M&     eval(const matrix_base& mat);
        static MATCL_MATREP_EXPORT M&           eval(matrix_base& mat);
        static MATCL_MATREP_EXPORT rvalue_holder<M>
                                                eval_move(matrix_base& mat);
        static MATCL_MATREP_EXPORT M&           get_ref(matrix_container_base*);
        static MATCL_MATREP_EXPORT void         destroy(matrix_container_base*);
    };

    template<class M>
    struct get_scalar_functor
    {
        static MATCL_MATREP_EXPORT M           eval_get(const matrix_base& mat);
        static MATCL_MATREP_EXPORT M&          eval_get_unique(matcl::Matrix& mat);
    };
    
    template<class T>
    struct constructor_helper
    {
        static MATCL_MATREP_EXPORT 
        bool    eval(const T& val,bool allow_conversions,matrix_base& mat);

        static MATCL_MATREP_EXPORT 
        bool    eval(T&& val,bool allow_conversions,matrix_base& mat);
    };

    class MATCL_MATREP_EXPORT matrix_base
    {
        public:
            using matrix_container      = details::matrix_container_base;
            using refcount_str          = refcount_str_default;

            union value
            {
                struct mat_info
                {
                    matrix_container*   m_mat_ptr;					
                    refcount_str*       m_refcount;
                }                       m_mat;
        
                Integer				    val_int;
                Real				    val_real;
                Float				    val_float;
                Real				    val_complex[2];
                Float				    val_fcomplex[2];
                char                    val_object[sizeof(Object)];
            };
            
            mutable value		        m_value;
            mutable mat_code	        m_type;

        public:
            Object&             get_object()        { return *reinterpret_cast<Object*>(m_value.val_object); };
            const Object&       get_object() const  { return *reinterpret_cast<const Object*>(m_value.val_object); };
            
            matrix_base(){}
            ~matrix_base();

            matrix_base(const matrix_base& other);
            matrix_base(matrix_base&& other);

            matrix_base&        operator=(const matrix_base& other);
            matrix_base&        operator=(matrix_base&& other);
            void                reset(matrix_base& other);

            static bool         is_same_matrix(const matrix_base& A1, const matrix_base& A2);
            bool                is_effective_unique_mat() const;

        private:
            void                increase_refcount() const;
            void                decrease_refcount() const;
            bool                check_effective_unique() const;

        protected:
            void                assign(const matrix_base& other) const;
            void                assign(matrix_base&& other) const;

            Object&             get_object_nc() const  { return *reinterpret_cast<Object*>(m_value.val_object); };            
    };
};

inline matcl::details::matrix_base::matrix_base(const matrix_base& other)
    :m_type(other.m_type),m_value(other.m_value)
{
    increase_refcount();
    if (m_type == mat_code::object_scalar)
        new (m_value.val_object) Object(other.get_object());
};

inline matcl::details::matrix_base::matrix_base(matrix_base&& other)
    :m_type(std::move(other.m_type)),m_value(std::move(other.m_value))

{
    other.m_type = mat_code::integer_scalar;
};

inline matcl::details::matrix_base& 
matcl::details::matrix_base::operator=(const matrix_base& other)
{
    assign(other);
    return *this;
};

inline void matcl::details::matrix_base::assign(const matrix_base& other) const
{
    other.increase_refcount();
    decrease_refcount();

    if (m_type == mat_code::object_scalar || other.m_type == mat_code::object_scalar)
    {
        if (m_type == mat_code::object_scalar && other.m_type == mat_code::object_scalar)
        {
            get_object_nc().reset(other.get_object());
            return;
        };
    
        if (m_type == mat_code::object_scalar)
        {
            get_object().~Object();
            m_type = other.m_type;
            m_value = other.m_value;
            return;
        }

        new(&get_object_nc()) Object(other.get_object());
        m_type = other.m_type;

        return;
    };

    m_type = other.m_type;
    m_value = other.m_value;
    return;
};

inline matcl::details::matrix_base&
matcl::details::matrix_base::operator=(matrix_base&& other)
{
    assign(std::move(other));
    return *this;
};

inline void matcl::details::matrix_base::assign(matrix_base&& other) const
{
    std::swap(m_value, other.m_value);
    std::swap(m_type, other.m_type);
};

inline matcl::details::matrix_base::~matrix_base()
{
    decrease_refcount();

    if (m_type == mat_code::object_scalar)
        get_object().~Object();
};

inline void matcl::details::matrix_base::reset(matrix_base& other)
{
    m_value         = other.m_value;
    m_type          = other.m_type;
    other.m_type    = mat_code::integer_scalar;
};

inline bool matcl::details::matrix_base::is_effective_unique_mat() const
{
    if (m_value.m_mat.m_refcount->is_unique())
        return true;

    return check_effective_unique();
};

template<class Mat>
inline matcl::details::rvalue_holder<Mat>::~rvalue_holder()
{ 
    get_functor<Mat>::destroy(m_container); 
}

template<class Mat>
inline matcl::details::rvalue_holder<Mat>::rvalue_holder(rvalue_holder&& other)
    :m_container(other.m_container)
{
    other.m_container = nullptr;
}

template<class Mat>
inline matcl::details::rvalue_holder<Mat>& 
    matcl::details::rvalue_holder<Mat>::operator=(rvalue_holder&& other)
{
    std::swap(m_container, other.m_container);
    return *this;
};

template<class Mat>
inline Mat& matcl::details::rvalue_holder<Mat>::get()
{ 
    return get_functor<Mat>::get_ref(m_container); 
}

};

