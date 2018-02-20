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

#include "matcl-matrep/details/fwd_decls.h"
#include "matcl-matrep/details/refcount.h"
#include "matcl-scalar/object.h"
#include "matcl-matrep/objects/details/type_info_object.h"
#include "matcl-matrep/matrix/struct_flag.h"
#include "matcl-core/IO/archive.h"

namespace matcl { namespace raw 
{

namespace md = matcl::details;

namespace details
{
    template<class Val>
    class MATCL_MATREP_EXPORT root_ptr
    {
        private:
            Val*        m_ptr;
            bool        m_foreign;

        public:
            root_ptr();

            void        set(Val* arr);
            void        set_foreign(Val* arr);
            void        destroy(Integer max_rows, Integer max_cols);
            ptrdiff_t   offset(const Val* ptr) const; //return ptr - root_ptr

            root_ptr(const root_ptr&);
            root_ptr&   operator=(const root_ptr&);
    };

    template<class val_type>
    class MATCL_MATREP_EXPORT mat_data
    {
        public:
            using type_info     = ti::ti_type<val_type>;
            using refcount_str  = matcl::details::refcount_str_default;

            mat_data(type_info ti);
            mat_data(type_info ti, Integer r, Integer c);

            //does not take ownership of arr
            mat_data(type_info ti, Integer r, Integer c, val_type* arr, Integer ld);

            ~mat_data();

            mat_data(const mat_data&);
            mat_data(mat_data&&);
            mat_data(const mat_data&&) = delete; 

            val_type*           ptr()                       { return m_ptr; }
            const val_type*     ptr() const                 { return m_ptr; }

            void                mark_unique(bool unique);
            bool                is_unique() const;
            refcount_str*       get_refstr() const          { return m_refcount; };

            type_info           get_type() const            { return m_ti; };

            const struct_flag&  get_struct() const          { return m_flag; };
            struct_flag&        get_struct()                { return m_flag; };

            void                reset_unique();
            void                reset_unique(Integer r, Integer c);
            void                destroy_data();

            mat_data&           assign_to_fresh(const mat_data&);
            mat_data&           assign_to_fresh(mat_data&&);

            void                serialize(oarchive_impl & ar, const unsigned int version) const;
            void                serialize(iarchive_impl & ar, const unsigned int version);          

        public:
            refcount_str*       m_refcount;
            Integer             m_rows;
            Integer             m_cols;
            Integer             m_max_rows;
            Integer             m_max_cols;
            Integer             m_ld;
            Integer             m_size;
            val_type*           m_ptr;
            root_ptr<val_type>  m_root_ptr;
            type_info           m_ti;  
            struct_flag         m_flag;
            bool                m_effective_unique;

        private:
            void                alloc(Integer s);
            void                take_foreign(val_type* arr);

        private:
            mat_data&           operator=(const mat_data&) = delete;
            mat_data&           operator=(mat_data&&) = delete;
    };
};

template<class value_type_, class struct_type_>
class MATCL_MATREP_EXPORT Matrix
{
    public:
        using value_type    = value_type_;
        using struct_type   = struct_type_;
};

template <class value_type_>
class MATCL_MATREP_EXPORT dense_matrix_base
{
    public:
        using value_type    = value_type_;
        using struct_type   = struct_dense;
        using refcount_str  = matcl::details::refcount_str_default;

    protected:
        using tinfo         = ti::ti_type<value_type>;

    public:
        dense_matrix_base(tinfo ti);
        dense_matrix_base(tinfo ti, Integer r, Integer c);
        //does not take ownership of arr
        dense_matrix_base(tinfo ti, Integer r, Integer c, value_type* arr, Integer ld);
        dense_matrix_base(const dense_matrix_base&);
        dense_matrix_base(dense_matrix_base&&);
        dense_matrix_base(const dense_matrix_base&&) = delete;

        //matrix must be unique
        dense_matrix_base&      reset_unique();
        dense_matrix_base&      reset_unique(Integer r, Integer c);

        const value_type*       ptr() const                     { return m_data.ptr(); };
        value_type*             ptr()                           { return m_data.ptr(); };

        void                    mark_unique(bool unique)        { m_data.mark_unique(unique); };        
        bool                    is_unique() const               { return m_data.is_unique(); };
        refcount_str*           get_refstr() const              { return m_data.get_refstr(); };

        const struct_flag&      get_struct() const              { return m_data.m_flag; };
        struct_flag&            get_struct()                    { return m_data.m_flag; };
        void                    set_struct(struct_flag f) const { m_data.m_flag.set(f); };
        void                    add_struct(struct_flag f) const { m_data.m_flag.add(f); };
        tinfo                   get_type() const;        

        Integer                 rows() const                    { return m_data.m_rows; }; 
        Integer                 cols() const                    { return m_data.m_cols; };
        Integer                 length() const                  { return (rows() == 0 || cols() == 0)? 0 : std::max(rows(), cols()); };
        Integer                 max_rows() const                { return m_data.m_max_rows; };
        Integer                 max_cols() const                { return m_data.m_max_cols; };
        Integer                 ld() const                      { return m_data.m_ld; };
        Integer                 size() const                    { return m_data.m_size; };
        Integer                 nnz() const;
        void                    read_from(const value_type* arr, Integer ld);

        void                    serialize(oarchive_impl & ar, const unsigned int version) const;
        void                    serialize(iarchive_impl & ar, const unsigned int version);    

    protected:
        details::mat_data<value_type>   m_data;

        dense_matrix_base&      assign_to_fresh(const dense_matrix_base&);
        dense_matrix_base&      assign_to_fresh(dense_matrix_base&&);

    private:
        dense_matrix_base&      operator=(const dense_matrix_base&) = delete;
        dense_matrix_base&      operator=(dense_matrix_base&&) = delete;
};

template<class V, class S>
bool is_real_matrix(const Matrix<V,S>&) 
{
    return md::is_float_real_scalar<V>::value || std::is_same<V,Integer>::value;
};

};};
