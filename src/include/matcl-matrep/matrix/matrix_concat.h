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
#include "matcl-matrep/details/mpl.h"
#include "matcl-matrep/details/isa.h"
#include "matcl-matrep/objects/details/type_info_object.h"

namespace matcl
{

#pragma warning(push)
#pragma warning(disable:4251)   //warning C4251: class needs to have dll-interface to be used by clients

// matrix horizontal concatenation helper, i.e. constructor of [A, B, ...]
// instances of this class cannot be shared between threads
class MATCL_MATREP_EXPORT mat_row
{
    private:
        using ref_ptr = std::shared_ptr<details::mat_cons_data>;

    private:
        Integer					m_rows;
        Integer					m_cols;
        Integer					m_nnz;

        ref_ptr					m_data;

    public:
        // construct empty matrix of Integer type
        mat_row();		

        // standard destructor
        ~mat_row();		

        // standard copy constructor and move constructor
        mat_row(const mat_row&);
        mat_row(mat_row&&);

        // standard assignment and move assignment;
        mat_row&				operator=(const mat_row&);
        mat_row&				operator=(mat_row&&);        

        // add scalar to the list of matrices
        template<class S>
        typename details::enable_if<details::is_scalar<S>::value,mat_row&>::type
                                add(S&&);

        // add scalar to the list of matrices; equivalent to add
        template<class S>
        typename details::enable_if<details::is_scalar<S>::value,mat_row&>::type
                                operator,(S&& v)            { return add(std::forward<S>(v)); };

        // add matrix to the list of matrices
        mat_row&				add(const Matrix&);
        mat_row&				add(Matrix&&);

        // add matrix to the list of matrices; equivalent to add
        mat_row&				operator,(const Matrix& v);
        mat_row&				operator,(Matrix&& v);

        // add concatenation helpers to the list of matrices
        mat_row&				add(const mat_row&);
        mat_row&				add(mat_row&&);
        mat_row&				add(const mat_col&);
        mat_row&				add(mat_col&&);

        // add concatenation helpers to the list of matrices; 
        // equivalent to add
        mat_row&				operator,(const mat_row& v);
        mat_row&				operator,(mat_row&& v);
        mat_row&				operator,(const mat_col& v);
        mat_row&				operator,(mat_col&& v);

        // convert to Matrix
        const Matrix			to_matrix() const;

        // make given instance of mat_row unique
        void					make_unique();

        // get type of stored elements
        ti::ti_object           get_type() const;

        // get value type of stored elements
        matcl::value_code       get_value_type() const;        

        friend mat_col;
        template<class M> friend struct details::sparse_matrix_constructor_row;
        template<class M> friend struct details::dense_matrix_constructor_row;
        template<class M> friend struct details::sparse_matrix_constructor_col;
        template<class M> friend struct details::dense_matrix_constructor_col;

    //internal use
    public:
        bool                    is_initialized() const;

    private:
        void					add_mat(const Matrix& mat);
        void					add_col(const mat_col&);
        void					add_row(const mat_row&);
        void					add_int(Integer val);
        void					add_real(Real val);
        void					add_float(Float val);
        void					add_complex(Real re,Real im);
        void					add_fcomplex(Float re, Float im);

        template<class S>
        mat_row&				add_scalar(const S& val);
        template<class S>
        mat_row&				add_scalar(typename std::decay<S>::type&& val);

        template<class S>
        void					add_scalar_impl(const S& val);

        Matrix                  build_matrix_inplace();
        Matrix                  build_matrix_new() const;
};

// matrix vertical concatenation helper, i.e. constructor of [A; B; ...]
// instances of this class cannot be shared between threads
class MATCL_MATREP_EXPORT mat_col
{
    private:
        using ref_ptr = std::shared_ptr<details::mat_cons_data>;

    private:
        Integer					m_rows;
        Integer					m_cols;
        Integer					m_nnz;

        ref_ptr					m_data;

    public:
        // construct empty matrix of Integer type
        mat_col();		

        // standard copy constructor and move constructor
        mat_col(const mat_col&);
        mat_col(mat_col&&);

        // standard assignment and move assignment;
        mat_col&				operator=(const mat_col&);
        mat_col&				operator=(mat_col&&);

        // standard destructor
        ~mat_col();

        // add scalar to the list of matrices
        template<class S>
        typename details::enable_if<details::is_scalar<S>::value,mat_col&>::type
                                add(S&&);

        // add scalar to the list of matrices; equivalent to add
        template<class S>
        typename details::enable_if<details::is_scalar<S>::value,mat_col&>::type
                                operator,(S&& v)            { return add(std::forward<S>(v)); };

        // add matrix to the list of matrices
        mat_col&				add(const Matrix&);
        mat_col&				add(Matrix&&);

        // add matrix to the list of matrices; equivalent to add
        mat_col&				operator,(const Matrix& v);
        mat_col&				operator,(Matrix&& v);

        // add concatenation helpers to the list of matrices
        mat_col&				add(const mat_row&);
        mat_col&				add(mat_row&&);
        mat_col&				add(const mat_col&);
        mat_col&				add(mat_col&&);

        // add concatenation helpers to the list of matrices;
        // equivalent to add
        mat_col&				operator,(const mat_row& v) { return add(v); };
        mat_col&				operator,(mat_row&& v)      { return add(std::move(v)); };
        mat_col&				operator,(const mat_col& v) { return add(v); };
        mat_col&				operator,(mat_col&& v)      { return add(std::move(v)); };    

        // convert to Matrix
        const Matrix			to_matrix() const;

        // make given instance of mat_row unique
        void					make_unique();

        // get type of stored elements
        ti::ti_object           get_type() const;

        // get value type of stored elements
        matcl::value_code       get_value_type() const;     

    //internal use
    public:
        bool                    is_initialized() const;

        friend mat_row;
        template<class M> friend struct details::sparse_matrix_constructor_row;
        template<class M> friend struct details::dense_matrix_constructor_row;
        template<class M> friend struct details::sparse_matrix_constructor_col;
        template<class M> friend struct details::dense_matrix_constructor_col;

    private:
        void					add_mat(const Matrix& mat);
        void					add_col(const mat_col&);
        void					add_row(const mat_row&);
        void					add_int(Integer val);
        void					add_real(Real val);
        void					add_float(Float val);
        void					add_complex(Real re, Real im);
        void					add_fcomplex(Float re, Float im);

        template<class S>
        mat_col&				add_scalar(const S& val);
        template<class S>
        mat_col&				add_scalar(typename std::decay<S>::type&& val);
        template<class S>
        void					add_scalar_impl(const S& val);

        Matrix                  build_matrix_inplace();
        Matrix                  build_matrix_new() const;
};

#pragma warning(pop)

};