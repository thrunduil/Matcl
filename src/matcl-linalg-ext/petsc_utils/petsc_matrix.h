/*
 *  This file is a part of Morfa Matrix Lib.
 *
 *  Copyright (c) Pawe³ Kowal 2011-2016
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

#if 0
TODO

#include "matcl-linalg/general/config_linalg.h"
#include "matcl-matrep/matrix/matrix.h"

#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include "matcl-linalg/petsc_utils/petsc_library.h"
#include "matcl-linalg/petsc_utils/petsc_objects.h"

/*
#include "matcl-linalg/linear_eq/ksp_solver.h"

#include <memory>

#include "matcl-linalg/petsc_utils/petsc_option.h"
*/

namespace matcl { namespace details
{

namespace mp = matcl::petsc;

class petsc_matrix
{
    private:
        using ptr_type      = std::shared_ptr<petsc_matrix>;
        using lock_helper   = matcl::petsc::petsc_lock_helper;

    public:
        virtual ~petsc_matrix() {};

        bool                all_finite() const      { return m_is_finite;}
        mp::smart_mat       get_Aksp() const        { return m_Aksp;}
        Integer             rows() const            { return m_rows; };
        Integer             cols() const            { return m_cols; };
        value_code          get_value_code() const  { return m_value_code;  };
        virtual bool        is_shell() const = 0;        

        virtual linear_operator get_linop() const = 0;
        virtual Matrix          get_matrix() const = 0;
        virtual Matrix          mmul_right(const Matrix& X, trans_type tA) const = 0;

        static ptr_type     create(const Matrix& A_in);
        static ptr_type     create(const linear_operator& A_in);

    protected:
        petsc_matrix(const Matrix& A_in);
        petsc_matrix(const linear_operator& A_in);

    private:
        petsc_matrix(const petsc_matrix& other) = delete;
        petsc_matrix& operator=(const petsc_matrix& other) = delete;

        lock_helper         m_lock;
        Integer             m_rows;
        Integer             m_cols;
        value_code          m_value_code;
        bool                m_is_finite;

    protected:
        mp::smart_mat       m_Aksp;
};

template <class S>
class petsc_matrix_str : public petsc_matrix
{};

template <>
class petsc_matrix_str<struct_dense> : public petsc_matrix
{
    private:
        raw::real_dense Araw;

    public:
        petsc_matrix_str(const Matrix& A_in);

        virtual ~petsc_matrix_str() {}

        virtual bool            is_shell() const override   { return false; }; 
        virtual linear_operator get_linop() const override  { return Matrix(Araw,false); };
        virtual Matrix          mmul_right(const Matrix& X, trans_type tA) const override;        
        virtual Matrix          get_matrix() const override { return Matrix(Araw,false); };
};

template <>
class petsc_matrix_str<struct_sparse> : public petsc_matrix
{
    private:
        using Mat           = raw::Matrix<Real,struct_sparse>;
        using Mat_I         = raw::Matrix<Integer,struct_dense>;

    private:
        matcl::Matrix       m_A;
        Mat                 m_At;
        Mat_I               m_column_indices;

    public:
        petsc_matrix_str(const Matrix& A_in);

        virtual ~petsc_matrix_str() {}        

        virtual bool            is_shell() const override   { return false; };
        virtual linear_operator get_linop() const override;
        virtual Matrix          get_matrix() const override;
        virtual Matrix          mmul_right(const Matrix& X, trans_type tA) const override;
};

class petsc_matrix_linop : public petsc_matrix
{
    private:
        linear_operator     m_Araw;

    public:
        petsc_matrix_linop(const linear_operator& A_in);

        virtual ~petsc_matrix_linop() {}        

        virtual bool            is_shell() const override   { return true; };
        virtual linear_operator get_linop() const override  { return m_Araw; };
        virtual Matrix          get_matrix() const override;
        virtual Matrix          mmul_right(const Matrix& X, trans_type tA) const override;
};
 
}};

#endif