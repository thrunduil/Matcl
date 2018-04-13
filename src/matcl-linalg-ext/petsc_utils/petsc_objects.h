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

#include "matcl-linalg/general/config_linalg.h"

#pragma warning(push)
#pragma warning(disable:4100)
#pragma warning(disable:4127)

#include "petsc.h"
#include "petscmat.h"

#pragma warning(pop)

#include "matcl-matrep/matrix/matrix.h"
#include "matcl-linalg/details/linalg_fwd.h"

#include "matcl-matrep/matrix/permvec.h"
#include "matcl-linalg/linear_eq/linsolve_object.h"
#include "matcl-internals/container/mat_s.h"
#include "matcl-internals/container/mat_d.h"
#include "matcl-internals/container/mat_b.h"

#include <memory>
#include <cassert>
#include <vector>

namespace matcl { namespace petsc
{

class smart_vec;

//-----------------------------------------------------------------
//              fully defined wrappers
//-----------------------------------------------------------------

/// class to wrap around a Petsc MatNullSpace
class petsc_null_mat
{
    private:
        using ptr_type      = std::shared_ptr<_p_MatNullSpace>;
        using this_ptr      = std::shared_ptr<petsc_null_mat>;
        using vec_vec       = std::vector<smart_vec>;
        using vec_raw_vec   = std::vector<Vec>;

    public:                
        ~petsc_null_mat();
        
                            operator MatNullSpace() const;
        static this_ptr     create(const Matrix& mat, bool with_const);

    private:
        ptr_type            m_impl;
        Matrix              m_data;
        vec_vec             m_smart_vectors;
        vec_raw_vec         m_vectors;

        static void         deleter(MatNullSpace& to_destroy);
};

/// class to wrap around a Petsc Mat, only matrix structure is stored
/// sparse matrices are not transposed
class petsc_sparse_struct
{
    private:
        using ptr_type      = std::shared_ptr<_p_Mat>;
        using this_ptr      = std::shared_ptr<petsc_sparse_struct>;
        using Mat_I         = raw::Matrix<Integer,struct_dense>;

    public:        
        ~petsc_sparse_struct();
        
                            operator Mat() const;
        bool                need_transpose() const;

        /// if const_access = true, then it is assumed that matrix will not
        /// be changed (and make_unique will not be called)
        static this_ptr     create(const Matrix& mat, bool const_access);

    private:
        ptr_type            m_impl;
        Matrix              m_data;
        bool                m_need_trans;
        Mat_I               m_column_indices;

        petsc_sparse_struct(const Matrix& mat, bool const_access);
        static void         deleter(Mat& to_destroy);
        void                initialize_sparse(bool const_access);
};

/// implement PCSHELL from linsolve_obj
class petsc_pcshell
{
    private:
        using ptr_type      = std::shared_ptr<_p_PC>;
        using this_ptr      = std::shared_ptr<petsc_pcshell>;

    public:                
        ~petsc_pcshell();
        
                            operator PC() const;
        static this_ptr     create(PC precond, const linsolve_obj& ls);
        static this_ptr     create(PC precond, const linear_operator& ls);

    private:
        ptr_type            m_impl;
        linsolve_obj        m_data_linsolve;
        linear_operator     m_data_linop;

        petsc_pcshell(const linsolve_obj& ls);
        petsc_pcshell(const linear_operator& ls);

        static void         deleter(PC& to_destroy);

        static PetscErrorCode   func_apply_linop(PC shell, Vec T,Vec AT);
        static PetscErrorCode   func_apply_linsolve(PC shell, Vec T,Vec AT);
        static PetscErrorCode   func_apply_trans_linop(PC shell, Vec T,Vec AT);
        static PetscErrorCode   func_apply_trans_linsolve(PC shell, Vec T,Vec AT);
        static PetscErrorCode   func_apply_trans_linsolve_impl(PC shell, Vec T,Vec AT, trans_type t);
        static PetscErrorCode   func_apply_trans_linop_impl(PC shell, Vec T,Vec AT, trans_type t);
};

/// class to wrap around a Petsc IS
class petsc_IS
{
    private:
        using ptr_type      = std::shared_ptr<_p_IS>;
        using this_ptr      = std::shared_ptr<petsc_IS>;
        using Mat_I         = raw::Matrix<Integer,struct_dense>;

    public:        
        ~petsc_IS();
        
                            operator IS() const;
        static this_ptr     create(const raw::integer_dense& mat, Integer max_val);

    private:
        ptr_type            m_impl;

        petsc_IS(const raw::integer_dense& mat, Integer max_val);
        static void         deleter(IS& to_destroy);
};

//-----------------------------------------------------------------
//                  low level wrappers
//-----------------------------------------------------------------

/// class to wrap around a Petsc C-object
///
/// handles the underlying object acording to RAII principle and should be exception safe
/// see petsc documentation for info on particular "Create" functions
 
class smart_vec
{
    private:
        std::shared_ptr<_p_Vec> m_v;

    public:
        operator Vec() const;

        static smart_vec    vec_create_seq_with_array(MPI_Comm comm, PetscInt bs, PetscInt n, 
                                               const PetscScalar arr[]);
        static smart_vec    vec_create_seq(MPI_Comm comm,PetscInt n);

    private:
        static void         deleter(Vec& to_destroy);
};

/// class to wrap around a Petsc C-object
class smart_mat
{
    private:
        using impl_type = std::shared_ptr<_p_Mat>;
        using linop_ptr = std::shared_ptr<linear_operator>;

    public:
        operator Mat() const;

        static smart_mat    create(MPI_Comm comm);
        static smart_mat    mat_create_seq_dense(MPI_Comm comm, PetscInt m, PetscInt n, PetscInt ld, 
                                PetscScalar *data);
        static smart_mat    mat_create_shell(MPI_Comm comm, const linear_operator& op);
        static smart_mat    mat_create_seq_aij_with_arrays(MPI_Comm comm, PetscInt m, PetscInt n, PetscInt *i,
                                PetscInt *j, PetscScalar *a);

        void                set_struct_flags(const Matrix& mat);

    private:
        static void             deleter(Mat& to_destroy);
        static PetscErrorCode   mat_mult(Mat shell, Vec T,Vec AT);
        static PetscErrorCode   mat_mult_trans(Mat shell, Vec T,Vec AT);
        static PetscErrorCode   mat_mult_ctrans(Mat shell, Vec T,Vec AT);
        static PetscErrorCode   mat_mult_trans_impl(Mat shell, Vec T,Vec AT, trans_type t);

    private:
        impl_type           m_impl;
        linop_ptr           m_linop;
};

/// class to wrap around a Petsc C-object
/// Handles the underlying object acording to RAII principle and should be exception safe
/// See Petsc documentation for info on particular "Create" functions    
class smart_KSP
{
    private:
        using ptr_type  = std::shared_ptr<smart_KSP>;
        using impl_type = std::shared_ptr<_p_KSP>;

    public:        
                            operator KSP() const;
        static smart_KSP    create(MPI_Comm comm);

    private:
        impl_type           m_impl;

        static void         destroy(KSP& to_destroy);
};

//-----------------------------------------------------------------
//                  converters
//-----------------------------------------------------------------

/// produce a column vector and copy contents of a Petsc Vector to it.
/// X: petsc vector to copy contents from
/// return created matrix from the petsc vector
matcl::Matrix petsc_to_mmlib(const Vec& X);

/// copy contents of an mmlib matrix to a created (size-conformant) Petsc vector
/// source: source column or row matrix
/// dest:  petsc vector with matching size
void mmlib_to_petsc(const matcl::Matrix& source, Vec& dest);

/// create Petsc MatNullSpace from mmlib Matrix
void mmlib_to_petsc_null(const matcl::Matrix& source, std::shared_ptr<petsc_null_mat>& dest, 
                         bool with_const);

/// create Petsc Mat containing only structure of the from mmlib Matrix
/// if const_access = true, then it is assumed that matrix will not
/// be changed (and make_unique will not be called)
void mmlib_to_sparse_struct(const matcl::Matrix& source, bool const_access,
                         std::shared_ptr<petsc_sparse_struct>& dest);

/// convert PETSC permutations to permvec
void IS_to_perm_vector(IS p, permvec& mp);

/// convert PETSC permutations to integer matrix
void IS_to_int_mat(IS p, Matrix& mp);

/// convert mmlib integer matrix to PETSC IS; max_val: maximum possible value
/// of indices in imat (for error handling)
void mmlib_to_IS(const raw::integer_dense& imat, Integer max_val, std::shared_ptr<petsc_IS>& is); 

/// create shell preconditioners from linear operator or linsolve object
void linsolve_to_pcshell(const linsolve_obj& ls, PC precod, std::shared_ptr<petsc_pcshell>& pc);
void linop_to_pcshell(const linear_operator& ls, PC precod, std::shared_ptr<petsc_pcshell>& pc);

/// convert PETSC Mat object to mmlib
Matrix petsc_mat_to_mmlib(Mat mat);

}}
