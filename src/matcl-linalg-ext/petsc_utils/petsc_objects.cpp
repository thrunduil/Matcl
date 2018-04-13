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

#if 0

TODO

#include "petsc_objects.h"
#include "matcl-matrep/matcl_matrep.h"
#include "matcl-linalg/general/mmlib_petsc_exception.h"
#include "matcl-linalg/special_matrices/matrix_functors.h"
#include "matcl-linalg/special_matrices/struct_flag_linalg.h"
#include "matcl-matrep/details/extract_type_switch.h"
#include "matcl-linalg/linear_eq/ksp_solver.h"

#pragma warning(push)
#pragma warning(disable: 4127)
#pragma warning(disable: 4100)

#include "petsc/../src/mat/impls/aij/seq/aij.h"
#include "petsc/../src/mat/impls/dense/seq/dense.h"

#pragma warning(pop)

//-----------------------------------------------------
//              utils
//-----------------------------------------------------
namespace matcl { namespace petsc { namespace details
{

struct vec_array_workspace
{
    public:
        vec_array_workspace(const Vec& vec_in);
        ~vec_array_workspace();

    public:
        PetscReal*  m_ptr;

    private:
        const Vec&  m_vec;
};


#pragma warning(push)
#pragma warning(disable : 4702) //unreachable code

vec_array_workspace::vec_array_workspace(const Vec& vec_in)
    : m_vec(vec_in), m_ptr(nullptr)
{
    try
    {
        VecGetArray(m_vec, &m_ptr);
    }
    catch(...)
    {
        if(m_ptr) 
            VecRestoreArray(m_vec, &m_ptr);

        // do not use rethrow here
        // exception propagation through PETSC is suspicious
        throw error::internal_petsc_error("VecGetArray failed");
    }
}

vec_array_workspace::~vec_array_workspace()
{
    try
    {
        VecRestoreArray(m_vec, &m_ptr);
    }
    catch(...) 
    {
        // destructors should not emit exceptions, and VecRestoreArray probably can. Better leak than crash
    } 
}

#pragma warning(pop)

}}};

namespace matcl { namespace petsc
{

//TODO: nonreal

//-----------------------------------------------------
//              conversions
//-----------------------------------------------------
Matrix petsc::petsc_to_mmlib(const Vec& X)
{
    PetscInt size;
    VecGetSize(X, &size);

    details::vec_array_workspace workspace(X);
    Matrix res = (size == 1 ? *(workspace.m_ptr) : make_real_dense(size, 1, workspace.m_ptr));
    return res;
};

void petsc::mmlib_to_petsc(const Matrix& source, Vec& dest)
{
    PetscInt size;
    VecGetSize(dest, &size);
    details::vec_array_workspace workspace(dest);

    dense_matrix<Real> sourceRep(source);
    const Real* ptr = sourceRep.ptr();

    if (sourceRep.cols() == 1)
    {
        if (size != sourceRep.rows())
        {            
            throw error::petsc_size_error(sourceRep.rows(), size);
        }
        for (Integer i = 0; i < sourceRep.rows(); i++)
        {
            workspace.m_ptr[i] = ptr[i];
        }
    }
    else if (sourceRep.rows() == 1)
    {
        if (size != sourceRep.cols())
        {            
            throw error::petsc_size_error(sourceRep.cols(), size);
        }

        for (Integer i = 0; i < sourceRep.cols(); i++)
        {
            workspace.m_ptr[i] = *ptr;
            ptr += sourceRep.ld();
        }
    }
    else 
        throw error::invalid_argument();
};

//-----------------------------------------------------
//              petsc objects
//-----------------------------------------------------
smart_vec smart_vec::vec_create_seq_with_array(MPI_Comm comm, PetscInt bs, PetscInt n,
                                                      const PetscScalar array[])
{
    smart_vec v;
    Vec help;
    
    PetscErrorCode ierr = ::VecCreateSeqWithArray(comm,bs,n,array,&help);

    if (ierr)
        throw error::unable_to_create_petsc_object("Vec");

    v.m_v.reset(help,deleter);
    return v;
}

smart_vec smart_vec::vec_create_seq(MPI_Comm comm,PetscInt n)
{
    smart_vec v;
    Vec help;
    PetscErrorCode ierr = ::VecCreateSeq(comm,n,&help);

    if (ierr)
        throw error::unable_to_create_petsc_object("Vec");

    v.m_v.reset(help,deleter);
    return v;
}

void smart_vec::deleter(Vec& to_destroy)
{
    try
    {
        VecDestroy(&to_destroy);
    }
    catch(...)
    {
        // cleanup must not leak exceptions, better leak memory
    } 
}

smart_mat smart_mat::create(MPI_Comm comm)
{
    smart_mat v;
    Mat help;
    PetscErrorCode ierr = ::MatCreate(comm,&help);

    if (ierr)
        throw error::unable_to_create_petsc_object("Mat");

    v.m_impl.reset(help,deleter);
    return v;
};

void smart_mat::deleter(Mat& to_destroy)
{
    try
    {
        MatDestroy(&to_destroy);
    }
    catch(...)
    {
        // cleanup must not leak exceptions, better leak memory
    }
}

smart_mat::operator Mat() const
{
    return m_impl.get();
};

smart_mat smart_mat::mat_create_seq_dense(MPI_Comm comm,PetscInt m,PetscInt n, PetscInt ld, PetscScalar *data)
{
    smart_mat v;
    Mat impl;
    PetscErrorCode ierr = ::MatCreateSeqDense(comm,m,n,data,&impl);

    if (ierr)
        throw error::unable_to_create_petsc_object("Mat seqdense");

    ::MatSeqDenseSetLDA(impl, ld);
    v.m_impl.reset(impl,deleter);
    return v;
}

void smart_mat::set_struct_flags(const Matrix& A_in)
{
    if (matcl::is_posdef(A_in.get_struct()))
        ::MatSetOption(m_impl.get(), MAT_SPD, PETSC_TRUE);

    bool is_real    = matcl::matrix_traits::is_real(A_in.get_value_code());
    bool is_her     = A_in.get_struct().is_hermitian(A_in.is_square(), is_real);
    bool is_sym     = A_in.get_struct().is_symmetric(A_in.is_square(), is_real);

    if (is_her && is_real || is_sym)
        ::MatSetOption(m_impl.get(), MAT_SYMMETRIC, PETSC_TRUE);
    else if (is_her)
        ::MatSetOption(m_impl.get(), MAT_HERMITIAN, PETSC_TRUE);
};

smart_mat smart_mat::mat_create_shell(MPI_Comm comm, const linear_operator& op)
{    
    smart_mat v;
    v.m_linop = linop_ptr(new linear_operator(op));

    linear_operator* ctx    = v.m_linop.get();

    Mat impl;
    PetscErrorCode ierr = ::MatCreateShell(comm, op.rows(), op.cols(), op.rows(), op.cols(), ctx, &impl);

    if (ierr)
        throw error::unable_to_create_petsc_object("Mat shell");

    v.m_impl.reset(impl,deleter);

    ::MatShellSetOperation(impl, MATOP_MULT, (void(*)(void)) smart_mat::mat_mult);
    ::MatShellSetOperation(impl, MATOP_MULT_TRANSPOSE, (void(*)(void)) smart_mat::mat_mult_trans);
    ::MatShellSetOperation(impl, MATOP_MULT_HERMITIAN_TRANSPOSE, (void(*)(void)) smart_mat::mat_mult_ctrans);
    
    return v;
}
PetscErrorCode smart_mat::mat_mult_trans_impl(Mat shell, Vec T,Vec AT, trans_type t)
{
    linear_operator* linop;
    ::MatShellGetContext(shell, (void**) &linop);

    PetscInt size_in;
    PetscInt size_out;
    ::VecGetSize(T, &size_in);
    ::VecGetSize(AT, &size_out);

    details::vec_array_workspace data_in(T);
    details::vec_array_workspace data_out(AT);

    Matrix mat_in   = make_dense_foreign(size_in, 1, data_in.m_ptr, size_in);
    Matrix mat_out  = make_dense_foreign(size_out, 1, data_out.m_ptr, size_out);

    linop->mmul_right(mat_in, t, mat_out);
    
    return 0;
};

PetscErrorCode smart_mat::mat_mult(Mat shell, Vec T,Vec AT)
{
    return mat_mult_trans_impl(shell, T, AT, trans_type::no_trans);
};
PetscErrorCode smart_mat::mat_mult_trans(Mat shell, Vec T,Vec AT)
{
    return mat_mult_trans_impl(shell, T, AT, trans_type::trans);
};
PetscErrorCode smart_mat::mat_mult_ctrans(Mat shell, Vec T,Vec AT)
{
    return mat_mult_trans_impl(shell, T, AT, trans_type::conj_trans);
};

smart_mat smart_mat::mat_create_seq_aij_with_arrays(MPI_Comm comm,PetscInt m,PetscInt n,
                                                    PetscInt *i,PetscInt *j,PetscScalar *a)
{
    smart_mat v;
    Mat impl;

    PetscErrorCode ierr;

    if (i == nullptr || j == nullptr || a == nullptr)
    {
        //it means that nnz = 0
        ierr = ::MatCreateSeqAIJFromTriple(PETSC_COMM_SELF, m, n, nullptr, nullptr, nullptr, &impl,
                                           0, PETSC_FALSE);
    }
    else
    {
        ierr = ::MatCreateSeqAIJWithArrays(comm,m,n,i,j,a,&impl);
    }

    if (ierr)
        throw error::unable_to_create_petsc_object("Mat seqaij");

    v.m_impl.reset(impl,deleter);
    return v;
}

smart_KSP smart_KSP::create(MPI_Comm comm)
{
    smart_KSP v;
    KSP help;
    PetscErrorCode ierr = ::KSPCreate(comm,&help);

    if (ierr)
        throw error::unable_to_create_petsc_object("KSP");

    v.m_impl.reset(help,destroy);
    return v;
}

smart_KSP::operator KSP() const
{
    return m_impl.get();
};

void smart_KSP::destroy(KSP& to_destroy)
{
    try
    {
        KSPDestroy(&to_destroy);
    }
    catch(...)
    {
        // cleanup must not leak exceptions, better leak memory
    }     
}

smart_vec::operator Vec() const
{ 
    return m_v.get();
};

//----------------------------------------------------------
//              petsc_null_mat
//----------------------------------------------------------
petsc_null_mat::~petsc_null_mat()
{};

std::shared_ptr<petsc_null_mat> petsc_null_mat::create(const Matrix& source, bool with_const)
{
    std::shared_ptr<petsc_null_mat> ret(new petsc_null_mat());
    ret->m_data         = source;

    using Mat_D         = raw::Matrix<Real, struct_dense>;
    const Mat_D& mat    = source.impl<Mat_D>();

    Integer M           = mat.rows();
    Integer N           = mat.cols();
    Integer ld          = mat.ld();

    const Real* ptr     = mat.ptr();

    ret->m_smart_vectors.reserve(N);
    ret->m_vectors.resize(N);

    for (Integer i = 0; i < N; ++i)
    {
        smart_vec tmp(smart_vec::vec_create_seq_with_array(PETSC_COMM_SELF, 1, M, ptr + i * ld));
        ret->m_smart_vectors.push_back(tmp);
        ret->m_vectors[i]   = tmp;
    };

    MatNullSpace impl;
    PetscErrorCode ierr;

    if (with_const == false)
        ierr = ::MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_FALSE, N, ret->m_vectors.data(), &impl);
    else
        ierr = ::MatNullSpaceCreate(PETSC_COMM_SELF, PETSC_TRUE, N, ret->m_vectors.data(), &impl);

    if (ierr)
        throw error::unable_to_create_petsc_object("MatNullSpace");

    ret->m_impl.reset(impl, deleter);

    return ret;
}

petsc_null_mat::operator MatNullSpace() const
{
    return m_impl.get();
};

void petsc_null_mat::deleter(MatNullSpace& to_destroy)
{
    try
    {
        MatNullSpaceDestroy(&to_destroy);
    }
    catch(...)
    {
        // cleanup must not leak exceptions, better leak memory
    }     
}

//----------------------------------------------------------
//              petsc_mat_struct
//----------------------------------------------------------
petsc_sparse_struct::~petsc_sparse_struct()
{};

petsc_sparse_struct::operator Mat() const
{
    return m_impl.get();
};

bool petsc_sparse_struct::need_transpose() const
{
    return m_need_trans;
};

void petsc_sparse_struct::deleter(Mat& to_destroy)
{
    try
    {
        MatDestroy(&to_destroy);
    }
    catch(...)
    {
        // cleanup must not leak exceptions, better leak memory
    }     
}

petsc_sparse_struct::petsc_sparse_struct(const Matrix& mat, bool const_access)
    :m_need_trans(false), m_data(sparse(mat)), m_column_indices(ti::ti_empty())
{
    initialize_sparse(const_access);
};

struct visit_sparse : public md::extract_type_switch<void, visit_sparse,true>
{
    using base_type = md::extract_type_switch<void, visit_sparse,true>;

    template<class V, class ... Arg>
    static void eval(const matcl::Matrix& h, const raw::Matrix<V,struct_sparse>& mat, 
                    Integer*& ptr_c, Integer*& ptr_r, matcl::Matrix& data)
    {
        data = h;

        using Mat_S = raw::Matrix<V, struct_sparse>;

        Mat_S tmp(mat);

        ptr_c   = tmp.rep().ptr_c();
        ptr_r   = tmp.rep().ptr_r();
    };

    template<class Mat, class ... Arg>
    static void eval(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };

    template<class Mat, class ... Arg>
    static void eval_scalar(const matcl::Matrix& h, const Mat&, Arg&& ... args)
    {
        return base_type::make<const Matrix&>(sparse(h), std::forward<Arg>(args)...);
    };
};

void petsc_sparse_struct::initialize_sparse(bool const_access)
{
    //implicit transposition (petsc is CSR)
    Integer M   = m_data.cols();
    Integer N   = m_data.rows();

    if (const_access == false)
        m_data.make_unique();

    Integer* ptr_c;
    Integer* ptr_r;

    visit_sparse::make<const Matrix&>(m_data, ptr_c, ptr_r, m_data);

    if (ptr_c == nullptr)
    {
        //it means that nnz = 0

        Mat impl;
        PetscErrorCode ierr;
        ierr = ::MatCreateSeqAIJFromTriple(PETSC_COMM_SELF, M, N, nullptr, nullptr, nullptr, &impl,
                                           0, PETSC_FALSE);

        if (ierr)
            throw error::unable_to_create_petsc_object("Mat seqaij");

        m_impl.reset(impl,deleter);
        m_need_trans = true;
        return;
    };

    Integer off     = ptr_c[0];
    Real* ptr_x     = nullptr;

    if (off == 0)
    {        
        Integer* ptr_i  = ptr_c;
        Integer* ptr_j  = ptr_r;

        Mat impl;
        PetscErrorCode ierr;
        ierr = ::MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, M, N, ptr_i, ptr_j, ptr_x, &impl);

        if (ierr)
            throw error::unable_to_create_petsc_object("Mat seqaij");

        m_impl.reset(impl,deleter);
        m_need_trans = true;
        return;
    };

    //Petsc requires that column indices starts from 0
    //offset can be nonzero if one creates column view of
    //sparse matrix and assign sym structure

    m_column_indices.assign_to_fresh(Mat_I(ti::ti_empty(), M, 1));

    Integer* ptr_new    = m_column_indices.ptr();

    for (Integer i = 0; i < M; ++i)
        ptr_new[i]  = ptr_c[i] - off;

    Integer* ptr_i  = ptr_new;
    Integer* ptr_j  = ptr_r;

    Mat impl;
    PetscErrorCode ierr;
    ierr = ::MatCreateSeqAIJWithArrays(PETSC_COMM_SELF, M, N, ptr_i, ptr_j, ptr_x, &impl);

    if (ierr)
        throw error::unable_to_create_petsc_object("Mat seqij");

    m_impl.reset(impl,deleter);
    m_need_trans = true;
    return;
};

petsc_sparse_struct::this_ptr petsc_sparse_struct::create(const Matrix& mat, bool const_access)
{
    std::shared_ptr<petsc_sparse_struct> ret(new petsc_sparse_struct(mat,const_access));
    return ret;
};

//----------------------------------------------------------
//              petsc_pcshell
//----------------------------------------------------------

petsc_pcshell::~petsc_pcshell()
{};

petsc_pcshell::operator PC() const
{
    return m_impl.get();
};
void petsc_pcshell::deleter(PC& to_destroy)
{
    try
    {
        PCDestroy(&to_destroy);
    }
    catch(...)
    {
        // cleanup must not leak exceptions, better leak memory
    }     
};

petsc_pcshell::petsc_pcshell(const linsolve_obj& ls)
    :m_data_linsolve(ls)
{};
petsc_pcshell::petsc_pcshell(const linear_operator& ls)
    :m_data_linop(ls)
{};

petsc_pcshell::this_ptr petsc_pcshell::create(PC pc, const linsolve_obj& ls)
{
    this_ptr ret(new petsc_pcshell(ls));

    linsolve_obj* ctx   = &ret->m_data_linsolve;

    ret->m_impl.reset(pc,deleter);
    PetscObjectReference((PetscObject)pc);

    ::PCSetType(pc,PCSHELL);
    ::PCShellSetContext(pc,ctx);
    ::PCShellSetApply(pc, func_apply_linsolve);
    ::PCShellSetApplyTranspose(pc, func_apply_trans_linsolve);

    return ret;
};
petsc_pcshell::this_ptr petsc_pcshell::create(PC pc, const linear_operator& ls)
{
    this_ptr ret(new petsc_pcshell(ls));

    linear_operator* ctx   = &ret->m_data_linop;

    ret->m_impl.reset(pc,deleter);
    PetscObjectReference((PetscObject)pc);

    ::PCSetType(pc,PCSHELL);
    ::PCShellSetContext(pc,ctx);
    ::PCShellSetApply(pc, func_apply_linop);
    ::PCShellSetApplyTranspose(pc, func_apply_trans_linop);

    return ret;
};

PetscErrorCode petsc_pcshell::func_apply_trans_linsolve_impl(PC shell, Vec T,Vec AT, trans_type t)
{
    linsolve_obj* linop;
    ::PCShellGetContext(shell, (void**) &linop);

    PetscInt size_in;
    PetscInt size_out;
    ::VecGetSize(T, &size_in);
    ::VecGetSize(AT, &size_out);

    details::vec_array_workspace data_in(T);
    details::vec_array_workspace data_out(AT);

    Matrix mat_in   = make_dense_foreign(size_in, 1, data_in.m_ptr, size_in);
    Matrix mat_out  = make_dense_foreign(size_out, 1, data_out.m_ptr, size_out);
    mat_out(colon()) = linop->solve(mat_in, t);
    
    return 0;
}
PetscErrorCode petsc_pcshell::func_apply_trans_linop_impl(PC shell, Vec T,Vec AT, trans_type t)
{
    linear_operator* linop;
    ::PCShellGetContext(shell, (void**) &linop);

    PetscInt size_in;
    PetscInt size_out;
    ::VecGetSize(T, &size_in);
    ::VecGetSize(AT, &size_out);

    details::vec_array_workspace data_in(T);
    details::vec_array_workspace data_out(AT);

    Matrix mat_in       = make_dense_foreign(size_in, 1, data_in.m_ptr, size_in);
    Matrix mat_out      = make_dense_foreign(size_out, 1, data_out.m_ptr, size_out);
    mat_out(colon())    = linop->mmul_right(mat_in, t);
    
    return 0;
}

PetscErrorCode petsc_pcshell::func_apply_linsolve(PC shell, Vec T,Vec AT)
{
    return func_apply_trans_linsolve_impl(shell, T, AT, trans_type::no_trans);
}
PetscErrorCode petsc_pcshell::func_apply_linop(PC shell, Vec T,Vec AT)
{
    return func_apply_trans_linop_impl(shell, T, AT, trans_type::no_trans);
}
PetscErrorCode petsc_pcshell::func_apply_trans_linsolve(PC shell, Vec T,Vec AT)
{
    return func_apply_trans_linsolve_impl(shell, T, AT, trans_type::conj_trans);
}
PetscErrorCode petsc_pcshell::func_apply_trans_linop(PC shell, Vec T,Vec AT)
{
    return func_apply_trans_linop_impl(shell, T, AT, trans_type::conj_trans);
}

//----------------------------------------------------------
//              petsc_IS
//----------------------------------------------------------
petsc_IS::~petsc_IS()
{};
        
petsc_IS::operator IS() const
{
    return m_impl.get();
};
void petsc_IS::deleter(IS& to_destroy)
{
    try
    {
        ISDestroy(&to_destroy);
    }
    catch(...)
    {
        // cleanup must not leak exceptions, better leak memory
    }     
}
petsc_IS::this_ptr petsc_IS::create(const raw::integer_dense& mat, Integer max_val)
{
    return this_ptr(new petsc_IS(mat, max_val));
}

petsc_IS::petsc_IS(const raw::integer_dense& mat, Integer max_val)
{
    IS ret;
    PetscCopyMode mode  = PETSC_COPY_VALUES;
    Integer N           = mat.length();
    Mat_I mat_ind       = mat.copy().make_explicit();
    Integer* ptr_ind    = mat_ind.ptr();

    // Petsc indices are 0-based
    for (Integer i = 0; i < N; ++i)
    {
        Integer tmp     = ptr_ind[i];

        if (tmp < 1 || tmp > max_val)
            throw error::invalid_single_index(tmp, max_val);

        ptr_ind[i]      = tmp - 1;
    };

    PetscErrorCode ierr = ::ISCreateGeneral(PETSC_COMM_SELF, N, mat_ind.ptr(), mode, &ret);

    if (ierr)
        throw error::unable_to_create_petsc_object("IS");

    m_impl.reset(ret, &deleter);
};


//----------------------------------------------------------
//              conversions
//----------------------------------------------------------
void petsc::mmlib_to_petsc_null(const matcl::Matrix& source, std::shared_ptr<petsc_null_mat>& dest, 
                                bool with_const)
{
    dest = petsc_null_mat::create(source, with_const);
};

void petsc::mmlib_to_sparse_struct(const matcl::Matrix& source, bool const_access, 
                                std::shared_ptr<petsc_sparse_struct>& dest)
{
    dest = petsc_sparse_struct::create(source, const_access);
};

void petsc::IS_to_perm_vector(IS p, permvec& mp)
{
    Matrix ip;
    IS_to_int_mat(p, ip);

    mp = permvec::from_matrix(ip);
}

void petsc::IS_to_int_mat(IS p, Matrix& mp)
{
    const PetscInt * indices;
    PetscInt size;

    ::ISGetIndices(p, &indices);
    ::ISGetSize(p, &size);
    ::ISRestoreIndices(p, &indices);

    Matrix ip       = make_integer_dense(size, 1, indices);
    Integer* arr    = ip.get_array_unique<Integer>();

    // change to 1-based
    for (Integer i = 0; i < size; ++i)
        arr[i]      += 1;

    mp = ip;
}

void petsc::mmlib_to_IS(const raw::integer_dense& imat, Integer max_val, std::shared_ptr<petsc_IS>& is)
{
    is = petsc_IS::create(imat, max_val);
}

void petsc::linsolve_to_pcshell(const linsolve_obj& ls, PC old_pc, std::shared_ptr<petsc_pcshell>& pc)
{
    pc = petsc_pcshell::create(old_pc, ls);
}

void petsc::linop_to_pcshell(const linear_operator& ls, PC old_pc, std::shared_ptr<petsc_pcshell>& pc)
{
    pc = petsc_pcshell::create(old_pc, ls);
}

Matrix petsc_sparse_mat_to_mmlib(Mat mat)
{
    ::Mat_SeqAIJ* aij = (::Mat_SeqAIJ*)mat->data;

    ::PetscInt m, n;
    ::MatGetSize(mat, &m, &n);

    const PetscInt* ptr_prow    = aij->i;   // pointer to beginning of each row
    const PetscInt* ptr_pcol    = aij->j;   // column values: j + i[k] - 1 is start of row k
    const Real* ptr_px          = aij->a;   // nonzero elements

    Integer nz      = ptr_prow[m];

    Matrix ret      = make_real_sparse_noinit(n, m, nz);

    using Mat_S     = raw::Matrix<Real, struct_sparse>;
    Mat_S mat_rep   = ret.get_impl_unique<Mat_S>();

    auto d          = mat_rep.rep();

    Integer* ptr_r  = d.ptr_r();
    Integer* ptr_c  = d.ptr_c();
    Real* ptr_x     = d.ptr_x();

    for (Integer i = 0; i <= m; ++i)
        ptr_c[i]    = ptr_prow[i];

    for (Integer i = 0; i < nz; ++i)
        ptr_r[i]    = ptr_pcol[i];

    for (Integer i = 0; i < nz; ++i)
        ptr_x[i]    = ptr_px[i];

    return ctrans(ret);
};

Matrix petsc_dense_mat_to_mmlib(Mat mat)
{
    ::Mat_SeqDense* aij  = (::Mat_SeqDense*)mat->data;

    ::PetscInt m, n;
    ::MatGetSize(mat, &m, &n);

    const Real* ptr_px  = aij->v;  
    Integer ld          = aij->lda;

    Real* ptr;
    Matrix ret          = make_real_dense_noinit(m, n, ptr);
    Integer ret_ld      = m;

    for (Integer j = 0; j < n; ++j)
    {
        for (Integer i = 0; i < m; ++i)
            ptr[i]      = ptr_px[i];

        ptr_px          += ld;
        ptr             += ret_ld;
    };

    return ret;
};

Matrix petsc::petsc_mat_to_mmlib(Mat mat)
{
    if (!mat)
        throw error::null_petsc_matrix();

    MatType mat_type;
    ::MatGetType(mat, &mat_type);

    std::string str_mat_type    = std::string(mat_type ? mat_type : "");

    if (str_mat_type == std::string(MATSEQAIJ))
        return petsc_sparse_mat_to_mmlib(mat);
    else if (str_mat_type == std::string(MATSEQDENSE))
        return petsc_dense_mat_to_mmlib(mat);
    else if (str_mat_type == std::string(MATSHELL))
        throw error::unable_convert_to_mmlib_shell();
    else
        throw error::unable_convert_to_mmlib_unsupported_format(str_mat_type);
            
    //MatGetSize(Mat,PetscInt*,PetscInt*);
};

}}

#if 0

/**
 * Class to wrap around a Petsc C-object
 *
 * Handles the underlying object acording to RAII principle and should be exception safe
 * See Petsc documentation for info on particular "Create" functions
 */    
class SmartSNES
{
public:
    typedef boost::shared_ptr<SmartSNES> Ptr;
    operator SNES() const {return v.get();};
    static SmartSNES SNESCreate(MPI_Comm comm)
    {
        SmartSNES v;
        SNES help;
        ::SNESCreate(comm,&help);
        v.v.reset(help,destroy);
        return v;
    }
private:
    boost::shared_ptr<_p_SNES> v;
    static void destroy(SNES& to_destroy)
    {
        try
        {
            SNESDestroy(&to_destroy);
        }
        catch(...)
        {
            assert(0 && "Petsc ?Destroy function failed - leak");
        } // cleanup must not leak exceptions, better leak memory
    }
};
#endif

#endif